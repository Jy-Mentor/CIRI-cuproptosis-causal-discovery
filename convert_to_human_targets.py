#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
将大鼠/小鼠DEG基因取并集并转换为人类同源基因
输出：human_targets_union.xlsx
"""

import pandas as pd
import requests
import time
from collections import defaultdict

# 读取Excel
input_file = r'c:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\temporal_degs_only.xlsx'
output_file = r'c:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\human_targets_union.xlsx'

print(f'读取: {input_file}')

# 读取所有sheet
sheets = ['3hr', '6hr', '12hr', '24hr', 'GSE61616_7d', 'GSE97537_24H']
all_genes = set()
sheet_genes = {}

for sheet in sheets:
    df = pd.read_excel(input_file, sheet_name=sheet)
    genes = set(df['gene_symbol'].dropna().tolist())
    sheet_genes[sheet] = genes
    all_genes.update(genes)
    print(f'  {sheet}: {len(genes)} 个基因')

print(f'\n并集基因总数: {len(all_genes)}')

# 分离小鼠/大鼠特有基因和可能的人类同源基因
# 小鼠基因通常以大小写混合，大鼠基因也是
# 我们需要使用BioMart或NCBI进行转换

# 首先尝试使用mygene.info API进行转换
print('\n正在转换基因...')

def batch_convert_genes(gene_list, species='mouse'):
    """使用mygene.info批量转换基因"""
    url = 'https://mygene.info/v3/query'
    
    # 分批处理，每批100个
    batch_size = 100
    all_results = {}
    
    for i in range(0, len(gene_list), batch_size):
        batch = gene_list[i:i+batch_size]
        
        # 构建查询
        query = ','.join(batch)
        params = {
            'q': query,
            'scopes': 'symbol',
            'species': species,
            'fields': 'symbol,homologene,genomic_pos',
            'dotfield': True
        }
        
        try:
            response = requests.post(url, params=params, timeout=30)
            if response.status_code == 200:
                results = response.json()
                for result in results:
                    if 'symbol' in result:
                        orig_gene = result['query']
                        # 获取人类同源基因
                        if 'homologene' in result and 'genes' in result['homologene']:
                            for gene_info in result['homologene']['genes']:
                                if len(gene_info) >= 3 and gene_info[1] == 9606:  # 人类taxid
                                    human_gene = gene_info[2]
                                    all_results[orig_gene] = human_gene
                                    break
                        else:
                            # 如果没有同源基因信息，保留原基因名（大写）
                            all_results[orig_gene] = orig_gene.upper()
            else:
                print(f'  API错误: {response.status_code}')
                
        except Exception as e:
            print(f'  请求错误: {e}')
        
        time.sleep(0.5)  # 避免请求过快
    
    return all_results

# 将基因分为小鼠和大鼠
# GSE104036是小鼠，GSE61616和GSE97537是大鼠
mouse_genes = set()
rat_genes = set()

for sheet, genes in sheet_genes.items():
    if 'GSE61616' in sheet or 'GSE97537' in sheet:
        rat_genes.update(genes)
    else:
        mouse_genes.update(genes)

print(f'\n小鼠基因数: {len(mouse_genes)}')
print(f'大鼠基因数: {len(rat_genes)}')

# 转换小鼠基因
print('\n转换小鼠基因...')
mouse_list = list(mouse_genes)
mouse_to_human = batch_convert_genes(mouse_list, 'mouse')

# 转换大鼠基因
print('转换大鼠基因...')
rat_list = list(rat_genes)
rat_to_human = batch_convert_genes(rat_list, 'rat')

# 合并所有人类基因
human_targets = set()
unmapped = []

# 处理小鼠基因
for gene in mouse_genes:
    if gene in mouse_to_human:
        human_targets.add(mouse_to_human[gene])
    else:
        unmapped.append(('mouse', gene))

# 处理大鼠基因
for gene in rat_genes:
    if gene in rat_to_human:
        human_targets.add(rat_to_human[gene])
    else:
        unmapped.append(('rat', gene))

print(f'\n转换完成:')
print(f'  成功映射: {len(human_targets)} 个人类基因')
print(f'  未映射: {len(unmapped)} 个')

# 创建输出DataFrame
human_df = pd.DataFrame({
    'human_gene_symbol': sorted(list(human_targets))
})

# 创建来源信息
source_info = []
for gene in sorted(list(human_targets)):
    sources = []
    # 检查来自哪些时序
    for sheet, genes in sheet_genes.items():
        # 需要反向查找原始基因
        for orig_gene in genes:
            if orig_gene in mouse_to_human and mouse_to_human[orig_gene] == gene:
                sources.append(sheet)
                break
            elif orig_gene in rat_to_human and rat_to_human[orig_gene] == gene:
                sources.append(sheet)
                break
    source_info.append(','.join(sources) if sources else 'unknown')

human_df['source_timepoints'] = source_info

# 保存到Excel
human_df.to_excel(output_file, index=False)
print(f'\n输出文件: {output_file}')
print(f'  包含 {len(human_df)} 个人类靶点基因')

# 同时保存未映射的基因
if unmapped:
    unmapped_df = pd.DataFrame(unmapped, columns=['species', 'gene_symbol'])
    unmapped_file = r'c:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\unmapped_genes.xlsx'
    unmapped_df.to_excel(unmapped_file, index=False)
    print(f'  未映射基因保存到: {unmapped_file}')

print('\n完成！')
