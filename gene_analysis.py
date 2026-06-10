#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
基因分析脚本：筛选hub基因，KEGG/GO富集分析，STRING分析，铜死亡相关分析
"""

import os
import sys
import pandas as pd
import requests
import json
from collections import Counter

# 基因列表
genes = [
    "TSPO", "AIF1", "CTSD", "TIMP1", "PTPN6", "KIF11", "CTSL", "PTPRC", "SREBF2", "CP",
    "MAN2B1", "GPX1", "TOP2A", "TBXAS1", "RENBP", "CTSB", "SERPINB1", "CASP8", "CTSC", "PTGS1",
    "COL1A1", "CNR2", "FDFT1", "HMOX1", "TGFB1", "IRF1", "B2M", "PTGR1", "STAT1", "ALDH9A1",
    "S100A6", "CCR5", "PLA2G4A", "ICAM1", "LYN", "PARP12", "GFAP", "PABPC1", "SQLE", "FABP5",
    "C3", "CTSS", "HPGDS", "IGKC", "FABP4", "TCN2", "GCH1", "SAT1", "GABRB3", "CYP51A1",
    "HTR2B", "HTR2A", "MAPKAPK2", "CCL2", "ALDH1A1", "CTSK", "CDK4", "CCND1", "ADORA1", "PCTP",
    "STAT3", "FAS", "IGFBP2", "NR1H3", "SPHK1", "MKNK2", "EGR1", "CNDP2", "RHOC", "IL10RA",
    "ITGA1", "ADRA1D", "SYN1", "ACTA1", "CHRM3", "XDH", "SERPINB10", "AOC3", "BST1", "IL6",
    "F3", "PTGES", "S1PR1", "GABRA3", "GALNT10", "LSS"
]

# 铜死亡相关基因（基于已有文献）
cuproptosis_related_genes = [
    "FDX1", "LIAS", "LIPT1", "DLAT", "PDHB", "DLD", "ATP7A", "ATP7B", "SLC31A1", "ATOX1",
    "MTF1", "GLS", "GLUD1", "GLUL", "CDKN2A", "TP53", "NFE2L2", "KEAP1", "HIF1A", "MT2A"
]


def get_gene_info(gene):
    """
    获取基因基本信息
    """
    try:
        url = f"https://rest.genenames.org/fetch/symbol/{gene}"
        headers = {"Accept": "application/json"}
        response = requests.get(url, headers=headers)
        if response.status_code == 200:
            data = response.json()
            if data.get('response', {}).get('docs'):
                return data['response']['docs'][0]
    except Exception as e:
        print(f"获取基因 {gene} 信息失败: {e}")
    return None


def analyze_cuproptosis_relation(genes):
    """
    分析基因与铜死亡的关系
    """
    cuproptosis_genes_in_list = [gene for gene in genes if gene in cuproptosis_related_genes]
    print(f"\n与铜死亡相关的基因 ({len(cuproptosis_genes_in_list)}):")
    for gene in cuproptosis_genes_in_list:
        print(f"  - {gene}")
    
    return cuproptosis_genes_in_list


def perform_string_analysis(genes):
    """
    执行STRING蛋白质相互作用分析
    """
    print("\n执行STRING蛋白质相互作用分析...")
    # 限制基因数量以符合STRING API要求
    if len(genes) > 100:
        genes_subset = genes[:100]
        print(f"基因数量超过100，使用前100个基因进行分析")
    else:
        genes_subset = genes
    
    try:
        url = "https://string-db.org/api/json/network"
        params = {
            "identifiers": "%0d".join(genes_subset),
            "species": 9606,  # 人类
            "required_score": 0.4
        }
        response = requests.get(url, params=params)
        if response.status_code == 200:
            data = response.json()
            interactions = len(data)
            print(f"STRING分析完成，找到 {interactions} 个蛋白质相互作用")
            return data
        else:
            print(f"STRING分析失败: {response.status_code}")
    except Exception as e:
        print(f"STRING分析错误: {e}")
    return None


def get_kegg_enrichment(genes):
    """
    获取KEGG通路富集分析
    """
    print("\n执行KEGG通路富集分析...")
    # 使用WebGestalt API进行富集分析
    try:
        url = "http://www.webgestalt.org/api/genest enrichment"
        data = {
            "organism": "hsapiens",
            "method": "ORA",
            "geneList": genes,
            "enrichmentDatabase": "KEGG_PATHWAY"
        }
        response = requests.post(url, json=data)
        if response.status_code == 200:
            results = response.json()
            print("KEGG富集分析完成")
            return results
        else:
            print(f"KEGG富集分析失败: {response.status_code}")
    except Exception as e:
        print(f"KEGG富集分析错误: {e}")
    return None


def get_go_enrichment(genes):
    """
    获取GO富集分析
    """
    print("\n执行GO富集分析...")
    # 使用WebGestalt API进行富集分析
    try:
        url = "http://www.webgestalt.org/api/genest enrichment"
        data = {
            "organism": "hsapiens",
            "method": "ORA",
            "geneList": genes,
            "enrichmentDatabase": "GO_Biological_Process"
        }
        response = requests.post(url, json=data)
        if response.status_code == 200:
            results = response.json()
            print("GO富集分析完成")
            return results
        else:
            print(f"GO富集分析失败: {response.status_code}")
    except Exception as e:
        print(f"GO富集分析错误: {e}")
    return None


def identify_hub_genes(genes, string_data):
    """
    基于STRING网络识别hub基因
    """
    print("\n识别hub基因...")
    if not string_data:
        print("无法识别hub基因，STRING数据不可用")
        return []
    
    # 计算每个基因的度数（连接数）
    degree_counter = Counter()
    for interaction in string_data:
        degree_counter[interaction.get('preferredName_A')] += 1
        degree_counter[interaction.get('preferredName_B')] += 1
    
    # 排序并选择前10个hub基因
    hub_genes = [gene for gene, _ in degree_counter.most_common(10)]
    print(f"识别到的hub基因: {hub_genes}")
    return hub_genes


def main():
    """
    主函数
    """
    print("基因分析开始...")
    print(f"总基因数: {len(genes)}")
    
    # 1. 分析与铜死亡的关系
    cuproptosis_genes = analyze_cuproptosis_relation(genes)
    
    # 2. 执行STRING分析
    string_data = perform_string_analysis(genes)
    
    # 3. 识别hub基因
    hub_genes = identify_hub_genes(genes, string_data)
    
    # 4. 执行KEGG富集分析
    kegg_results = get_kegg_enrichment(genes)
    
    # 5. 执行GO富集分析
    go_results = get_go_enrichment(genes)
    
    # 6. 保存分析结果
    results = {
        "total_genes": len(genes),
        "cuproptosis_related_genes": cuproptosis_genes,
        "hub_genes": hub_genes,
        "string_interactions_count": len(string_data) if string_data else 0
    }
    
    # 保存为JSON
    output_file_json = "C:\\Users\\Jy-Mentor-7\\Desktop\\你信我 这不是重蹈覆辙\\基因分析结果.json"
    with open(output_file_json, 'w', encoding='utf-8') as f:
        json.dump(results, f, ensure_ascii=False, indent=2)
    
    # 保存为Excel
    output_file_excel = "C:\\Users\\Jy-Mentor-7\\Desktop\\你信我 这不是重蹈覆辙\\基因分析结果.xlsx"
    
    # 创建数据框
    data = []
    for gene in genes:
        is_cuproptosis = "是" if gene in cuproptosis_related_genes else "否"
        is_hub = "是" if gene in hub_genes else "否"
        data.append([gene, is_cuproptosis, is_hub])
    
    df = pd.DataFrame(data, columns=['基因', '是否与铜死亡相关', '是否为hub基因'])
    
    # 创建Excel写入器
    with pd.ExcelWriter(output_file_excel, engine='openpyxl') as writer:
        # 写入基因列表
        df.to_excel(writer, sheet_name='基因列表', index=False)
        
        # 写入分析摘要
        summary_data = [
            ['总基因数', len(genes)],
            ['与铜死亡相关基因数', len(cuproptosis_related_genes)],
            ['Hub基因数', len(hub_genes)],
            ['蛋白质相互作用数', len(string_data) if string_data else 0]
        ]
        summary_df = pd.DataFrame(summary_data, columns=['项目', '数值'])
        summary_df.to_excel(writer, sheet_name='分析摘要', index=False)
        
        # 写入铜死亡相关基因
        if cuproptosis_related_genes:
            cuproptosis_df = pd.DataFrame(cuproptosis_related_genes, columns=['铜死亡相关基因'])
            cuproptosis_df.to_excel(writer, sheet_name='铜死亡相关基因', index=False)
        
        # 写入Hub基因
        if hub_genes:
            hub_df = pd.DataFrame(hub_genes, columns=['Hub基因'])
            hub_df.to_excel(writer, sheet_name='Hub基因', index=False)
    
    print(f"\n分析结果已保存到:")
    print(f"JSON文件: {output_file_json}")
    print(f"Excel文件: {output_file_excel}")
    print("\n分析完成！")


if __name__ == "__main__":
    main()