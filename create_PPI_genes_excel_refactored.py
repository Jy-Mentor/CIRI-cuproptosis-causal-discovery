#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PPI网络输入基因Excel报告生成器（重构版）
使用ExcelReportBuilder工具类，消除重复代码

原文件: create_PPI_genes_excel.py
重构日期: 2025-01-24
"""

import os
import sys

# 添加utils目录到路径
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'utils'))

from excel_report_builder import ExcelReportBuilder


def read_gene_mapping(ppi_file: str, mapping_file: str) -> tuple:
    """
    读取PPI基因并映射到人类基因
    
    Args:
        ppi_file: PPI输入基因文件路径
        mapping_file: 基因映射库文件路径
        
    Returns:
        (mapped_results, unmapped_genes, df_unique) 元组
    """
    # 读取PPI输入基因（大鼠基因）
    with open(ppi_file, 'r') as f:
        rat_genes = [line.strip().upper() for line in f if line.strip()]
    
    print(f"PPI输入基因数（大鼠）: {len(rat_genes)}")
    
    # 读取映射库，建立大鼠→人类映射
    mapping_lines = open(mapping_file, 'r').readlines()
    header_line = None
    for i, line in enumerate(mapping_lines):
        if line.startswith('RAT_GENE_SYMBOL'):
            header_line = i
            break
    
    # 解析映射文件
    rat_to_human = {}
    for line in mapping_lines[header_line+1:]:
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            rat_gene = parts[0].strip().upper()
            human_ortholog = parts[1].strip().upper()
            
            if rat_gene and human_ortholog and human_ortholog != 'N/A':
                # 处理多个人类同源基因（用|分隔）
                human_genes = [g.strip() for g in human_ortholog.split('|') if g.strip()]
                if rat_gene not in rat_to_human:
                    rat_to_human[rat_gene] = []
                rat_to_human[rat_gene].extend(human_genes)
    
    # 大鼠基因映射到人类基因
    results = []
    mapped_count = 0
    unmapped_genes = []
    
    for rat_gene in rat_genes:
        if rat_gene in rat_to_human:
            human_genes = list(set(rat_to_human[rat_gene]))  # 去重
            for human_gene in human_genes:
                results.append({
                    'Rat_Gene': rat_gene,
                    'Human_Gene': human_gene,
                    'Mapping_Status': 'Mapped'
                })
            mapped_count += 1
        else:
            results.append({
                'Rat_Gene': rat_gene,
                'Human_Gene': 'N/A',
                'Mapping_Status': 'Unmapped'
            })
            unmapped_genes.append(rat_gene)
    
    # 去重（保留唯一的Human_Gene）
    df_unique = [r for r in results if r['Mapping_Status'] == 'Mapped']
    seen = set()
    df_unique_clean = []
    for r in df_unique:
        if r['Human_Gene'] not in seen:
            seen.add(r['Human_Gene'])
            df_unique_clean.append(r)
    
    print(f"成功映射基因数: {mapped_count}")
    print(f"未映射基因数: {len(unmapped_genes)}")
    print(f"去重后人类基因数: {len(df_unique_clean)}")
    
    return df_unique_clean, unmapped_genes, results


def create_ppi_genes_excel():
    """
    创建PPI基因映射Excel报告（重构版）
    使用ExcelReportBuilder工具类简化代码
    """
    # 文件路径配置
    ppi_file = "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/GSE61616_cluster_ssGSEA_PPI_results_v2/String_PPI_input_genes.txt"
    mapping_file = "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/大鼠 小鼠 人类映射库.txt"
    output_file = "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/PPI_Input_Genes_Human.xlsx"
    
    # 读取和映射基因数据
    df_unique, unmapped_genes, all_results = read_gene_mapping(ppi_file, mapping_file)
    
    # 创建报告构建器（一行代码完成所有初始化）
    builder = ExcelReportBuilder(output_file, theme='default')
    
    # ========== Sheet 1: PPI输入基因（已映射） ==========
    ws1 = builder.create_sheet("PPI输入基因 (Human)")
    
    # 添加标题（简洁的API）
    builder.add_title(ws1, "PPI网络输入基因 (转换为人类基因符号)", 
                     row=1, col=1, merge_range='A1:D1')
    
    # 添加统计信息（使用add_table简化数据写入）
    stats_data = [
        [f"原始大鼠基因数: {len([r for r in all_results]) // 2}"],  # 估算
        [f"成功映射: {len(df_unique)}"],
        [f"未映射: {len(unmapped_genes)}"],
        [f"去重后人类基因数: {len(df_unique)}"]
    ]
    builder.add_table(ws1, stats_data, start_row=3, start_col=1)
    
    # 添加表头和数据
    headers = ['序号', '大鼠基因', '人类基因', '映射状态']
    data_rows = []
    for idx, row in enumerate(df_unique, 1):
        data_rows.append([idx, row['Rat_Gene'], row['Human_Gene'], row['Mapping_Status']])
    
    # 合并表头和数据，使用add_table一次性写入
    all_data = [headers] + data_rows
    builder.add_table(ws1, all_data, start_row=8, has_header=True)
    
    # 设置列宽（批量设置）
    builder.set_column_widths(ws1, {
        'A': 10,
        'B': 15,
        'C': 15,
        'D': 12
    })
    
    # ========== Sheet 2: 未映射基因（如果存在） ==========
    if unmapped_genes:
        ws2 = builder.create_sheet("未映射基因")
        builder.add_title(ws2, "未映射的大鼠基因", row=1, col=1)
        
        ws2['A3'] = "以下基因在映射库中未找到对应的人类同源基因:"
        
        # 添加未映射基因列表
        unmapped_data = [[gene] for gene in unmapped_genes]
        builder.add_table(ws2, unmapped_data, start_row=4)
        
        builder.set_column_width(ws2, 'A', 20)
    
    # 保存文件
    builder.save()
    
    print(f"\n✅ Excel文件已保存: {output_file}")
    print(f"包含 {len(df_unique)} 个人类基因，每行一个基因")
    print(f"工作表列表: {builder.get_sheetnames()}")


if __name__ == "__main__":
    create_ppi_genes_excel()
