#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
从 deg_temporal_summary.xlsx 提取各时序DEG基因列表
输出：temporal_degs_only.xlsx (仅基因名，无统计值)
"""

import pandas as pd
import os

# 输入文件
input_file = r'c:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\deg_temporal_summary.xlsx'
output_file = r'c:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\temporal_degs_only.xlsx'

print(f'读取: {input_file}')

# 读取各个时序的DEG sheet
sheets_to_extract = {
    'GSE104036_all_DEGs': 'All_Timepoints',
    'GSE104036_3hr_DEGs': '3hr',
    'GSE104036_6hr_DEGs': '6hr',
    'GSE104036_12hr_DEGs': '12hr',
    'GSE104036_24hr_DEGs': '24hr',
    'GSE61616_DEGs': 'GSE61616_7d',
    'GSE97537_DEGs': 'GSE97537_24H',
}

writer = pd.ExcelWriter(output_file, engine='openpyxl')

for sheet_name, output_name in sheets_to_extract.items():
    try:
        df = pd.read_excel(input_file, sheet_name=sheet_name)
        
        # 只保留 gene_symbol 列，并筛选 significant=True 的
        if 'gene_symbol' in df.columns and 'significant' in df.columns:
            sig_df = df[df['significant'] == True][['gene_symbol']].copy()
        else:
            sig_df = df[['gene_symbol']].copy() if 'gene_symbol' in df.columns else df.iloc[:, [0]].copy()
        
        # 去重并排序
        sig_df = sig_df.drop_duplicates().sort_values('gene_symbol').reset_index(drop=True)
        
        # 写入Excel
        sig_df.to_excel(writer, sheet_name=output_name, index=False)
        print(f'  {output_name}: {len(sig_df)} 个DEG')
        
    except Exception as e:
        print(f'  警告: 无法读取 {sheet_name}: {e}')

writer.close()
print(f'\n输出完成: {output_file}')
