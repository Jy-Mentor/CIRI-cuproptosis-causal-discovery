#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 汇总所有 MR 分析结果 CSV 到一个 Excel 文件

import pandas as pd
import os
from pathlib import Path

# 配置
RESULTS_DIR = Path('c:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/mr_results_138genes_integrated')
OUTPUT_FILE = RESULTS_DIR / 'MR_分析结果汇总.xlsx'

# 获取所有 CSV 文件
csv_files = list(RESULTS_DIR.glob('*.csv'))
print(f'找到 {len(csv_files)} 个 CSV 文件')

# 创建 Excel writer
with pd.ExcelWriter(OUTPUT_FILE, engine='openpyxl') as writer:
    # 分类处理文件
    summary_files = []
    harmonised_files = []
    mr_result_files = []
    other_files = []
    
    for csv_file in csv_files:
        filename = csv_file.stem
        if 'summary' in filename or 'detailed' in filename or 'fdr' in filename:
            summary_files.append(csv_file)
        elif 'harmonised' in filename:
            harmonised_files.append(csv_file)
        elif '_mr_results' in filename:
            mr_result_files.append(csv_file)
        else:
            other_files.append(csv_file)
    
    # 1. 写入汇总结果（放在第一个 sheet）
    print("\n写入汇总结果...")
    for csv_file in summary_files:
        sheet_name = csv_file.stem.replace('_', ' ')[:31]  # Excel sheet name 长度限制
        try:
            df = pd.read_csv(csv_file, encoding='utf-8')
            df.to_excel(writer, sheet_name=sheet_name, index=False)
            print(f"  ✓ {sheet_name}: {len(df)} 行")
        except Exception as e:
            print(f"  ✗ {csv_file.name}: {e}")
    
    # 2. 写入显著基因的详细数据（harmonised + mr_results）
    print("\n写入显著基因详细数据...")
    significant_genes = ['ADRB1', 'SREBF1', 'ACADVL', 'PABPC1', 'PTPRJ', 'RHOC', 'AIF1', 'CNDP2', 'HSD17B4']
    
    for gene in significant_genes:
        # 查找该基因的 harmonised 和 mr_results 文件
        harmonised_file = RESULTS_DIR / f'{gene}_harmonised.csv'
        mr_file = RESULTS_DIR / f'{gene}_mr_results.csv'
        
        if harmonised_file.exists():
            try:
                df = pd.read_csv(harmonised_file, encoding='utf-8')
                sheet_name = f'{gene}_harmonised'[:31]
                df.to_excel(writer, sheet_name=sheet_name, index=False)
                print(f"  ✓ {sheet_name}: {len(df)} 行")
            except Exception as e:
                print(f"  ✗ {harmonised_file.name}: {e}")
        
        if mr_file.exists():
            try:
                df = pd.read_csv(mr_file, encoding='utf-8')
                sheet_name = f'{gene}_mr_results'[:31]
                df.to_excel(writer, sheet_name=sheet_name, index=False)
                print(f"  ✓ {sheet_name}: {len(df)} 行")
            except Exception as e:
                print(f"  ✗ {mr_file.name}: {e}")
    
    # 3. 写入铜死亡相关基因数据
    print("\n写入铜死亡相关基因数据...")
    cuproptosis_genes = ['NFKB1', 'FDX1', 'ATP7B', 'ATOX1']
    
    for gene in cuproptosis_genes:
        if gene in significant_genes:
            continue  # 已经添加过了
        
        harmonised_file = RESULTS_DIR / f'{gene}_harmonised.csv'
        mr_file = RESULTS_DIR / f'{gene}_mr_results.csv'
        
        if harmonised_file.exists():
            try:
                df = pd.read_csv(harmonised_file, encoding='utf-8')
                sheet_name = f'{gene}_harmonised'[:31]
                df.to_excel(writer, sheet_name=sheet_name, index=False)
                print(f"  ✓ {sheet_name}: {len(df)} 行")
            except Exception as e:
                print(f"  ✗ {harmonised_file.name}: {e}")
        
        if mr_file.exists():
            try:
                df = pd.read_csv(mr_file, encoding='utf-8')
                sheet_name = f'{gene}_mr_results'[:31]
                df.to_excel(writer, sheet_name=sheet_name, index=False)
                print(f"  ✓ {sheet_name}: {len(df)} 行")
            except Exception as e:
                print(f"  ✗ {mr_file.name}: {e}")
    
    # 4. 写入跳过的基因列表
    skipped_file = RESULTS_DIR / 'skipped_genes.csv'
    if skipped_file.exists():
        try:
            df = pd.read_csv(skipped_file, encoding='utf-8')
            df.to_excel(writer, sheet_name='Skipped_Genes', index=False)
            print(f"\n  ✓ Skipped_Genes: {len(df)} 行")
        except Exception as e:
            print(f"  ✗ skipped_genes.csv: {e}")

print(f"\n完成！Excel 文件已保存到：{OUTPUT_FILE}")
print(f"\nExcel 文件包含:")
print(f"  - 汇总结果 sheets")
print(f"  - 显著基因详细数据 sheets")
print(f"  - 铜死亡相关基因数据 sheets")
print(f"  - 跳过的基因列表")
