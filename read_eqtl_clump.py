#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
读取eQTLgen clumping数据，了解数据结构
"""

import pandas as pd
import os

# 读取其中一个clumping文件
eqtl_file = r"D:\EQTL\clump\eQTLgen_allgene_p_5e-8_kb_10000_r2_0.001.xlsx"

print("读取eQTLgen clumping文件...")
print(f"文件: {eqtl_file}")

try:
    df = pd.read_excel(eqtl_file)
    print(f"\n数据维度: {df.shape}")
    print(f"\n列名: {list(df.columns)}")
    print(f"\n前5行:")
    print(df.head())
    
    # 检查基因名格式
    print(f"\n基因名示例 (Gene列前10个):")
    if 'Gene' in df.columns:
        print(df['Gene'].head(10).tolist())
    elif 'gene' in df.columns:
        print(df['gene'].head(10).tolist())
    elif 'Phenotype' in df.columns:
        print(df['Phenotype'].head(10).tolist())
    
    # 检查SNP相关列
    print(f"\nSNP列信息:")
    if 'SNP' in df.columns:
        print(f"SNP列存在，示例: {df['SNP'].head(5).tolist()}")
    elif 'rsid' in df.columns:
        print(f"rsid列存在，示例: {df['rsid'].head(5).tolist()}")
        
except Exception as e:
    print(f"错误: {e}")
