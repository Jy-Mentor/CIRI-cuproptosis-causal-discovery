#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Microglia亚群内NFKB1-FDX1相关性分析
从已保存的AnnData对象中读取数据进行分析
"""

import os
import pandas as pd
import scanpy as sc
from scipy import stats

# 设置路径
WORK_DIR = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙"
OUTPUT_DIR = os.path.join(WORK_DIR, "GSE174574_SCISSOR_Results")

# 读取AnnData对象
print("读取AnnData对象...")
adata = sc.read_h5ad(os.path.join(OUTPUT_DIR, "GSE174574_SCISSOR_final.h5ad"))
print(f"总细胞数: {adata.n_obs}")

# 检查NFKB1和FDX1评分是否存在
if 'NFKB1_Score' not in adata.obs.columns or 'FDX1_Score' not in adata.obs.columns:
    print("错误: NFKB1_Score 或 FDX1_Score 不存在于数据中")
    exit(1)

# ============================================
# Microglia亚群内NFKB1-FDX1相关性分析
# ============================================
print("\n" + "="*60)
print("Microglia亚群内NFKB1-FDX1相关性分析")
print("="*60)

# 提取Microglia亚群
microglia_mask = adata.obs['CellType'] == 'Microglia'
microglia_data = adata[microglia_mask, :]

print(f"\nMicroglia亚群细胞数: {microglia_data.n_obs}")
print(f"  MCAO组: {sum(microglia_data.obs['group'] == 'MCAO')}")
print(f"  Sham组: {sum(microglia_data.obs['group'] == 'Sham')}")

# 在Microglia内分别计算MCAO和Sham组的相关性
microglia_cor_results = []

for group in ['MCAO', 'Sham']:
    group_mask = microglia_data.obs['group'] == group
    sub_data = microglia_data[group_mask, :]
    
    if sub_data.n_obs > 10:
        nfkb1_scores = sub_data.obs['NFKB1_Score']
        fdx1_scores = sub_data.obs['FDX1_Score']
        
        rho, pval = stats.spearmanr(nfkb1_scores, fdx1_scores)
        
        print(f"\n{group}组 Microglia:")
        print(f"  细胞数: {sub_data.n_obs}")
        print(f"  Spearman rho = {rho:.4f}")
        print(f"  P-value = {pval:.2e}")
        
        microglia_cor_results.append({
            'CellType': 'Microglia',
            'Group': group,
            'Correlation': rho,
            'P_value': pval,
            'N': sub_data.n_obs
        })

# 全细胞MCAO组相关性
mcmo_mask = adata.obs['group'] == 'MCAO'
nfkb1_scores_all = adata.obs.loc[mcmo_mask, 'NFKB1_Score']
fdx1_scores_all = adata.obs.loc[mcmo_mask, 'FDX1_Score']
rho_all, pval_all = stats.spearmanr(nfkb1_scores_all, fdx1_scores_all)

print(f"\n全细胞MCA0组:")
print(f"  细胞数: {sum(mcmo_mask)}")
print(f"  Spearman rho = {rho_all:.4f}")
print(f"  P-value = {pval_all:.2e}")

# 保存结果
cor_results_all = [
    {'CellType': 'All_Cells', 'Group': 'MCAO', 'Correlation': rho_all, 'P_value': pval_all, 'N': sum(mcmo_mask)}
] + microglia_cor_results

cor_df = pd.DataFrame(cor_results_all)
cor_df.to_csv(os.path.join(OUTPUT_DIR, "03_nfkb1_fdx1_correlation.csv"), index=False)

print("\n" + "="*60)
print("相关性分析结果汇总:")
print("="*60)
print(cor_df.to_string(index=False))

# 结果解读
print("\n" + "="*60)
print("结果解读:")
print("="*60)

microglia_mcao_rho = cor_df[(cor_df['CellType'] == 'Microglia') & (cor_df['Group'] == 'MCAO')]['Correlation'].values
if len(microglia_mcao_rho) > 0:
    microglia_rho = microglia_mcao_rho[0]
    if microglia_rho > 0.1:
        print(f"✅ Microglia内NFKB1-FDX1相关性 (rho={microglia_rho:.3f}) > 0.1")
        print("   → 比全细胞更强，支持亚群特异性共表达")
    elif microglia_rho > 0.05:
        print(f"⚠️ Microglia内NFKB1-FDX1相关性 (rho={microglia_rho:.3f}) ~ 0.05-0.1")
        print("   → 中等程度相关，可能存在间接调控")
    else:
        print(f"❌ Microglia内NFKB1-FDX1相关性 (rho={microglia_rho:.3f}) < 0.05")
        print("   → 即使在靶细胞内，二者也是弱相关/间接调控")

print("\n结果已保存到:", os.path.join(OUTPUT_DIR, "03_nfkb1_fdx1_correlation.csv"))
