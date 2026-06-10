#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Microglia虚拟敲除分析
基于基因表达相关性预测敲除效应
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
import seaborn as sns

# 设置路径
WORK_DIR = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙"
OUTPUT_DIR = os.path.join(WORK_DIR, "GSE174574_SCISSOR_Results")
KNOCKDOWN_DIR = os.path.join(OUTPUT_DIR, "virtual_knockdown")
os.makedirs(KNOCKDOWN_DIR, exist_ok=True)

print("="*60)
print("Microglia虚拟敲除分析")
print("="*60)

# 1. 读取AnnData对象
print("\n读取AnnData对象...")
adata = sc.read_h5ad(os.path.join(OUTPUT_DIR, "GSE174574_SCISSOR_final.h5ad"))
print(f"总细胞数: {adata.n_obs}")

# 2. 提取Microglia
print("\n提取Microglia亚群...")
micro = adata[adata.obs['CellType'] == 'Microglia', :].copy()
print(f"Microglia细胞数: {micro.n_obs}")
print(f"  MCAO组: {sum(micro.obs['group'] == 'MCAO')}")
print(f"  Sham组: {sum(micro.obs['group'] == 'Sham')}")

# 分别分析MCAO和Sham组
for group in ['MCAO', 'Sham']:
    print(f"\n{'='*60}")
    print(f"{group}组Microglia虚拟敲除分析")
    print(f"{'='*60}")
    
    # 提取该组的Microglia
    group_mask = micro.obs['group'] == group
    micro_group = micro[group_mask, :].copy()
    
    if micro_group.n_obs < 100:
        print(f"细胞数不足 ({micro_group.n_obs})，跳过")
        continue
    
    print(f"细胞数: {micro_group.n_obs}")
    
    # 3. 目标基因（精简为6个核心）
    targets = ['Nfkb1', 'Fdx1', 'Stat3', 'Hspa5', 'Hmox1', 'Gpx4']
    targets = [g for g in targets if g in micro_group.var_names]
    print(f"\n可用目标基因: {targets}")
    
    if len(targets) < 2:
        print("目标基因不足，跳过")
        continue
    
    # 4. 筛选Top高表达细胞（基于目标基因的表达）
    print("\n筛选Top高表达细胞...")
    X_targets = micro_group[:, targets].X.toarray() if hasattr(micro_group[:, targets].X, 'toarray') else np.array(micro_group[:, targets].X)
    
    # 计算每个细胞的覆盖率和平均表达
    coverage = (X_targets > 0).mean(axis=1)
    mean_expr = X_targets.mean(axis=1)
    
    # 综合评分
    score = 0.6 * (coverage / (coverage.max() + 1e-10)) + 0.4 * (mean_expr / (mean_expr.max() + 1e-10))
    
    # 选择Top 4000细胞（或所有细胞如果少于4000）
    n_select = min(4000, len(score))
    top_idx = np.argsort(score)[-n_select:]
    micro_sub = micro_group[top_idx, :].copy()
    
    print(f"选择Top {n_select}细胞进行分析")
    
    # 5. 计算Spearman相关性矩阵
    print("\n计算基因间Spearman相关性...")
    
    # 获取表达量最高的4000个基因（不分层扩大基因池）
    print("  选择表达量最高的4000个基因...")
    mean_expression = np.array(micro_sub.X.mean(axis=0)).flatten()
    top_genes_idx = np.argsort(mean_expression)[-4000:]
    genes_for_corr = micro_sub.var_names[top_genes_idx].tolist()
    
    # 确保目标基因包含在内
    genes_for_corr = list(set(genes_for_corr + targets))
    
    print(f"用于相关性分析的基因数: {len(genes_for_corr)}")
    
    # 提取表达矩阵
    X_full = micro_sub[:, genes_for_corr].X.toarray() if hasattr(micro_sub[:, genes_for_corr].X, 'toarray') else np.array(micro_sub[:, genes_for_corr].X)
    
    # 计算相关性矩阵（基因 x 基因）- 使用pandas优化
    print("  使用pandas计算Spearman相关性...")
    df_expr = pd.DataFrame(X_full, columns=genes_for_corr)
    corr_df = df_expr.corr(method='spearman')
    
    # 处理NaN值
    corr_df = corr_df.fillna(0)
    
    # 6. 虚拟敲除预测
    print("\n虚拟敲除预测...")
    
    # 方法: 基于相关性的线性预测
    # 假设敲除某基因后，其他基因的变化与该基因的相关性成正比
    # pred_logFC = -correlation * knockdown_factor
    knockdown_factor = 2.0  # 假设敲除效率导致2倍变化基准
    
    pred_logFC = pd.DataFrame(index=genes_for_corr, columns=targets)
    
    for target in targets:
        if target in corr_df.columns:
            # 预测敲除该基因后其他基因的变化
            # 负相关意味着敲除后该基因上调，正相关意味着下调
            pred_logFC[target] = -corr_df[target] * knockdown_factor
            pred_logFC.loc[target, target] = -knockdown_factor  # 自身完全敲除
    
    pred_logFC = pred_logFC.astype(float)
    
    # 7. 关键输出
    print("\n" + "="*40)
    print("核心结果")
    print("="*40)
    
    if 'Nfkb1' in targets and 'Fdx1' in targets:
        nfkb1_to_fdx1 = pred_logFC.loc['Fdx1', 'Nfkb1']
        fdx1_to_nfkb1 = pred_logFC.loc['Nfkb1', 'Fdx1']
        
        print(f"\nNFKB1敲除 → FDX1 predicted_logFC: {nfkb1_to_fdx1:.4f}")
        print(f"FDX1敲除 → NFKB1 predicted_logFC: {fdx1_to_nfkb1:.4f}")
        
        # 判断调控方向
        if nfkb1_to_fdx1 < -0.1:
            print("  → NFKB1正向调控FDX1 (敲除NFKB1导致FDX1下调)")
        elif nfkb1_to_fdx1 > 0.1:
            print("  → NFKB1负向调控FDX1 (敲除NFKB1导致FDX1上调)")
        else:
            print("  → NFKB1对FDX1调控效应较弱")
    
    # 8. DEG统计（|logFC| > 0.25）
    print("\n" + "="*40)
    print("各基因敲除后差异基因统计 (|logFC| > 0.25)")
    print("="*40)
    
    deg_stats = {}
    for target in targets:
        if target in pred_logFC.columns:
            up = (pred_logFC[target] > 0.25).sum()
            down = (pred_logFC[target] < -0.25).sum()
            deg_stats[target] = {'up': up, 'down': down, 'total': up + down}
    
    deg_df = pd.DataFrame(deg_stats).T
    print(deg_df)
    
    # 9. 共同下游基因分析
    if 'Nfkb1' in targets and 'Fdx1' in targets:
        print("\n" + "="*40)
        print("NFKB1与FDX1共同下游基因分析")
        print("="*40)
        
        nfkb1_degs = set(pred_logFC['Nfkb1'][abs(pred_logFC['Nfkb1']) > 0.25].index)
        fdx1_degs = set(pred_logFC['Fdx1'][abs(pred_logFC['Fdx1']) > 0.25].index)
        common_degs = nfkb1_degs & fdx1_degs
        
        print(f"NFKB1敲除差异基因: {len(nfkb1_degs)}个")
        print(f"FDX1敲除差异基因: {len(fdx1_degs)}个")
        print(f"共同下游基因: {len(common_degs)}个")
        
        if len(common_degs) > 0:
            print(f"\n共同下游基因列表 (前20个):")
            for i, gene in enumerate(list(common_degs)[:20]):
                nfkb1_fc = pred_logFC.loc[gene, 'Nfkb1']
                fdx1_fc = pred_logFC.loc[gene, 'Fdx1']
                print(f"  {gene}: NFKB1_logFC={nfkb1_fc:.3f}, FDX1_logFC={fdx1_fc:.3f}")
        
        # 计算Jaccard相似性
        jaccard = len(common_degs) / len(nfkb1_degs | fdx1_degs) if len(nfkb1_degs | fdx1_degs) > 0 else 0
        print(f"\nJaccard相似性指数: {jaccard:.4f}")
    
    # 10. 保存结果
    print(f"\n保存{group}组结果...")
    
    # 保存预测logFC矩阵
    pred_logFC.to_csv(os.path.join(KNOCKDOWN_DIR, f"predicted_logFC_{group}.csv"))
    
    # 保存DEG统计
    deg_df.to_csv(os.path.join(KNOCKDOWN_DIR, f"knockdown_stats_{group}.csv"))
    
    # 保存共同下游基因
    if 'Nfkb1' in targets and 'Fdx1' in targets:
        pd.Series(list(common_degs)).to_csv(
            os.path.join(KNOCKDOWN_DIR, f"common_genes_{group}.csv"), 
            index=False, 
            header=False
        )
        
        # 保存相关性矩阵（目标基因部分）
        corr_df.loc[targets, targets].to_csv(
            os.path.join(KNOCKDOWN_DIR, f"target_correlation_{group}.csv")
        )
    
    # 11. 可视化
    print(f"\n生成{group}组可视化...")
    
    # 热图：目标基因的相关性
    plt.figure(figsize=(8, 6))
    sns.heatmap(corr_df.loc[targets, targets], 
                annot=True, 
                cmap='RdBu_r', 
                center=0,
                vmin=-1, 
                vmax=1,
                square=True)
    plt.title(f'{group} Microglia: Target Gene Correlation')
    plt.tight_layout()
    plt.savefig(os.path.join(KNOCKDOWN_DIR, f"target_correlation_{group}.png"), dpi=300)
    plt.close()
    
    # 热图：预测的敲除效应（前50个变化最大的基因）
    top_genes_by_change = pred_logFC.abs().sum(axis=1).sort_values(ascending=False).head(50).index
    plt.figure(figsize=(10, 12))
    sns.heatmap(pred_logFC.loc[top_genes_by_change, targets], 
                cmap='RdBu_r', 
                center=0,
                vmin=-2, 
                vmax=2)
    plt.title(f'{group} Microglia: Predicted Knockdown Effects (Top 50 genes)')
    plt.tight_layout()
    plt.savefig(os.path.join(KNOCKDOWN_DIR, f"knockdown_heatmap_{group}.png"), dpi=300)
    plt.close()

print("\n" + "="*60)
print("Microglia虚拟敲除分析完成!")
print("="*60)
print(f"\n结果保存在: {KNOCKDOWN_DIR}")
