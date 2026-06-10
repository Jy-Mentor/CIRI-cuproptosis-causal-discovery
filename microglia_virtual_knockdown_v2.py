#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Microglia虚拟敲除分析 V2 - 优化版
- 使用pandas corr(method='spearman')优化计算效率
- 支持全量Microglia分析或分块计算
- 内存监控
"""

import os
import sys
import time
import psutil
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

# 内存监控函数
def get_memory_usage():
    """获取当前内存使用(MB)"""
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / 1024 / 1024

# 设置路径
WORK_DIR = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙"
OUTPUT_DIR = os.path.join(WORK_DIR, "GSE174574_SCISSOR_Results")
KNOCKDOWN_DIR = os.path.join(OUTPUT_DIR, "virtual_knockdown_v2")
os.makedirs(KNOCKDOWN_DIR, exist_ok=True)

print("="*70)
print("Microglia虚拟敲除分析 V2 - 优化版")
print("="*70)
print(f"初始内存使用: {get_memory_usage():.1f} MB")

# 1. 读取AnnData对象
print("\n读取AnnData对象...")
start_time = time.time()
adata = sc.read_h5ad(os.path.join(OUTPUT_DIR, "GSE174574_SCISSOR_final.h5ad"))
print(f"总细胞数: {adata.n_obs}")
print(f"读取耗时: {time.time() - start_time:.1f}s")
print(f"当前内存: {get_memory_usage():.1f} MB")

# 2. 提取Microglia
print("\n提取Microglia亚群...")
micro = adata[adata.obs['CellType'] == 'Microglia', :].copy()
print(f"Microglia细胞数: {micro.n_obs}")
print(f"  MCAO组: {sum(micro.obs['group'] == 'MCAO')}")
print(f"  Sham组: {sum(micro.obs['group'] == 'Sham')}")
print(f"当前内存: {get_memory_usage():.1f} MB")

# 目标基因
targets = ['Nfkb1', 'Fdx1', 'Stat3', 'Hspa5', 'Hmox1', 'Gpx4']
targets = [g for g in targets if g in micro.var_names]
print(f"\n可用目标基因: {targets}")

if len(targets) < 2:
    print("目标基因不足，退出")
    sys.exit(1)

# 对每个组进行分析
for group in ['MCAO', 'Sham']:
    print(f"\n{'='*70}")
    print(f"{group}组Microglia虚拟敲除分析")
    print(f"{'='*70}")
    
    # 提取该组的Microglia
    group_mask = micro.obs['group'] == group
    micro_group = micro[group_mask, :].copy()
    
    if micro_group.n_obs < 100:
        print(f"细胞数不足 ({micro_group.n_obs})，跳过")
        continue
    
    print(f"细胞数: {micro_group.n_obs}")
    print(f"当前内存: {get_memory_usage():.1f} MB")
    
    # 筛选高表达细胞（全部使用，不降采样）
    print("\n筛选高表达细胞（使用全部可用细胞）...")
    X_targets = micro_group[:, targets].X.toarray() if hasattr(micro_group[:, targets].X, 'toarray') else np.array(micro_group[:, targets].X)
    
    coverage = (X_targets > 0).mean(axis=1)
    mean_expr = X_targets.mean(axis=1)
    score = 0.6 * (coverage / (coverage.max() + 1e-10)) + 0.4 * (mean_expr / (mean_expr.max() + 1e-10))
    
    # 使用全部细胞（或Top 8000如果超过）
    n_select = min(8000, len(score)) if len(score) > 8000 else len(score)
    if n_select < len(score):
        print(f"细胞数超过8000，选择Top {n_select}高表达细胞")
        top_idx = np.argsort(score)[-n_select:]
        micro_sub = micro_group[top_idx, :].copy()
    else:
        print(f"使用全部 {n_select} 细胞")
        micro_sub = micro_group.copy()
    
    print(f"分析细胞数: {micro_sub.n_obs}")
    print(f"当前内存: {get_memory_usage():.1f} MB")
    
    # 3. 获取基因池（扩大至6000个高表达基因）
    print("\n选择表达量最高的6000个基因...")
    mean_expression = np.array(micro_sub.X.mean(axis=0)).flatten()
    top_genes_idx = np.argsort(mean_expression)[-6000:]
    genes_for_corr = micro_sub.var_names[top_genes_idx].tolist()
    
    # 确保目标基因包含在内
    genes_for_corr = list(set(genes_for_corr + targets))
    
    print(f"用于相关性分析的基因数: {len(genes_for_corr)}")
    print(f"当前内存: {get_memory_usage():.1f} MB")
    
    # 4. 使用pandas计算Spearman相关性（优化版本）
    print("\n使用pandas计算Spearman相关性...")
    start_time = time.time()
    
    # 提取表达矩阵
    X_full = micro_sub[:, genes_for_corr].X.toarray() if hasattr(micro_sub[:, genes_for_corr].X, 'toarray') else np.array(micro_sub[:, genes_for_corr].X)
    
    # 使用pandas corr（优化版本）
    df_expr = pd.DataFrame(X_full, columns=genes_for_corr)
    corr_df = df_expr.corr(method='spearman')
    corr_df = corr_df.fillna(0)
    
    print(f"相关性计算耗时: {time.time() - start_time:.1f}s")
    print(f"相关性矩阵维度: {corr_df.shape}")
    print(f"当前内存: {get_memory_usage():.1f} MB")
    
    # 5. 虚拟敲除预测
    print("\n虚拟敲除预测...")
    knockdown_factor = 2.0
    pred_logFC = pd.DataFrame(index=genes_for_corr, columns=targets)
    
    for target in targets:
        if target in corr_df.columns:
            pred_logFC[target] = -corr_df[target] * knockdown_factor
            pred_logFC.loc[target, target] = -knockdown_factor
    
    pred_logFC = pred_logFC.astype(float)
    
    # 6. 核心结果
    print("\n" + "="*50)
    print("核心结果")
    print("="*50)
    
    if 'Nfkb1' in targets and 'Fdx1' in targets:
        nfkb1_to_fdx1 = pred_logFC.loc['Fdx1', 'Nfkb1']
        fdx1_to_nfkb1 = pred_logFC.loc['Nfkb1', 'Fdx1']
        
        print(f"\nNFKB1敲除 → FDX1 predicted_logFC: {nfkb1_to_fdx1:.4f}")
        print(f"FDX1敲除 → NFKB1 predicted_logFC: {fdx1_to_nfkb1:.4f}")
        
        if abs(nfkb1_to_fdx1) > 0.1:
            direction = "正向" if nfkb1_to_fdx1 < 0 else "负向"
            print(f"  → NFKB1{direction}调控FDX1")
        else:
            print(f"  → NFKB1对FDX1调控效应较弱")
    
    # 7. DEG统计
    print("\n各基因敲除后差异基因统计 (|logFC| > 0.25):")
    deg_stats = {}
    for target in targets:
        if target in pred_logFC.columns:
            up = (pred_logFC[target] > 0.25).sum()
            down = (pred_logFC[target] < -0.25).sum()
            deg_stats[target] = {'up': up, 'down': down, 'total': up + down}
    
    deg_df = pd.DataFrame(deg_stats).T
    print(deg_df)
    
    # 8. 共同下游基因分析
    if 'Nfkb1' in targets and 'Fdx1' in targets:
        print("\n" + "="*50)
        print("NFKB1与FDX1共同下游基因分析")
        print("="*50)
        
        nfkb1_degs = set(pred_logFC['Nfkb1'][abs(pred_logFC['Nfkb1']) > 0.25].index)
        fdx1_degs = set(pred_logFC['Fdx1'][abs(pred_logFC['Fdx1']) > 0.25].index)
        common_degs = nfkb1_degs & fdx1_degs
        
        print(f"NFKB1敲除差异基因: {len(nfkb1_degs)}个")
        print(f"FDX1敲除差异基因: {len(fdx1_degs)}个")
        print(f"共同下游基因: {len(common_degs)}个")
        
        if len(common_degs) > 0:
            print(f"\n共同下游基因列表:")
            for gene in list(common_degs)[:20]:
                nfkb1_fc = pred_logFC.loc[gene, 'Nfkb1']
                fdx1_fc = pred_logFC.loc[gene, 'Fdx1']
                print(f"  {gene}: NFKB1_logFC={nfkb1_fc:.3f}, FDX1_logFC={fdx1_fc:.3f}")
        
        jaccard = len(common_degs) / len(nfkb1_degs | fdx1_degs) if len(nfkb1_degs | fdx1_degs) > 0 else 0
        print(f"\nJaccard相似性指数: {jaccard:.4f}")
    
    # 9. 保存结果
    print(f"\n保存{group}组结果...")
    pred_logFC.to_csv(os.path.join(KNOCKDOWN_DIR, f"predicted_logFC_{group}_v2.csv"))
    deg_df.to_csv(os.path.join(KNOCKDOWN_DIR, f"knockdown_stats_{group}_v2.csv"))
    corr_df.loc[targets, targets].to_csv(os.path.join(KNOCKDOWN_DIR, f"target_correlation_{group}_v2.csv"))
    
    if 'Nfkb1' in targets and 'Fdx1' in targets:
        pd.Series(list(common_degs)).to_csv(
            os.path.join(KNOCKDOWN_DIR, f"common_genes_{group}_v2.csv"),
            index=False, header=False
        )
    
    # 10. 可视化
    print(f"\n生成{group}组可视化...")
    
    # 相关性热图
    plt.figure(figsize=(8, 6))
    sns.heatmap(corr_df.loc[targets, targets], annot=True, cmap='RdBu_r', center=0, vmin=-1, vmax=1, square=True)
    plt.title(f'{group} Microglia: Target Gene Correlation (V2)')
    plt.tight_layout()
    plt.savefig(os.path.join(KNOCKDOWN_DIR, f"target_correlation_{group}_v2.png"), dpi=300)
    plt.close()
    
    # 敲除效应热图
    top_genes = pred_logFC.abs().sum(axis=1).sort_values(ascending=False).head(50).index
    plt.figure(figsize=(10, 12))
    sns.heatmap(pred_logFC.loc[top_genes, targets], cmap='RdBu_r', center=0, vmin=-2, vmax=2)
    plt.title(f'{group} Microglia: Predicted Knockdown Effects (Top 50)')
    plt.tight_layout()
    plt.savefig(os.path.join(KNOCKDOWN_DIR, f"knockdown_heatmap_{group}_v2.png"), dpi=300)
    plt.close()

print("\n" + "="*70)
print("Microglia虚拟敲除分析 V2 完成!")
print("="*70)
print(f"结果保存在: {KNOCKDOWN_DIR}")
print(f"最终内存使用: {get_memory_usage():.1f} MB")
