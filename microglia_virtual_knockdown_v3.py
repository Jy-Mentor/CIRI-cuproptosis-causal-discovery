#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Microglia虚拟敲除分析 V3 - 关键修复版
修复: 按单个靶基因筛选细胞，而非综合评分
逻辑: Nfkb1敲除只选Nfkb1阳性且高表达的细胞
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
import gc

# 内存监控函数
def get_memory_usage():
    """获取当前内存使用(MB)"""
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / 1024 / 1024

# 设置路径
WORK_DIR = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙"
OUTPUT_DIR = os.path.join(WORK_DIR, "GSE174574_SCISSOR_Results")
KNOCKDOWN_DIR = os.path.join(OUTPUT_DIR, "virtual_knockdown_v3")
os.makedirs(KNOCKDOWN_DIR, exist_ok=True)

print("="*70)
print("Microglia虚拟敲除分析 V3 - 关键修复版")
print("修复: 按单个靶基因筛选细胞")
print("="*70)
print(f"初始内存使用: {get_memory_usage():.1f} MB")

# 1. 读取AnnData对象
print("\n读取AnnData对象...")
start_time = time.time()
adata = sc.read_h5ad(os.path.join(OUTPUT_DIR, "GSE174574_SCISSOR_v3.h5ad"))
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
    
    print(f"总细胞数: {micro_group.n_obs}")
    print(f"当前内存: {get_memory_usage():.1f} MB")
    
    # 获取基因池（6000个高表达基因）- 在全组上计算
    print("\n选择表达量最高的6000个基因作为分析池...")
    mean_expression = np.array(micro_group.X.mean(axis=0)).flatten()
    top_genes_idx = np.argsort(mean_expression)[-6000:]
    genes_for_corr = micro_group.var_names[top_genes_idx].tolist()
    genes_for_corr = list(set(genes_for_corr + targets))
    print(f"用于相关性分析的基因数: {len(genes_for_corr)}")
    
    # 存储每个靶基因的结果
    all_pred_logFC = {}
    all_target_corrs = {}
    
    # 🔑 关键修复: 对每个靶基因单独循环
    for target in targets:
        print(f"\n{'='*50}")
        print(f"目标基因: {target} 虚拟敲除")
        print(f"{'='*50}")
        
        # 1. 提取当前目标基因的表达量
        target_expr = micro_group[:, target].X.toarray().flatten()
        
        # 2. 只保留表达量 > 0 的细胞（关键修复！）
        positive_mask = target_expr > 0
        n_positive = positive_mask.sum()
        print(f"  {target} 阳性细胞: {n_positive} / {len(target_expr)}")
        
        if n_positive < 100:
            print(f"  ⚠️ 阳性细胞不足，跳过")
            continue
        
        # 3. 在阳性细胞中，按该基因表达量排序，取 Top 50%（至少500个）
        positive_expr = target_expr[positive_mask]
        positive_cells_idx = np.where(positive_mask)[0]
        
        n_select = max(500, int(n_positive * 0.5))  # 至少500个，最多50%
        n_select = min(n_select, n_positive)  # 不超过阳性细胞总数
        
        # 按表达量排序，取Top
        top_idx_in_positive = np.argsort(positive_expr)[-n_select:]
        selected_cell_idx = positive_cells_idx[top_idx_in_positive]
        
        # 创建子集
        micro_sub = micro_group[selected_cell_idx, :].copy()
        print(f"  筛选后用于分析的细胞: {micro_sub.n_obs} (Top {n_select} {target}高表达细胞)")
        print(f"  {target}表达范围: {positive_expr[top_idx_in_positive].min():.2f} - {positive_expr[top_idx_in_positive].max():.2f}")
        
        # 4. 在这个细胞子集上，计算 target 与基因池的 Spearman 相关性
        print(f"  计算{target}与基因池的相关性...")
        start_time = time.time()
        
        # 提取表达矩阵
        X_full = micro_sub[:, genes_for_corr].X.toarray() if hasattr(micro_sub[:, genes_for_corr].X, 'toarray') else np.array(micro_sub[:, genes_for_corr].X)
        df_expr = pd.DataFrame(X_full, columns=genes_for_corr)
        
        # 计算相关性
        target_correlations = df_expr.corrwith(df_expr[target], method='spearman')
        target_correlations = target_correlations.fillna(0)
        
        print(f"  相关性计算耗时: {time.time() - start_time:.1f}s")
        
        # 保存该target的相关性结果
        all_target_corrs[target] = target_correlations
        
        # 5. 虚拟敲除预测（仅针对当前target）
        knockdown_factor = 2.0
        pred_logFC_single = -target_correlations * knockdown_factor
        pred_logFC_single.loc[target] = -knockdown_factor  # 自身敲除为-2
        
        all_pred_logFC[target] = pred_logFC_single
        
        # 6. 打印该target的核心结果
        # 检查该target敲除对FDX1的影响（如果是Nfkb1）
        if target == 'Nfkb1' and 'Fdx1' in genes_for_corr:
            fdx1_logfc = pred_logFC_single.loc['Fdx1']
            print(f"\n  💡 Nfkb1敲除 → Fdx1 predicted_logFC: {fdx1_logfc:.4f}")
            if abs(fdx1_logfc) > 0.1:
                direction = "下调" if fdx1_logfc < 0 else "上调"
                print(f"     提示: Nfkb1可能{direction}调控Fdx1")
        
        # 检查该target敲除对NFKB1的影响（如果是Fdx1）
        if target == 'Fdx1' and 'Nfkb1' in genes_for_corr:
            nfkb1_logfc = pred_logFC_single.loc['Nfkb1']
            print(f"\n  💡 Fdx1敲除 → Nfkb1 predicted_logFC: {nfkb1_logfc:.4f}")
            if abs(nfkb1_logfc) > 0.1:
                direction = "下调" if nfkb1_logfc < 0 else "上调"
                print(f"     提示: Fdx1可能{direction}调控Nfkb1")
        
        # DEG统计
        up = (pred_logFC_single > 0.25).sum()
        down = (pred_logFC_single < -0.25).sum()
        print(f"  预测DEGs: Up={up}, Down={down}, Total={up+down}")
        
        # 清理内存
        del micro_sub, X_full, df_expr
        gc.collect()
    
    # 合并所有target的预测结果
    print(f"\n{'='*50}")
    print(f"{group}组汇总结果")
    print(f"{'='*50}")
    
    if all_pred_logFC:
        pred_logFC_df = pd.DataFrame(all_pred_logFC)
        
        # 保存完整预测矩阵
        pred_logFC_df.to_csv(os.path.join(KNOCKDOWN_DIR, f"predicted_logFC_{group}_v3.csv"))
        
        # 保存target间相关性
        target_corr_matrix = pd.DataFrame({t: all_target_corrs[t] for t in all_target_corrs.keys() if t in all_target_corrs}).loc[targets, targets]
        target_corr_matrix.to_csv(os.path.join(KNOCKDOWN_DIR, f"target_correlation_{group}_v3.csv"))
        
        # NFKB1-FDX1互作分析
        if 'Nfkb1' in all_pred_logFC and 'Fdx1' in all_pred_logFC:
            print("\n【NFKB1-FDX1互作分析】")
            
            nfkb1_to_fdx1 = pred_logFC_df.loc['Fdx1', 'Nfkb1']
            fdx1_to_nfkb1 = pred_logFC_df.loc['Nfkb1', 'Fdx1']
            
            print(f"  Nfkb1敲除 → Fdx1 predicted_logFC: {nfkb1_to_fdx1:.4f}")
            print(f"  Fdx1敲除 → Nfkb1 predicted_logFC: {fdx1_to_nfkb1:.4f}")
            
            # 共同下游基因
            nfkb1_degs = set(pred_logFC_df['Nfkb1'][abs(pred_logFC_df['Nfkb1']) > 0.25].index)
            fdx1_degs = set(pred_logFC_df['Fdx1'][abs(pred_logFC_df['Fdx1']) > 0.25].index)
            common_degs = nfkb1_degs & fdx1_degs
            
            print(f"\n  NFKB1敲除差异基因: {len(nfkb1_degs)}个")
            print(f"  FDX1敲除差异基因: {len(fdx1_degs)}个")
            print(f"  共同下游基因: {len(common_degs)}个")
            
            if len(common_degs) > 0:
                print(f"\n  共同下游基因示例 (Top 10):")
                for gene in list(common_degs)[:10]:
                    nfkb1_fc = pred_logFC_df.loc[gene, 'Nfkb1']
                    fdx1_fc = pred_logFC_df.loc[gene, 'Fdx1']
                    print(f"    {gene}: NFKB1_logFC={nfkb1_fc:.3f}, FDX1_logFC={fdx1_fc:.3f}")
                
                # 保存共同下游基因
                pd.Series(list(common_degs)).to_csv(
                    os.path.join(KNOCKDOWN_DIR, f"common_genes_{group}_v3.csv"),
                    index=False, header=False
                )
            
            jaccard = len(common_degs) / len(nfkb1_degs | fdx1_degs) if len(nfkb1_degs | fdx1_degs) > 0 else 0
            print(f"\n  Jaccard相似性指数: {jaccard:.4f}")
        
        # 保存DEG统计
        deg_stats = {}
        for target in all_pred_logFC.keys():
            up = (pred_logFC_df[target] > 0.25).sum()
            down = (pred_logFC_df[target] < -0.25).sum()
            deg_stats[target] = {'up': up, 'down': down, 'total': up + down}
        deg_df = pd.DataFrame(deg_stats).T
        deg_df.to_csv(os.path.join(KNOCKDOWN_DIR, f"knockdown_stats_{group}_v3.csv"))
        print(f"\n各基因敲除后差异基因统计:")
        print(deg_df)
        
        # 可视化
        print(f"\n生成{group}组可视化...")
        
        # Target间相关性热图
        if len(all_target_corrs) > 1:
            plt.figure(figsize=(8, 6))
            corr_for_plot = pd.DataFrame({t: all_target_corrs[t] for t in all_target_corrs.keys()}).loc[targets, targets]
            sns.heatmap(corr_for_plot, annot=True, cmap='RdBu_r', center=0, vmin=-1, vmax=1, square=True)
            plt.title(f'{group} Microglia: Target Gene Correlation (V3)')
            plt.tight_layout()
            plt.savefig(os.path.join(KNOCKDOWN_DIR, f"target_correlation_{group}_v3.png"), dpi=300)
            plt.close()
        
        # 敲除效应热图 (Top 50 genes)
        top_genes = pred_logFC_df.abs().sum(axis=1).sort_values(ascending=False).head(50).index
        plt.figure(figsize=(10, 12))
        sns.heatmap(pred_logFC_df.loc[top_genes, list(all_pred_logFC.keys())], 
                   cmap='RdBu_r', center=0, vmin=-2, vmax=2)
        plt.title(f'{group} Microglia: Predicted Knockdown Effects (Top 50)')
        plt.tight_layout()
        plt.savefig(os.path.join(KNOCKDOWN_DIR, f"knockdown_heatmap_{group}_v3.png"), dpi=300)
        plt.close()
    else:
        print("⚠️ 没有成功计算的target")
    
    # 清理内存
    gc.collect()
    print(f"\n当前内存: {get_memory_usage():.1f} MB")

print("\n" + "="*70)
print("Microglia虚拟敲除分析 V3 完成!")
print("="*70)
print(f"结果保存在: {KNOCKDOWN_DIR}")
print(f"最终内存使用: {get_memory_usage():.1f} MB")
