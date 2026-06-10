#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GSE174574 SCISSOR-like 单细胞分析流程 V3 - 正确内存优化版
关键修复: 在scale之前将adata子集化到HVG (3000基因)
"""

import os
import sys
import time
import gc
import psutil
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# 内存监控
class MemoryMonitor:
    def __init__(self):
        self.initial_mem = self.get_mem()
        self.peak_mem = self.initial_mem
        
    def get_mem(self):
        process = psutil.Process(os.getpid())
        return process.memory_info().rss / 1024 / 1024
    
    def checkpoint(self, label):
        current = self.get_mem()
        self.peak_mem = max(self.peak_mem, current)
        delta = current - self.initial_mem
        print(f"[内存] {label}: {current:.1f} MB (Δ{delta:+.1f} MB)")
        return current

sc.settings.verbosity = 2

# ============================================
# 参数配置
# ============================================
WORK_DIR = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙"
OUTPUT_DIR = os.path.join(WORK_DIR, "GSE174574_SCISSOR_Results")
DATA_DIR = r"C:\Users\Jy-Mentor-7\Desktop\虚拟敲除"
os.makedirs(OUTPUT_DIR, exist_ok=True)

MAX_GENES_FOR_MODULE = 3000
HUB_GENES_HUMAN = ["NFKB1", "FDX1", "HSPA5", "HMOX1", "STAT3", 
                   "HIF1A", "TNF", "IL6", "GPX4", "DLAT"]
MCAO_SIGNATURE = ["Il6", "Tnf", "Nfkb1", "Ccl2", "Icam1", "Vcam1", 
                  "Sele", "Ptgs2", "Mmp9", "Hif1a", "Stat3", "Rela",
                  "Hmox1", "Sod2", "Gpx4", "Cat", "Nqo1", "Hspa5"]
SHAM_SIGNATURE = ["Bdnf", "Ngf", "Nt3", "Gria1", "Grin1", "Syn1",
                  "Syp", "Snap25", "Vamp2", "Stx1a", "Cplx1", "Rab3a",
                  "Camk2a", "Creb1", "Arc", "Fos", "Egr1", "Nr4a1"]
MAPPING_FILE = os.path.join(WORK_DIR, "大创", "大鼠 小鼠 人类映射库.txt")

print("="*70)
print("GSE174574 SCISSOR-like 分析 V3 - 正确内存优化版")
print("="*70)

monitor = MemoryMonitor()
monitor.checkpoint("脚本启动")

# ============================================
# 阶段0: 数据读取
# ============================================
print("\n阶段0: 数据读取")
print("="*70)

SAMPLES = {
    'sham1': (os.path.join(DATA_DIR, "GSM5319987_sham1_barcodes.tsv"),
              os.path.join(DATA_DIR, "GSM5319987_sham1_genes.tsv"),
              os.path.join(DATA_DIR, "GSM5319987_sham1_matrix.mtx"), 'Sham'),
    'sham2': (os.path.join(DATA_DIR, "GSM5319988_sham2_barcodes.tsv"),
              os.path.join(DATA_DIR, "GSM5319988_sham2_genes.tsv"),
              os.path.join(DATA_DIR, "GSM5319988_sham2_matrix.mtx"), 'Sham'),
    'sham3': (os.path.join(DATA_DIR, "GSM5319989_sham3_barcodes.tsv"),
              os.path.join(DATA_DIR, "GSM5319989_sham3_genes.tsv"),
              os.path.join(DATA_DIR, "GSM5319989_sham3_matrix.mtx"), 'Sham'),
    'mcao1': (os.path.join(DATA_DIR, "GSM5319990_MCAO1_barcodes.tsv"),
              os.path.join(DATA_DIR, "GSM5319990_MCAO1_genes.tsv"),
              os.path.join(DATA_DIR, "GSM5319990_MCAO1_matrix.mtx"), 'MCAO'),
    'mcao2': (os.path.join(DATA_DIR, "GSM5319991_MCAO2_barcodes.tsv"),
              os.path.join(DATA_DIR, "GSM5319991_MCAO2_genes.tsv"),
              os.path.join(DATA_DIR, "GSM5319991_MCAO2_matrix.mtx"), 'MCAO'),
    'mcao3': (os.path.join(DATA_DIR, "GSM5319992_MCAO3_barcodes.tsv"),
              os.path.join(DATA_DIR, "GSM5319992_MCAO3_genes.tsv"),
              os.path.join(DATA_DIR, "GSM5319992_MCAO3_matrix.mtx"), 'MCAO')
}

print("\n读取样本数据...")
adata_list = []

for sample_name, (barcode_file, feature_file, matrix_file, group) in SAMPLES.items():
    print(f"  读取 {sample_name}...", end=' ')
    
    # 读取基因名
    gene_mapping = {}
    with open(feature_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            gene_mapping[parts[0]] = parts[1] if len(parts) >= 2 else parts[0]
    genes = list(gene_mapping.keys())
    gene_symbols = list(gene_mapping.values())
    
    # 读取barcodes
    with open(barcode_file, 'r') as f:
        barcodes = [line.strip() for line in f]
    
    # 读取MTX矩阵
    with open(matrix_file, 'r') as f:
        lines = f.readlines()
    
    data_start = 0
    for i, line in enumerate(lines):
        if line.startswith('%%'):
            data_start = i + 1
        elif line.strip() and not line.startswith('%'):
            dims = line.strip().split()
            n_rows, n_cols, n_entries = int(dims[0]), int(dims[1]), int(dims[2])
            data_start = i + 1
            break
    
    row_indices, col_indices, data_values = [], [], []
    for line in lines[data_start:]:
        if line.strip():
            parts = line.strip().split()
            row_indices.append(int(parts[0]) - 1)
            col_indices.append(int(parts[1]) - 1)
            data_values.append(int(parts[2]))
    
    from scipy.sparse import coo_matrix
    matrix = coo_matrix((data_values, (row_indices, col_indices)), 
                        shape=(n_rows, n_cols)).T.tocsr()
    
    adata = sc.AnnData(X=matrix)
    adata.var_names = gene_symbols
    adata.var['ensembl_id'] = genes
    adata.obs_names = [f"{sample_name}_{bc}" for bc in barcodes]
    adata.obs['sample'] = sample_name
    adata.obs['group'] = group
    
    # 去重
    if len(adata.var_names) != len(set(adata.var_names)):
        unique_genes = ~adata.var_names.duplicated()
        adata = adata[:, unique_genes].copy()
    
    print(f"{adata.n_obs}×{adata.n_vars}")
    adata_list.append(adata)

print("\n合并所有样本...")
adata = sc.concat(adata_list, label='sample')
print(f"合并后: {adata.n_obs} 细胞 × {adata.n_vars} 基因")
monitor.checkpoint("数据读取完成")

# ============================================
# 阶段1: 质控
# ============================================
print("\n阶段1: 质控")
print("="*70)

adata.var['mt'] = adata.var_names.str.startswith('mt-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

cells_before = adata.n_obs
adata = adata[adata.obs.n_genes_by_counts > 200, :]
adata = adata[adata.obs.n_genes_by_counts < 7500, :]
adata = adata[adata.obs.pct_counts_mt < 10, :]

print(f"质控后: {adata.n_obs} / {cells_before} 细胞")
print(f"  Sham组: {sum(adata.obs['group'] == 'Sham')} 细胞")
print(f"  MCAO组: {sum(adata.obs['group'] == 'MCAO')} 细胞")
monitor.checkpoint("质控完成")

# ============================================
# 阶段2: 标准化和高变基因选择
# ============================================
print("\n阶段2: 标准化和高变基因选择")
print("="*70)

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=MAX_GENES_FOR_MODULE)
monitor.checkpoint("标准化完成")

# 强制包含关键Hub基因
print("\n强制包含关键Hub基因到HVG池...")
essential_hub_genes = ['Nfkb1', 'Fdx1', 'Stat3', 'Hspa5', 'Hmox1', 'Gpx4', 'Hif1a', 'Dlat']
essential_in_data = [g for g in essential_hub_genes if g in adata.var_names]
print(f"关键Hub基因: {essential_in_data}")

for gene in essential_in_data:
    if gene in adata.var_names:
        adata.var.loc[gene, 'highly_variable'] = True

# 如果超过限制，调整
n_hvg = adata.var['highly_variable'].sum()
if n_hvg > MAX_GENES_FOR_MODULE:
    non_essential = adata.var_names[adata.var['highly_variable'] & ~adata.var_names.isin(essential_in_data)]
    keep = list(essential_in_data) + list(non_essential[:MAX_GENES_FOR_MODULE - len(essential_in_data)])
    adata.var['highly_variable'] = adata.var_names.isin(keep)

print(f"最终HVG数: {adata.var['highly_variable'].sum()}")
monitor.checkpoint("Hub基因加入HVG")

# ============================================
# 阶段3: 🔑关键修复 - adata子集化到HVG
# ============================================
print("\n阶段3: 将adata子集化到HVG (关键内存优化)")
print("="*70)
print(f"子集化前: {adata.n_obs} × {adata.n_vars}")

# 🔑 关键修复：adata本身子集化
adata = adata[:, adata.var['highly_variable']].copy()
print(f"子集化后: {adata.n_obs} × {adata.n_vars}")

# 强制垃圾回收
gc.collect()
monitor.checkpoint("adata子集化完成")

# ============================================
# 阶段4: 缩放和降维 (现在在3000基因上安全执行)
# ============================================
print("\n阶段4: 缩放和降维 (安全)")
print("="*70)

sc.pp.scale(adata)  # 现在安全：57k × 3k ≈ 1.4GB if dense
sc.tl.pca(adata, svd_solver='arpack', n_comps=30)
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=20)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.8)
adata.obs['Cluster'] = 'Cluster_' + adata.obs['leiden'].astype(str)
monitor.checkpoint("降维完成")

# ============================================
# 阶段5: 细胞类型注释
# ============================================
print("\n阶段5: 细胞类型注释")
print("="*70)

marker_dict = {
    'Microglia': ['Cx3cr1', 'Tmem119', 'P2ry12', 'C1qa'],
    'Astrocytes': ['Gfap', 'S100b', 'Aqp4', 'Aldh1l1'],
    'Neuron': ['Rbfox3', 'Snap25', 'Syn1', 'Dlg4'],
    'Oligodendrocytes': ['Mbp', 'Mog', 'Plp1', 'Olig2'],
    'Endothelial': ['Pecam1', 'Cldn5', 'Flt1', 'Kdr'],
    'OPC': ['Pdgfra', 'Cspg4', 'Sox10'],
    'Pericytes': ['Rgs5', 'Pdgfrb', 'Acta2']
}

marker_genes_available = {ct: [g for g in genes if g in adata.var_names] 
                          for ct, genes in marker_dict.items() 
                          if any(g in adata.var_names for g in genes)}

# 计算marker scores
for cell_type, genes in marker_genes_available.items():
    sc.tl.score_genes(adata, genes, score_name=f'{cell_type}_score')

# 分配细胞类型
cluster_to_celltype = {}
for cluster in adata.obs['leiden'].unique():
    cluster_mask = adata.obs['leiden'] == cluster
    best_celltype = 'Unknown'
    best_score = -float('inf')
    
    for cell_type in marker_genes_available.keys():
        mean_score = adata.obs.loc[cluster_mask, f'{cell_type}_score'].mean()
        if mean_score > best_score:
            best_score = mean_score
            best_celltype = cell_type
    
    cluster_to_celltype[cluster] = best_celltype

adata.obs['CellType'] = adata.obs['leiden'].astype(str).map(cluster_to_celltype)

print(f"\n细胞类型分布:")
print(adata.obs['CellType'].value_counts())
monitor.checkpoint("细胞类型注释完成")

# ============================================
# 阶段6: Hub基因映射
# ============================================
print("\n阶段6: Hub基因映射")
print("="*70)

hub_genes_in_data = []
if os.path.exists(MAPPING_FILE):
    with open(MAPPING_FILE, 'r') as f:
        data_lines = [l for l in f.readlines() if not l.startswith('#')]
    
    import csv
    reader = csv.DictReader(data_lines, delimiter='\t')
    
    human_to_mouse = {}
    for row in reader:
        hg = row['HUMAN_ORTHOLOG_SYMBOL'].upper().strip()
        mg = row['MOUSE_ORTHOLOG_SYMBOL'].upper().strip()
        if hg and mg and hg != 'N/A' and mg != 'N/A':
            if hg not in human_to_mouse:
                human_to_mouse[hg] = []
            human_to_mouse[hg].extend([g.strip() for g in mg.split('|')])
    
    hub_genes_mouse = list(set(sum([human_to_mouse.get(hg, []) for hg in HUB_GENES_HUMAN], [])))
    
    adata_genes_upper = [g.upper() for g in adata.var_names.tolist()]
    for g in hub_genes_mouse:
        if g.upper() in adata_genes_upper:
            hub_genes_in_data.append(adata.var_names.tolist()[adata_genes_upper.index(g.upper())])
    
    print(f"Hub基因在数据中: {hub_genes_in_data}")
monitor.checkpoint("Hub基因映射完成")

# ============================================
# 阶段7: SCISSOR-like评分
# ============================================
print("\n阶段7: SCISSOR-like评分")
print("="*70)

mcmo_in_data = [g for g in MCAO_SIGNATURE if g in adata.var_names]
sham_in_data = [g for g in SHAM_SIGNATURE if g in adata.var_names]

print(f"MCAO特征基因: {len(mcmo_in_data)}个")
print(f"Sham特征基因: {len(sham_in_data)}个")

if mcmo_in_data:
    sc.tl.score_genes(adata, mcmo_in_data, score_name='MCAO_Score')
if sham_in_data:
    sc.tl.score_genes(adata, sham_in_data, score_name='Sham_Score')

if 'MCAO_Score' in adata.obs.columns and 'Sham_Score' in adata.obs.columns:
    adata.obs['Net_Score'] = adata.obs['MCAO_Score'] - adata.obs['Sham_Score']
    print("✅ Net_Score计算完成")

monitor.checkpoint("SCISSOR评分完成")

# 差异检验
print("\n按细胞类型进行差异检验...")
net_score_results = []

for cell_type in adata.obs['CellType'].unique():
    mcmo_mask = (adata.obs['CellType'] == cell_type) & (adata.obs['group'] == 'MCAO')
    sham_mask = (adata.obs['CellType'] == cell_type) & (adata.obs['group'] == 'Sham')
    
    mcmo_scores = adata.obs.loc[mcmo_mask, 'Net_Score']
    sham_scores = adata.obs.loc[sham_mask, 'Net_Score']
    
    if len(mcmo_scores) > 3 and len(sham_scores) > 3:
        statistic, pvalue = stats.mannwhitneyu(mcmo_scores, sham_scores, alternative='two-sided')
        n1, n2 = len(mcmo_scores), len(sham_scores)
        z = stats.norm.ppf(pvalue / 2) if pvalue > 0 else 0
        effect_size = z / np.sqrt(n1 + n2)
        
        net_score_results.append({
            'CellType': cell_type,
            'MCAO_Median': mcmo_scores.median(),
            'Sham_Median': sham_scores.median(),
            'P_value': pvalue,
            'Effect_Size': effect_size,
            'MCAO_N': n1,
            'Sham_N': n2
        })

if net_score_results:
    results_df = pd.DataFrame(net_score_results)
    from statsmodels.stats.multitest import multipletests
    _, p_adj, _, _ = multipletests(results_df['P_value'], method='bonferroni')
    results_df['P_adj'] = p_adj
    results_df = results_df.sort_values('P_value')
    results_df.to_csv(os.path.join(OUTPUT_DIR, "02_net_score_by_celltype_v3.csv"), index=False)
    print("\n差异检验结果:")
    print(results_df.to_string())

monitor.checkpoint("差异检验完成")

# ============================================
# 阶段8: Hub模块评分
# ============================================
print("\n阶段8: Hub模块评分")
print("="*70)

if hub_genes_in_data:
    sc.tl.score_genes(adata, hub_genes_in_data, score_name='Hub_Module_Score')
    
    nfkb1_genes = [g for g in hub_genes_in_data if 'NFKB' in g.upper()]
    fdx1_genes = [g for g in hub_genes_in_data if 'FDX' in g.upper()]
    
    if nfkb1_genes:
        sc.tl.score_genes(adata, nfkb1_genes, score_name='NFKB1_Score')
    if fdx1_genes:
        sc.tl.score_genes(adata, fdx1_genes, score_name='FDX1_Score')
    
    if 'NFKB1_Score' in adata.obs.columns and 'FDX1_Score' in adata.obs.columns:
        mcmo_mask = adata.obs['group'] == 'MCAO'
        rho, pval = stats.spearmanr(
            adata.obs.loc[mcmo_mask, 'NFKB1_Score'],
            adata.obs.loc[mcmo_mask, 'FDX1_Score']
        )
        print(f"\nMCAO组NFKB1-FDX1相关性: rho={rho:.4f}, p={pval:.2e}")
        
        cor_df = pd.DataFrame({
            'CellType': ['All_Cells'],
            'Group': ['MCAO'],
            'Correlation': [rho],
            'P_value': [pval],
            'N': [sum(mcmo_mask)]
        })
        cor_df.to_csv(os.path.join(OUTPUT_DIR, "03_nfkb1_fdx1_correlation_v3.csv"), index=False)

monitor.checkpoint("Hub模块评分完成")

# ============================================
# 保存结果
# ============================================
print("\n保存结果...")
adata.write(os.path.join(OUTPUT_DIR, "GSE174574_SCISSOR_v3.h5ad"))

print("\n" + "="*70)
print("全量分析 V3 完成!")
print("="*70)
print(f"总细胞数: {adata.n_obs}")
print(f"使用基因数: {adata.n_vars}")
print(f"峰值内存: {monitor.peak_mem:.1f} MB ({monitor.peak_mem/1024:.1f} GB)")
monitor.checkpoint("完成")
