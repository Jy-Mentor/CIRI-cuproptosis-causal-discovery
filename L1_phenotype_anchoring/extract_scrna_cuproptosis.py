#!/usr/bin/env python3
"""
从 scRNA-seq .h5ad 文件中提取铜死亡基因表达

权威方法参考：
  - Wolf et al. 2018, "SCANPY: large-scale single-cell gene expression data analysis"
    Genome Biology, DOI: 10.1186/s13059-017-1382-0
  - Luecken & Theis 2019, "Current best practices in single-cell RNA-seq analysis: a tutorial"
    Molecular Systems Biology, DOI: 10.15252/msb.20188746
  - Scanpy 官方文档: https://scanpy.readthedocs.io/

标准流程：
  1. 从 .raw 提取原始整数计数
  2. sc.pp.normalize_total(adata, target_sum=1e4)  — 文库大小归一化
  3. sc.pp.log1p(adata)                              — 对数变换 ln(1 + x)
  4. 按组计算 log1p 均值，差值除以 ln(2) 得 log2FC
  5. Wilcoxon 秩和检验计算 p 值
"""

import scanpy as sc
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
import warnings
warnings.filterwarnings('ignore')

cuproptosis_genes = [
    "Fdx1", "Lias", "Dld", "Dlat", "Dlst", "Pdha1", "Pdhb", "Gls",
    "Gcsh", "Lipt1", "Lipt2", "Cdkn2a", "Nfe2l2", "Nlrp3",
    "Slc31a1", "Slc31a2", "Slc11a2", "Steap3", "Atp7a", "Atp7b",
    "Atox1", "Ccs", "Cox17", "Cox11", "Sco1", "Sco2",
    "Mt1a", "Mt2a", "Alb", "Cp", "Sod1", "Sod3",
    "Commd1", "Mtf1"
]

h5ad_path = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\results\cuproptosis_singlecell\sc_adata_cuproptosis.h5ad"
print(f"Loading {h5ad_path}...")
adata = sc.read_h5ad(h5ad_path)

if adata.raw is None:
    raise ValueError("adata.raw is None — 无法获取原始计数")

print(f"Shape: {adata.shape}")
print(f"Raw shape: {adata.raw.shape}")
print(f"Raw X min/max: {adata.raw.X.min():.0f}/{adata.raw.X.max():.0f}")
print(f"Condition distribution: {adata.obs['condition'].value_counts().to_dict()}")

# 从 .raw 提取原始计数（整数矩阵），仅保留目标基因
raw_adata = adata.raw.to_adata()

gene_name_map = {}
for gene in cuproptosis_genes:
    matches = [g for g in raw_adata.var_names if g.upper() == gene.upper()]
    if matches:
        gene_name_map[gene] = matches[0]

found_genes = list(gene_name_map.values())
print(f"\nFound {len(found_genes)}/{len(cuproptosis_genes)} cuproptosis genes in raw")

# 提取原始计数子集
adata_raw_subset = raw_adata[:, found_genes].copy()
adata_raw_subset.obs['condition'] = adata.obs['condition'].values

print(f"\nSubset shape: {adata_raw_subset.shape}")
print(f"Raw counts range: [{adata_raw_subset.X.min():.0f}, {adata_raw_subset.X.max():.0f}]")

# ============================================================
# 权威标准化流程 (Wolf et al. 2018; Luecken & Theis 2019)
# ============================================================
adata_norm = adata_raw_subset.copy()

# Step 1: 文库大小归一化 — 每个细胞归一化到 10,000 总计数
sc.pp.normalize_total(adata_norm, target_sum=1e4)

# Step 2: 对数变换 — ln(1 + x)
sc.pp.log1p(adata_norm)

print(f"\nNormalized data range: [{adata_norm.X.min():.4f}, {adata_norm.X.max():.4f}]")

# ============================================================
# 计算每组均值和 log2FC
# ============================================================
mcao_mask = adata_norm.obs['condition'].str.lower() == 'mcao'
sham_mask = adata_norm.obs['condition'].str.lower() == 'sham'

print(f"\nMCAO cells: {mcao_mask.sum()}, Sham cells: {sham_mask.sum()}")

results = []
for gene_name, raw_gene_name in gene_name_map.items():
    gene_idx = list(adata_norm.var_names).index(raw_gene_name)
    expr_values = adata_norm.X[:, gene_idx].toarray().flatten()

    mcao_expr = expr_values[mcao_mask.values]
    sham_expr = expr_values[sham_mask.values]

    mcao_mean = float(np.mean(mcao_expr))
    sham_mean = float(np.mean(sham_expr))

    # log2FC = (ln(1+x_mcao) - ln(1+x_sham)) / ln(2)
    log2fc_val = (mcao_mean - sham_mean) / np.log(2)

    # Wilcoxon 秩和检验（标准单细胞差异分析方法）
    try:
        if len(mcao_expr) > 0 and len(sham_expr) > 0:
            stat, p_val = mannwhitneyu(mcao_expr, sham_expr, alternative='two-sided')
        else:
            stat, p_val = np.nan, np.nan
    except Exception:
        stat, p_val = np.nan, np.nan

    results.append({
        'gene': gene_name,
        'sham_mean_log1p': round(sham_mean, 6),
        'mcao_mean_log1p': round(mcao_mean, 6),
        'log2FC': round(log2fc_val, 6),
        'wilcoxon_stat': round(stat, 2) if not np.isnan(stat) else np.nan,
        'p_value': p_val,
        'p_value_log10': -np.log10(max(p_val, 1e-300)) if not np.isnan(p_val) else np.nan,
        'pct_mcao': round(float(np.mean(mcao_expr > 0) * 100), 1),
        'pct_sham': round(float(np.mean(sham_expr > 0) * 100), 1),
        'found': True
    })
    p_str = f"{p_val:.2e}" if not np.isnan(p_val) else "NA"
    print(f"  {gene_name:10s} | Sham: {sham_mean:.4f} | MCAO: {mcao_mean:.4f} | log2FC: {log2fc_val:+.4f} | p: {p_str}")

# 补充未找到的基因
found_set = set(gene_name_map.keys())
for gene in cuproptosis_genes:
    if gene not in found_set:
        results.append({
            'gene': gene,
            'sham_mean_log1p': np.nan,
            'mcao_mean_log1p': np.nan,
            'log2FC': np.nan,
            'wilcoxon_stat': np.nan,
            'p_value': np.nan,
            'p_value_log10': np.nan,
            'pct_mcao': np.nan,
            'pct_sham': np.nan,
            'found': False
        })

df = pd.DataFrame(results)
output_path = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\CIRI-cuproptosis-causal-discovery\results\L1_phenotype_anchoring\scRNA_cuproptosis_all_genes.csv"
df.to_csv(output_path, index=False)

print(f"\n{'='*60}")
print(f"Saved to: {output_path}")
print(f"Found {df['found'].sum()}/{len(cuproptosis_genes)} genes")
print(f"Standard workflow applied: normalize_total(1e4) → log1p → log2FC(Wolf2018, LueckenTheis2019)")
print(f"{'='*60}")