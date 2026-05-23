#!/usr/bin/env python3
"""
从 .h5ad 提取各细胞类型铜死亡基因表达数据
用于 R 绘图：
  图2: 细胞类型特异性铜死亡差异表达谱（上调/下调分群）
  图3: 各亚群差异小提琴图
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
output_dir = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\CIRI-cuproptosis-causal-discovery\results\L1_phenotype_anchoring"

print(f"Loading {h5ad_path}...")
adata = sc.read_h5ad(h5ad_path)

if adata.raw is None:
    raise ValueError("adata.raw is None")
raw_adata = adata.raw.to_adata()

gene_name_map = {}
for gene in cuproptosis_genes:
    matches = [g for g in raw_adata.var_names if g.upper() == gene.upper()]
    if matches:
        gene_name_map[gene] = matches[0]

found_genes = list(gene_name_map.values())
print(f"Found {len(found_genes)}/{len(cuproptosis_genes)} cuproptosis genes")

adata_raw_subset = raw_adata[:, found_genes].copy()
adata_raw_subset.obs['condition'] = adata.obs['condition'].values
adata_raw_subset.obs['cell_type'] = adata.obs['cell_type'].values

sc.pp.normalize_total(adata_raw_subset, target_sum=1e4)
sc.pp.log1p(adata_raw_subset)

cell_types = adata_raw_subset.obs['cell_type'].unique().tolist()
print(f"Cell types: {cell_types}")

results_rows = []
violin_data_rows = []

for ct in cell_types:
    ct_mask = adata_raw_subset.obs['cell_type'] == ct
    ct_adata = adata_raw_subset[ct_mask.values, :]
    n_cells = ct_adata.n_obs

    mcao_mask = ct_adata.obs['condition'].str.lower() == 'mcao'
    sham_mask = ct_adata.obs['condition'].str.lower() == 'sham'

    n_mcao = mcao_mask.sum()
    n_sham = sham_mask.sum()

    if n_mcao < 3 or n_sham < 3:
        print(f"  {ct}: 跳过（MCAO={n_mcao}, Sham={n_sham}）")
        continue

    for gene_name, raw_gene_name in gene_name_map.items():
        gene_idx = list(ct_adata.var_names).index(raw_gene_name)
        expr_values = ct_adata.X[:, gene_idx].toarray().flatten()

        mcao_expr = expr_values[mcao_mask.values]
        sham_expr = expr_values[sham_mask.values]

        mcao_mean = float(np.mean(mcao_expr))
        sham_mean = float(np.mean(sham_expr))
        log2fc_val = (mcao_mean - sham_mean) / np.log(2)

        try:
            stat, p_val = mannwhitneyu(mcao_expr, sham_expr, alternative='two-sided')
        except Exception:
            stat, p_val = np.nan, np.nan

        results_rows.append({
            'cell_type': ct,
            'gene': gene_name,
            'n_cells': n_cells,
            'n_mcao': n_mcao,
            'n_sham': n_sham,
            'sham_mean': round(sham_mean, 6),
            'mcao_mean': round(mcao_mean, 6),
            'log2FC': round(log2fc_val, 6),
            'wilcoxon_stat': round(stat, 2) if not np.isnan(stat) else np.nan,
            'p_value': p_val,
            'pct_mcao': round(float(np.mean(mcao_expr > 0) * 100), 1),
            'pct_sham': round(float(np.mean(sham_expr > 0) * 100), 1)
        })

        for cell_idx, val in enumerate(expr_values):
            violin_data_rows.append({
                'cell_type': ct,
                'gene': gene_name,
                'condition': 'MCAO' if ct_adata.obs['condition'].values[cell_idx].lower() == 'mcao' else 'Sham',
                'expression': val
            })

    print(f"  {ct}: {n_cells} cells (MCAO={n_mcao}, Sham={n_sham})")

df_results = pd.DataFrame(results_rows)
df_results['p_adjust'] = np.nan
for ct in df_results['cell_type'].unique():
    ct_mask = df_results['cell_type'] == ct
    pvals = df_results.loc[ct_mask, 'p_value'].values
    df_results.loc[ct_mask, 'p_adjust'] = pd.Series(
        np.where(~np.isnan(pvals), 
                 pd.Series(pvals).fillna(1).apply(
                     lambda x: min(x * sum(ct_mask) * 2 / (sum(~np.isnan(pvals)) or 1), 1.0) if not np.isnan(x) else np.nan
                 ),
                 np.nan)
    ).values

df_violin = pd.DataFrame(violin_data_rows)

df_results.to_csv(f"{output_dir}/celltype_cuproptosis_DEGs.csv", index=False)
df_violin.to_csv(f"{output_dir}/celltype_cuproptosis_violin_data.csv", index=False)

print(f"\nResults saved:")
print(f"  {output_dir}/celltype_cuproptosis_DEGs.csv ({len(df_results)} rows)")
print(f"  {output_dir}/celltype_cuproptosis_violin_data.csv ({len(df_violin)} rows)")
print(f"\nSummary per cell type:")
for ct in df_results['cell_type'].unique():
    sub = df_results[df_results['cell_type'] == ct]
    n_up = sum((sub['log2FC'] > 0.25) & (sub['p_adjust'] < 0.05))
    n_down = sum((sub['log2FC'] < -0.25) & (sub['p_adjust'] < 0.05))
    print(f"  {ct}: up={n_up}, down={n_down}")