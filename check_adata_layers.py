#!/usr/bin/env python3
import scanpy as sc

h5ad_path = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\results\cuproptosis_singlecell\sc_adata_cuproptosis.h5ad"
adata = sc.read_h5ad(h5ad_path)

print("Layers:", list(adata.layers.keys()))
print("X min/max:", adata.X.min(), adata.X.max())
print("X mean:", adata.X.mean())

if 'counts' in adata.layers:
    print("counts min/max:", adata.layers['counts'].min(), adata.layers['counts'].max())
if 'raw' in adata.layers:
    print("raw min/max:", adata.layers['raw'].min(), adata.layers['raw'].max())

# 检查一个基因的分布
import numpy as np
gene_expr = adata[:, 'Fdx1'].X.toarray().flatten()
print(f"\nFdx1 expression stats:")
print(f"  min={gene_expr.min():.4f}, max={gene_expr.max():.4f}")
print(f"  mean={gene_expr.mean():.4f}, median={np.median(gene_expr):.4f}")
print(f"  sham mean={gene_expr[adata.obs['condition'] == 'Sham'].mean():.4f}")
print(f"  mcao mean={gene_expr[adata.obs['condition'] == 'MCAO'].mean():.4f}")
