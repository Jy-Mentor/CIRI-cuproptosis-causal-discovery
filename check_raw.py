#!/usr/bin/env python3
"""
从 scRNA-seq .h5ad 文件中提取铜死亡基因表达
遵循权威方法：
- Scanpy 官方流程 (Wolf et al. 2018, Genome Biology)
- Luecken & Theis 2019 (Mol Syst Biol) 
- 使用 .raw 原始计数 → 重新标准归一化 → 计算 log2FC
"""

import scanpy as sc
import pandas as pd
import numpy as np

# 铜死亡基因列表 (小鼠基因名，小写)
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

print(f"Shape: {adata.shape}")
print(f"Has raw: {adata.raw is not None}")
if adata.raw is not None:
    print(f"Raw shape: {adata.raw.shape}")
    print(f"Raw X type: {type(adata.raw.X)}")
    print(f"Raw X min/max: {adata.raw.X.min():.2f}/{adata.raw.X.max():.2f}")
else:
    print("No raw attribute found!")

print(f"Layers: {list(adata.layers.keys())}")
print(f"X min/max: {adata.X.min():.4f}/{adata.X.max():.4f}")
print(f"Condition: {adata.obs['condition'].value_counts().to_dict()}")

# 检查原始计数中铜死亡基因的存在
if adata.raw is not None:
    raw_genes = [g for g in adata.raw.var_names]
    for gene in cuproptosis_genes:
        matches = [g for g in raw_genes if g.upper() == gene.upper()]
        if matches:
            print(f"  {gene} found in raw: {matches[0]}")
        else:
            print(f"  {gene} NOT found in raw")
    
    # 检查Fdx1在raw中的表达
    if 'Fdx1' in adata.raw.var_names:
        fdx1_raw = adata.raw[:, 'Fdx1'].X.toarray().flatten()
        print(f"\nFdx1 raw counts: min={fdx1_raw.min()}, max={fdx1_raw.max()}, mean={fdx1_raw.mean():.4f}")
        print(f"  Fdx1 raw > 0: {(fdx1_raw > 0).mean()*100:.1f}%")
else:
    # 如果没有raw，从adata.var_names中搜索 (不区分大小写)
    var_genes = [g.upper() for g in adata.var_names]
    for gene in cuproptosis_genes:
        if gene.upper() in var_genes:
            idx = var_genes.index(gene.upper())
            print(f"  {gene} found in var at position {idx}: '{adata.var_names[idx]}'")