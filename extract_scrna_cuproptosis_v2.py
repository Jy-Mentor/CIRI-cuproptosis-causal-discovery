#!/usr/bin/env python3
"""
从 scRNA-seq .h5ad 文件中提取所有铜死亡基因的表达均值
使用标准化后的数据计算差异（不取log2，直接计算标准化均值差异）
"""

import scanpy as sc
import pandas as pd
import numpy as np

# 铜死亡基因列表
cuproptosis_genes = [
    "Fdx1", "Lias", "Dld", "Dlat", "Dlst", "Pdha1", "Pdhb", "Gls",
    "Gcsh", "Lipt1", "Lipt2", "Cdkn2a", "Nfe2l2", "Nlrp3",
    "Slc31a1", "Slc31a2", "Slc11a2", "Steap3", "Atp7a", "Atp7b",
    "Atox1", "Ccs", "Cox17", "Cox11", "Sco1", "Sco2",
    "Mt1a", "Mt2a", "Alb", "Cp", "Sod1", "Sod3",
    "Commd1", "Mtf1"
]

# 加载 scRNA-seq 数据
h5ad_path = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\results\cuproptosis_singlecell\sc_adata_cuproptosis.h5ad"
print(f"Loading {h5ad_path}...")
adata = sc.read_h5ad(h5ad_path)

print(f"Data shape: {adata.shape}")

# 计算每个组的平均表达（标准化数据）
results = []
for gene in cuproptosis_genes:
    gene_upper = gene.upper()
    matching_genes = [g for g in adata.var_names if g.upper() == gene_upper]
    
    if not matching_genes:
        print(f"  Warning: {gene} not found")
        results.append({
            'gene': gene,
            'sham_mean': np.nan,
            'mcao_mean': np.nan,
            'diff': np.nan,
            'log2FC': np.nan,
            'found': False
        })
        continue
    
    gene_name = matching_genes[0]
    expr = adata[:, gene_name].X.toarray().flatten()
    
    sham_mask = adata.obs['condition'].str.lower() == 'sham'
    mcao_mask = adata.obs['condition'].str.lower() == 'mcao'
    
    sham_mean = np.mean(expr[sham_mask])
    mcao_mean = np.mean(expr[mcao_mask])
    
    # 标准化数据的差异 = mcao - sham（不是log2FC）
    diff = mcao_mean - sham_mean
    
    # 同时计算伪log2FC（用于与Bulk比较）
    # 将标准化数据平移到正值后计算
    min_val = min(sham_mean, mcao_mean)
    offset = abs(min_val) + 0.01 if min_val < 0 else 0.01
    pseudo_log2fc = np.log2((mcao_mean + offset) / (sham_mean + offset))
    
    results.append({
        'gene': gene,
        'sham_mean': sham_mean,
        'mcao_mean': mcao_mean,
        'diff': diff,
        'log2FC': pseudo_log2fc,
        'found': True
    })
    
    print(f"  {gene}: sham={sham_mean:.4f}, mcao={mcao_mean:.4f}, diff={diff:.4f}, log2FC={pseudo_log2fc:.4f}")

# 保存结果
df = pd.DataFrame(results)
output_path = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\CIRI-cuproptosis-causal-discovery\results\L1_phenotype_anchoring\scRNA_cuproptosis_all_genes.csv"
df.to_csv(output_path, index=False)
print(f"\nSaved to {output_path}")
print(f"Found {df['found'].sum()}/{len(cuproptosis_genes)} genes")
