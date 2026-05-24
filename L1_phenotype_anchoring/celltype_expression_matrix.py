#!/usr/bin/env python3
"""
权威方法：细胞类型特异性基因筛选 (Top ≤200)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
方法1 (主): scanpy rank_genes_groups (Wilcoxon) → |log2FC|排序 → Top200
方法2 (辅): Pseudobulk + DESeq2 (需样本信息)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
参考:
  Squair et al., 2021, Nature Communications
  Crowell et al., 2020, Genome Biology  
  Lee & Han, 2024, Bioinformatics
  sc-best-practices.org (Theis Lab)
"""

import scanpy as sc
import pandas as pd
import numpy as np
import os
import warnings
warnings.filterwarnings('ignore')
np.random.seed(42)

H5AD_PATH = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\results\cuproptosis_singlecell\sc_adata_cuproptosis.h5ad"
RAW_10X_DIR = r"D:\反向网络药理学\L1 数据集\RNA-seq\GSE174574_extracted"
OUTPUT_DIR = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\CIRI-cuproptosis-causal-discovery\results\L1_phenotype_anchoring"
D_DISK_TABLE = r"D:\反向网络药理学\L1 数据集\L1结果（只允许放表格和图，且用时间标注）\表格\GSE174574_scRNA_24h"
D_DISK_RNASEQ = r"D:\反向网络药理学\L1 数据集\RNA-seq"

TOP_K = 200
MIN_CELLS_PER_GROUP = 10
LOG2FC_THRESHOLD = 0.1
PADJ_THRESHOLD = 0.05

EXCLUDE_TYPES = ['Unknown']

print("=" * 70)
print("权威方法: 细胞类型特异性基因筛选 (Top ≤200)")
print("=" * 70)
print()
print("参考:")
print("  [1] Squair et al., 2021, Nat Commun — Pseudobulk是金标准")
print("  [2] Crowell et al., 2020, Genome Biol — 单细胞DE需伪批量")
print("  [3] Lee & Han, 2024, Bioinformatics — 数学证明伪批量等价GLMM")
print("  [4] sc-best-practices.org (Theis Lab) — pseudobulk+DESeq2")
print("  [5] 哈佛HBC Training — 专用Pseudobulk DESeq2教程")
print()

# ===========================================================================
# 方法1: scanpy rank_genes_groups (Wilcoxon) 标准化筛选
# ===========================================================================
print("=" * 70)
print("方法1: Wilcoxon rank-sum test per cell type, Top200 ranking")
print("=" * 70)

print(f"\n[1/4] 加载数据...")
adata = sc.read_h5ad(H5AD_PATH)
print(f"  数据: {adata.n_obs} cells x {adata.n_vars} genes")

if adata.raw is not None:
    print("  使用 .raw 层（原始计数→标准化+log1p）")
    work_adata = adata.raw.to_adata()
    sc.pp.normalize_total(work_adata, target_sum=1e4)
    sc.pp.log1p(work_adata)
    work_adata.obs['condition'] = adata.obs['condition'].values
    work_adata.obs['cell_type'] = adata.obs['cell_type'].values
else:
    work_adata = adata.copy()

cell_types_all = work_adata.obs['cell_type'].value_counts()
cell_types = [ct for ct in cell_types_all.index if ct not in EXCLUDE_TYPES]
print(f"  细胞类型: {cell_types}")

print(f"\n[2/4] 对每种细胞类型: Wilcoxon rank_genes_groups...")
de_results = {}
for ct in cell_types:
    ct_mask = (work_adata.obs['cell_type'] == ct).values
    ct_adata = work_adata[ct_mask, :].copy()

    n_mcao = (ct_adata.obs['condition'].str.lower().str.strip() == 'mcao').sum()
    n_sham = (ct_adata.obs['condition'].str.lower().str.strip() == 'sham').sum()
    n_total = ct_adata.n_obs

    if n_mcao < MIN_CELLS_PER_GROUP or n_sham < MIN_CELLS_PER_GROUP:
        print(f"  {ct}: ⚠ 跳过 (MCAO={n_mcao}, Sham={n_sham})")
        continue

    print(f"  {ct}: {n_total} cells (MCAO={n_mcao}, Sham={n_sham})", end=" ")
    try:
        sc.tl.rank_genes_groups(
            ct_adata, groupby='condition',
            reference='Sham', method='wilcoxon',
            n_genes=TOP_K + 50
        )
        de_df = sc.get.rank_genes_groups_df(ct_adata, group='MCAO')
        de_df['cell_type'] = ct
        de_df['abs_log2fc'] = de_df['logfoldchanges'].abs()
        de_df = de_df.sort_values('abs_log2fc', ascending=False).head(TOP_K)
        de_df = de_df.reset_index(drop=True)
        de_results[ct] = {
            'de_df': de_df,
            'n_cells': n_total, 'n_mcao': n_mcao, 'n_sham': n_sham
        }
        n_sig = ((de_df['logfoldchanges'].abs() > LOG2FC_THRESHOLD) &
                 (de_df['pvals_adj'] < PADJ_THRESHOLD)).sum()
        print(f"→ {len(de_df)} genes (|log2FC|>0.1 & padj<0.05: {n_sig})")
    except Exception as e:
        print(f"  FAILED: {e}")
        continue

print(f"\n  完成 {len(de_results)} / {len(cell_types)} 细胞类型")

print(f"\n[3/4] 构建表达式矩阵...")

all_selected_genes = set()
for ct, res in de_results.items():
    all_selected_genes.update(res['de_df']['names'].tolist())

print(f"  基因池: {len(all_selected_genes)} unique genes")

expr_matrix_data = {}
for ct, res in de_results.items():
    ct_mask = (work_adata.obs['cell_type'] == ct).values
    ct_adata = work_adata[ct_mask, :]

    selected_genes = res['de_df']['names'].tolist()
    found_genes = [g for g in selected_genes if g in ct_adata.var_names]

    if not found_genes:
        continue

    ct_subset = ct_adata[:, found_genes]
    X = ct_subset.X
    if hasattr(X, 'toarray'):
        X = X.toarray()

    mcao_mask = (ct_adata.obs['condition'].str.lower().str.strip() == 'mcao').values
    sham_mask = (ct_adata.obs['condition'].str.lower().str.strip() == 'sham').values

    for i, gene in enumerate(found_genes):
        if gene not in expr_matrix_data:
            expr_matrix_data[gene] = {}
        col = X[:, i]
        expr_matrix_data[gene][f"{ct}_MCAO_mean"] = round(float(col[mcao_mask].mean()), 6)
        expr_matrix_data[gene][f"{ct}_Sham_mean"] = round(float(col[sham_mask].mean()), 6)
        expr_matrix_data[gene][f"{ct}_log2FC"] = round(
            float(col[mcao_mask].mean() - col[sham_mask].mean()), 6
        )
        expr_matrix_data[gene][f"{ct}_MCAO_pct"] = round(
            float((col[mcao_mask] > 0).mean() * 100), 2
        )
        expr_matrix_data[gene][f"{ct}_Sham_pct"] = round(
            float((col[sham_mask] > 0).mean() * 100), 2
        )

ct_order_final = [ct for ct in ['Neuron', 'OPC', 'Oligodendrocyte', 'Astrocyte',
                                   'Microglia', 'Pericyte', 'Endothelial', 'Ependymal']
                  if ct in de_results]

print(f"\n[4/4] 写入Excel...")

xlsx_path = os.path.join(OUTPUT_DIR, "L1_CellType_ExpressionMatrix.xlsx")

with pd.ExcelWriter(xlsx_path, engine='openpyxl') as writer:
    summary_rows = []
    for ct in ct_order_final:
        res = de_results[ct]
        n_sig = ((res['de_df']['logfoldchanges'].abs() > LOG2FC_THRESHOLD) &
                 (res['de_df']['pvals_adj'] < PADJ_THRESHOLD)).sum()
        summary_rows.append({
            'cell_type': ct,
            'n_cells': res['n_cells'],
            'n_mcao': res['n_mcao'],
            'n_sham': res['n_sham'],
            'n_genes_selected': len(res['de_df']),
            'n_sig': n_sig,
            'method': 'Wilcoxon rank-sum (scanpy)',
            'ranking': '|log2FC| descending'
        })
    pd.DataFrame(summary_rows).to_excel(writer, sheet_name='Summary', index=False)

    for ct in ct_order_final:
        res = de_results[ct]
        df_out = res['de_df'].copy()
        df_out['rank'] = range(1, len(df_out) + 1)
        df_out = df_out[['rank', 'names', 'logfoldchanges', 'pvals', 'pvals_adj',
                          'scores', 'cell_type']]
        df_out.columns = ['rank', 'gene', 'log2FC', 'p_value', 'p_adjust',
                           'wilcoxon_score', 'cell_type']
        sheet_name = f"{ct}_Top{TOP_K}"[:31]
        df_out.to_excel(writer, sheet_name=sheet_name, index=False)

    df_matrix = pd.DataFrame.from_dict(expr_matrix_data, orient='index')
    df_matrix.index.name = 'gene'
    df_matrix = df_matrix.sort_index()
    df_matrix.to_excel(writer, sheet_name='Expression_Matrix')

    log2fc_cols = [c for c in df_matrix.columns if '_log2FC' in c]
    mcao_cols = [c for c in df_matrix.columns if '_MCAO_mean' in c]
    sham_cols = [c for c in df_matrix.columns if '_Sham_mean' in c]
    mcao_pct_cols = [c for c in df_matrix.columns if '_MCAO_pct' in c]
    sham_pct_cols = [c for c in df_matrix.columns if '_Sham_pct' in c]

    if log2fc_cols:
        sub = df_matrix[log2fc_cols].copy()
        sub.columns = [c.replace('_log2FC', '') for c in sub.columns]
        sub.to_excel(writer, sheet_name='log2FC_Matrix')
    if mcao_cols:
        sub = df_matrix[mcao_cols].copy()
        sub.columns = [c.replace('_MCAO_mean', '') for c in sub.columns]
        sub.to_excel(writer, sheet_name='MCAO_MeanExpr')
    if sham_cols:
        sub = df_matrix[sham_cols].copy()
        sub.columns = [c.replace('_Sham_mean', '') for c in sub.columns]
        sub.to_excel(writer, sheet_name='Sham_MeanExpr')
    if mcao_pct_cols:
        sub = df_matrix[mcao_pct_cols].copy()
        sub.columns = [c.replace('_MCAO_pct', '') for c in sub.columns]
        sub.to_excel(writer, sheet_name='MCAO_PctExpr')
    if sham_pct_cols:
        sub = df_matrix[sham_pct_cols].copy()
        sub.columns = [c.replace('_Sham_pct', '') for c in sub.columns]
        sub.to_excel(writer, sheet_name='Sham_PctExpr')

print(f"\n{'='*70}")
print("输出完成")
print("=" * 70)
for ct in ct_order_final:
    res = de_results[ct]
    top3_genes = res['de_df']['names'].head(3).tolist()
    top3_lfc = res['de_df']['logfoldchanges'].head(3).tolist()
    print(f"\n  {ct} (n={res['n_cells']}, Top {len(res['de_df'])}):")
    for g, lfc in zip(top3_genes, top3_lfc):
        print(f"    {g:15s}  log2FC={lfc:+.4f}")

print(f"\n✅ 主输出: {xlsx_path}")

for dst_dir in [D_DISK_TABLE, D_DISK_RNASEQ]:
    try:
        import shutil
        dst_path = os.path.join(dst_dir, os.path.basename(xlsx_path))
        shutil.copy2(xlsx_path, dst_path)
        print(f"✅ 已同步: {dst_path}")
    except Exception as e:
        print(f"⚠ 同步失败 {dst_dir}: {e}")

print(f"\n{'='*70}")
print("全部分析完成")
print("=" * 70)
print(f"""
方法与局限性说明:
  ✓ 方法: scanpy rank_genes_groups (Wilcoxon) → |log2FC|排序 → Top {TOP_K}
  ✓ 参考: scanpy官方文档、sc-best-practices.org
  ⚠ 局限性: 单细胞Wilcoxon存在伪重复偏差(Squair 2021)
    - 极限情况下p值可能膨胀10^3-10^4倍
    - 本方法仅用于基因排序(Ranking)，非假设检验
    - 理想方法Pseudobulk+DESeq2需6个独立样本ID，但此h5ad仅有batch=0/1
""")

# ===========================================================================
# 方法2: Pseudobulk + DESeq2 验证 (使用原始10X数据)
# ===========================================================================
print("\n" + "=" * 70)
print("方法2: Pseudobulk + DESeq2 (需要样本信息, 仅显示可行性)")
print("=" * 70)
print("  由于h5ad中batch=0(Sham)/1(MCAO)丢失了6个独立样本信息,")
print("  无法执行真正的Pseudobulk DESeq2.")
print("  若需要，可从原始10X数据(C_10X_RAW_DIR)重建样本注释.")
print(f"  原始10X数据路径: {RAW_10X_DIR}")