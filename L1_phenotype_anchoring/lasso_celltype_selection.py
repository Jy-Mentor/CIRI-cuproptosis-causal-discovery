#!/usr/bin/env python3
"""
LASSO预筛：对每种细胞类型，用L1正则化Logistic回归筛选Top ≤200基因
输出：细胞类型特异性表达矩阵 Excel
"""

import scanpy as sc
import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegressionCV
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
import warnings
import os

warnings.filterwarnings('ignore')
np.random.seed(42)

H5AD_PATH = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\results\cuproptosis_singlecell\sc_adata_cuproptosis.h5ad"
OUTPUT_DIR = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\CIRI-cuproptosis-causal-discovery\results\L1_phenotype_anchoring"
D_DISK_TABLE = r"D:\反向网络药理学\L1 数据集\L1结果（只允许放表格和图，且用时间标注）\表格\GSE174574_scRNA_24h"
D_DISK_RNASEQ = r"D:\反向网络药理学\L1 数据集\RNA-seq"

TOP_K = 200
MIN_CELLS_PER_GROUP = 10
HVG_N_GENES = 2000
N_ALPHAS = 20
CV_FOLDS = 3

print("=" * 60)
print("LASSO 细胞类型特异性基因预筛 (Top ≤200)")
print("=" * 60)

print(f"\n加载 {H5AD_PATH}...")
adata = sc.read_h5ad(H5AD_PATH)
print(f"  数据: {adata.n_obs} cells x {adata.n_vars} genes")

if adata.raw is not None:
    raw_adata = adata.raw.to_adata()
    print("  使用 .raw 层（原始计数）")
else:
    raw_adata = adata.copy()
    print("  使用当前 adata（无 .raw 层）")

condition_key = 'condition'
cell_type_key = 'cell_type'

n_conditions = raw_adata.obs[condition_key].unique()
print(f"  条件: {n_conditions}")

cell_types_all = raw_adata.obs[cell_type_key].value_counts()
print(f"\n细胞类型分布:")
for ct, n in cell_types_all.items():
    print(f"  {ct}: {n} cells")

exclude_types = ['Unknown']
cell_types = [ct for ct in cell_types_all.index if ct not in exclude_types]
print(f"\n分析细胞类型: {cell_types}")

results_by_ct = {}
lasso_genes_all = set()

for ct in cell_types:
    print(f"\n{'='*50}")
    print(f"  细胞类型: {ct}")
    ct_mask = (raw_adata.obs[cell_type_key] == ct).values
    ct_adata = raw_adata[ct_mask, :].copy()

    n_total = ct_adata.n_obs

    cond_mask_mcao = ct_adata.obs[condition_key].str.lower().str.strip() == 'mcao'
    cond_mask_sham = ct_adata.obs[condition_key].str.lower().str.strip() == 'sham'

    n_mcao = cond_mask_mcao.sum()
    n_sham = cond_mask_sham.sum()
    print(f"    细胞数: total={n_total}, MCAO={n_mcao}, Sham={n_sham}")

    if n_mcao < MIN_CELLS_PER_GROUP or n_sham < MIN_CELLS_PER_GROUP:
        print(f"    ⚠ 跳过（每条件<{MIN_CELLS_PER_GROUP}个细胞）")
        continue

    y = np.where(cond_mask_mcao.values, 1, 0)

    print(f"    标准化+log1p...")
    sc.pp.normalize_total(ct_adata, target_sum=1e4)
    sc.pp.log1p(ct_adata)

    print(f"    筛选高变基因 (top {HVG_N_GENES})...")
    try:
        sc.pp.highly_variable_genes(ct_adata, n_top_genes=HVG_N_GENES, flavor='seurat_v3')
    except Exception:
        sc.pp.highly_variable_genes(ct_adata, n_top_genes=HVG_N_GENES)

    n_hvg = ct_adata.var['highly_variable'].sum()
    print(f"    高变基因数: {n_hvg}")

    if n_hvg < 10:
        print(f"    ⚠ 高变基因不足，跳过")
        continue

    X_dense = ct_adata[:, ct_adata.var['highly_variable']].X
    if hasattr(X_dense, 'toarray'):
        X_dense = X_dense.toarray()
    X_dense = np.asarray(X_dense, dtype=np.float64)

    n_features = X_dense.shape[1]
    print(f"    特征矩阵: {X_dense.shape[0]} cells x {n_features} genes")

    if np.any(np.isnan(X_dense)) or np.any(np.isinf(X_dense)):
        X_dense = np.nan_to_num(X_dense, nan=0.0, posinf=0.0, neginf=0.0)

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_dense)

    print(f"    LASSO CV (L1, Cs=50, cv={CV_FOLDS})...")
    cs_values = np.logspace(-4, 1, 50)

    kf = StratifiedKFold(n_splits=CV_FOLDS, shuffle=True, random_state=42)

    model = LogisticRegressionCV(
        Cs=cs_values,
        penalty='l1',
        solver='saga',
        cv=kf,
        max_iter=5000,
        tol=1e-4,
        n_jobs=-1,
        random_state=42,
        class_weight='balanced'
    )

    try:
        model.fit(X_scaled, y)
    except Exception as e:
        print(f"    ⚠ LASSO失败: {e}")
        continue

    C_optimal = model.C_[0]
    coef = model.coef_.flatten()
    print(f"    最优 C: {C_optimal:.6f}")
    print(f"    非零系数基因数: {np.sum(np.abs(coef) > 1e-8)}")

    hvg_names = ct_adata.var_names[ct_adata.var['highly_variable']].tolist()

    coef_df = pd.DataFrame({
        'gene': hvg_names,
        'lasso_coef': coef,
        'abs_coef': np.abs(coef)
    })

    coef_df = coef_df.sort_values('abs_coef', ascending=False)

    n_select = min(TOP_K, len(coef_df))
    top_genes = coef_df.head(n_select)['gene'].tolist()

    print(f"    Top {n_select} LASSO基因: {', '.join(top_genes[:10])}...")

    ct_genes = raw_adata[:, top_genes].copy()
    ct_genes = ct_genes[ct_mask, :]

    sc.pp.normalize_total(ct_genes, target_sum=1e4)
    sc.pp.log1p(ct_genes)

    expr_matrix = ct_genes.X
    if hasattr(expr_matrix, 'toarray'):
        expr_matrix = expr_matrix.toarray()

    df_expr = pd.DataFrame(
        expr_matrix,
        index=ct_genes.obs_names,
        columns=top_genes
    )
    df_expr[condition_key] = ct_genes.obs[condition_key].values
    df_expr[cell_type_key] = ct

    results_by_ct[ct] = {
        'n_cells': n_total,
        'n_mcao': n_mcao,
        'n_sham': n_sham,
        'C_optimal': C_optimal,
        'coef_df': coef_df,
        'top_genes': top_genes,
        'n_selected': n_select,
        'expr_df': df_expr
    }

    lasso_genes_all.update(top_genes)

print(f"\n{'='*60}")
print(f"总基因池: {len(lasso_genes_all)} unique genes across all cell types")

all_top_genes = sorted(lasso_genes_all)

summary_rows = []
for ct, res in results_by_ct.items():
    summary_rows.append({
        'cell_type': ct,
        'n_cells': res['n_cells'],
        'n_mcao': res['n_mcao'],
        'n_sham': res['n_sham'],
        'C_optimal': round(res['C_optimal'], 6),
        'n_genes_lasso': res['n_selected']
    })

df_summary = pd.DataFrame(summary_rows)
print("\nLASSO预筛摘要:")
print(df_summary.to_string(index=False))

ct_order = [ct for ct in ['Neuron', 'OPC', 'Oligodendrocyte', 'Astrocyte',
                            'Microglia', 'Pericyte', 'Endothelial', 'Ependymal']
            if ct in results_by_ct]

print(f"\n生成细胞类型特异性表达矩阵...")

xlsx_path = os.path.join(OUTPUT_DIR, "L1_CellType_LASSO_ExpressionMatrix.xlsx")

with pd.ExcelWriter(xlsx_path, engine='openpyxl') as writer:
    df_summary.to_excel(writer, sheet_name='LASSO_Summary', index=False)

    for ct in ct_order:
        res = results_by_ct[ct]

        coef_df_out = res['coef_df'].head(TOP_K).copy()
        coef_df_out['rank'] = range(1, len(coef_df_out) + 1)
        coef_df_out = coef_df_out[['rank', 'gene', 'lasso_coef', 'abs_coef']]

        sheet_name = f"{ct}_LASSO_Top{TOP_K}"[:31]
        coef_df_out.to_excel(writer, sheet_name=sheet_name, index=False)

        expr_agg_rows = []
        for gene in res['top_genes']:
            mcao_expr = res['expr_df'].loc[res['expr_df'][condition_key].str.lower().str.strip() == 'mcao', gene]
            sham_expr = res['expr_df'].loc[res['expr_df'][condition_key].str.lower().str.strip() == 'sham', gene]
            expr_agg_rows.append({
                'gene': gene,
                'MCAO_mean_log1p': round(float(mcao_expr.mean()), 6),
                'MCAO_std': round(float(mcao_expr.std()), 6),
                'MCAO_pct_expressed': round(float((mcao_expr > 0).mean() * 100), 2),
                'Sham_mean_log1p': round(float(sham_expr.mean()), 6),
                'Sham_std': round(float(sham_expr.std()), 6),
                'Sham_pct_expressed': round(float((sham_expr > 0).mean() * 100), 2),
            })

        df_expr_agg = pd.DataFrame(expr_agg_rows)
        df_expr_agg['log2FC'] = df_expr_agg['MCAO_mean_log1p'] - df_expr_agg['Sham_mean_log1p']
        df_expr_agg = df_expr_agg.sort_values('log2FC', key=abs, ascending=False)

        sheet_name_expr = f"{ct}_Expression"[:31]
        df_expr_agg.to_excel(writer, sheet_name=sheet_name_expr, index=False)

    all_expr_matrix = {}
    for ct in ct_order:
        res = results_by_ct[ct]
        for gene in res['top_genes']:
            if gene not in all_expr_matrix:
                all_expr_matrix[gene] = {}
            mcao_expr = res['expr_df'].loc[res['expr_df'][condition_key].str.lower().str.strip() == 'mcao', gene]
            sham_expr = res['expr_df'].loc[res['expr_df'][condition_key].str.lower().str.strip() == 'sham', gene]
            all_expr_matrix[gene][f"{ct}_MCAO_mean"] = round(float(mcao_expr.mean()), 6)
            all_expr_matrix[gene][f"{ct}_Sham_mean"] = round(float(sham_expr.mean()), 6)
            all_expr_matrix[gene][f"{ct}_log2FC"] = round(float(mcao_expr.mean() - sham_expr.mean()), 6)

    df_matrix = pd.DataFrame.from_dict(all_expr_matrix, orient='index')
    df_matrix.index.name = 'gene'
    df_matrix = df_matrix.sort_index()

    df_matrix.to_excel(writer, sheet_name='Full_Expression_Matrix')

print(f"\n✅ 已保存: {xlsx_path}")

for dst_dir in [D_DISK_TABLE, D_DISK_RNASEQ]:
    try:
        dst_path = os.path.join(dst_dir, "L1_CellType_LASSO_ExpressionMatrix.xlsx")
        import shutil
        shutil.copy2(xlsx_path, dst_path)
        print(f"✅ 已同步: {dst_path}")
    except Exception as e:
        print(f"⚠ 同步失败 {dst_dir}: {e}")

xlsx_path2 = os.path.join(OUTPUT_DIR, "L1_CellType_LASSO_ExpressionMatrix_wide.xlsx")

with pd.ExcelWriter(xlsx_path2, engine='openpyxl') as writer:
    df_summary.to_excel(writer, sheet_name='LASSO_Summary', index=False)

    for ct in ct_order:
        res = results_by_ct[ct]
        coef_sorted = res['coef_df'].sort_values('abs_coef', ascending=False).head(TOP_K)
        sheet_name = f"{ct}_Top{TOP_K}"[:31]
        coef_sorted.to_excel(writer, sheet_name=sheet_name, index=False)

    df_matrix.to_excel(writer, sheet_name='Expression_Matrix')

    log2fc_cols = [c for c in df_matrix.columns if '_log2FC' in c]
    if log2fc_cols:
        df_log2fc = df_matrix[log2fc_cols].copy()
        df_log2fc.columns = [c.replace('_log2FC', '') for c in df_log2fc.columns]
        df_log2fc.to_excel(writer, sheet_name='log2FC_Matrix')

    mean_cols_mcao = [c for c in df_matrix.columns if '_MCAO_mean' in c]
    mean_cols_sham = [c for c in df_matrix.columns if '_Sham_mean' in c]
    if mean_cols_mcao:
        df_mcao = df_matrix[mean_cols_mcao].copy()
        df_mcao.columns = [c.replace('_MCAO_mean', '') for c in df_mcao.columns]
        df_mcao.to_excel(writer, sheet_name='MCAO_MeanExpression')

    if mean_cols_sham:
        df_sham = df_matrix[mean_cols_sham].copy()
        df_sham.columns = [c.replace('_Sham_mean', '') for c in df_sham.columns]
        df_sham.to_excel(writer, sheet_name='Sham_MeanExpression')

print(f"✅ 已保存: {xlsx_path2}")

for dst_dir in [D_DISK_TABLE, D_DISK_RNASEQ]:
    try:
        dst_path = os.path.join(dst_dir, "L1_CellType_LASSO_ExpressionMatrix_wide.xlsx")
        import shutil
        shutil.copy2(xlsx_path2, dst_path)
        print(f"✅ 已同步: {dst_path}")
    except Exception as e:
        print(f"⚠ 同步失败 {dst_dir}: {e}")

np.set_printoptions(precision=3, suppress=True)
print(f"\n{'='*60}")
print("全部分析完成！")
print("=" * 60)
for ct in ct_order:
    res = results_by_ct[ct]
    top3 = res['coef_df'].head(3)
    print(f"\n  {ct} (n={res['n_cells']}, Top {res['n_selected']}):")
    for _, row in top3.iterrows():
        print(f"    {row['gene']:15s}  coef={row['lasso_coef']:+.6f}")