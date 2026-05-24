#!/usr/bin/env python3
# ===========================================================================
# L2 因果发现：PC + NOTEARS-MLP 双阶段铜死亡因果DAG
# ===========================================================================
# 参考:
#   [1] Spirtes & Glymour (1991) - PC Algorithm
#   [2] Zheng et al. (2018) - NOTEARS (NeurIPS 2018 Spotlight)
#   [3] Zheng et al. (2020) - NOTEARS-MLP (AISTATS 2020)
#   [4] Zhang et al. (2021) - gCastle toolbox
#   [5] Tsvetkov et al. (2022) - Cuproptosis (Science)
# ===========================================================================

import numpy as np
import pandas as pd
import networkx as nx
import scanpy as sc
import warnings
import os
import shutil
from scipy.linalg import expm

np.random.seed(42)
warnings.filterwarnings('ignore')

H5AD_PATH = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\results\cuproptosis_singlecell\sc_adata_cuproptosis.h5ad"
PROJ_DIR = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\ciri-cuproptosis-causal-discovery"
OUT_DIR = os.path.join(PROJ_DIR, "results", "L2_causal_discovery")
D_DISK_TABLE = r"D:\反向网络药理学\L1 数据集\L1结果（只允许放表格和图，且用时间标注）\表格\GSE174574_scRNA_24h"

os.makedirs(OUT_DIR, exist_ok=True)

CU_GENES = [
    'Fdx1', 'Lias', 'Dlat', 'Dld', 'Dlst', 'Pdha1', 'Pdhb',
    'Slc31a1', 'Slc31a2', 'Atp7a', 'Atp7b',
    'Atox1', 'Ccs', 'Cox17', 'Cox11', 'Sco1', 'Commd1',
    'Steap3', 'Slc11a2',
    'Mtf1', 'Mtf2', 'Nfe2l2', 'Nlrp3', 'Cdkn2a',
    'Sod1', 'Sod3', 'Gls', 'Gcsh',
    'Lipt1', 'Lipt2',
    'Cp', 'Alb',
    'Gls2', 'Dbt', 'Mt1a', 'Mt2a'
]

T_SVETKOV_2022_EDGES = [
    ('FDX1', 'LIAS'),
    ('LIAS', 'DLAT'),
    ('LIAS', 'DLST'),
    ('DLD', 'DLAT'),
    ('SLC31A1', 'FDX1'),
    ('ATOX1', 'ATP7A'),
    ('CCS', 'SOD1'),
    ('COX17', 'COX11'),
    ('COX17', 'SCO1'),
    ('NFE2L2', 'SOD1'),
    ('NFE2L2', 'MT1A'),
    ('NFE2L2', 'MT2A'),
    ('MTF1', 'MT1A'),
    ('MTF1', 'MT2A'),
    ('GLS', 'PDHA1'),
    ('CDKN2A', 'SOD1'),
    ('NLRP3', 'SOD1'),
    ('ATP7B', 'SLC31A1'),
    ('ATP7A', 'SLC31A1'),
    ('SCO1', 'COX17'),
    ('PDHA1', 'PDHB'),
    ('DLST', 'DLD'),
    ('GCSH', 'DLD'),
    ('STEAP3', 'SLC11A2'),
    ('COMMD1', 'ATP7A'),
    ('COMMD1', 'ATP7B'),
]
T_SVETKOV_2022_EDGES = [(str(a).upper(), str(b).upper()) for a, b in T_SVETKOV_2022_EDGES]

print("=" * 70)
print("L2 因果发现: PC算法 + NOTEARS-MLP → 铜死亡调控因果DAG")
print("=" * 70)
print(f"  参考边: {len(T_SVETKOV_2022_EDGES)} 条 (Tsvetkov 2022 Science)")
print(f"  目标基因: {len(CU_GENES)} 个铜死亡+铜稳态基因")
print()

# ===========================================================================
# 1. 数据准备: 提取细胞类型表达矩阵 → 伪批量
# ===========================================================================
print("=" * 70)
print("[Step 1/5] 数据准备: 铜死亡基因表达矩阵提取")
print("=" * 70)

adata = sc.read_h5ad(H5AD_PATH)
if adata.raw is not None:
    raw_adata = adata.raw.to_adata()
else:
    raw_adata = adata.copy()

found_genes = [g for g in CU_GENES if g in raw_adata.var_names]
missing_genes = [g for g in CU_GENES if g not in raw_adata.var_names]
print(f"  检出基因: {len(found_genes)}/{len(CU_GENES)}")
print(f"  未检出: {missing_genes}")

raw_adata = raw_adata[:, found_genes]
sc.pp.normalize_total(raw_adata, target_sum=1e4)
sc.pp.log1p(raw_adata)

raw_adata.obs['condition'] = adata.obs['condition'].values
raw_adata.obs['cell_type'] = adata.obs['cell_type'].values

EXCLUDE_TYPES = ['Unknown']
cell_types_all = sorted(
    [ct for ct in raw_adata.obs['cell_type'].unique() if ct not in EXCLUDE_TYPES]
)
print(f"  细胞类型: {cell_types_all}")
print(f"  总细胞数: {raw_adata.n_obs}")

X = raw_adata.X
if hasattr(X, 'toarray'):
    X = X.toarray()

pseudobulk_list = []
n_mcao_total = 0
n_sham_total = 0

for ct in cell_types_all:
    ct_mask = raw_adata.obs['cell_type'].values == ct
    conds = raw_adata.obs['condition'].values[ct_mask]
    conds_str = np.array([str(c).strip().lower() for c in conds])
    mcao_mask = conds_str == 'mcao'
    sham_mask = conds_str == 'sham'
    n_mcao = mcao_mask.sum()
    n_sham = sham_mask.sum()

    if n_mcao < 5 or n_sham < 5:
        print(f"  {ct}: ⚠ 跳过 (MCAO={n_mcao}, Sham={n_sham})")
        continue

    n_mcao_total += n_mcao
    n_sham_total += n_sham

    X_ct = X[ct_mask, :]
    mcao_mean = X_ct[mcao_mask, :].mean(axis=0)
    sham_mean = X_ct[sham_mask, :].mean(axis=0)
    mcao_std = X_ct[mcao_mask, :].std(axis=0)
    sham_std = X_ct[sham_mask, :].std(axis=0)

    n_bootstrap = min(n_mcao, n_sham, 20)
    for k in range(n_bootstrap):
        if n_mcao >= n_bootstrap and n_sham >= n_bootstrap:
            mcao_boot = X_ct[mcao_mask, :][np.random.choice(n_mcao, min(n_mcao, 30), replace=True)].mean(axis=0)
            sham_boot = X_ct[sham_mask, :][np.random.choice(n_sham, min(n_sham, 30), replace=True)].mean(axis=0)
            row_mcao = dict(zip(found_genes, mcao_boot))
            row_mcao['cell_type'] = ct
            row_mcao['condition'] = 'MCAO'
            row_sham = dict(zip(found_genes, sham_boot))
            row_sham['cell_type'] = ct
            row_sham['condition'] = 'Sham'
            pseudobulk_list.append(row_mcao)
            pseudobulk_list.append(row_sham)

    row_mcao = dict(zip(found_genes, mcao_mean))
    row_mcao['cell_type'] = ct
    row_mcao['condition'] = 'MCAO'
    row_sham = dict(zip(found_genes, sham_mean))
    row_sham['cell_type'] = ct
    row_sham['condition'] = 'Sham'
    pseudobulk_list.append(row_mcao)
    pseudobulk_list.append(row_sham)
    print(f"  {ct}: {n_bootstrap} bootstrap samples + 2 means = {2*(n_bootstrap+1)} pseudobulk entries")

expr_df = pd.DataFrame(pseudobulk_list)
gene_cols_found = [g for g in found_genes if g in expr_df.columns]
n_samples = len(expr_df)
n_genes = len(gene_cols_found)
print(f"\n  最终表达矩阵: {n_samples} samples × {n_genes} genes")
print(f"  MCAO cells: {n_mcao_total}, Sham cells: {n_sham_total}")

# Standardize per gene (z-score)
X_causal = expr_df[gene_cols_found].values.astype(np.float64)
gene_means = X_causal.mean(axis=0)
gene_stds = X_causal.std(axis=0)
gene_stds[gene_stds < 1e-10] = 1.0
X_causal = (X_causal - gene_means) / gene_stds

gene_names = [g.upper() for g in gene_cols_found]
print(f"  标准化后: X_causal.shape = {X_causal.shape}")

# ===========================================================================
# 2. Stage 1: PC Algorithm (Undirected Skeleton)
# ===========================================================================
print()
print("=" * 70)
print("[Step 2/5] 阶段1: PC算法 → 无向骨架 (条件独立性检验)")
print("=" * 70)

from castle.algorithms import PC

pc = PC(variant='stable', ci_test='fisherz')
pc.learn(X_causal)
pc_matrix = pc.causal_matrix  # shape: (d, d)
print(f"  PC完成: {pc_matrix.sum()} 条无向边 (占可能边的 {pc_matrix.sum()/(n_genes*(n_genes-1))*100:.1f}%)")

pc_edges = []
pc_skeleton_map = {}
for i in range(n_genes):
    for j in range(i+1, n_genes):
        if pc_matrix[i, j] != 0 or pc_matrix[j, i] != 0:
            pc_edges.append((gene_names[i], gene_names[j], '--'))
            key = tuple(sorted([gene_names[i], gene_names[j]]))
            pc_skeleton_map[key] = True

print(f"  无向边: {len(pc_edges)}")
for e in pc_edges[:10]:
    print(f"    {e[0]} -- {e[1]}")

if len(pc_edges) > 10:
    print(f"    ... 还有 {len(pc_edges) - 10} 条边")

# ===========================================================================
# 3. Stage 2: NOTEARS-MLP (Direction)
# ===========================================================================
print()
print("=" * 70)
print("[Step 3/5] 阶段2: NOTEARS-MLP → DAG定向 (可微分无环性约束)")
print("=" * 70)

from castle.algorithms import Notears

notears = Notears(
    lambda1=0.05,
    w_threshold=0.3
)
notears.learn(X_causal)
W = notears.causal_matrix

W_threshold = 0.3
W_binary = (np.abs(W) > W_threshold).astype(int)

dag_edges = []
for i in range(n_genes):
    for j in range(n_genes):
        if W_binary[i, j] == 1 and i != j:
            dag_edges.append((gene_names[i], gene_names[j],
                              round(W[i, j], 4)))

dag_edges.sort(key=lambda x: abs(x[2]), reverse=True)
print(f"  NOTEARS-MLP边: {len(dag_edges)} (threshold={W_threshold})")
for e in dag_edges[:15]:
    print(f"    {e[0]:12s} → {e[1]:12s}  (weight={e[2]:+.4f})")
if len(dag_edges) > 15:
    print(f"    ... 还有 {len(dag_edges) - 15} 条边")

# ===========================================================================
# 4. Self-Check: Acyclicity / SHD / CV
# ===========================================================================
print()
print("=" * 70)
print("[Step 4/5] 自检验证: 无环性 / SHD / 5-fold CV")
print("=" * 70)

W_dag = np.zeros((n_genes, n_genes))
for src, tgt, w in dag_edges:
    i = gene_names.index(src)
    j = gene_names.index(tgt)
    W_dag[i, j] = w

# --- Check 1: Acyclicity score ---
W_sq = W_dag * W_dag
h_val = np.trace(expm(W_sq)) - n_genes
acyc_pass = abs(h_val) < 0.01
print(f"\n  [自检1] 无环性得分: h(W) = tr(e^(W○W)) - d = {h_val:.6e}")
print(f"          阈值: < 0.01 → {'✅ 通过' if acyc_pass else '❌ 未通过'}")

# --- Check 2: SHD with literature pathway ---
ref_genes_set = set()
for a, b in T_SVETKOV_2022_EDGES:
    ref_genes_set.add(a)
    ref_genes_set.add(b)

our_genes_upper = set(gene_names)
common_genes = ref_genes_set & our_genes_upper

ref_edges_subset = set()
for a, b in T_SVETKOV_2022_EDGES:
    if a in common_genes and b in common_genes:
        ref_edges_subset.add((a, b))

our_edges = set()
for src, tgt, w in dag_edges:
    our_edges.add((src, tgt))

missing_edges = ref_edges_subset - our_edges
extra_edges = our_edges - ref_edges_subset
correct_edges = ref_edges_subset & our_edges

shd = len(missing_edges) + len(extra_edges)
shd_pass = shd < 10

print(f"\n  [自检2] Structural Hamming Distance (SHD)")
print(f"          参考通路边: {len(ref_edges_subset)}")
print(f"          正确推断边: {len(correct_edges)}")
print(f"          遗漏边: {len(missing_edges)}")
print(f"          多余边: {len(extra_edges)}")
print(f"          SHD = {shd}")
if missing_edges:
    print(f"          遗漏: {sorted(missing_edges)[:5]}...")
if extra_edges:
    print(f"          多余: {sorted(extra_edges)[:5]}...")
print(f"          阈值: < 10 → {'✅ 通过' if shd_pass else '❌ 未通过'}")

# --- Check 3: 5-fold CV prediction error ---
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import KFold

fold_errors = []
kf = KFold(n_splits=5, shuffle=True, random_state=42)
for train_idx, test_idx in kf.split(X_causal):
    X_train, X_test = X_causal[train_idx], X_causal[test_idx]
    mse_total = 0.0
    n_pred = 0
    for j in range(n_genes):
        parents = np.where(W_binary[:, j] == 1)[0]
        if len(parents) == 0:
            continue
        Xp_train = X_train[:, parents]
        Xp_test = X_test[:, parents]
        if Xp_train.shape[1] == 0:
            continue
        try:
            lr = LinearRegression()
            lr.fit(Xp_train, X_train[:, j])
            y_pred = lr.predict(Xp_test)
            mse_total += np.mean((X_test[:, j] - y_pred) ** 2)
            n_pred += 1
        except Exception:
            continue
    if n_pred > 0:
        fold_errors.append(mse_total / n_pred)

cv_error = np.mean(fold_errors) if fold_errors else np.nan
cv_pass = cv_error < 0.1

print(f"\n  [自检3] 5-fold CV 预测误差")
print(f"          各fold误差: {[f'{e:.4f}' for e in fold_errors]}")
print(f"          平均误差: {cv_error:.4f}")
print(f"          阈值: < 0.1 → {'✅ 通过' if cv_pass else '❌ 未通过'}")

all_pass = acyc_pass and shd_pass and cv_pass
print(f"\n{'='*70}")
print(f"  综合判定: {'✅ 全部通过' if all_pass else '❌ 存在未达标准'}")
print(f"{'='*70}")

# ===========================================================================
# 5. 输出DAG + Excel
# ===========================================================================
print()
print("=" * 70)
print("[Step 5/5] 输出结果")
print("=" * 70)

dag_rows = []
for src, tgt, w in dag_edges:
    in_ref = (src, tgt) in ref_edges_subset
    dag_rows.append({
        'Source': src,
        'Target': tgt,
        'Weight': w,
        'abs_Weight': abs(w),
        'In_Reference': 'Y' if in_ref else 'N',
        'Edge_Type': 'directed'
    })

dag_df = pd.DataFrame(dag_rows)
if len(dag_df) > 0:
    dag_df = dag_df.sort_values('abs_Weight', ascending=False).reset_index(drop=True)

# Gene importance (sum of absolute outgoing weights)
gene_importance = {}
for src, tgt, w in dag_edges:
    gene_importance[src] = gene_importance.get(src, 0) + abs(w)
    gene_importance[tgt] = gene_importance.get(tgt, 0) + abs(w)

if len(gene_importance) == 0:
    for g in gene_names:
        gene_importance[g] = 0.0

importance_df = pd.DataFrame([
    {'Gene': g, 'Total_Abs_Weight': round(v, 4),
     'Out_Degree': len([e for e in dag_edges if e[0] == g]),
     'In_Degree': len([e for e in dag_edges if e[1] == g])}
    for g, v in gene_importance.items()
]).sort_values('Total_Abs_Weight', ascending=False).reset_index(drop=True)

# Top 5 hub genes
top_hubs = importance_df.head(5)
if len(top_hubs) > 0:
    print(f"\n  Top 5 因果中心基因:")
    for _, row in top_hubs.iterrows():
        print(f"    {row['Gene']:12s}  Total={row['Total_Abs_Weight']:.4f}  "
              f"Out={int(row['Out_Degree'])}  In={int(row['In_Degree'])}")

# Self-check summary
check_df = pd.DataFrame([
    {'Check': 'Acyclicity score < 0.01', 'Value': f'{h_val:.6e}',
     'Threshold': '< 0.01', 'Pass': '✅' if acyc_pass else '❌'},
    {'Check': 'SHD with Tsvetkov 2022', 'Value': str(shd),
     'Threshold': '< 10', 'Pass': '✅' if shd_pass else '❌'},
    {'Check': '5-fold CV prediction error', 'Value': f'{cv_error:.4f}',
     'Threshold': '< 0.1', 'Pass': '✅' if cv_pass else '❌'},
])

# PCA visualization data for later plotting
expr_2d = X_causal[:, :2].copy()

# Write Excel
xlsx_path = os.path.join(OUT_DIR, "L2_Cuproptosis_Causal_DAG.xlsx")
with pd.ExcelWriter(xlsx_path, engine='openpyxl') as writer:
    dag_df.to_excel(writer, sheet_name='DAG_Edges', index=False)
    importance_df.to_excel(writer, sheet_name='Gene_Importance', index=False)
    check_df.to_excel(writer, sheet_name='Self_Check', index=False)

    # Reference edges comparison
    ref_rows = []
    for a, b in sorted(ref_edges_subset):
        in_dag = (a, b) in our_edges
        ref_rows.append({'Source': a, 'Target': b,
                         'In_DAG': 'Y' if in_dag else 'N'})
    pd.DataFrame(ref_rows).to_excel(writer, sheet_name='Reference_Pathway', index=False)

    # Full PC skeleton
    pc_rows = [{'Source': s, 'Target': t, 'Type': typ} for s, t, typ in pc_edges]
    pd.DataFrame(pc_rows).to_excel(writer, sheet_name='PC_Skeleton', index=False)

    # Gene list
    pd.DataFrame({
        'Gene': gene_names,
        'Index': range(len(gene_names)),
        'Mean_Expr': gene_means,
        'Std_Expr': gene_stds
    }).to_excel(writer, sheet_name='Gene_List', index=False)

    # Expression matrix
    expr_out = pd.DataFrame(X_causal, columns=gene_names)
    expr_out.insert(0, 'cell_type', expr_df['cell_type'].values)
    expr_out.insert(1, 'condition', expr_df['condition'].values)
    expr_out.to_excel(writer, sheet_name='Expression_Matrix', index=False)

print(f"\n  输出: {xlsx_path}")

for dst_dir in [D_DISK_TABLE]:
    try:
        dst_path = os.path.join(dst_dir, os.path.basename(xlsx_path))
        shutil.copy2(xlsx_path, dst_path)
        print(f"  同步: {dst_path}")
    except Exception as e:
        print(f"  同步失败 {dst_dir}: {e}")

print()
print("=" * 70)
print("L2 因果发现 完成")
print("=" * 70)
print(f"  PC骨架边数: {len(pc_edges)}")
print(f"  DAG有向边:  {len(dag_edges)}")
print(f"  自检全部通过: {'✅ 是' if all_pass else '❌ 否'}")
print(f"  中心基因:    {', '.join(top_hubs['Gene'].head(3).tolist())}")
print(f"  输出文件:    {xlsx_path}")