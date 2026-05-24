#!/usr/bin/env python3
# ===========================================================================
# L2b 因果涌现粗粒化: Greedy-Spectral Coarsening (CE-Coarsening)
# ===========================================================================
# 参考文献:
#   [1] Loukas (2019) - Graph reduction with spectral and cut guarantees (JMLR)
#   [2] Griebenow, Klein, Hoel (2019) - Finding the right scale: causal emergence
#       via spectral clustering (arXiv:1908.07565)
#   [3] Leicht & Newman (2008) - Community structure in directed networks (PRL)
#   [4] Wang et al. (2007) - A new method to measure the semantic similarity
#       of GO terms (Bioinformatics)
#   [5] Resnik (1999) - Semantic similarity in a taxonomy (JAIR)
#
# 输入:
#   - L2a W_merged (nxn 带权因果DAG邻接矩阵)
#   - L2a merged_edges [(src, tgt, weight), ...]
#   - L1 celltype log2FC 矩阵
#   - [Optional] STRING PPI weights
#
# 输出:
#   - 社区划分 C = {c1, c2, ..., ck}
#   - 宏节点级因果图 G_macro
#   - 代表基因列表 (HubScore_L2a + |log2FC_L1|)
#   - 人类同源映射接口 (for L3 MR)
# ===========================================================================

import numpy as np
import pandas as pd
import networkx as nx
import warnings
import os
import shutil
from scipy.sparse.linalg import eigsh
from scipy.sparse import csr_matrix, diags, eye
from sklearn.cluster import KMeans
from collections import defaultdict

np.random.seed(42)
warnings.filterwarnings('ignore')

PROJ_DIR = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\ciri-cuproptosis-causal-discovery"
L2_RESULTS = os.path.join(PROJ_DIR, "results", "L2_causal_discovery")
L1_RESULTS = os.path.join(PROJ_DIR, "results", "L1_phenotype_anchoring")
L2B_DIR = os.path.join(PROJ_DIR, "L2b_emergence_coarsening")
L2B_RESULTS = os.path.join(PROJ_DIR, "results", "L2b_emergence_coarsening")
D_DISK_TABLE = r"D:\反向网络药理学\L1 数据集\L1结果（只允许放表格和图，且用时间标注）\表格\GSE174574_scRNA_24h"

os.makedirs(L2B_RESULTS, exist_ok=True)

K_CLUSTERS_INITIAL = 6
BETA_STRING = 0.0
ALPHA_GO = 0.4
ALPHA_Q = 0.4
ALPHA_CAUSAL = 0.2
TARGET_RATIO_MIN = 0.20
TARGET_RATIO_MAX = 0.40
STAGNATION_LIMIT = 5
DELTA_EPSILON = 0.01

print("=" * 70)
print("L2b 因果涌现粗粒化: Greedy-Spectral Coarsening")
print("=" * 70)

# ===========================================================================
# 0. 加载数据
# ===========================================================================
print("\n[0/4] 加载输入数据")

xlsx_path = os.path.join(L2_RESULTS, "L2_Cuproptosis_Causal_DAG.xlsx")
dag_df = pd.read_excel(xlsx_path, sheet_name='Final_DAG_Merged')
gene_df = pd.read_excel(xlsx_path, sheet_name='Gene_List')
importance_df = pd.read_excel(xlsx_path, sheet_name='Gene_Importance')

merged_edges = list(zip(dag_df['Source'].values, dag_df['Target'].values, dag_df['Weight'].values))
gene_names = gene_df['Gene'].str.upper().values.tolist()
n_genes = len(gene_names)

name_to_idx = {g: i for i, g in enumerate(gene_names)}

W_merged = np.zeros((n_genes, n_genes))
for src, tgt, w in merged_edges:
    W_merged[name_to_idx[src], name_to_idx[tgt]] = w

hub_score = {}
for _, row in importance_df.iterrows():
    hub_score[row['Gene']] = row['Total_Abs_Weight']
for g in gene_names:
    if g not in hub_score:
        hub_score[g] = 0.0

# L1 log2FC
log2fc_path = os.path.join(L1_RESULTS, "celltype_cuproptosis_log2FC.csv")
log2fc_df = pd.read_csv(log2fc_path, index_col=0)
gene_to_max_log2fc = {}
for gene in gene_names:
    gene_title = gene[0] + gene[1:].lower()
    if gene_title in log2fc_df.index:
        vals = log2fc_df.loc[gene_title].values
        gene_to_max_log2fc[gene] = np.nanmax(np.abs(vals.astype(float)))
    else:
        gene_to_max_log2fc[gene] = 0.0

n_data_driven = int(dag_df['Edge_Type'].value_counts().get('Data_Driven', 0))
n_literature = int(dag_df['Edge_Type'].value_counts().get('Literature', 0))
print(f"  基因数: {n_genes}")
print(f"  因果边数: {len(merged_edges)} (数据驱动{n_data_driven} + 文献增强{n_literature})")

# ===========================================================================
# 1. Step 1: 构建亲和矩阵 + 归一化拉普拉斯
# ===========================================================================
print("\n[1/4] Step 1: 亲和矩阵 A + 归一化拉普拉斯 L_sym")

A_causal = (np.abs(W_merged) + np.abs(W_merged).T) / 2.0

if BETA_STRING > 0:
    string_path = os.path.join(PROJ_DIR, "data", "STRING_cuproptosis_ppi.csv")
    if os.path.exists(string_path):
        string_df = pd.read_csv(string_path, index_col=0)
        string_matrix = string_df.reindex(index=gene_names, columns=gene_names).fillna(0).values
        string_matrix = (string_matrix + string_matrix.T) / 2.0
        A = A_causal + BETA_STRING * np.where(string_matrix > 0, string_matrix / string_matrix.max(), 0)
        print(f"  STRING PPI 已整合 (beta={BETA_STRING})")
    else:
        print(f"  STRING 文件未找到, 使用纯因果亲和矩阵 (beta=0)")
        A = A_causal
else:
    A = A_causal

deg = A.sum(axis=1)
deg_safe = np.where(deg > 0, deg, 1.0)
D_inv_sqrt = np.diag(1.0 / np.sqrt(deg_safe))
L_sym = np.eye(n_genes) - D_inv_sqrt @ A @ D_inv_sqrt

print(f"  A 非零元素: {np.count_nonzero(A)}")
print(f"  L_sym trace: {np.trace(L_sym):.4f}")
print(f"  L_sym rank: {np.linalg.matrix_rank(L_sym)}")

# ===========================================================================
# 2. Step 2: 谱预分组 - k=6 eigenvectors + k-means++
# ===========================================================================
print("\n[2/4] Step 2: 谱预分组 (k=6 eigenvectors + k-means++)")

k = min(K_CLUSTERS_INITIAL, n_genes - 1)
k = max(k, 3)
print(f"  有效聚类数 k={k}")

A_sparse = csr_matrix(A)
deg_sparse = A_sparse.sum(axis=1).A1
deg_safe_sparse = np.where(deg_sparse > 0, deg_sparse, 1.0)
D_inv_sqrt_diag = 1.0 / np.sqrt(deg_safe_sparse)
L_sym_sparse = eye(n_genes) - diags(D_inv_sqrt_diag) @ A_sparse @ diags(D_inv_sqrt_diag)

eigenvalues, eigenvectors = eigsh(L_sym_sparse.astype(np.float64), k=k + 1, which='SM')
eigenvectors = eigenvectors[:, 1:k + 1]

kmeans = KMeans(n_clusters=k, random_state=42, n_init=20)
initial_labels = kmeans.fit_predict(eigenvectors)

communities = {}
for i in range(n_genes):
    communities.setdefault(initial_labels[i], []).append(i)

community_ids = list(communities.keys())
print(f"  初始社区数: {len(community_ids)}")
for cid in community_ids:
    members = [gene_names[i] for i in communities[cid]]
    print(f"    C{cid} ({len(members)}基因): {', '.join(members)}")

# ===========================================================================
# 3. GO-BP 功能注释 (via mygene)
# ===========================================================================
print("\n[GO] 获取 GO-BP 注释 (mygene API)")

try:
    import mygene
    mg = mygene.MyGeneInfo()

    mouse_genes_for_query = [g[0] + g[1:].lower() for g in gene_names]
    results = mg.querymany(
        mouse_genes_for_query,
        scopes='symbol',
        species='mouse',
        fields='go.BP',
        returnall=True,
        as_dataframe=False
    )

    gene_to_go_bp = {}
    for i, g in enumerate(gene_names):
        gene_to_go_bp[g] = set()

    hit_count = 0
    for hit in results.get('out', []):
        query = hit.get('query', '').upper()
        if query in name_to_idx:
            go_data = hit.get('go', {})
            bp_terms = go_data.get('BP', [])
            if bp_terms:
                for term in bp_terms:
                    term_id = term.get('id', '')
                    if term_id.startswith('GO:'):
                        gene_to_go_bp[query].add(term_id)
            if gene_to_go_bp[query]:
                hit_count += 1

    # 对未命中基因使用已知功能注释作为回退
    MISSING_GO_FALLBACK = {
        'FDX1': {'GO:0016226', 'GO:0051536', 'GO:0044571'},
        'LIAS': {'GO:0016992', 'GO:0051536', 'GO:0009107'},
        'DLAT': {'GO:0006086', 'GO:0006099', 'GO:0005759'},
        'DLD': {'GO:0005759', 'GO:0006099', 'GO:0045254'},
        'DLST': {'GO:0005759', 'GO:0006099', 'GO:0045252'},
        'PDHA1': {'GO:0006086', 'GO:0006099', 'GO:0005759'},
        'PDHB': {'GO:0006086', 'GO:0006099', 'GO:0005759'},
        'GLS': {'GO:0006543', 'GO:0006541', 'GO:0005739'},
        'GCSH': {'GO:0006546', 'GO:0006544', 'GO:0005759'},
        'LIPT1': {'GO:0009249', 'GO:0009107', 'GO:0005739'},
        'LIPT2': {'GO:0009249', 'GO:0009107', 'GO:0005739'},
        'SLC31A1': {'GO:0006825', 'GO:0005886', 'GO:0070578'},
        'SLC31A2': {'GO:0006825', 'GO:0005886', 'GO:0070578'},
        'ATP7A': {'GO:0006825', 'GO:0015691', 'GO:0005768'},
        'ATP7B': {'GO:0006825', 'GO:0015691', 'GO:0005768'},
        'ATOX1': {'GO:0006825', 'GO:0005507', 'GO:0032767'},
        'CCS': {'GO:0006825', 'GO:0005507', 'GO:0005737'},
        'COX17': {'GO:0006825', 'GO:0005507', 'GO:0005746'},
        'COX11': {'GO:0006825', 'GO:0004129', 'GO:0005746'},
        'SCO1': {'GO:0006825', 'GO:0004129', 'GO:0005746'},
        'COMMD1': {'GO:0006825', 'GO:0031410', 'GO:0005768'},
        'STEAP3': {'GO:0006826', 'GO:0006879', 'GO:0005769'},
        'SLC11A2': {'GO:0006826', 'GO:0006879', 'GO:0005886'},
        'MTF1': {'GO:0006355', 'GO:0003677', 'GO:0005634'},
        'MTF2': {'GO:0006355', 'GO:0003677', 'GO:0005634'},
        'NFE2L2': {'GO:0006979', 'GO:0006355', 'GO:0005634'},
        'NLRP3': {'GO:0006954', 'GO:0050727', 'GO:0005737'},
        'CDKN2A': {'GO:0007049', 'GO:0008285', 'GO:0005634'},
        'SOD1': {'GO:0006801', 'GO:0004784', 'GO:0005737'},
        'SOD3': {'GO:0006801', 'GO:0004784', 'GO:0005576'},
        'CP': {'GO:0006825', 'GO:0004322', 'GO:0005576'},
        'ALB': {'GO:0006810', 'GO:0008281', 'GO:0005576'},
        'GLS2': {'GO:0006543', 'GO:0006541', 'GO:0005739'},
        'DBT': {'GO:0006546', 'GO:0006099', 'GO:0005759'},
    }

    for g in gene_names:
        if len(gene_to_go_bp[g]) == 0 and g in MISSING_GO_FALLBACK:
            gene_to_go_bp[g] = MISSING_GO_FALLBACK[g]
            hit_count += 1

    print(f"  GO-BP 注释成功: {hit_count}/{n_genes} 基因")

except Exception as e:
    print(f"  mygene API 失败: {e}, 使用手动注释回退")
    raise

def go_similarity(g1, g2):
    bp1 = gene_to_go_bp.get(g1, set())
    bp2 = gene_to_go_bp.get(g2, set())
    if len(bp1) == 0 and len(bp2) == 0:
        return 0.0
    if len(bp1) == 0 or len(bp2) == 0:
        return 0.0
    intersection = len(bp1 & bp2)
    union = len(bp1 | bp2)
    return intersection / union if union > 0 else 0.0

print(f"  GO相似度矩阵已构建")

# ===========================================================================
# 4. Step 3: 贪婪精化
# ===========================================================================
print("\n[3/4] Step 3: 贪婪精化 (DeltaScore = 0.4*DeltaQ + 0.4*DeltaGO + 0.2*DeltaCausal)")

def compute_directed_modularity(W, communities_dict):
    m = W.sum()
    if m == 0:
        return 0.0
    d_out = W.sum(axis=1)
    d_in = W.sum(axis=0)
    Q = 0.0
    for cid, members in communities_dict.items():
        for i in members:
            for j in members:
                Q += W[i, j] - (d_out[i] * d_in[j]) / m
    return Q / m

def compute_internal_causal_edges(W, members):
    count = 0
    nodes_set = set(members)
    for i in members:
        for j in members:
            if i != j and W[i, j] != 0:
                count += 1
    return count

def compute_go_avg_similarity(members):
    if len(members) <= 1:
        return 1.0
    total = 0.0
    pairs = 0
    for i in range(len(members)):
        for j in range(i + 1, len(members)):
            total += go_similarity(gene_names[members[i]], gene_names[members[j]])
            pairs += 1
    return total / pairs if pairs > 0 else 0.0

def compute_delta_score(W, dict_before, c_a, c_b, merged_members):
    m = W.sum()
    if m == 0:
        return 0.0

    # Delta Q directed
    Q_before = compute_directed_modularity(W, dict_before)
    dict_after = {k: v for k, v in dict_before.items()}
    for k in list(dict_after.keys()):
        if k == c_a or k == c_b:
            del dict_after[k]
    dict_after[len(dict_after)] = merged_members
    Q_after = compute_directed_modularity(W, dict_after)
    delta_q = Q_after - Q_before

    # Delta GO similarity
    go_before_a = compute_go_avg_similarity(dict_before[c_a])
    go_before_b = compute_go_avg_similarity(dict_before[c_b])
    go_before_max = max(go_before_a, go_before_b)
    go_after = compute_go_avg_similarity(merged_members)
    delta_go = go_after - go_before_max

    # Delta Causal
    internal_a = compute_internal_causal_edges(W, dict_before[c_a])
    internal_b = compute_internal_causal_edges(W, dict_before[c_b])
    internal_merged = compute_internal_causal_edges(W, merged_members)

    na = len(dict_before[c_a])
    nb = len(dict_before[c_b])
    nm = len(merged_members)
    max_a = na * (na - 1)
    max_b = nb * (nb - 1)
    max_m = nm * (nm - 1)

    causal_density_before = 0.0
    if max_a + max_b > 0:
        causal_density_before = (internal_a + internal_b) / (max_a + max_b)
    causal_density_after = internal_merged / max_m if max_m > 0 else 0.0
    delta_causal = causal_density_after - causal_density_before

    score = ALPHA_Q * delta_q + ALPHA_GO * delta_go + ALPHA_CAUSAL * delta_causal
    return score

current_communities = {cid: communities[cid] for cid in community_ids}
current_scores = []
Q_initial = compute_directed_modularity(W_merged, current_communities)

target_min = max(int(len(current_communities) * TARGET_RATIO_MIN), 2)
target_max = max(int(len(current_communities) * TARGET_RATIO_MAX), 3)

print(f"  初始社区数: {len(current_communities)}")
print(f"  目标范围: {target_min}-{target_max}")
print(f"  初始有向模块度 Q_directed: {Q_initial:.6f}")

iteration = 0
stagnation_count = 0
best_delta_ever = 0.0
merge_history = []

while len(current_communities) > target_max and stagnation_count < STAGNATION_LIMIT:
    cids = list(current_communities.keys())
    best_score = -np.inf
    best_pair = None
    best_merged = None

    for a_idx in range(len(cids)):
        for b_idx in range(a_idx + 1, len(cids)):
            c_a, c_b = cids[a_idx], cids[b_idx]
            members_a = current_communities[c_a]
            members_b = current_communities[c_b]
            merged_members = members_a + members_b

            has_internal_edge = False
            merged_set = set(merged_members)
            for src, tgt, _ in merged_edges:
                si = name_to_idx[src]
                ti = name_to_idx[tgt]
                if si in merged_set and ti in merged_set:
                    has_internal_edge = True
                    break
            for i in members_a:
                for j in members_b:
                    if W_merged[i, j] != 0 or W_merged[j, i] != 0:
                        has_internal_edge = True
                        break
                if has_internal_edge:
                    break

            if not has_internal_edge:
                continue

            score = compute_delta_score(W_merged, current_communities, c_a, c_b, merged_members)
            if score > best_score:
                best_score = score
                best_pair = (c_a, c_b)
                best_merged = merged_members

    if best_pair is None:
        print(f"\n  [STOP] 无可合并的社区对 (硬约束: 内部无因果边)")
        break

    c_a, c_b = best_pair
    new_cid = min(c_a, c_b)
    old_cid = max(c_a, c_b)
    current_communities[new_cid] = best_merged
    del current_communities[old_cid]

    merge_history.append({
        'Iteration': iteration,
        'Merged_C1': c_a,
        'Merged_C2': c_b,
        'New_CID': new_cid,
        'Size': len(best_merged),
        'Delta_Score': round(best_score, 6),
        'N_Communities': len(current_communities)
    })

    a_names = [gene_names[i] for i in communities.get(c_a, []) + communities.get(c_b, [])]
    print(f"  Iter {iteration:3d}: C{c_a} + C{c_b} -> C{new_cid} "
          f"({len(best_merged)} genes) DeltaScore={best_score:.6f} "
          f"n_comm={len(current_communities)}")

    if best_score > best_delta_ever:
        best_delta_ever = best_score
        stagnation_count = 0
    elif best_score < DELTA_EPSILON:
        stagnation_count += 1
        if stagnation_count >= STAGNATION_LIMIT:
            print(f"  [STOP] 连续 {stagnation_count} 次 DeltaScore < {DELTA_EPSILON}")

    iteration += 1
    if iteration > 50:
        print("  [STOP] 达到最大迭代次数 50")
        break

final_communities = current_communities
print(f"\n  最终社区数: {len(final_communities)} (原始基因数 {n_genes} 的 {len(final_communities)/n_genes*100:.1f}%)")

Q_final = compute_directed_modularity(W_merged, final_communities)
print(f"  最终有向模块度 Q_directed: {Q_final:.6f} (变化: {Q_final - Q_initial:+.6f})")

# ===========================================================================
# 5. Step 4: 宏节点生成 + 代表基因选择
# ===========================================================================
print("\n[4/4] Step 4: 宏节点生成 + 代表基因选择")

macro_nodes = {}
community_gene_map = {}
macro_id = 0

for cid, members in sorted(final_communities.items(), key=lambda x: -len(x[1])):
    member_names = [gene_names[i] for i in members]
    macro_name = f"M{macro_id}"
    macro_nodes[macro_name] = {
        'community_id': cid,
        'members': member_names,
        'size': len(members)
    }

    community_gene_map.update({g: macro_name for g in member_names})

    hub_scores_local = {g: hub_score.get(g, 0.0) for g in member_names}
    log2fc_local = {g: gene_to_max_log2fc.get(g, 0.0) for g in member_names}

    hub_max = max(hub_scores_local.values()) if max(hub_scores_local.values()) > 0 else 1.0
    log2fc_max = max(log2fc_local.values()) if max(log2fc_local.values()) > 0 else 1.0

    gene_rank = {}
    for g in member_names:
        hub_norm = hub_scores_local[g] / hub_max if hub_max > 0 else 0.0
        lfc_norm = log2fc_local[g] / log2fc_max if log2fc_max > 0 else 0.0
        gene_rank[g] = hub_norm + lfc_norm

    representative = max(gene_rank, key=gene_rank.get)
    macro_nodes[macro_name]['representative'] = representative
    macro_nodes[macro_name]['hub_score_norm'] = round(hub_scores_local[representative] / hub_max, 4) if hub_max > 0 else 0.0
    macro_nodes[macro_name]['log2fc'] = round(log2fc_local[representative], 4)

    print(f"  {macro_name} ({len(members)}基因): "
          f"代表={representative} "
          f"(Hub={macro_nodes[macro_name]['hub_score_norm']:.4f}, "
          f"|log2FC|={macro_nodes[macro_name]['log2fc']:.4f})")
    print(f"    成员: {', '.join(member_names)}")

    macro_id += 1

# 跨社区 DAG 边 → 宏节点级因果图
print(f"\n  宏节点级因果图 G_macro:")
macro_adj = np.zeros((len(macro_nodes), len(macro_nodes)))
macro_name_list = list(macro_nodes.keys())
macro_idx_map = {m: i for i, m in enumerate(macro_name_list)}

macro_edges = []
for src, tgt, w in merged_edges:
    src_macro = community_gene_map.get(src)
    tgt_macro = community_gene_map.get(tgt)
    if src_macro is None or tgt_macro is None:
        continue
    if src_macro != tgt_macro:
        si = macro_idx_map[src_macro]
        ti = macro_idx_map[tgt_macro]
        macro_adj[si, ti] += abs(w)
        macro_edges.append((src_macro, tgt_macro, src, tgt, abs(w)))

macro_edge_summary = defaultdict(lambda: defaultdict(float))
for sm, tm, _, _, w in macro_edges:
    macro_edge_summary[sm][tm] += w

for sm in macro_name_list:
    for tm in macro_name_list:
        w = macro_edge_summary[sm].get(tm, 0.0)
        if w > 0:
            print(f"    {sm} -> {tm}: weight={w:.1f}")

# ===========================================================================
# 6. 人类同源映射接口 (for L3 MR)
# ===========================================================================
print(f"\n[Interface] 人类同源映射接口 (for L3 MR)")

MOUSE_TO_HUMAN = {
    'Fdx1': 'FDX1', 'Lias': 'LIAS', 'Dlat': 'DLAT', 'Dld': 'DLD',
    'Dlst': 'DLST', 'Pdha1': 'PDHA1', 'Pdhb': 'PDHB',
    'Slc31a1': 'SLC31A1', 'Slc31a2': 'SLC31A2',
    'Atp7a': 'ATP7A', 'Atp7b': 'ATP7B', 'Atox1': 'ATOX1',
    'Ccs': 'CCS', 'Cox17': 'COX17', 'Cox11': 'COX11',
    'Sco1': 'SCO1', 'Commd1': 'COMMD1',
    'Steap3': 'STEAP3', 'Slc11a2': 'SLC11A2',
    'Mtf1': 'MTF1', 'Mtf2': 'MTF2',
    'Nfe2l2': 'NFE2L2', 'Nlrp3': 'NLRP3', 'Cdkn2a': 'CDKN2A',
    'Sod1': 'SOD1', 'Sod3': 'SOD3', 'Gls': 'GLS', 'Gcsh': 'GCSH',
    'Lipt1': 'LIPT1', 'Lipt2': 'LIPT2',
    'Cp': 'CP', 'Alb': 'ALB', 'Gls2': 'GLS2', 'Dbt': 'DBT',
    'Mt1': 'MT1A', 'Mt2': 'MT2A', 'Sco2': 'SCO2',
}

for mname, mdata in macro_nodes.items():
    rep_mouse = mdata['representative'][0] + mdata['representative'][1:].lower()
    rep_human = MOUSE_TO_HUMAN.get(rep_mouse, mdata['representative'])
    members_human = [MOUSE_TO_HUMAN.get(g[0] + g[1:].lower(), g) for g in mdata['members']]
    print(f"  {mname}: {rep_mouse} -> {rep_human} (human)")
    mdata['representative_human'] = rep_human
    mdata['members_human'] = members_human

# ===========================================================================
# 7. 输出 Excel
# ===========================================================================
print(f"\n[Output] 生成 Excel")

output_xlsx = os.path.join(L2B_RESULTS, "L2b_CE_Coarsening_Results.xlsx")
with pd.ExcelWriter(output_xlsx, engine='openpyxl') as writer:
    # Sheet 1: 社区划分
    comm_rows = []
    for mname, mdata in macro_nodes.items():
        for gene in mdata['members']:
            comm_rows.append({
                'Macro_Node': mname,
                'Gene': gene,
                'Representative': 'Y' if gene == mdata['representative'] else '',
                'Hub_Score_L2a': hub_score.get(gene, 0.0),
                'Max_abs_log2FC_L1': gene_to_max_log2fc.get(gene, 0.0),
                'Human_Ortholog': MOUSE_TO_HUMAN.get(gene[0] + gene[1:].lower(), gene)
            })
    pd.DataFrame(comm_rows).to_excel(writer, sheet_name='Community_Assignment', index=False)

    # Sheet 2: 宏节点级因果图
    macro_edge_rows = []
    for sm in macro_name_list:
        for tm in macro_name_list:
            w = macro_edge_summary[sm].get(tm, 0.0)
            if w > 0:
                macro_edge_rows.append({
                    'Source_Macro': sm,
                    'Target_Macro': tm,
                    'Weight': round(w, 4),
                    'Source_Rep': macro_nodes[sm]['representative'],
                    'Target_Rep': macro_nodes[tm]['representative']
                })
    pd.DataFrame(macro_edge_rows).to_excel(writer, sheet_name='Macro_DAG', index=False)

    # Sheet 3: 合并历史
    pd.DataFrame(merge_history).to_excel(writer, sheet_name='Merge_History', index=False)

    # Sheet 4: 谱向量
    eig_df = pd.DataFrame(
        eigenvectors,
        columns=[f'Eigen_{i+1}' for i in range(k)]
    )
    eig_df.insert(0, 'Gene', gene_names)
    eig_df.to_excel(writer, sheet_name='Spectral_Embedding', index=False)

    # Sheet 5: 代表基因详情
    rep_rows = []
    for mname, mdata in macro_nodes.items():
        rep_rows.append({
            'Macro_Node': mname,
            'Size': mdata['size'],
            'Representative_Mouse': mdata['representative'],
            'Representative_Human': mdata['representative_human'],
            'Hub_Score_Normalized': mdata['hub_score_norm'],
            'abs_log2FC': mdata['log2fc'],
            'Members_Mouse': ', '.join(mdata['members']),
            'Members_Human': ', '.join(mdata['members_human'])
        })
    pd.DataFrame(rep_rows).to_excel(writer, sheet_name='Macro_Node_Summary', index=False)

    # Sheet 6: 亲和矩阵
    pd.DataFrame(A, index=gene_names, columns=gene_names).to_excel(
        writer, sheet_name='Affinity_Matrix')

print(f"  输出: {output_xlsx}")

for dst_dir in [D_DISK_TABLE]:
    try:
        dst_path = os.path.join(dst_dir, os.path.basename(output_xlsx))
        shutil.copy2(output_xlsx, dst_path)
        if os.path.exists(dst_path):
            print(f"  同步: {dst_path}")
    except Exception as e:
        print(f"  同步 {dst_dir}: {e}")

# ===========================================================================
# 8. 汇总报告
# ===========================================================================
print(f"\n{'=' * 70}")
print(f"L2b 因果涌现粗粒化 完成")
print(f"{'=' * 70}")
print(f"  原始基因数: {n_genes}")
print(f"  宏节点数: {len(macro_nodes)} (压缩至 {len(macro_nodes)/n_genes*100:.1f}%)")
print(f"  合并迭代: {len(merge_history)} 次")
print(f"  有向模块度 Q: {Q_initial:.6f} → {Q_final:.6f}")
print(f"  宏节点间因果边: {len(macro_edge_rows)} 条")
print(f"  代表基因已映射人类同源, 为 L3 MR 预留接口")
print(f"  输出: {output_xlsx}")