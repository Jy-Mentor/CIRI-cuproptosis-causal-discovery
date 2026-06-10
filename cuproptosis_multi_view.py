"""
铜死亡基因 × 共识桥接基因：分模块检验 + 多视角分析 (v2 — 高效版)
======================================================
分析1: 铜死亡分功能模块超几何检验
分析2: 全量数据中铜死亡基因的预测分数分布
分析3: PPI网络邻近度 — 多源BFS优化
"""
import os
import sys
import time
import warnings
import numpy as np
import pandas as pd
import networkx as nx
from collections import deque
from scipy.stats import hypergeom, mannwhitneyu
warnings.filterwarnings("ignore")

# ============ 路径配置 ============
GAT_DATA_DIR = r"C:\Users\Jy-Mentor-7\Desktop\GAT"
BASE_DIR = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙"
RF_DIR = r"C:\Users\Jy-Mentor-7\Desktop\随机森林"

PPI_FILE = os.path.join(GAT_DATA_DIR, "ppi_subgraph.csv")
LGBM_FILE = os.path.join(RF_DIR, "bridge_predictions_v8_final.csv")
GAT_FILE = os.path.join(BASE_DIR, "all_bridge_genes.csv")
CONSENSUS_FILE = os.path.join(BASE_DIR, "consensus_unknown_genes.csv")

OUT_ANALYSIS1 = os.path.join(BASE_DIR, "cuproptosis_module_enrichment.csv")
OUT_ANALYSIS2 = os.path.join(BASE_DIR, "cuproptosis_full_scores.csv")
OUT_ANALYSIS3 = os.path.join(BASE_DIR, "cuproptosis_ppi_proximity.csv")

# ============ 铜死亡基因列表及功能模块 ============
CUPROPTOSIS_MODULES = {
    "Core_Ketoacid_Dehydrogenase": [
        "FDX1", "LIAS", "LIPT1", "LIPT2",
        "DLD", "DLAT", "DLST", "DBT", "GCSH",
        "PDHA1", "PDHB",
    ],
    "Copper_Transport_Homeostasis": [
        "SLC31A1", "ATP7A", "ATP7B", "SLC11A2",
        "ATOX1", "CCS", "COX17", "SCO1", "SCO2",
        "COX11", "MT1A", "MT2A", "MTF1",
    ],
    "Antioxidant_ROS_Defense": [
        "SOD1", "SOD3", "CP", "TYR",
        "LOX", "LOXL1", "LOXL2", "LOXL3", "LOXL4",
        "ALB", "HAMP",
    ],
    "Regulatory_Signaling": [
        "GLS", "CDKN2A", "NLRP3", "NFE2L2",
    ],
}
ALL_CUPROPTOSIS = [g for genes in CUPROPTOSIS_MODULES.values() for g in genes]

# ============ 多源BFS（邻接表，比 NetworkX 快10倍+）============
def build_adjacency_from_csv(ppi_path):
    """从CSV文件构建邻接表 (无向图)"""
    t0 = time.time()
    df = pd.read_csv(ppi_path, sep="\t")
    col_a = df.columns[0]
    col_b = df.columns[1]
    adj = {}
    for _, row in df.iterrows():
        a = str(row[col_a]).upper()
        b = str(row[col_b]).upper()
        adj.setdefault(a, set()).add(b)
        adj.setdefault(b, set()).add(a)
    print(f"  PPI图: {len(adj):,} 节点, {len(df):,} 边 ({time.time()-t0:.1f}s)")
    return adj


def multi_source_bfs(adj, sources):
    """多源BFS — 同时从所有sources出发，一次遍历得到所有节点到最近source的距离"""
    dist = {}
    queue = deque()
    for s in sources:
        if s in adj:
            dist[s] = 0
            queue.append(s)
    while queue:
        node = queue.popleft()
        d = dist[node] + 1
        for nb in adj.get(node, set()):
            if nb not in dist:
                dist[nb] = d
                queue.append(nb)
    return dist


def compute_proximity_fast(adj, cupro_genes, bridge_genes):
    """高效邻近度计算: 多源BFS from bridge genes"""
    t0 = time.time()
    bridge_in_graph = [b for b in bridge_genes if b in adj]
    cupr_in_graph = [c for c in cupro_genes if c in adj]
    print(f"  桥基因在图中: {len(bridge_in_graph)}/{len(bridge_genes)}")
    print(f"  铜死亡基因在图中: {len(cupr_in_graph)}/{len(cupro_genes)}")

    # 多源BFS: 从所有桥基因同时出发
    bridge_dist = multi_source_bfs(adj, bridge_in_graph)
    print(f"  多源BFS完成: {len(bridge_dist)} 可达节点 ({time.time()-t0:.1f}s)")

    # 查询铜死亡基因距离
    results = []
    dists = []
    for gene in cupro_genes:
        if gene in bridge_dist:
            results.append({
                "gene": gene, "in_graph": True,
                "dist_to_nearest_bridge": bridge_dist[gene],
            })
            dists.append(bridge_dist[gene])
        elif gene in adj:
            results.append({
                "gene": gene, "in_graph": True,
                "dist_to_nearest_bridge": None,
            })
        else:
            results.append({
                "gene": gene, "in_graph": False,
                "dist_to_nearest_bridge": None,
            })

    obs_mean = np.mean(dists) if dists else float("inf")

    # 置换检验 (100次)
    np.random.seed(42)
    all_nodes = list(adj.keys())
    n_cupr = len(cupr_in_graph)
    rand_means = []
    for _ in range(100):
        rand_genes = np.random.choice(all_nodes, size=n_cupr, replace=False)
        rand_dists = [bridge_dist.get(rg, None) for rg in rand_genes]
        rand_dists = [d for d in rand_dists if d is not None]
        if rand_dists:
            rand_means.append(np.mean(rand_dists))

    rand_mean = np.mean(rand_means)
    rand_std = np.std(rand_means)
    z_score = (obs_mean - rand_mean) / rand_std if rand_std > 0 else 0
    p_proximity = np.mean([1 if rm <= obs_mean else 0 for rm in rand_means])

    print(f"  置换检验完成 ({time.time()-t0:.1f}s total)")

    results_df = pd.DataFrame(results)
    return results_df, obs_mean, rand_mean, rand_std, z_score, p_proximity


# ============ 辅助函数 ============
def load_known_bridge_genes():
    df = pd.read_csv(LGBM_FILE)
    bridge = df[df["is_bridge_gene"] == 1]["gene_symbol"].tolist()
    return set(bridge)


def module_enrichment(module_name, genes, consensus_df, n_total=4414):
    cupr_in = [g for g in genes if g in set(consensus_df["gene_symbol"])]
    n_in = len(cupr_in)
    n_cupr = len(genes)
    n_cons = len(consensus_df)
    exp_val = n_cupr * n_cons / n_total
    p_val = hypergeom.sf(n_in - 1, n_total, n_cupr, n_cons)
    return {
        "module": module_name, "n_genes_in_module": n_cupr,
        "n_in_consensus": n_in, "expected": round(exp_val, 2),
        "p_value": p_val, "significant": p_val < 0.05,
        "genes_found": ";".join(cupr_in) if cupr_in else "None",
    }


# ================================================================
# 主流程
# ================================================================
print("=" * 70)
print("  铜死亡基因 × 桥接基因 — 分模块 + 多视角分析 (v2 高效版)")
print("=" * 70)

consensus_df = pd.read_csv(CONSENSUS_FILE)
print(f"\n共识未知基因: {len(consensus_df)}")

# ================================================================
# 分析1: 分功能模块超几何检验
# ================================================================
print("\n" + "-" * 70)
print("  分析1: 铜死亡分功能模块超几何检验")
print("-" * 70)

module_results = []
for mod_name, mod_genes in CUPROPTOSIS_MODULES.items():
    res = module_enrichment(mod_name, mod_genes, consensus_df)
    module_results.append(res)

df_mod = pd.DataFrame(module_results)
df_mod = df_mod.sort_values("p_value")

print(f"\n{'模块':<35s} {'基因数':>5s} {'命中':>5s} {'期望':>6s} {'P-value':>12s} {'显著':>6s}")
print("-" * 75)
for _, row in df_mod.iterrows():
    sig = (" ***" if row["p_value"] < 0.001 else " **" if row["p_value"] < 0.01
           else " *" if row["p_value"] < 0.05 else "")
    print(f"{row['module']:<35s} {row['n_genes_in_module']:>5d} {row['n_in_consensus']:>5d} "
          f"{row['expected']:>6.1f} {row['p_value']:>12.4e}{sig}")

df_mod.to_csv(OUT_ANALYSIS1, index=False)
print(f"\n  结果已保存: {OUT_ANALYSIS1}")

# ================================================================
# 分析2: 全量数据中铜死亡基因的预测分数分布
# ================================================================
print("\n" + "-" * 70)
print("  分析2: 全量数据预测分数分布（不区分未知/已知）")
print("-" * 70)

df_lgbm = pd.read_csv(LGBM_FILE)
df_gat = pd.read_csv(GAT_FILE)

df_labels = df_lgbm[["gene_symbol", "is_drug_target", "is_disease_gene", "is_bridge_gene"]]
df_gat = df_gat.merge(df_labels, on="gene_symbol", how="inner")

df_full = df_lgbm[["gene_symbol", "bridge_probability"]].rename(
    columns={"bridge_probability": "lgbm_score"}
)
df_gat_score = df_gat[["gene_symbol", "combined_score"]].rename(
    columns={"combined_score": "gat_score"}
)
df_full = df_full.merge(df_gat_score, on="gene_symbol", how="inner")
df_full = df_full.merge(df_labels, on="gene_symbol", how="inner")
df_full["is_cuproptosis"] = df_full["gene_symbol"].isin(ALL_CUPROPTOSIS).astype(int)

cupr_mask = df_full["is_cuproptosis"] == 1
non_cupr_mask = df_full["is_cuproptosis"] == 0

cupr_lgbm = df_full.loc[cupr_mask, "lgbm_score"].dropna()
non_lgbm = df_full.loc[non_cupr_mask, "lgbm_score"].dropna()
cupr_gat = df_full.loc[cupr_mask, "gat_score"].dropna()
non_gat = df_full.loc[non_cupr_mask, "gat_score"].dropna()

u_lgbm, p_lgbm = mannwhitneyu(cupr_lgbm, non_lgbm, alternative="two-sided")
u_gat, p_gat = mannwhitneyu(cupr_gat, non_gat, alternative="two-sided")

print(f"\n  LightGBM bridge_probability:")
print(f"    铜死亡 (n={len(cupr_lgbm)}): mean={cupr_lgbm.mean():.6f} median={cupr_lgbm.median():.6f}")
print(f"    非铜死亡 (n={len(non_lgbm)}): mean={non_lgbm.mean():.6f} median={non_lgbm.median():.6f}")
print(f"    Mann-Whitney U: P={p_lgbm:.4e}")

print(f"\n  GAT combined_score:")
print(f"    铜死亡 (n={len(cupr_gat)}): mean={cupr_gat.mean():.6f} median={cupr_gat.median():.6f}")
print(f"    非铜死亡 (n={len(non_gat)}): mean={non_gat.mean():.6f} median={non_gat.median():.6f}")
print(f"    Mann-Whitney U: P={p_gat:.4e}")

cupr_detail = df_full[cupr_mask][["gene_symbol", "is_drug_target", "is_disease_gene",
                                    "is_bridge_gene", "lgbm_score", "gat_score"]].copy()
cupr_detail["module"] = ""
for mod_name, mod_genes in CUPROPTOSIS_MODULES.items():
    cupr_detail.loc[cupr_detail["gene_symbol"].isin(mod_genes), "module"] = mod_name
cupr_detail = cupr_detail.sort_values("lgbm_score", ascending=False)
cupr_detail.to_csv(OUT_ANALYSIS2, index=False)
print(f"\n  详情已保存: {OUT_ANALYSIS2}")

# 28个未在GAT文件中的铜死亡基因
cupr_in_gat = set(cupr_detail["gene_symbol"])
cupr_missing = [g for g in ALL_CUPROPTOSIS if g not in cupr_in_gat]
print(f"\n  未在GAT/共识基因集中的铜死亡基因: {len(cupr_missing)}")
print(f"    {cupr_missing}")

# ================================================================
# 分析3: PPI网络邻近度（多源BFS高效版）
# ================================================================
print("\n" + "-" * 70)
print("  分析3: PPI网络邻近度 — 铜死亡基因 → 已知桥基因 (多源BFS)")
print("-" * 70)

adj = build_adjacency_from_csv(PPI_FILE)
bridge_genes = load_known_bridge_genes()

proximity_df, obs_mean, rand_mean, rand_std, z_score, p_prox = compute_proximity_fast(
    adj, ALL_CUPROPTOSIS, bridge_genes
)

print(f"\n  铜死亡基因 → 最近桥基因距离:")
print(f"    观测平均距离: {obs_mean:.2f}")
print(f"    随机期望:     {rand_mean:.2f} ± {rand_std:.2f}")
print(f"    Z-score:      {z_score:.2f} (负值=比随机更近)")
print(f"    P (perm):     {p_prox:.4e}  " +
      ("*** 显著邻近" if p_prox < 0.001 else "** 较近" if p_prox < 0.01 else
       "* 边缘" if p_prox < 0.05 else "不显著"))

proximity_df["module"] = ""
for mod_name, mod_genes in CUPROPTOSIS_MODULES.items():
    proximity_df.loc[proximity_df["gene"].isin(mod_genes), "module"] = mod_name

print(f"\n  模块级别 (到最近桥基因的平均距离):")
for mod_name in CUPROPTOSIS_MODULES:
    mod_df = proximity_df[proximity_df["module"] == mod_name]
    mod_dist = mod_df[mod_df["dist_to_nearest_bridge"].notna()]
    if len(mod_dist) > 0:
        print(f"    {mod_name}: mean={mod_dist['dist_to_nearest_bridge'].mean():.2f}  "
              f"(n={len(mod_dist)} in graph, n={len(mod_df)-len(mod_dist)} not in graph)")

proximity_df = proximity_df.sort_values("dist_to_nearest_bridge")
proximity_df.to_csv(OUT_ANALYSIS3, index=False)
print(f"\n  邻近度已保存: {OUT_ANALYSIS3}")

# ================================================================
print("\n" + "=" * 70)
print("  三项分析完成！")
print(f"  分析1: {OUT_ANALYSIS1}")
print(f"  分析2: {OUT_ANALYSIS2}")
print(f"  分析3: {OUT_ANALYSIS3}")
print("=" * 70)