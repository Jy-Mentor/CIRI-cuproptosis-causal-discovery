#!/usr/bin/env python3
"""
铜死亡基因子网络隔离分析（权威范式）

流程：
1. 提取 GAT 预测 Top200 候选靶点 + 15 个铜死亡基因（CRGs）。
2. 用 STRING 构建这 215 个基因的诱导子网络（induced subgraph），combined_score >= 0.4。
3. 在该子网络内重新计算度中心性、Betweenness、PageRank。
4. 分析 CRGs 在子网络中的拓扑地位，以及它们与 Top200 靶点的直接/间接连接路径。

依据：Frontiers in Cellular Neuroscience 2024 — 铜死亡基因在全基因组 PPI 中天然为低度节点，
     禁止在全基因组 6000+ 节点网络中要求铜死亡基因进入 Top200。
"""

import json
import logging
from pathlib import Path
from typing import Dict, List, Set, Tuple

import numpy as np
import pandas as pd
import networkx as nx
from scipy import stats

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
logger = logging.getLogger(__name__)

CUPROPTOSIS_GENES = {
    "FDX1", "LIAS", "LIPT1", "DLAT", "PDHB", "PDHX",
    "SLC31A1", "ATP7A", "ATP7B", "ATOX1", "NFE2L2",
    "HIF1A", "MTOR", "NFKB1", "GPX4",
}


def build_induced_subgraph(
    top_genes: List[str],
    cupro_genes: Set[str],
    edge_df: pd.DataFrame,
    min_weight: float = 0.4,
) -> nx.Graph:
    """
    构建诱导子网络（induced subgraph）

    Args:
        top_genes: Top200 候选靶点基因列表
        cupro_genes: 铜死亡基因集合
        edge_df: 全基因组边 DataFrame（Source, Target, Weight）
        min_weight: 最小边权重阈值（combined_score >= 0.4）

    Returns:
        G: 诱导子网络（无向图）
    """
    selected_genes = set(top_genes) | cupro_genes

    # 过滤边：两端都在 selected_genes 中，且权重 >= min_weight
    sub_edges = edge_df[
        (edge_df["Source"].isin(selected_genes)) &
        (edge_df["Target"].isin(selected_genes)) &
        (edge_df["Weight"] >= min_weight)
    ].copy()

    G = nx.Graph()
    for _, row in sub_edges.iterrows():
        G.add_edge(row["Source"], row["Target"], weight=row["Weight"])

    # 确保所有 selected_genes 都在图中（孤立节点也保留）
    for g in selected_genes:
        if g not in G:
            G.add_node(g)

    logger.info(f"诱导子网络: {G.number_of_nodes()} 节点, {G.number_of_edges()} 边")
    return G


def compute_subgraph_topology(G: nx.Graph, cupro_genes: Set[str]) -> pd.DataFrame:
    """
    在子网络内重新计算拓扑指标

    Returns:
        df: 每个基因的拓扑指标 DataFrame
    """
    genes = sorted(G.nodes())

    # 度中心性（子网络内）
    degree_dict = dict(G.degree())
    # Betweenness（子网络内）
    betweenness_dict = nx.betweenness_centrality(G, weight="weight")
    # PageRank（子网络内）
    pagerank_dict = nx.pagerank(G, weight="weight")
    # 聚类系数（子网络内）
    clustering_dict = nx.clustering(G)

    records = []
    for g in genes:
        records.append({
            "GeneSymbol": g,
            "Degree": degree_dict.get(g, 0),
            "Betweenness": betweenness_dict.get(g, 0.0),
            "PageRank": pagerank_dict.get(g, 0.0),
            "ClusteringCoefficient": clustering_dict.get(g, 0.0),
            "is_cuproptosis": 1 if g in cupro_genes else 0,
        })

    df = pd.DataFrame(records)

    # 子网络内的度中心性归一化（Min-Max）
    d_min = df["Degree"].min()
    d_max = df["Degree"].max()
    if d_max > d_min:
        df["Degree_norm"] = (df["Degree"] - d_min) / (d_max - d_min)
    else:
        df["Degree_norm"] = 0.0

    return df


def analyze_cuproptosis_connectivity(
    G: nx.Graph,
    cupro_genes: Set[str],
    top_genes: List[str],
) -> pd.DataFrame:
    """
    分析铜死亡基因与 Top200 靶点的连接路径

    Returns:
        df: 每个铜死亡基因的连通性分析结果
    """
    records = []
    for cg in sorted(cupro_genes):
        if cg not in G:
            continue

        # 直接连接的 Top200 靶点
        direct_neighbors = set(G.neighbors(cg)) & set(top_genes)
        n_direct = len(direct_neighbors)

        # 间接连接（2-hop）的 Top200 靶点
        indirect_neighbors = set()
        for neighbor in G.neighbors(cg):
            for nn in G.neighbors(neighbor):
                if nn in top_genes and nn != cg:
                    indirect_neighbors.add(nn)
        n_indirect = len(indirect_neighbors - direct_neighbors)

        # 平均边权重（直接连接）
        avg_weight = 0.0
        if n_direct > 0:
            weights = [G[cg][n]["weight"] for n in direct_neighbors if "weight" in G[cg][n]]
            if weights:
                avg_weight = np.mean(weights)

        # 到 Top200 靶点的最短路径长度（平均）
        path_lengths = []
        for tg in top_genes:
            if tg in G and nx.has_path(G, cg, tg):
                path_lengths.append(nx.shortest_path_length(G, cg, tg))
        avg_path_length = np.mean(path_lengths) if path_lengths else np.nan

        records.append({
            "Cuproptosis_Gene": cg,
            "Direct_Top200_Connections": n_direct,
            "Indirect_Top200_Connections": n_indirect,
            "Avg_Edge_Weight": avg_weight,
            "Avg_Shortest_Path_to_Top200": avg_path_length,
            "Directly_Connected_Top200": ", ".join(sorted(direct_neighbors)) if direct_neighbors else "None",
        })

    return pd.DataFrame(records)


def compare_cuproptosis_vs_others(topo_df: pd.DataFrame) -> Dict:
    """
    比较铜死亡基因与非铜死亡基因在子网络中的拓扑地位

    Returns:
        stats_dict: 统计比较结果
    """
    cupro_df = topo_df[topo_df["is_cuproptosis"] == 1]
    other_df = topo_df[topo_df["is_cuproptosis"] == 0]

    results = {}
    for col in ["Degree", "Betweenness", "PageRank", "ClusteringCoefficient"]:
        cupro_vals = cupro_df[col].dropna()
        other_vals = other_df[col].dropna()

        if len(cupro_vals) > 0 and len(other_vals) > 0:
            # Mann-Whitney U 检验（非参数）
            statistic, pvalue = stats.mannwhitneyu(cupro_vals, other_vals, alternative="two-sided")
            results[col] = {
                "cupro_median": float(cupro_vals.median()),
                "other_median": float(other_vals.median()),
                "mannwhitney_u": float(statistic),
                "pvalue": float(pvalue),
                "significant": pvalue < 0.05,
            }
        else:
            results[col] = {
                "cupro_median": float(cupro_vals.median()) if len(cupro_vals) > 0 else np.nan,
                "other_median": float(other_vals.median()) if len(other_vals) > 0 else np.nan,
                "mannwhitney_u": np.nan,
                "pvalue": np.nan,
                "significant": False,
            }

    return results


def main():
    processed_dir = Path("./processed")
    results_dir = Path("./results")
    output_dir = Path("./subgraph_analysis")
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1. 加载 Top200 候选靶点
    top200_path = results_dir / "top_targets_200.csv"
    if not top200_path.exists():
        logger.error(f"Top200 文件不存在: {top200_path}")
        return

    top200_df = pd.read_csv(top200_path)
    top_genes = top200_df["GeneSymbol"].tolist()
    logger.info(f"加载 Top200 候选靶点: {len(top_genes)} 个")

    # 2. 加载全基因组边
    edge_df = pd.read_csv(processed_dir / "edge_index.csv")
    logger.info(f"加载全基因组边: {len(edge_df)} 条")

    # 3. 构建诱导子网络
    G = build_induced_subgraph(top_genes, CUPROPTOSIS_GENES, edge_df, min_weight=0.4)

    # 4. 计算子网络拓扑指标
    topo_df = compute_subgraph_topology(G, CUPROPTOSIS_GENES)
    topo_df.to_csv(output_dir / "subgraph_topology.csv", index=False)
    logger.info(f"子网络拓扑指标已保存: {output_dir / 'subgraph_topology.csv'}")

    # 5. 铜死亡基因连通性分析
    connectivity_df = analyze_cuproptosis_connectivity(G, CUPROPTOSIS_GENES, top_genes)
    connectivity_df.to_csv(output_dir / "cuproptosis_connectivity.csv", index=False)
    logger.info(f"铜死亡基因连通性分析已保存: {output_dir / 'cuproptosis_connectivity.csv'}")

    # 6. 统计比较
    comparison = compare_cuproptosis_vs_others(topo_df)
    with open(output_dir / "topology_comparison.json", "w", encoding="utf-8") as f:
        # 将 numpy bool 转换为 Python bool
        comparison_serializable = {}
        for col, stat in comparison.items():
            comparison_serializable[col] = {k: bool(v) if isinstance(v, (np.bool_, bool)) else float(v) if isinstance(v, (np.floating, float)) else int(v) if isinstance(v, (np.integer, int)) else v for k, v in stat.items()}
        json.dump(comparison_serializable, f, indent=2, ensure_ascii=False)
    logger.info(f"拓扑统计比较已保存: {output_dir / 'topology_comparison.json'}")

    # 7. 打印摘要
    logger.info("=" * 60)
    logger.info("铜死亡基因子网络隔离分析摘要")
    logger.info("=" * 60)

    cupro_topo = topo_df[topo_df["is_cuproptosis"] == 1].sort_values("Degree", ascending=False)
    logger.info(f"铜死亡基因在子网络中的度中心性排名:")
    for _, row in cupro_topo.iterrows():
        logger.info(
            f"  {row['GeneSymbol']:8s} 度={row['Degree']:2.0f} "
            f"Betweenness={row['Betweenness']:.4f} "
            f"PageRank={row['PageRank']:.6f}"
        )

    logger.info("")
    logger.info(f"铜死亡基因与 Top200 靶点的直接连接数:")
    for _, row in connectivity_df.iterrows():
        logger.info(
            f"  {row['Cuproptosis_Gene']:8s} 直接={row['Direct_Top200_Connections']:2d} "
            f"间接={row['Indirect_Top200_Connections']:2d} "
            f"平均路径={row['Avg_Shortest_Path_to_Top200']:.2f}"
        )

    logger.info("")
    logger.info("拓扑统计比较 (Mann-Whitney U):")
    for col, stat in comparison.items():
        sig = "*" if stat["significant"] else ""
        logger.info(
            f"  {col:20s}: 铜死亡中位数={stat['cupro_median']:.4f}, "
            f"其他中位数={stat['other_median']:.4f}, p={stat['pvalue']:.4f}{sig}"
        )

    logger.info("=" * 60)
    logger.info("分析完成")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
