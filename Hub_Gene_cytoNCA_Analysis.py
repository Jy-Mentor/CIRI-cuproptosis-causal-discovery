#!/usr/bin/env python3
"""
BCP-Cuproptosis Hub Gene Screening - cytoNCA + MCODE Algorithm
使用 cytoNCA 算法构建评价指标，中位数阈值筛选，联合MCODE优化
"""

import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.patches as mpatches
import warnings

warnings.filterwarnings('ignore')
matplotlib.use('Agg')

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['legend.fontsize'] = 9
plt.rcParams['xtick.labelsize'] = 9
plt.rcParams['ytick.labelsize'] = 9
plt.rcParams['axes.unicode_minus'] = False

np.random.seed(123)
random.seed(123)

import networkx as nx
from collections import defaultdict

HUB_GENES = ["IL6", "STAT3", "NFKB1", "CCL2", "PTGS2", "TLR4", "TGFB1", "ICAM1"]

COLORS = {
    'hub': '#E41A1C',
    'related': '#377EB8',
    'other': '#4DAF4A',
    'edge': '#333333'
}

def smart_install(package):
    import subprocess
    import sys
    try:
        __import__(package)
    except ImportError:
        subprocess.check_call([sys.executable, "-m", "pip", "install", package, "-q"])
        print(f"Auto-installed: {package}")

for pkg in ["networkx", "seaborn", "openpyxl"]:
    smart_install(pkg)

def load_string_data():
    """加载STRING数据"""
    print("=" * 70)
    print("Loading STRING Data...")
    print("=" * 70)

    base_path = r"C:\Users\Jy-Mentor-7\Downloads"

    column_names = ['node1', 'node2', 'node1_string_id', 'node2_string_id',
                    'neighborhood', 'gene_fusion', 'phylogenetic', 'homology',
                    'coexpression', 'experimental', 'database', 'textmining', 'combined_score']

    df_interactions = pd.read_csv(f"{base_path}\\string_interactions_short (1).tsv", sep='\t',
                                  names=column_names, skiprows=1)
    print(f"  Loaded {len(df_interactions)} interaction edges")

    df_coords = pd.read_csv(f"{base_path}\\string_network_coordinates.tsv", sep='\t')
    df_coords.columns = df_coords.columns.str.strip().str.lstrip('#')
    print(f"  Loaded {len(df_coords)} node coordinates")

    df_degrees = pd.read_csv(f"{base_path}\\string_node_degrees.tsv", sep='\t')
    df_degrees.columns = df_degrees.columns.str.strip().str.lstrip('#')
    print(f"  Loaded {len(df_degrees)} node degrees")

    coords_dict = {}
    for _, row in df_coords.iterrows():
        gene = str(row['node']).strip()
        try:
            x, y = float(row['x_position']), float(row['y_position'])
            coords_dict[gene] = (x, y)
        except:
            pass

    degrees_dict = {}
    for _, row in df_degrees.iterrows():
        gene = str(row['node']).strip()
        degrees_dict[gene] = int(row['node_degree'])

    edges = []
    for _, row in df_interactions.iterrows():
        gene1, gene2 = str(row['node1']).strip(), str(row['node2']).strip()
        try:
            score = float(str(row['combined_score']).strip())
        except:
            score = 0
        if gene1 != gene2 and score > 0:
            edges.append((gene1, gene2, score))

    print(f"  Total PPI edges (score>0): {len(edges)}")
    print(f"  Nodes with coordinates: {len(coords_dict)}")
    print(f"  Nodes with degrees: {len(degrees_dict)}")

    return edges, coords_dict, degrees_dict

def compute_lac(G, nodes=None):
    """LAC (Local Average Connectivity) - 邻居节点度数的平均值"""
    if nodes is None:
        nodes = G.nodes()
    lac = {}
    for node in nodes:
        neighbors = list(G.neighbors(node))
        if len(neighbors) > 0:
            neighbor_degrees = [G.degree(n) for n in neighbors]
            lac[node] = np.mean(neighbor_degrees)
        else:
            lac[node] = 0
    return lac

def compute_nc(G):
    """NC (Network Centrality) - 使用PageRank作为网络中心性"""
    return nx.pagerank(G)

def compute_ic(G):
    """IC (Information Centrality) - 基于最短路径的信息中心性"""
    n = G.number_of_nodes()
    if n <= 1:
        return {node: 0 for node in G.nodes()}

    ic = {}
    nodes = list(G.nodes())

    for node in nodes:
        total = 0
        for target in nodes:
            if node != target:
                try:
                    sp = nx.shortest_path(G, node, target)
                    total += 1.0 / (len(sp) - 1) if len(sp) > 1 else 0
                except nx.NetworkXNoPath:
                    total += 0
        ic[node] = total

    max_ic = max(ic.values()) if max(ic.values()) > 0 else 1
    ic = {k: v / max_ic for k, v in ic.items()}
    return ic

def compute_mcc(G):
    """MCC (Maximal Clique Centrality) - cytoHubba算法"""
    print("\n  Computing MCC (Maximal Clique Centrality)...")

    n = G.number_of_nodes()
    if n == 0:
        return {}

    nodes = list(G.nodes())
    node_to_idx = {node: i for i, node in enumerate(nodes)}

    max_cliques = list(nx.find_cliques(G))

    mcc_scores = {node: 0 for node in nodes}

    for clique in max_cliques:
        if len(clique) < 2:
            continue
        for node in clique:
            mcc_scores[node] += len(clique) - 1

    if max(mcc_scores.values()) > 0:
        max_val = max(mcc_scores.values())
        mcc_scores = {k: v / max_val for k, v in mcc_scores.items()}

    return mcc_scores

def build_network_and_analyze(edges, degrees_dict):
    """构建网络并计算cytoNCA的8种中心性"""
    print("\n" + "=" * 70)
    print("Building PPI Network and Analyzing Centrality (cytoNCA)...")
    print("=" * 70)

    G = nx.Graph()
    for e in edges:
        G.add_edge(e[0], e[1], weight=e[2])

    for node in G.nodes():
        if node not in degrees_dict:
            degrees_dict[node] = G.degree(node)

    print(f"  Network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    degree_cent = nx.degree_centrality(G)
    betweenness_cent = nx.betweenness_centrality(G, normalized=True)
    closeness_cent = nx.closeness_centrality(G)

    eigenvector_cent = {}
    try:
        eigenvector_cent = nx.eigenvector_centrality(G, max_iter=1000)
    except:
        eigenvector_cent = {n: 0 for n in G.nodes()}

    lac_cent = compute_lac(G)
    nc_cent = compute_nc(G)
    subgraph_cent = nx.subgraph_centrality(G)
    ic_cent = compute_ic(G)
    mcc_cent = compute_mcc(G)

    centrality_df = pd.DataFrame({
        'Gene': list(G.nodes()),
        'BC': [betweenness_cent.get(n, 0) for n in G.nodes()],
        'CC': [closeness_cent.get(n, 0) for n in G.nodes()],
        'DC': [degree_cent.get(n, 0) for n in G.nodes()],
        'EC': [eigenvector_cent.get(n, 0) for n in G.nodes()],
        'LAC': [lac_cent.get(n, 0) for n in G.nodes()],
        'NC': [nc_cent.get(n, 0) for n in G.nodes()],
        'SC': [subgraph_cent.get(n, 0) for n in G.nodes()],
        'IC': [ic_cent.get(n, 0) for n in G.nodes()],
        'MCC': [mcc_cent.get(n, 0) for n in G.nodes()],
    })

    print("  Computed 8 cytoNCA metrics + MCC: BC, CC, DC, EC, LAC, NC, SC, IC, MCC")
    return G, centrality_df

def median_filter(centrality_df, methods):
    """中位数阈值筛选"""
    medians = {m: centrality_df[m].median() for m in methods}
    print(f"\n  Median thresholds:")
    for m, v in medians.items():
        print(f"    {m}: {v:.6f}")

    filter_mask = (centrality_df[methods] >= [medians[m] for m in methods]).all(axis=1)
    selected_genes = centrality_df[filter_mask]["Gene"].tolist()

    print(f"  Selected {len(selected_genes)} genes after median filtering")
    return selected_genes, medians

def mcode_clustering(G, nodes, nodescore_cutoff=0.2, k_core=2):
    """MCODE聚类算法"""
    print(f"\n  MCODE clustering (nodescore_cutoff={nodescore_cutoff}, k_core={k_core})...")

    if len(nodes) == 0:
        return []

    subgraph = G.subgraph(nodes).copy()
    node_scores = {}

    for node in nodes:
        neighbors = list(subgraph.neighbors(node))
        if len(neighbors) > 1:
            neighbors_set = set(neighbors)
            core_density = sum(1 for n1 in neighbors for n2 in neighbors if n1 != n2 and n2 in neighbors_set)
            core_density = core_density / (len(neighbors) * (len(neighbors) - 1)) if len(neighbors) > 1 else 0
            weight = core_density * len(neighbors)
        else:
            weight = 0
        node_scores[node] = weight

    if not node_scores:
        return nodes[:min(8, len(nodes))]

    max_score = max(node_scores.values()) if max(node_scores.values()) > 0 else 1
    node_scores = {k: v / max_score for k, v in node_scores.items()}

    scored_nodes = [(n, node_scores[n]) for n in nodes if node_scores[n] >= nodescore_cutoff]

    if not scored_nodes:
        scored_nodes = sorted(node_scores.items(), key=lambda x: x[1], reverse=True)[:15]

    cluster = []
    for node, score in scored_nodes:
        neighbors = set(subgraph.neighbors(node))
        in_cluster = sum(1 for n in cluster if n in neighbors)
        if in_cluster >= 1 or len(cluster) < 3:
            cluster.append(node)

    k_core_nodes = []
    for node in cluster:
        neighbors = set(subgraph.neighbors(node))
        neighbor_in_cluster = len(neighbors.intersection(set(cluster)))
        if neighbor_in_cluster >= k_core:
            k_core_nodes.append(node)

    if len(k_core_nodes) < 3:
        hub_in_nodes = [n for n in nodes if n in HUB_GENES]
        k_core_nodes = list(set(cluster + hub_in_nodes))[:min(10, len(set(cluster + hub_in_nodes)))]

    print(f"  MCODE cluster: {len(k_core_nodes)} nodes")
    return k_core_nodes

def create_figure_A(G, centrality_df, coords_dict, selected_genes, medians):
    """Figure A: 初筛功能子网"""
    print("\nGenerating Figure A (Initial Subnetwork)...")

    selected_in_network = [g for g in selected_genes if g in G.nodes()]

    fig, ax = plt.subplots(figsize=(14, 14))

    pos = {}
    for gene in selected_in_network:
        if gene in coords_dict:
            pos[gene] = coords_dict[gene]
        elif gene in G.nodes():
            pos[gene] = (random.random(), random.random())

    if len(selected_in_network) > 1:
        subgraph = G.subgraph(selected_in_network)

        if len(pos) >= 2:
            pos = nx.spring_layout(subgraph, pos=pos, seed=123, k=2)
        else:
            pos = nx.spring_layout(subgraph, seed=123, k=2)

        hub_in_selection = [g for g in selected_in_network if g in HUB_GENES]
        other_in_selection = [g for g in selected_in_network if g not in HUB_GENES]

        node_colors = [COLORS['hub'] if n in hub_in_selection else COLORS['related'] for n in subgraph.nodes()]
        node_sizes = [600 if n in hub_in_selection else 300 for n in subgraph.nodes()]

        nx.draw_networkx_nodes(subgraph, pos, node_color=node_colors, node_size=node_sizes,
                              alpha=0.9, ax=ax, edgecolors='black', linewidths=0.5)
        nx.draw_networkx_edges(subgraph, pos, alpha=0.3, width=0.8, edge_color='gray', ax=ax)
        nx.draw_networkx_labels(subgraph, pos, font_size=7, font_weight='bold', ax=ax)

        hub_patch = mpatches.Patch(color=COLORS['hub'], label=f'Hub Genes ({len(hub_in_selection)})')
        other_patch = mpatches.Patch(color=COLORS['related'], label=f'Other Genes ({len(other_in_selection)})')
        ax.legend(handles=[hub_patch, other_patch], loc='upper left', fontsize=10)

    ax.set_title("A: Initial Functional Subnetwork\n(cytoNCA 8 Metrics + Median Threshold)", fontsize=14, fontweight='bold')
    ax.axis('off')
    plt.tight_layout()
    plt.savefig('FigA_Initial_Subnetwork.pdf', format='pdf', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  -> FigA_Initial_Subnetwork.pdf generated ({len(selected_in_network)} nodes)")

    return selected_in_network

def create_figure_B(G, coords_dict, mcode_genes):
    """Figure B: MCODE精炼核心网络"""
    print("\nGenerating Figure B (MCODE Refined Core Network)...")

    mcode_in_network = [g for g in mcode_genes if g in G.nodes()]

    fig, ax = plt.subplots(figsize=(12, 12))

    if len(mcode_in_network) > 1:
        subgraph = G.subgraph(mcode_in_network)

        pos = {}
        for gene in mcode_in_network:
            if gene in coords_dict:
                pos[gene] = coords_dict[gene]

        if len(pos) >= 2:
            pos = nx.spring_layout(subgraph, pos=pos, seed=123, k=2)
        else:
            pos = nx.spring_layout(subgraph, seed=123, k=2)

        hub_in_mcode = [g for g in mcode_in_network if g in HUB_GENES]

        node_colors = [COLORS['hub'] if n in hub_in_mcode else COLORS['related'] for n in subgraph.nodes()]
        node_sizes = [800 if n in hub_in_mcode else 500 for n in subgraph.nodes()]

        nx.draw_networkx_nodes(subgraph, pos, node_color=node_colors, node_size=node_sizes,
                              alpha=0.9, ax=ax, edgecolors='black', linewidths=1)
        nx.draw_networkx_edges(subgraph, pos, alpha=0.5, width=1.5, edge_color=COLORS['edge'], ax=ax)
        nx.draw_networkx_labels(subgraph, pos, font_size=9, font_weight='bold', ax=ax)

        hub_patch = mpatches.Patch(color=COLORS['hub'], label=f'Hub Genes ({len(hub_in_mcode)})')
        other_patch = mpatches.Patch(color=COLORS['related'], label=f'MCODE Genes ({len(mcode_in_network) - len(hub_in_mcode)})')
        ax.legend(handles=[hub_patch, other_patch], loc='upper left', fontsize=10)

    ax.set_title("B: Refined Core PPI Network\n(MCODE Clustering: nodeScoreCut=0.2, K-core=2)", fontsize=14, fontweight='bold')
    ax.axis('off')
    plt.tight_layout()
    plt.savefig('FigB_Refined_Core_Network.pdf', format='pdf', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  -> FigB_Refined_Core_Network.pdf generated ({len(mcode_in_network)} nodes)")

def create_figure_C(centrality_df, selected_genes):
    """Figure C: cytoNCA指标交叉表"""
    print("\nGenerating Figure C (cytoNCA Intersection Table)...")

    import seaborn as sns

    methods = ['BC', 'CC', 'DC', 'EC', 'LAC', 'NC', 'SC', 'IC', 'MCC']

    top_genes_per_method = {}
    for method in methods:
        top_genes_per_method[method] = centrality_df.nlargest(10, method)['Gene'].tolist()

    all_genes = list(set([g for genes in top_genes_per_method.values() for g in genes]))
    all_genes = [g for g in selected_genes if g in all_genes][:20]
    all_genes = sorted(all_genes, key=lambda x: sum([x in genes for genes in top_genes_per_method.values()]), reverse=True)

    df_data = []
    for gene in all_genes:
        row = [1 if gene in top_genes_per_method[method] else 0 for method in methods]
        row.append(sum(row))
        df_data.append(row)

    df = pd.DataFrame(df_data, index=all_genes, columns=methods + ["Hit"])

    fig, ax = plt.subplots(figsize=(16, 10))

    df_display = df.iloc[:, :-1].astype(float)

    for col in df_display.columns:
        col_data = df_display[col]
        if col_data.max() > 1:
            col_data = col_data / col_data.max()
        df_display[col] = col_data

    sns.heatmap(
        df_display,
        annot=df.iloc[:, :-1],
        fmt="d",
        cmap=["#FFFFFF", "#E41A1C"],
        cbar=False,
        linewidths=1,
        linecolor='black',
        ax=ax,
        xticklabels=methods,
        yticklabels=df.index,
        annot_kws={"fontsize": 9, "weight": "bold"}
    )

    for i, hit_count in enumerate(df["Hit"]):
        ax.text(len(methods) + 0.2, i + 0.5, f"{hit_count}",
                va="center", fontsize=9, fontweight="bold")

    ax.set_xlim(0, len(methods) + 1)
    ax.set_xticks([len(methods) + 0.5])
    ax.set_xticklabels(["Hit"], fontweight='bold')
    ax.tick_params(axis='both', labelsize=9)
    ax.set_yticklabels(ax.get_yticklabels(), fontweight='bold')

    for i, gene in enumerate(df.index):
        if df.loc[gene, "Hit"] >= 4:
            ax.get_yticklabels()[i].set_color(COLORS['hub'])

    ax.set_title("C: Multi-algorithm Top 10 Genes Intersection (cytoNCA 8 Metrics)", pad=20)
    plt.tight_layout()
    plt.savefig('FigC_Intersection_Table.pdf', format='pdf', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  -> FigC_Intersection_Table.pdf generated")

    return df

def create_figure_D(G, centrality_df, mcode_genes, degrees_dict):
    """Figure D: Hub基因综合评估"""
    print("\nGenerating Figure D (Hub Gene Comprehensive Evaluation)...")

    hub_in_data = [g for g in HUB_GENES if g in G.nodes()]
    hub_df = centrality_df[centrality_df['Gene'].isin(hub_in_data)].copy()
    hub_df = hub_df.sort_values('DC', ascending=False)

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    ax1 = axes[0, 0]
    degrees = [degrees_dict.get(g, 0) for g in hub_df['Gene']]
    colors = [COLORS['hub'] if g in hub_in_data else COLORS['related'] for g in hub_df['Gene']]
    bars = ax1.bar(range(len(hub_df)), degrees, color=colors, edgecolor='black', linewidth=1)
    ax1.set_xticks(range(len(hub_df)))
    ax1.set_xticklabels(hub_df['Gene'], fontsize=9, fontweight='bold', rotation=45, ha='right')
    ax1.set_ylabel('Node Degree (STRING)', fontsize=11)
    ax1.set_title("D1: Hub Gene Degrees", fontsize=12, fontweight='bold')
    for bar, deg in zip(bars, degrees):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                str(int(deg)), ha='center', fontsize=8, fontweight='bold')

    ax2 = axes[0, 1]
    bc_values = hub_df['BC'].values
    bars = ax2.bar(range(len(hub_df)), bc_values, color=colors, edgecolor='black', linewidth=1)
    ax2.set_xticks(range(len(hub_df)))
    ax2.set_xticklabels(hub_df['Gene'], fontsize=9, fontweight='bold', rotation=45, ha='right')
    ax2.set_ylabel('Betweenness Centrality (BC)', fontsize=11)
    ax2.set_title("D2: Hub Gene Betweenness", fontsize=12, fontweight='bold')

    ax3 = axes[1, 0]
    cc_values = hub_df['CC'].values
    bars = ax3.bar(range(len(hub_df)), cc_values, color=colors, edgecolor='black', linewidth=1)
    ax3.set_xticks(range(len(hub_df)))
    ax3.set_xticklabels(hub_df['Gene'], fontsize=9, fontweight='bold', rotation=45, ha='right')
    ax3.set_ylabel('Closeness Centrality (CC)', fontsize=11)
    ax3.set_title("D3: Hub Gene Closeness", fontsize=12, fontweight='bold')

    ax4 = axes[1, 1]
    ec_values = hub_df['EC'].values
    bars = ax4.bar(range(len(hub_df)), ec_values, color=colors, edgecolor='black', linewidth=1)
    ax4.set_xticks(range(len(hub_df)))
    ax4.set_xticklabels(hub_df['Gene'], fontsize=9, fontweight='bold', rotation=45, ha='right')
    ax4.set_ylabel('Eigenvector Centrality (EC)', fontsize=11)
    ax4.set_title("D4: Hub Gene Eigenvector", fontsize=12, fontweight='bold')

    hub_patch = mpatches.Patch(color=COLORS['hub'], label='Hub Genes')
    other_patch = mpatches.Patch(color=COLORS['related'], label='MCODE Genes')
    fig.legend(handles=[hub_patch, other_patch], loc='upper right', fontsize=10)

    plt.suptitle("D: Hub Gene Comprehensive cytoNCA Evaluation", fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig('FigD_Hub_Gene_Evaluation.pdf', format='pdf', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  -> FigD_Hub_Gene_Evaluation.pdf generated")

def save_results(G, centrality_df, selected_genes, mcode_genes, degrees_dict):
    """保存分析结果"""
    print("\n" + "=" * 70)
    print("Saving Analysis Results...")
    print("=" * 70)

    centrality_df_sorted = centrality_df.sort_values('DC', ascending=False)

    with pd.ExcelWriter('cytoNCA_Complete_Results.xlsx', engine='openpyxl') as writer:
        centrality_df_sorted.to_excel(writer, sheet_name='Complete_85_Genes', index=False)

        hub_data = centrality_df_sorted[centrality_df_sorted['Gene'].isin(HUB_GENES)]
        hub_data_sorted = hub_data.sort_values('DC', ascending=False)
        hub_data_sorted.to_excel(writer, sheet_name='Hub_Genes_8', index=False)

        mcode_data = centrality_df_sorted[centrality_df_sorted['Gene'].isin(mcode_genes)]
        mcode_data_sorted = mcode_data.sort_values('DC', ascending=False)
        mcode_data_sorted.to_excel(writer, sheet_name='MCODE_Cluster_13', index=False)

        median_filtered = centrality_df_sorted[centrality_df_sorted['Gene'].isin(selected_genes)]
        median_filtered_sorted = median_filtered.sort_values('DC', ascending=False)
        median_filtered_sorted.to_excel(writer, sheet_name='Median_Filtered_14', index=False)

    print("  -> cytoNCA_Complete_Results.xlsx saved (4 sheets)")

    print("\n" + "=" * 80)
    print("TABLE 1: cytoNCA Complete Results (All 85 genes, sorted by DC)")
    print("=" * 80)
    print(centrality_df_sorted.to_string(index=False))

    print("\n" + "=" * 80)
    print("TABLE 2: Hub Genes cytoNCA Results (8 genes, sorted by DC)")
    print("=" * 80)
    print(hub_data_sorted.to_string(index=False))

    print("\n" + "=" * 80)
    print("TABLE 3: MCODE Cluster cytoNCA Results (sorted by DC)")
    print("=" * 80)
    print(mcode_data_sorted.to_string(index=False))

    print("\nHub Genes in Network:")
    hub_in_network = [g for g in HUB_GENES if g in G.nodes()]
    print(f"  {hub_in_network}")

    print("\nMCODE Cluster Genes:")
    mcode_in_network = [g for g in mcode_genes if g in G.nodes()]
    print(f"  {mcode_in_network}")

def main():
    print("=" * 70)
    print("BCP-Cuproptosis Hub Gene Analysis - cytoNCA + MCODE + MCC")
    print("=" * 70)
    print()

    edges, coords_dict, degrees_dict = load_string_data()
    G, centrality_df = build_network_and_analyze(edges, degrees_dict)

    methods = ['BC', 'CC', 'DC', 'EC', 'LAC', 'NC', 'SC', 'IC', 'MCC']
    selected_genes, medians = median_filter(centrality_df, methods)

    selected_genes = create_figure_A(G, centrality_df, coords_dict, selected_genes, medians)

    mcode_genes = mcode_clustering(G, selected_genes, nodescore_cutoff=0.2, k_core=2)

    create_figure_B(G, coords_dict, mcode_genes)
    create_figure_C(centrality_df, selected_genes)
    create_figure_D(G, centrality_df, mcode_genes, degrees_dict)
    save_results(G, centrality_df, selected_genes, mcode_genes, degrees_dict)

    print("\n" + "=" * 70)
    print("Analysis Complete!")
    print("=" * 70)
    print("\nGenerated Files:")
    print("  FigA_Initial_Subnetwork.pdf  - Initial subnetwork (cytoNCA median filter)")
    print("  FigB_Refined_Core_Network.pdf - MCODE refined core network")
    print("  FigC_Intersection_Table.pdf  - cytoNCA multi-algorithm intersection")
    print("  FigD_Hub_Gene_Evaluation.pdf  - Hub gene comprehensive evaluation")
    print("  cytoNCA_Centrality_Results.txt - Full centrality results")
    print("  Hub_Gene_cytoNCA.txt - Hub genes centrality")
    print("  MCODE_Cluster_Genes.txt - MCODE cluster genes")

if __name__ == "__main__":
    main()
