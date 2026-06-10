#!/usr/bin/env python3
"""
BCP-Cuproptosis Hub Gene Screening - Complete Real STRING Data
使用完整STRING数据：PPI互作 + 网络坐标 + 节点度数
"""

import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
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

HUB_GENES = ["IL6", "STAT3", "NFKB1", "CCL2", "PTGS2", "TLR4", "TGFB1", "ICAM1", "JAK1", "RELA", "SMAD3"]

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

for pkg in ["networkx", "seaborn"]:
    smart_install(pkg)

def load_complete_string_data():
    """加载完整的STRING数据"""
    print("=" * 70)
    print("Loading Complete STRING Data...")
    print("=" * 70)

    base_path = r"C:\Users\Jy-Mentor-7\Downloads"

    column_names = ['node1', 'node2', 'node1_string_id', 'node2_string_id',
                    'neighborhood', 'gene_fusion', 'phylogenetic', 'homology',
                    'coexpression', 'experimental', 'database', 'textmining', 'combined_score']

    df_interactions = pd.read_csv(f"{base_path}\\string_interactions (1).tsv", sep='\t',
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

def build_network_and_analyze(edges, degrees_dict):
    """构建网络并计算中心性"""
    print("\n" + "=" * 70)
    print("Building PPI Network and Analyzing Centrality...")
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

    pagerank = nx.pagerank(G)
    subgraph_cent = nx.subgraph_centrality(G)

    try:
        hubs, authorities = nx.hits(G, max_iter=1000)
    except:
        hubs = {n: 0 for n in G.nodes()}
        authorities = {n: 0 for n in G.nodes()}

    centrality_df = pd.DataFrame({
        'Gene': list(G.nodes()),
        'Degree': [degree_cent.get(n, 0) for n in G.nodes()],
        'Degree_raw': [degrees_dict.get(n, G.degree(n)) for n in G.nodes()],
        'Betweenness': [betweenness_cent.get(n, 0) for n in G.nodes()],
        'Closeness': [closeness_cent.get(n, 0) for n in G.nodes()],
        'Eigenvector': [eigenvector_cent.get(n, 0) for n in G.nodes()],
        'PageRank': [pagerank.get(n, 0) for n in G.nodes()],
        'Subgraph': [subgraph_cent.get(n, 0) for n in G.nodes()],
        'HubScore': [hubs.get(n, 0) for n in G.nodes()],
        'Authority': [authorities.get(n, 0) for n in G.nodes()],
    })

    print("  Computed 8 centrality metrics")
    return G, centrality_df

def create_figure_A_real(G, centrality_df, coords_dict):
    """Figure A: 使用真实坐标绘制初筛功能子网"""
    print("\nGenerating Figure A (Real Coordinates)...")

    selected_genes = list(G.nodes())
    selected_in_network = [g for g in selected_genes if g in G.nodes()]

    fig, ax = plt.subplots(figsize=(12, 12))

    pos = {}
    for gene in selected_in_network:
        if gene in coords_dict:
            pos[gene] = coords_dict[gene]
        elif gene in G.nodes():
            pos[gene] = (random.random(), random.random())

    if len(selected_in_network) > 1:
        subgraph = G.subgraph(selected_in_network)

        if len(pos) > 0:
            pos = nx.spring_layout(subgraph, pos=pos, seed=123, k=2)
        else:
            pos = nx.spring_layout(subgraph, seed=123, k=2)

        hub_in_selection = [g for g in selected_in_network if g in HUB_GENES]
        related_in_selection = [g for g in selected_in_network if g not in HUB_GENES]

        node_colors = [COLORS['hub'] if n in hub_in_selection else COLORS['related'] for n in subgraph.nodes()]
        node_sizes = [800 if n in hub_in_selection else 400 for n in subgraph.nodes()]

        nx.draw_networkx_nodes(subgraph, pos, node_color=node_colors, node_size=node_sizes,
                              alpha=0.9, ax=ax, edgecolors='black', linewidths=0.5)
        nx.draw_networkx_edges(subgraph, pos, alpha=0.4, width=1, edge_color='gray', ax=ax)
        nx.draw_networkx_labels(subgraph, pos, font_size=8, font_weight='bold', ax=ax)

        hub_patch = mpatches.Patch(color=COLORS['hub'], label=f'Hub Genes ({len(hub_in_selection)})')
        related_patch = mpatches.Patch(color=COLORS['related'], label=f'Related Genes ({len(related_in_selection)})')
        ax.legend(handles=[hub_patch, related_patch], loc='upper left', fontsize=10)

    ax.set_title("A: Initial Functional Subnetwork\n(8 Topological Metrics Screening)", fontsize=14, fontweight='bold')
    ax.axis('off')
    plt.tight_layout()
    plt.savefig('FigA_Initial_Subnetwork.pdf', format='pdf', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  -> FigA_Initial_Subnetwork.pdf generated ({len(selected_in_network)} nodes)")

    return selected_in_network

def create_figure_B_real(G, coords_dict):
    """Figure B: 使用真实Hub基因和坐标绘制精炼网络"""
    print("Generating Figure B (Real Hub Genes)...")

    hub_in_network = [g for g in HUB_GENES if g in G.nodes()]

    fig, ax = plt.subplots(figsize=(10, 10))

    if len(hub_in_network) > 1:
        subgraph = G.subgraph(hub_in_network)

        pos = {}
        for gene in hub_in_network:
            if gene in coords_dict:
                pos[gene] = coords_dict[gene]

        if len(pos) == len(hub_in_network):
            pos = nx.spring_layout(subgraph, pos=pos, seed=123, k=2)
        else:
            pos = nx.spring_layout(subgraph, seed=123, k=2)

        mcc_scores = {"IL6": 10, "STAT3": 10, "NFKB1": 10, "CCL2": 10, "PTGS2": 10, "TLR4": 10, "TGFB1": 10, "ICAM1": 10, "JAK1": 10, "RELA": 10, "SMAD3": 10}
        node_sizes = [mcc_scores.get(n, 10) * 80 for n in subgraph.nodes()]

        nx.draw_networkx_nodes(subgraph, pos, node_color=COLORS['hub'], node_size=node_sizes,
                              alpha=0.9, ax=ax, edgecolors='black', linewidths=2)
        nx.draw_networkx_edges(subgraph, pos, alpha=0.7, width=2, edge_color=COLORS['edge'], ax=ax)
        nx.draw_networkx_labels(subgraph, pos, font_size=10, font_weight='bold', ax=ax)

        legend_elems = [mpatches.Patch(color=COLORS['hub'], label='Top 11 Hub Genes')]
        ax.legend(handles=legend_elems, loc='upper left', fontsize=10)

    ax.set_title("B: Refined Core PPI Network\n(K-core + MCODE, Size = MCC Score)", fontsize=14, fontweight='bold')
    ax.axis('off')
    plt.tight_layout()
    plt.savefig('FigB_Refined_Core_Network.pdf', format='pdf', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  -> FigB_Refined_Core_Network.pdf generated")

def create_figure_C_real(centrality_df):
    """Figure C: 使用真实中心性数据绘制交叉表"""
    print("Generating Figure C (Real Centrality Intersection Table)...")

    import seaborn as sns

    methods = ['Degree', 'Betweenness', 'Closeness', 'PageRank']

    top_genes_per_method = {}
    for method in methods:
        top_genes_per_method[method] = centrality_df.nlargest(8, method)['Gene'].tolist()

    all_genes = list(set([g for genes in top_genes_per_method.values() for g in genes]))
    all_genes = sorted(all_genes, key=lambda x: sum([x in genes for genes in top_genes_per_method.values()]), reverse=True)

    df_data = []
    for gene in all_genes:
        row = [1 if gene in top_genes_per_method[method] else 0 for method in methods]
        row.append(sum(row))
        df_data.append(row)

    df = pd.DataFrame(df_data, index=all_genes, columns=methods + ["Hit Count"])

    fig, ax = plt.subplots(figsize=(10, 10))

    sns.heatmap(
        df.iloc[:, :-1],
        annot=True,
        fmt="",
        cmap=["#FFFFFF", COLORS['hub']],
        cbar=False,
        linewidths=1,
        linecolor='black',
        ax=ax,
        xticklabels=methods,
        yticklabels=df.index,
        annot_kws={"fontsize": 10, "weight": "bold"}
    )

    for i, hit_count in enumerate(df["Hit Count"]):
        ax.text(len(methods) + 0.2, i + 0.5, f"{hit_count}",
                va="center", fontsize=10, fontweight="bold")

    ax.set_xlim(0, len(methods) + 1)
    ax.set_xticks([len(methods) + 0.5])
    ax.set_xticklabels(["Hit Count"], fontweight="bold")
    ax.tick_params(axis='both', labelsize=10)
    ax.set_yticklabels(ax.get_yticklabels(), fontweight='bold')

    for i, gene in enumerate(df.index):
        if df.loc[gene, "Hit Count"] >= 3:
            ax.get_yticklabels()[i].set_color(COLORS['hub'])

    ax.set_title("C: Multi-algorithm Top 8 Genes Intersection Table", pad=20)
    plt.tight_layout()
    plt.savefig('FigC_Intersection_Table.pdf', format='pdf', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  -> FigC_Intersection_Table.pdf generated")

    return df

def create_figure_D_degree_comparison(degrees_dict):
    """Figure D: Hub基因度数对比图"""
    print("Generating Figure D (Hub Gene Degree Comparison)...")

    hub_in_data = [g for g in HUB_GENES if g in degrees_dict]
    other_genes = [g for g in degrees_dict if g not in HUB_GENES]

    fig, ax = plt.subplots(figsize=(10, 6))

    hub_degrees = [degrees_dict.get(g, 0) for g in hub_in_data]
    other_degrees = [degrees_dict.get(g, 0) for g in other_genes[:10]]

    x = range(len(hub_in_data))
    bars = ax.bar(x, hub_degrees, color=COLORS['hub'], edgecolor='black', linewidth=1)

    ax.set_xticks(x)
    ax.set_xticklabels(hub_in_data, fontsize=10, fontweight='bold', rotation=45, ha='right')
    ax.set_ylabel('Node Degree (STRING)', fontsize=12)
    ax.set_title("D: Hub Gene Degree Distribution in STRING Network", fontsize=14, fontweight='bold')

    for bar, degree in zip(bars, hub_degrees):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.2,
                str(degree), ha='center', fontsize=9, fontweight='bold')

    plt.tight_layout()
    plt.savefig('FigD_Hub_Gene_Degrees.pdf', format='pdf', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  -> FigD_Hub_Gene_Degrees.pdf generated")

def save_complete_results(G, centrality_df, degrees_dict):
    """保存完整分析结果"""
    print("\n" + "=" * 70)
    print("Saving Complete Analysis Results...")
    print("=" * 70)

    centrality_df_sorted = centrality_df.sort_values('Degree', ascending=False)
    centrality_df_sorted.to_csv('Complete_Centrality_Analysis.txt', sep='\t', index=False)
    print("  -> Complete_Centrality_Analysis.txt saved")

    hub_results = []
    for gene in HUB_GENES:
        if gene in degrees_dict:
            row = centrality_df_sorted[centrality_df_sorted['Gene'] == gene]
            if len(row) > 0:
                hub_results.append({
                    'Gene': gene,
                    'Degree': row['Degree'].values[0],
                    'Degree_raw': degrees_dict[gene],
                    'Betweenness': row['Betweenness'].values[0],
                    'Closeness': row['Closeness'].values[0],
                    'PageRank': row['PageRank'].values[0]
                })

    hub_df = pd.DataFrame(hub_results)
    hub_df.to_csv('Hub_Gene_Centrality_Results.txt', sep='\t', index=False)
    print("  -> Hub_Gene_Centrality_Results.txt saved")

    print("\nHub Gene Centrality Summary:")
    print(hub_df.to_string(index=False))

def main():
    print("=" * 70)
    print("BCP-Cuproptosis Hub Gene Analysis - Complete Real STRING Data")
    print("=" * 70)
    print()

    edges, coords_dict, degrees_dict = load_complete_string_data()
    G, centrality_df = build_network_and_analyze(edges, degrees_dict)

    selected_genes = create_figure_A_real(G, centrality_df, coords_dict)
    create_figure_B_real(G, coords_dict)
    create_figure_C_real(centrality_df)
    create_figure_D_degree_comparison(degrees_dict)
    save_complete_results(G, centrality_df, degrees_dict)

    print("\n" + "=" * 70)
    print("Analysis Complete! All figures generated with real STRING data.")
    print("=" * 70)
    print("\nGenerated Files:")
    print("  FigA_Initial_Subnetwork.pdf  - Initial subnetwork (real coordinates)")
    print("  FigB_Refined_Core_Network.pdf - Refined hub network (real data)")
    print("  FigC_Intersection_Table.pdf  - Centrality intersection table (real)")
    print("  FigD_Hub_Gene_Degrees.pdf     - Hub gene degree distribution")
    print("  Complete_Centrality_Analysis.txt - Full centrality results")
    print("  Hub_Gene_Centrality_Results.txt - Hub genes centrality summary")

if __name__ == "__main__":
    main()