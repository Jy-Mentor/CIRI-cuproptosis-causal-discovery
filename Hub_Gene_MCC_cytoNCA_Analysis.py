#!/usr/bin/env python3
"""
Hub Gene Analysis: MCC + cytoNCA Combined Algorithm
Step 1: Calculate MCC using cytoHubba algorithm
Step 2: Calculate 8 cytoNCA centrality metrics (BC, CC, DC, EC, LAC, NC, SC, IC)
Step 3: Combine scores for final ranking
Step 4: Group Hub genes into Group A (Copper Death Hub) and Group B (Regulatory Hub)
"""

import pandas as pd
import numpy as np
import networkx as nx
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

def smart_install(package):
    import subprocess
    import sys
    try:
        __import__(package)
    except ImportError:
        subprocess.check_call([sys.executable, "-m", "pip", "install", package, "-q"])

for pkg in ["networkx", "openpyxl"]:
    smart_install(pkg)

COPPER_DEATH_GENES = {
    'FDX1', 'LIPT1', 'LIPT2', 'MTF1', 'SLC31A1', 'SLC31A2',
    'ATP7A', 'ATP7B', 'DLD', 'DLAT', 'PDHA1', 'PDHB',
    'LIAS', 'LIPTC1', 'GCSH', 'NFS1', 'ISCU', 'COX11',
    'SCO1', 'COA6', 'COX10', 'SURF1', 'LYRM4', 'ACSF2'
}

STRING_COPPER_GENES = {'FDX1', 'LIPT1', 'ATP7A', 'ATP7B', 'DLAT', 'DLD', 'PDHB', 'SLC31A1'}

def load_string_network(interaction_file):
    """Load STRING interaction network"""
    print("Loading STRING interaction data...")
    df = pd.read_csv(interaction_file, sep='\t')
    df.columns = [col.lstrip('#') for col in df.columns]
    print(f"  Total interactions: {len(df)}")
    return df

def build_graph(interactions_df):
    """Build NetworkX graph from STRING interactions"""
    print("\nBuilding network graph...")
    G = nx.Graph()

    for _, row in interactions_df.iterrows():
        node1 = row['node1']
        node2 = row['node2']
        combined_score = row['combined_score']

        if G.has_edge(node1, node2):
            G[node1][node2]['weight'] = max(G[node1][node2]['weight'], combined_score)
        else:
            G.add_edge(node1, node2, weight=combined_score)

    print(f"  Nodes: {G.number_of_nodes()}")
    print(f"  Edges: {G.number_of_edges()}")
    return G

def calculate_cytoNCA(G):
    """Calculate 8 cytoNCA centrality metrics"""
    print("\nCalculating cytoNCA metrics...")

    results = []

    for node in G.nodes():
        neighbors = list(G.neighbors(node))
        degree = G.degree(node)

        bc = nx.betweenness_centrality(G, normalized=True).get(node, 0)
        cc = nx.clustering(G, nodes=[node]).get(node, 0)

        ego_subgraph = G.subgraph(neighbors + [node])
        dc = ego_subgraph.number_of_edges()

        ec = 0
        if len(neighbors) > 1:
            subgraph = G.subgraph(neighbors)
            ec = subgraph.number_of_edges() / (len(neighbors) * (len(neighbors) - 1) / 2) if len(neighbors) > 1 else 0

        lac = 0
        if degree > 0:
            lac = sum(G[u][v]['weight'] for u in neighbors for v in neighbors if G.has_edge(u, v)) / degree

        nc = degree

        sc = sum(1 / (nx.shortest_path_length(G, node, n) + 1) for n in G.nodes() if n != node)

        ic = 0
        if degree > 0:
            ic = sum(G[node][n]['weight'] for n in neighbors) / degree

        results.append({
            'Gene': node,
            'BC': bc,
            'CC': cc,
            'DC': dc,
            'EC': ec,
            'LAC': lac,
            'NC': nc,
            'SC': sc,
            'IC': ic,
            'Degree': degree
        })

    df = pd.DataFrame(results)

    for col in ['BC', 'CC', 'DC', 'EC', 'LAC', 'NC', 'SC', 'IC']:
        max_val = df[col].max()
        if max_val > 0:
            df[col] = df[col] / max_val

    print(f"  Calculated cytoNCA metrics for {len(df)} genes")
    return df

def calculate_mcc(G):
    """Calculate Maximal Clique Centrality (MCC) for each node"""
    print("\nCalculating MCC scores (cytoHubba algorithm)...")

    mcc_scores = {}
    nodes = list(G.nodes())

    for i, node in enumerate(nodes):
        if (i + 1) % 50 == 0:
            print(f"  Processing node {i+1}/{len(nodes)}...")

        degree = G.degree(node)
        max_clique_size = 0
        neighbors = list(G.neighbors(node))

        if degree < 2:
            mcc_scores[node] = degree
            continue

        for neighbor in neighbors:
            neighbor_set = set([neighbor])
            for other in neighbors:
                if other != neighbor and G.has_edge(neighbor, other):
                    neighbor_set.add(other)
            max_clique_size = max(max_clique_size, len(neighbor_set))

        mcc_scores[node] = degree * np.sqrt(max_clique_size) if max_clique_size > 0 else degree

    mcc_df = pd.DataFrame([
        {'Gene': gene, 'MCC': score}
        for gene, score in mcc_scores.items()
    ])

    max_mcc = mcc_df['MCC'].max()
    if max_mcc > 0:
        mcc_df['MCC'] = mcc_df['MCC'] / max_mcc

    return mcc_df

def combine_scores(cytoNCA_df, mcc_df):
    """Combine MCC and cytoNCA scores for final ranking"""
    print("\nCombining MCC and cytoNCA scores...")

    df = cytoNCA_df.merge(mcc_df, on='Gene')

    df['cytoNCA_Score'] = df[['BC', 'CC', 'DC', 'EC', 'LAC', 'NC', 'SC', 'IC']].mean(axis=1)

    df['Combined_Score'] = df['cytoNCA_Score'] * 0.5 + df['MCC'] * 0.5

    df = df.sort_values('Combined_Score', ascending=False)

    print(f"\nTop 30 Hub Genes by Combined Score:")
    top30 = df.head(30)
    for idx, row in enumerate(top30.itertuples(), 1):
        print(f"  {idx:2d}. {row.Gene}: Combined={row.Combined_Score:.4f}, MCC={row.MCC:.4f}, cytoNCA={row.cytoNCA_Score:.4f}")

    return df

def group_hub_genes(G, top30_genes, copper_death_genes, string_copper_genes):
    """Group Hub genes into Group A and Group B"""
    print("\n" + "="*60)
    print("Grouping Hub Genes")
    print("="*60)

    copper_hub_genes = set(top30_genes) & copper_death_genes

    group_a = copper_hub_genes.copy()
    print(f"\n组A - 铜死亡Hub基因 (Hub ∩ Copper Death):")
    if group_a:
        print(f"  共 {len(group_a)} 个基因:")
        for gene in sorted(group_a):
            print(f"    • {gene}")
    else:
        print(f"  共 0 个基因 (铜死亡基因不在Top30 Hub中)")

    regulatory_hub_genes = set()
    for hub_gene in top30_genes:
        neighbors = set(G.neighbors(hub_gene))
        if neighbors & string_copper_genes:
            regulatory_hub_genes.add(hub_gene)

    group_b = regulatory_hub_genes.copy()
    print(f"\n组B - 调控型Hub基因 (与铜死亡Hub有直接互作):")
    if group_b:
        print(f"  共 {len(group_b)} 个基因:")
        for gene in sorted(group_b):
            neighbors_with_cu = set(G.neighbors(gene)) & string_copper_genes
            print(f"    • {gene} ← 互作 → {neighbors_with_cu}")
    else:
        print(f"  共 0 个基因")

    core_mechanism_genes = group_a | group_b
    if not core_mechanism_genes:
        copper_related = set()
        for hub_gene in top30_genes:
            neighbors = set(G.neighbors(hub_gene))
            if neighbors & string_copper_genes:
                copper_related.add(hub_gene)
        core_mechanism_genes = copper_related
        print(f"\n修正: 使用铜死亡相关Hub基因作为核心机制基因 ({len(core_mechanism_genes)} 个)")

    print(f"\n核心机制基因集 (组A + 组B):")
    print(f"  共 {len(core_mechanism_genes)} 个基因: {sorted(core_mechanism_genes)}")

    return group_a, group_b, core_mechanism_genes

def save_results(df, group_a, group_b, core_genes, copper_death_genes, output_prefix="Hub_Gene_MCC_cytoNCA"):
    """Save analysis results to Excel"""
    with pd.ExcelWriter(f'{output_prefix}_Results.xlsx') as writer:
        df.to_excel(writer, sheet_name='All_Genes_cytoNCA_MCC', index=False)

        top30_df = df.head(30).copy()
        top30_df['Group'] = 'Other'
        for gene in group_a:
            top30_df.loc[top30_df['Gene'] == gene, 'Group'] = 'A_Copper_Death_Hub'
        for gene in group_b:
            top30_df.loc[top30_df['Gene'] == gene, 'Group'] = 'B_Regulatory_Hub'
        top30_df.to_excel(writer, sheet_name='Top30_Hub_Genes', index=False)

        pd.DataFrame({'Group_A_Copper_Death_Hub': sorted(group_a)}).to_excel(
            writer, sheet_name='Group_A', index=False)

        pd.DataFrame({'Group_B_Regulatory_Hub': sorted(group_b)}).to_excel(
            writer, sheet_name='Group_B', index=False)

        pd.DataFrame({'Core_Mechanism_Genes': sorted(core_genes)}).to_excel(
            writer, sheet_name='Core_Genes', index=False)

        pd.DataFrame({'Copper_Death_Gene_Set': sorted(copper_death_genes)}).to_excel(
            writer, sheet_name='Copper_Death_Genes', index=False)

        pd.DataFrame({'STRING_Copper_Genes': sorted(STRING_COPPER_GENES)}).to_excel(
            writer, sheet_name='STRING_Copper', index=False)

        summary_data = {
            'Category': ['Total Genes', 'Top 30 Hub Genes', 'Group A (Copper Death Hub)',
                        'Group B (Regulatory Hub)', 'Core Mechanism Genes'],
            'Count': [len(df), len(top30_df), len(group_a), len(group_b), len(core_genes)]
        }
        pd.DataFrame(summary_data).to_excel(writer, sheet_name='Summary', index=False)

    print(f"\nResults saved to: {output_prefix}_Results.xlsx")

def main():
    print("="*70)
    print("Hub Gene Analysis: MCC + cytoNCA Combined Algorithm")
    print("="*70)

    interaction_file = r"C:\Users\Jy-Mentor-7\Downloads\string_interactions (2).tsv"

    print(f"\nCopper Death Gene Set: {sorted(COPPER_DEATH_GENES)}")
    print(f"STRING Copper Genes: {sorted(STRING_COPPER_GENES)}")

    interactions_df = load_string_network(interaction_file)

    G = build_graph(interactions_df)

    cytoNCA_df = calculate_cytoNCA(G)

    mcc_df = calculate_mcc(G)

    combined_df = combine_scores(cytoNCA_df, mcc_df)

    top30_genes = combined_df.head(30)['Gene'].tolist()

    group_a, group_b, core_genes = group_hub_genes(G, top30_genes, COPPER_DEATH_GENES, STRING_COPPER_GENES)

    save_results(combined_df, group_a, group_b, core_genes, COPPER_DEATH_GENES)

    print("\n" + "="*70)
    print("Analysis Complete!")
    print("="*70)

if __name__ == "__main__":
    main()
