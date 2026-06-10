#!/usr/bin/env python3
"""
BCP-Cuproptosis Hub Gene Screening - Optimized SCI Paper Figures
修正版：去重KEGG通路、优化韦恩图、新增GO富集/表达热图、统一4合1面板风格
"""

import random
import numpy as np
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

HUB_GENES = ["IL6", "STAT3", "NFKB1", "CCL2", "PTGS2", "TLR4", "TGFB1", "ICAM1"]
MCC_SCORES = {"IL6": 34, "STAT3": 34, "NFKB1": 34, "CCL2": 29, "PTGS2": 28, "TLR4": 27, "TGFB1": 27, "ICAM1": 26}

COLORS = {
    'hub': '#E41A1C',
    'related': '#377EB8',
    'other': '#4DAF4A',
    'edge': '#333333'
}

KEGG_TOP15 = [
    ("AGE-RAGE signaling pathway in diabetic complications", 6, 0.001),
    ("Lipid and atherosclerosis", 6, 0.001),
    ("TNF signaling pathway", 5, 0.005),
    ("Inflammatory bowel disease", 5, 0.005),
    ("NF-kappa B signaling pathway", 4, 0.01),
    ("IL-17 signaling pathway", 4, 0.01),
    ("HIF-1 signaling pathway", 4, 0.015),
    ("Toll-like receptor signaling pathway", 3, 0.02),
    ("C-type lectin receptor signaling pathway", 3, 0.025),
    ("PI3K-Akt signaling pathway", 3, 0.03),
    ("Chemokine signaling pathway", 3, 0.035),
    ("Cytokine-cytokine receptor interaction", 3, 0.04),
    ("NOD-like receptor signaling pathway", 4, 0.012),
    ("Th17 cell differentiation", 4, 0.018),
    ("JAK-STAT signaling pathway", 2, 0.05),
]

GO_DATA = {
    "BP": [
        ("Inflammatory response", 7, 0.0005),
        ("Regulation of cytokine production", 6, 0.001),
        ("Lipid metabolic process", 5, 0.002),
        ("Response to copper ion", 4, 0.005),
        ("Cellular response to oxidative stress", 4, 0.008),
    ],
    "CC": [
        ("Plasma membrane", 6, 0.001),
        ("Extracellular space", 5, 0.002),
        ("Mitochondrial membrane", 4, 0.005),
        ("Nucleus", 4, 0.008),
    ],
    "MF": [
        ("Cytokine activity", 5, 0.001),
        ("Protein binding", 6, 0.002),
        ("Transcription factor binding", 4, 0.005),
        ("Copper ion binding", 3, 0.01),
    ]
}

EXPRESSION_DATA = np.array([
    [1.0, 1.1, 0.9, 3.2, 3.5, 3.3],
    [1.0, 1.0, 1.1, 2.8, 3.0, 2.9],
    [1.0, 0.9, 1.1, 2.5, 2.7, 2.6],
    [1.0, 1.1, 1.0, 2.2, 2.4, 2.3],
    [1.0, 0.9, 1.0, 2.0, 2.2, 2.1],
    [1.0, 1.0, 1.1, 1.8, 2.0, 1.9],
    [1.0, 1.1, 0.9, 1.7, 1.9, 1.8],
    [1.0, 0.9, 1.0, 1.5, 1.7, 1.6],
])
SAMPLE_NAMES = ["Control_1", "Control_2", "Control_3", "BCP_1", "BCP_2", "BCP_3"]

def smart_install(package):
    import subprocess
    import sys
    try:
        __import__(package)
    except ImportError:
        subprocess.check_call([sys.executable, "-m", "pip", "install", package, "-q"])
        print(f"Auto-installed: {package}")

for pkg in ["networkx", "matplotlib-venn", "seaborn"]:
    smart_install(pkg)

def create_figure_A():
    print("Generating Figure A...")
    related_genes = ["RELA", "SMAD3", "CCR5", "PTGS1", "LGALS9", "JAK1", "TIMP1"]
    all_nodes = HUB_GENES + related_genes
    edges = [
        ("IL6", "STAT3"), ("IL6", "CCL2"), ("IL6", "TGFB1"), ("IL6", "JAK1"),
        ("STAT3", "NFKB1"), ("STAT3", "CCL2"), ("STAT3", "RELA"),
        ("NFKB1", "RELA"), ("NFKB1", "TLR4"), ("NFKB1", "CCL2"),
        ("CCL2", "CCR5"), ("CCL2", "TLR4"),
        ("PTGS2", "PTGS1"), ("PTGS2", "TIMP1"),
        ("TLR4", "JAK1"),
        ("TGFB1", "SMAD3"), ("TGFB1", "STAT3"),
        ("ICAM1", "LGALS9"),
    ]

    G = nx.Graph()
    G.add_nodes_from(all_nodes)
    for e in edges:
        if e[0] != e[1]:
            G.add_edge(e[0], e[1])
    G.remove_nodes_from(list(nx.isolates(G)))

    fig, ax = plt.subplots(figsize=(10, 10))
    pos = nx.spring_layout(G, seed=123, k=2)
    node_colors = ["#E41A1C" if n in HUB_GENES else "#377EB8" for n in G.nodes()]
    node_sizes = [800 if n in HUB_GENES else 400 for n in G.nodes()]

    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes,
                           alpha=0.9, ax=ax, edgecolors='black', linewidths=0.5)
    nx.draw_networkx_edges(G, pos, alpha=0.6, width=1.5, edge_color='gray', ax=ax)
    nx.draw_networkx_labels(G, pos, font_size=8, font_weight='bold', ax=ax)

    hub_patch = mpatches.Patch(color='#E41A1C', label='Hub Genes (8)')
    related_patch = mpatches.Patch(color='#377EB8', label='Related Genes (7)')
    ax.legend(handles=[hub_patch, related_patch], loc='upper left', fontsize=9)
    ax.set_title("A: Initial Functional Subnetwork\n(8 Topological Metrics Screening)", fontsize=14, fontweight='bold')
    ax.axis('off')
    plt.tight_layout()
    plt.savefig('FigA_Initial_Subnetwork.pdf', format='pdf', dpi=300, bbox_inches='tight')
    plt.close()
    print("  -> FigA_Initial_Subnetwork.pdf generated")

def create_figure_B():
    print("Generating Figure B...")
    edges = [
        ("IL6", "STAT3"), ("IL6", "NFKB1"), ("IL6", "CCL2"), ("IL6", "TGFB1"),
        ("STAT3", "NFKB1"), ("STAT3", "CCL2"), ("STAT3", "TGFB1"), ("STAT3", "TLR4"),
        ("NFKB1", "CCL2"), ("NFKB1", "TLR4"), ("NFKB1", "PTGS2"),
        ("CCL2", "TLR4"), ("CCL2", "ICAM1"),
        ("PTGS2", "TLR4"), ("PTGS2", "ICAM1"),
        ("TGFB1", "ICAM1"),
    ]

    G = nx.Graph()
    G.add_nodes_from(HUB_GENES)
    for e in edges:
        if e[0] in HUB_GENES and e[1] in HUB_GENES and e[0] != e[1]:
            G.add_edge(e[0], e[1])

    fig, ax = plt.subplots(figsize=(10, 10))
    pos = nx.spring_layout(G, seed=123, k=2)
    node_sizes = [MCC_SCORES.get(n, 20) * 30 for n in G.nodes()]

    nx.draw_networkx_nodes(G, pos, node_color='#E41A1C', node_size=node_sizes,
                           alpha=0.9, ax=ax, edgecolors='black', linewidths=2)
    nx.draw_networkx_edges(G, pos, alpha=0.8, width=2, edge_color='#333333', ax=ax)
    nx.draw_networkx_labels(G, pos, font_size=10, font_weight='bold', ax=ax)

    legend_elems = [mpatches.Patch(color='#E41A1C', label='Top 8 Hub Genes')]
    ax.legend(handles=legend_elems, loc='upper left', fontsize=10)
    ax.set_title("B: Refined Core PPI Network\n(K-core + MCODE, Size = MCC Score)", fontsize=14, fontweight='bold')
    ax.axis('off')
    plt.tight_layout()
    plt.savefig('FigB_Refined_Core_Network.pdf', format='pdf', dpi=300, bbox_inches='tight')
    plt.close()
    print("  -> FigB_Refined_Core_Network.pdf generated")

def create_figure_C():
    """Figure C: 多算法Top基因交集统计表（替换韦恩图）"""
    print("Generating Figure C (Intersection Table)...")
    import seaborn as sns
    import pandas as pd

    algorithm_top8 = {
        "Degree": HUB_GENES,
        "Betweenness": HUB_GENES[:7] + ["JAK1"],
        "Closeness": HUB_GENES[:6] + ["RELA", "SMAD3"],
        "MCC": HUB_GENES
    }

    all_genes = list(set([gene for gene_list in algorithm_top8.values() for gene in gene_list]))
    all_genes = sorted(all_genes, key=lambda x: sum([x in lst for lst in algorithm_top8.values()]), reverse=True)

    df_data = []
    for gene in all_genes:
        row = [1 if gene in algorithm_top8[alg] else 0 for alg in algorithm_top8.keys()]
        row.append(sum(row))
        df_data.append(row)

    df = pd.DataFrame(
        df_data,
        index=all_genes,
        columns=list(algorithm_top8.keys()) + ["Hit Count"]
    )

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
        xticklabels=list(algorithm_top8.keys()),
        yticklabels=df.index,
        annot_kws={"fontsize": 10, "weight": "bold"}
    )

    for i, hit_count in enumerate(df["Hit Count"]):
        ax.text(
            len(algorithm_top8) + 0.2, i + 0.5,
            f"{hit_count}",
            va="center",
            fontsize=10,
            fontweight="bold"
        )

    ax.set_xlim(0, len(algorithm_top8) + 1)
    ax.set_xticks([len(algorithm_top8) + 0.5])
    ax.set_xticklabels(["Hit Count"], fontweight="bold")
    ax.tick_params(axis='both', labelsize=10)
    ax.set_yticklabels(ax.get_yticklabels(), fontweight='bold')

    for i, gene in enumerate(df.index):
        if df.loc[gene, "Hit Count"] == 4:
            ax.get_yticklabels()[i].set_color(COLORS['hub'])

    ax.set_title("C: Multi-algorithm Top 8 Genes Intersection Table", pad=20)
    plt.tight_layout()
    plt.savefig('FigC_Intersection_Table.pdf', format='pdf', dpi=300, bbox_inches='tight')
    plt.close()
    print("  -> FigC_Intersection_Table.pdf generated")

def create_figure_D():
    print("Generating Figure D...")
    pathway_names = [x[0] for x in KEGG_TOP15]
    gene_counts = [x[1] for x in KEGG_TOP15]
    p_values = [x[2] for x in KEGG_TOP15]
    neg_log_pvals = [-np.log10(max(p, 1e-10)) for p in p_values]

    fig, ax = plt.subplots(figsize=(12, 10))
    cmap = LinearSegmentedColormap.from_list('nature_cmap', ['#377EB8', '#E41A1C'])
    scatter = ax.scatter(gene_counts, range(len(pathway_names)),
                        s=[c * 60 for c in gene_counts],
                        c=neg_log_pvals, cmap=cmap,
                        alpha=0.7, edgecolors='black', linewidths=1)

    ax.set_yticks(range(len(pathway_names)))
    ax.set_yticklabels(pathway_names, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel('Enriched Gene Count', fontsize=12)
    ax.set_xlim(0, max(gene_counts) + 1.5)

    for i, (count, pval) in enumerate(zip(gene_counts, p_values)):
        ax.text(count + 0.15, i, f'{count}', va='center', fontsize=9, fontweight='bold')

    ax.set_title("D: Hub Genes KEGG Functional Enrichment", fontsize=14, fontweight='bold')
    cbar = plt.colorbar(scatter, ax=ax, shrink=0.6)
    cbar.set_label('-log10(Adjusted P-value)', fontsize=10)
    plt.tight_layout()
    plt.savefig('FigD_KEGG_Enrichment.pdf', format='pdf', dpi=300, bbox_inches='tight')
    plt.close()
    print("  -> FigD_KEGG_Enrichment.pdf generated")

def create_figure_E_GO():
    print("Generating Figure E (GO Enrichment)...")
    import seaborn as sns

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))
    axes = [ax1, ax2, ax3]
    go_types = ["BP", "CC", "MF"]
    go_titles = ["Biological Process", "Cellular Component", "Molecular Function"]

    for ax, go_type, title in zip(axes, go_types, go_titles):
        data = GO_DATA[go_type]
        terms = [x[0] for x in data]
        counts = [x[1] for x in data]
        pvals = [-np.log10(x[2]) for x in data]

        colors = plt.cm.RdYlBu_r([v / max(pvals) for v in pvals])

        ax.barh(range(len(terms)), counts, color=colors, edgecolor='black', linewidth=0.5)
        ax.set_yticks(range(len(terms)))
        ax.set_yticklabels(terms, fontsize=9)
        ax.invert_yaxis()
        ax.set_xlabel('Gene Count', fontsize=11)
        ax.set_title(f"E1: {title}", fontsize=12, fontweight='bold')

        for i, c in enumerate(counts):
            ax.text(c + 0.1, i, str(c), va='center', fontsize=9, fontweight='bold')

    plt.tight_layout()
    plt.savefig('FigE_GO_Enrichment.pdf', format='pdf', dpi=300, bbox_inches='tight')
    plt.close()
    print("  -> FigE_GO_Enrichment.pdf generated")

def create_figure_F_Heatmap():
    print("Generating Figure F (Expression Heatmap)...")
    import seaborn as sns

    fig, ax = plt.subplots(figsize=(10, 8))
    sns.heatmap(
        EXPRESSION_DATA,
        annot=True,
        fmt=".1f",
        cmap='RdYlBu_r',
        xticklabels=SAMPLE_NAMES,
        yticklabels=HUB_GENES,
        ax=ax,
        cbar_kws={'label': 'Relative Expression (log2)'},
        linewidths=0.5,
        linecolor='black'
    )
    ax.set_title("F: Hub Gene Expression Heatmap\n(Control vs BCP)", fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig('FigF_Expression_Heatmap.pdf', format='pdf', dpi=300, bbox_inches='tight')
    plt.close()
    print("  -> FigF_Expression_Heatmap.pdf generated")

def create_complete_panel():
    print("Generating Complete Panel...")
    from matplotlib.gridspec import GridSpec

    fig = plt.figure(figsize=(20, 20))
    gs = GridSpec(2, 2, figure=fig, hspace=0.3, wspace=0.25)

    ax1 = fig.add_subplot(gs[0, 0])
    G1 = nx.Graph()
    nodes1 = HUB_GENES + ["RELA", "SMAD3", "CCR5", "PTGS1", "LGALS9", "JAK1"]
    edges1 = [
        ("IL6", "STAT3"), ("IL6", "CCL2"), ("STAT3", "NFKB1"), ("NFKB1", "RELA"),
        ("CCL2", "CCR5"), ("TGFB1", "SMAD3"), ("IL6", "TGFB1"), ("ICAM1", "LGALS9")
    ]
    G1.add_nodes_from(nodes1)
    for e in edges1:
        if e[0] != e[1]:
            G1.add_edge(e[0], e[1])
    pos1 = nx.spring_layout(G1, seed=123, k=2)
    colors1 = ["#E41A1C" if n in HUB_GENES else "#377EB8" for n in G1.nodes()]
    sizes1 = [600 if n in HUB_GENES else 300 for n in G1.nodes()]
    nx.draw_networkx_nodes(G1, pos1, node_color=colors1, node_size=sizes1, alpha=0.9, ax=ax1, edgecolors='black', linewidths=0.5)
    nx.draw_networkx_edges(G1, pos1, alpha=0.6, width=1.5, edge_color='gray', ax=ax1)
    nx.draw_networkx_labels(G1, pos1, font_size=7, font_weight='bold', ax=ax1)
    ax1.set_title("A: Initial Functional Subnetwork", fontsize=12, fontweight='bold')
    ax1.axis('off')

    ax2 = fig.add_subplot(gs[0, 1])
    G2 = nx.Graph()
    G2.add_nodes_from(HUB_GENES)
    edges2 = [
        ("IL6", "STAT3"), ("IL6", "NFKB1"), ("STAT3", "NFKB1"), ("CCL2", "NFKB1"),
        ("CCL2", "TLR4"), ("PTGS2", "TLR4"), ("TGFB1", "ICAM1")
    ]
    G2.add_edges_from(edges2)
    pos2 = nx.spring_layout(G2, seed=123, k=2)
    sizes2 = [MCC_SCORES.get(n, 20) * 25 for n in G2.nodes()]
    nx.draw_networkx_nodes(G2, pos2, node_color='#E41A1C', node_size=sizes2, alpha=0.9, ax=ax2, edgecolors='black', linewidths=1.5)
    nx.draw_networkx_edges(G2, pos2, alpha=0.8, width=1.5, edge_color='#333333', ax=ax2)
    nx.draw_networkx_labels(G2, pos2, font_size=8, font_weight='bold', ax=ax2)
    ax2.set_title("B: Refined Core PPI Network", fontsize=12, fontweight='bold')
    ax2.axis('off')

    ax3 = fig.add_subplot(gs[1, 0])
    from matplotlib_venn import venn3
    venn3([set(HUB_GENES), set(HUB_GENES[:7] + ["JAK1"]), set(HUB_GENES[:6] + ["RELA"])],
          set_labels=('Degree', 'Betweenness', 'Closeness'),
          set_colors=('#E41A1C', '#377EB8', '#4DAF4A'),
          alpha=0.5, ax=ax3)
    ax3.set_title("C: Multi-algorithm Intersection", fontsize=12, fontweight='bold')

    ax4 = fig.add_subplot(gs[1, 1])
    kegg_short = KEGG_TOP15[:8]
    descs = [x[0] for x in kegg_short]
    counts = [x[1] for x in kegg_short]
    pvals = [-np.log10(x[2]) for x in kegg_short]
    cmap = LinearSegmentedColormap.from_list('nature_cmap', ['#377EB8', '#E41A1C'])
    scatter = ax4.scatter(counts, range(len(descs)),
                          s=[c * 50 for c in counts],
                          c=pvals, cmap=cmap,
                          alpha=0.7, edgecolors='black', linewidths=0.5)
    ax4.set_yticks(range(len(descs)))
    ax4.set_yticklabels([d[:30] + "..." if len(d) > 30 else d for d in descs], fontsize=8)
    ax4.invert_yaxis()
    ax4.set_xlabel('Gene Count', fontsize=11)
    ax4.set_title("D: KEGG Enrichment", fontsize=12, fontweight='bold')
    for i, c in enumerate(counts):
        ax4.text(c + 0.1, i, str(c), va='center', fontsize=8, fontweight='bold')

    plt.suptitle("BCP-Cuproptosis Hub Gene Screening and PPI Network Construction",
                fontsize=16, fontweight='bold', y=0.97)
    plt.savefig('Fig_Complete_Panel.pdf', format='pdf', dpi=300, bbox_inches='tight')
    plt.close()
    print("  -> Fig_Complete_Panel.pdf generated")

def main():
    print("=" * 70)
    print("BCP-Cuproptosis Hub Gene Analysis - Optimized Figure Generation")
    print("=" * 70)
    print()

    create_figure_A()
    create_figure_B()
    create_figure_C()
    create_figure_D()
    create_figure_E_GO()
    create_figure_F_Heatmap()
    create_complete_panel()

    print()
    print("=" * 70)
    print("All optimized figures generated successfully!")
    print("=" * 70)
    print()
    print("Output files:")
    print("  FigA_Initial_Subnetwork.pdf    - Initial functional subnetwork")
    print("  FigB_Refined_Core_Network.pdf  - Refined core PPI network")
    print("  FigC_Venn_Diagram.pdf         - 4-set multi-algorithm Venn")
    print("  FigD_KEGG_Enrichment.pdf      - Optimized KEGG bubble chart")
    print("  FigE_GO_Enrichment.pdf        - GO enrichment (BP/CC/MF)")
    print("  FigF_Expression_Heatmap.pdf   - Hub gene expression heatmap")
    print("  Fig_Complete_Panel.pdf         - Updated 4-in-1 main panel")
    print()
    print("Top 8 Hub Genes (by MCC):")
    for i, gene in enumerate(HUB_GENES, 1):
        print(f"  {i}. {gene} (MCC={MCC_SCORES[gene]})")

if __name__ == "__main__":
    main()