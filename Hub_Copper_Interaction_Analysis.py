#!/usr/bin/env python3
"""
Hub Gene & Copper Death Gene Interaction Analysis
找出Hub基因与铜死亡基因之间的直接连接
"""

import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore')

def smart_install(package):
    import subprocess
    import sys
    try:
        __import__(package)
    except ImportError:
        subprocess.check_call([sys.executable, "-m", "pip", "install", package, "-q"])

for pkg in ["openpyxl"]:
    smart_install(pkg)

COPPER_DEATH_GENES = {
    'FDX1', 'DLAT', 'DLD', 'LIPT1', 'PDHX', 'PDHB', 'SLC31A1',
    'ATP7B', 'ATP7A', 'ATOX1', 'COMMD1', 'MT2A', 'NFKB1', 'ATF4', 'TLR4'
}

HUB_GENES_29 = {
    'IL6', 'STAT3', 'PPARG', 'TGFB1', 'CCL2', 'PTGS2',
    'CCND1', 'STAT1', 'ICAM1', 'PTPRC', 'RELA', 'CASP8', 'CXCR4',
    'NOTCH1', 'MAPK1', 'MDM2', 'HSPA5', 'PARP1', 'JAK1', 'CREBBP',
    'MMP2', 'SREBF1', 'CDC42', 'STAT5A', 'NFE2L2', 'HMOX1', 'IGF1R'
}

def main():
    print("="*70)
    print("Hub Gene & Copper Death Gene Interaction Analysis")
    print("="*70)

    print(f"\nHub Genes (29): {sorted(HUB_GENES_29)}")
    print(f"\nCopper Death Genes: {sorted(COPPER_DEATH_GENES)}")

    interaction_file = r"C:\Users\Jy-Mentor-7\Downloads\string_interactions (3).tsv"

    print("\nLoading STRING interaction data...")
    df = pd.read_csv(interaction_file, sep='\t')
    df.columns = [col.lstrip('#') for col in df.columns]
    print(f"  Total interactions: {len(df)}")

    all_genes = set(df['node1'].unique()) | set(df['node2'].unique())
    print(f"  Total unique genes: {len(all_genes)}")

    hub_in_network = HUB_GENES_29 & all_genes
    copper_in_network = COPPER_DEATH_GENES & all_genes

    print(f"\nHub genes in network: {len(hub_in_network)}")
    print(f"Copper death genes in network: {len(copper_in_network)}")

    print(f"\nHub genes NOT in network: {HUB_GENES_29 - all_genes}")
    print(f"Copper genes NOT in network: {COPPER_DEATH_GENES - all_genes}")

    print("\n" + "="*70)
    print("Filtering: Hub genes <-> Copper Death genes interactions")
    print("="*70)

    hub_copper_edges = df[
        ((df['node1'].isin(hub_in_network)) & (df['node2'].isin(copper_in_network))) |
        ((df['node1'].isin(copper_in_network)) & (df['node2'].isin(hub_in_network)))
    ].copy()

    print(f"\nFound {len(hub_copper_edges)} direct interactions")

    if len(hub_copper_edges) > 0:
        hub_copper_edges['Hub_Gene'] = hub_copper_edges.apply(
            lambda x: x['node1'] if x['node1'] in hub_in_network else x['node2'], axis=1
        )
        hub_copper_edges['Copper_Gene'] = hub_copper_edges.apply(
            lambda x: x['node2'] if x['node1'] in hub_in_network else x['node1'], axis=1
        )

        hub_copper_edges = hub_copper_edges[['Hub_Gene', 'Copper_Gene', 'combined_score',
                                             'experimentally_determined_interaction',
                                             'database_annotated', 'automated_textmining']]

        hub_copper_edges = hub_copper_edges.sort_values(['Hub_Gene', 'combined_score'], ascending=[True, False])

        print("\nHub-Copper Death Gene Interactions:")
        print("-" * 50)
        for hub in sorted(hub_copper_edges['Hub_Gene'].unique()):
            copper_list = hub_copper_edges[hub_copper_edges['Hub_Gene'] == hub]['Copper_Gene'].tolist()
            print(f"  {hub} → {copper_list}")

        output_file = "Hub_Copper_Interactions.xlsx"
        with pd.ExcelWriter(output_file) as writer:
            hub_copper_edges.to_excel(writer, sheet_name='Hub_Copper_Edges', index=False)

            pd.DataFrame({'Hub_Genes': sorted(hub_in_network)}).to_excel(
                writer, sheet_name='Hub_Genes_in_Network', index=False)

            pd.DataFrame({'Copper_Genes': sorted(copper_in_network)}).to_excel(
                writer, sheet_name='Copper_Genes_in_Network', index=False)

        print(f"\nResults saved to: {output_file}")

    else:
        print("\nNo direct interactions found between Hub genes and Copper Death genes.")

    print("\n" + "="*70)
    print("Analysis Complete!")
    print("="*70)

if __name__ == "__main__":
    main()
