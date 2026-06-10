#!/usr/bin/env python3
"""检查铜死亡基因在网络中的存在情况"""

import pandas as pd
import numpy as np
import networkx as nx
import warnings
warnings.filterwarnings('ignore')

COPPER_DEATH_GENES = {
    'FDX1', 'LIPT1', 'LIPT2', 'MTF1', 'SLC31A1', 'SLC31A2',
    'ATP7A', 'ATP7B', 'DLD', 'DLAT', 'PDHA1', 'PDHB',
    'LIAS', 'LIPTC1', 'GCSH', 'NFS1', 'ISCU', 'COX11',
    'SCO1', 'COA6', 'COX10', 'SURF1', 'LYRM4', 'ACSF2'
}

print("Loading STRING interaction data...")
df = pd.read_csv(r"C:\Users\Jy-Mentor-7\Downloads\string_interactions (2).tsv", sep='\t')
df.columns = [col.lstrip('#') for col in df.columns]

print(f"Total interactions: {len(df)}")

all_genes = set(df['node1'].unique()) | set(df['node2'].unique())
print(f"Total unique genes in network: {len(all_genes)}")

print("\n" + "="*60)
print("Checking Copper Death Genes in Network")
print("="*60)

found_copper_genes = all_genes & COPPER_DEATH_GENES
not_found_copper_genes = COPPER_DEATH_GENES - all_genes

print(f"\nFound in network ({len(found_copper_genes)}):")
for gene in sorted(found_copper_genes):
    neighbors = set(df[(df['node1']==gene) | (df['node2']==gene)]['node1'].tolist() +
                   df[(df['node1']==gene) | (df['node2']==gene)]['node2'].tolist())
    neighbors.discard(gene)
    print(f"  • {gene}: {len(neighbors)} neighbors")

print(f"\nNOT found in network ({len(not_found_copper_genes)}):")
for gene in sorted(not_found_copper_genes):
    print(f"  ✗ {gene}")

hub_genes = {'IL6', 'STAT3', 'NFKB1', 'TGFB1', 'PPARG', 'CCL2', 'TLR4', 'PTGS2',
             'ICAM1', 'PTPRC', 'CCND1', 'STAT1', 'MAPK1', 'CXCR4', 'CASP8',
             'NOTCH1', 'MDM2', 'RELA', 'HSPA5', 'CREBBP', 'PARP1', 'SREBF1',
             'MMP2', 'CDC42', 'HMOX1', 'STAT5A', 'JAK1', 'IRF1', 'NFE2L2', 'FASN'}

print("\n" + "="*60)
print("Hub Genes in Network")
print("="*60)
print(f"Hub genes found: {len(hub_genes & all_genes)}")

print("\nHub genes that interact with Copper Death genes:")
for hub in sorted(hub_genes & all_genes):
    neighbors = set(df[(df['node1']==hub) | (df['node2']==hub)]['node1'].tolist() +
                   df[(df['node1']==hub) | (df['node2']==hub)]['node2'].tolist())
    neighbors.discard(hub)
    cu_neighbors = neighbors & COPPER_DEATH_GENES
    if cu_neighbors:
        print(f"  • {hub} → {cu_neighbors}")
