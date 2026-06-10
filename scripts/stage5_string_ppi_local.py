# -*- coding: utf-8 -*-
"""
Stage 5: STRING PPI Network from Local Files (v9.0)
====================================================

Uses locally downloaded STRING files instead of API:
- string_interactions_short.tsv: PPI edges with combined_score
- string_node_degrees.tsv: Node degree information
- enrichment.all.tsv: GO enrichment results

Input:
  - Local STRING files
  - stage4/seed_pool_genes.txt
  
Output:
  - string_ppi.tsv: PPI network edges
  - ppi_topology.json: Network topology metrics
"""

import os
import sys
import json
import logging
from datetime import datetime
from typing import Dict, Set, Optional

import numpy as np
import pandas as pd
import networkx as nx

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import RESULTS_DIR

STAGE_DIR = os.path.join(RESULTS_DIR, "stage5_string_ppi")
os.makedirs(STAGE_DIR, exist_ok=True)

logger = logging.getLogger("stage5_local")
logger.setLevel(logging.INFO)

file_handler = logging.FileHandler(os.path.join(STAGE_DIR, "stage5_local.log"), encoding="utf-8", mode="w")
file_handler.setFormatter(logging.Formatter("%(asctime)s | %(levelname)s | %(message)s"))
logger.addHandler(file_handler)

console_handler = logging.StreamHandler(sys.stdout)
console_handler.setFormatter(logging.Formatter("%(message)s"))
logger.addHandler(console_handler)

# Local file paths
LOCAL_FILES = {
    'interactions': r"C:\Users\Jy-Mentor-7\Downloads\string_interactions_short (7).tsv",
    'degrees': r"C:\Users\Jy-Mentor-7\Downloads\string_node_degrees (3).tsv",
    'enrichment': r"C:\Users\Jy-Mentor-7\Downloads\enrichment.all (1).tsv",
}


def load_local_interactions() -> Optional[pd.DataFrame]:
    """Load PPI interactions from local STRING file"""
    file_path = LOCAL_FILES['interactions']
    
    if not os.path.exists(file_path):
        logger.error(f"Local file not found: {file_path}")
        return None
    
    logger.info(f"Loading local STRING interactions: {file_path}")
    
    # Read TSV (header starts with #, skip that character)
    df = pd.read_csv(file_path, sep='\t')
    
    # Fix column names (remove # from first column)
    df.columns = [col.lstrip('#') for col in df.columns]
    
    logger.info(f"Loaded {len(df)} interactions")
    logger.info(f"Columns: {list(df.columns)}")
    
    # Rename columns to standard format
    if 'node1' in df.columns and 'node2' in df.columns:
        df = df.rename(columns={
            'node1': 'preferredName_A',
            'node2': 'preferredName_B',
            'combined_score': 'score'
        })
    
    # Normalize score to 0-1 if needed
    if df['score'].max() > 1:
        logger.info("Normalizing combined_score to 0-1 range")
        df['score'] = df['score'] / 1000.0
    
    logger.info(f"Score range: [{df['score'].min():.3f}, {df['score'].max():.3f}]")
    
    return df


def load_local_degrees() -> Optional[pd.DataFrame]:
    """Load node degrees from local STRING file"""
    file_path = LOCAL_FILES['degrees']
    
    if not os.path.exists(file_path):
        logger.error(f"Local file not found: {file_path}")
        return None
    
    logger.info(f"Loading local STRING degrees: {file_path}")
    
    df = pd.read_csv(file_path, sep='\t')
    
    # Fix column names (remove # from first column)
    df.columns = [col.lstrip('#') for col in df.columns]
    
    logger.info(f"Loaded {len(df)} nodes with degrees")
    logger.info(f"Columns: {list(df.columns)}")
    
    return df


def build_seed_genes() -> Set[str]:
    """Build seed gene set from Stage4"""
    seed_file = os.path.join(RESULTS_DIR, "stage4_seed_wgcna", "seed_pool_genes.txt")
    
    if not os.path.exists(seed_file):
        logger.error(f"Seed file not found: {seed_file}")
        return set()
    
    with open(seed_file, 'r', encoding='utf-8') as f:
        seed_genes = set(line.strip().upper() for line in f if line.strip())
    
    logger.info(f"Loaded {len(seed_genes)} seed genes from Stage4")
    return seed_genes


def filter_ppi_by_seed(ppi_df: pd.DataFrame, seed_genes: Set[str], 
                        min_score: float = 0.4) -> pd.DataFrame:
    """Filter PPI network to include only seed genes and their neighbors"""
    logger.info("Filtering PPI network...")
    
    # Filter by confidence score
    filtered = ppi_df[ppi_df['score'] >= min_score].copy()
    logger.info(f"After score filter (>{min_score}): {len(filtered)} edges")
    
    # Get all genes in network
    all_genes = set(filtered['preferredName_A'].str.upper()) | \
                set(filtered['preferredName_B'].str.upper())
    
    # Find seed genes present in network
    seed_in_network = seed_genes & all_genes
    logger.info(f"Seed genes in network: {len(seed_in_network)}/{len(seed_genes)}")
    
    # Keep edges where at least one node is a seed gene
    mask = (
        filtered['preferredName_A'].str.upper().isin(seed_genes) |
        filtered['preferredName_B'].str.upper().isin(seed_genes)
    )
    
    filtered = filtered[mask].copy()
    logger.info(f"After seed filtering: {len(filtered)} edges")
    
    return filtered


def compute_topology(ppi_df: pd.DataFrame, degrees_df: Optional[pd.DataFrame] = None) -> Dict:
    """Compute network topology metrics"""
    logger.info("Computing network topology...")
    
    G = nx.Graph()
    
    for _, row in ppi_df.iterrows():
        g1 = row["preferredName_A"].upper()
        g2 = row["preferredName_B"].upper()
        w = row.get("score", 1.0)
        G.add_edge(g1, g2, weight=w)
    
    logger.info(f"Network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    
    # Compute centrality metrics
    degree_dict = dict(G.degree())
    betweenness_dict = nx.betweenness_centrality(G, weight='weight')
    closeness_dict = nx.closeness_centrality(G)
    
    # Add local degree info if available
    local_degrees = {}
    if degrees_df is not None:
        for _, row in degrees_df.iterrows():
            gene = str(row['node']).upper()
            local_degrees[gene] = row['node_degree']
    
    topology = {}
    for node in G.nodes():
        topology[node] = {
            'degree': degree_dict.get(node, 0),
            'betweenness': betweenness_dict.get(node, 0.0),
            'closeness': closeness_dict.get(node, 0.0),
            'local_degree': local_degrees.get(node, degree_dict.get(node, 0)),
        }
    
    return topology


def save_results(ppi_df: pd.DataFrame, topology: Dict):
    """Save PPI and topology results"""
    # Save PPI edges
    ppi_file = os.path.join(STAGE_DIR, "string_ppi.tsv")
    ppi_df.to_csv(ppi_file, sep='\t', index=False)
    logger.info(f"PPI saved: {ppi_file}")
    
    # Save topology
    topo_file = os.path.join(STAGE_DIR, "ppi_topology.json")
    with open(topo_file, 'w', encoding='utf-8') as f:
        json.dump(topology, f, indent=2)
    logger.info(f"Topology saved: {topo_file}")
    
    # Save statistics
    stats = {
        'n_nodes': len(topology),
        'n_edges': len(ppi_df),
        'min_score': float(ppi_df['score'].min()),
        'max_score': float(ppi_df['score'].max()),
        'mean_score': float(ppi_df['score'].mean()),
        'source': 'local_string_files',
    }
    
    stats_file = os.path.join(STAGE_DIR, "ppi_stats.json")
    with open(stats_file, 'w', encoding='utf-8') as f:
        json.dump(stats, f, indent=2)
    
    logger.info(f"Statistics: {stats}")


def main():
    """Main workflow using local STRING files"""
    logger.info("=" * 60)
    logger.info("Stage 5: STRING PPI from Local Files (v9.0)")
    logger.info("=" * 60)
    logger.info(f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Load local files
    ppi_df = load_local_interactions()
    if ppi_df is None:
        logger.error("Failed to load local interactions")
        return
    
    degrees_df = load_local_degrees()
    
    # Load seed genes
    seed_genes = build_seed_genes()
    if not seed_genes:
        logger.error("No seed genes loaded")
        return
    
    # Filter by seed genes
    ppi_filtered = filter_ppi_by_seed(ppi_df, seed_genes)
    
    # Compute topology
    topology = compute_topology(ppi_filtered, degrees_df)
    
    # Save results
    save_results(ppi_filtered, topology)
    
    logger.info("=" * 60)
    logger.info("Stage 5 completed successfully (using local files)")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
