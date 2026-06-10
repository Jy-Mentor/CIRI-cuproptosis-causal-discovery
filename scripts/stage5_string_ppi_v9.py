# -*- coding: utf-8 -*-
"""
Stage 5: STRING PPI Network Construction (v9.0)
===============================================

FIX in v9.0:
- Use POST method for STRING API (correct batch handling)
- Unified score normalization (0-1 range)
- Proper error handling and retry logic
- Species validation

Reference:
  - STRING API documentation: https://string-db.org/help/api/
  - POST method for multiple identifiers

Input:
  - stage4/seed_pool_genes.txt
  
Output:
  - string_ppi.tsv: PPI network edges
  - ppi_topology.json: Network topology metrics
  - stage5.log: Run log
"""

import os
import sys
import json
import logging
import time
from datetime import datetime
from typing import List, Dict, Optional, Set

import numpy as np
import pandas as pd
import networkx as nx

# v9.0: Use requests for POST method
import requests

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import RESULTS_DIR

STAGE_DIR = os.path.join(RESULTS_DIR, "stage5_string_ppi")
os.makedirs(STAGE_DIR, exist_ok=True)

logger = logging.getLogger("stage5_v9")
logger.setLevel(logging.INFO)

file_handler = logging.FileHandler(os.path.join(STAGE_DIR, "stage5_v9.log"), encoding="utf-8", mode="w")
file_handler.setFormatter(logging.Formatter("%(asctime)s | %(levelname)s | %(message)s"))
logger.addHandler(file_handler)

console_handler = logging.StreamHandler(sys.stdout)
console_handler.setFormatter(logging.Formatter("%(message)s"))
logger.addHandler(console_handler)

# v9.0: Configuration
STRING_API_URL = "https://string-db.org/api/tsv/interactions"
STRING_NETWORK_URL = "https://string-db.org/api/tsv/network"
STRING_SPECIES = 9606  # Human
STRING_CONFIDENCE = 0.4  # Minimum confidence score
MAX_RETRIES = 5
BATCH_SIZE = 500  # v9.0: Larger batch size with POST
RATE_LIMIT_DELAY = 1.0  # Seconds between batches


def fetch_ppi_from_string_post(gene_list: List[str]) -> Optional[pd.DataFrame]:
    """
    v9.0: Fetch PPI from STRING using POST method
    
    POST method correctly handles multiple identifiers
    """
    logger.info(f"Fetching PPI from STRING (POST method, {len(gene_list)} genes)...")
    
    all_results = []
    n_batches = (len(gene_list) + BATCH_SIZE - 1) // BATCH_SIZE
    
    for batch_idx in range(n_batches):
        batch_genes = gene_list[batch_idx * BATCH_SIZE:(batch_idx + 1) * BATCH_SIZE]
        logger.info(f"  Batch {batch_idx + 1}/{n_batches}: {len(batch_genes)} genes")
        
        # v9.0: Use POST method with proper payload
        payload = {
            'identifiers': '\r\n'.join(batch_genes),  # STRING expects newline-separated
            'species': STRING_SPECIES,
            'limit': 1000000,
            'caller_identity': 'cir_bcp_pipeline_v9'
        }
        
        success = False
        for attempt in range(MAX_RETRIES):
            try:
                response = requests.post(STRING_API_URL, data=payload, timeout=60)
                
                if response.status_code == 200:
                    # Parse TSV response
                    from io import StringIO
                    df = pd.read_csv(StringIO(response.text), sep='\t')
                    
                    if len(df) > 0:
                        # v9.0: Normalize score to 0-1 range
                        if df['score'].max() > 1:
                            df['score'] = df['score'] / 1000.0
                        
                        all_results.append(df)
                        logger.info(f"    Retrieved {len(df)} interactions")
                        success = True
                        break
                    else:
                        logger.warning(f"    Empty response")
                        success = True  # Empty is not a failure
                        break
                else:
                    logger.warning(f"    HTTP {response.status_code}: {response.text[:100]}")
                    
            except Exception as e:
                wait_time = 2 ** attempt
                logger.warning(f"    Attempt {attempt + 1} failed: {e}")
                if attempt < MAX_RETRIES - 1:
                    logger.info(f"    Waiting {wait_time}s before retry...")
                    time.sleep(wait_time)
        
        if not success:
            logger.error(f"  Batch {batch_idx + 1} failed after {MAX_RETRIES} attempts")
        
        # Rate limiting
        if batch_idx < n_batches - 1:
            time.sleep(RATE_LIMIT_DELAY)
    
    if not all_results:
        logger.error("No PPI data retrieved")
        return None
    
    # Merge all batches
    merged_df = pd.concat(all_results, ignore_index=True)
    
    # Remove duplicates
    merged_df = merged_df.drop_duplicates(subset=['preferredName_A', 'preferredName_B'])
    
    # Filter by confidence
    merged_df = merged_df[merged_df['score'] >= STRING_CONFIDENCE]
    
    logger.info(f"Total PPI interactions: {len(merged_df)}")
    logger.info(f"Score range: [{merged_df['score'].min():.3f}, {merged_df['score'].max():.3f}]")
    
    return merged_df


def load_local_ppi() -> Optional[pd.DataFrame]:
    """Load locally cached PPI network"""
    ppi_file = os.path.join(STAGE_DIR, "string_ppi.tsv")
    
    if os.path.exists(ppi_file):
        logger.info(f"Loading local PPI: {ppi_file}")
        df = pd.read_csv(ppi_file, sep='\t')
        
        # v9.0: Ensure score normalization
        if 'score' in df.columns:
            if df['score'].max() > 1:
                logger.info("Normalizing scores to 0-1 range")
                df['score'] = df['score'] / 1000.0
        
        logger.info(f"Loaded {len(df)} interactions")
        return df
    
    return None


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


def compute_topology(ppi_df: pd.DataFrame) -> tuple:
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
    
    topology = {}
    for node in G.nodes():
        topology[node] = {
            'degree': degree_dict.get(node, 0),
            'betweenness': betweenness_dict.get(node, 0.0),
            'closeness': closeness_dict.get(node, 0.0),
        }
    
    return G, topology


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
    }
    
    stats_file = os.path.join(STAGE_DIR, "ppi_stats.json")
    with open(stats_file, 'w', encoding='utf-8') as f:
        json.dump(stats, f, indent=2)
    
    logger.info(f"Statistics: {stats}")


def main():
    """Main workflow"""
    logger.info("=" * 60)
    logger.info("Stage 5: STRING PPI Network Construction (v9.0)")
    logger.info("=" * 60)
    logger.info(f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Species: {STRING_SPECIES} (Homo sapiens)")
    logger.info(f"Confidence threshold: {STRING_CONFIDENCE}")
    logger.info(f"API method: POST (v9.0)")
    
    # Load seed genes
    seed_genes = build_seed_genes()
    if not seed_genes:
        logger.error("No seed genes loaded")
        return
    
    # Try to fetch from STRING
    ppi_df = fetch_ppi_from_string_post(list(seed_genes))
    
    # Fallback to local cache
    if ppi_df is None:
        logger.info("Falling back to local cache...")
        ppi_df = load_local_ppi()
    
    if ppi_df is None or len(ppi_df) == 0:
        logger.error("Unable to obtain PPI data")
        return
    
    # Compute topology
    G, topology = compute_topology(ppi_df)
    
    # Save results
    save_results(ppi_df, topology)
    
    logger.info("=" * 60)
    logger.info("Stage 5 completed successfully")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
