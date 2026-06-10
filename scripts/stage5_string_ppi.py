# -*- coding: utf-8 -*-
"""
Stage 5: STRING PPI Network Construction & Topology Analysis
============================================================

FIX: [P0-3][species mismatch + over-filtering + retry]
  - Changed species from 10116 (rat) to 9606 (human) since seed pool uses human symbols
  - Filter logic: "at least one end in seed pool" (1-hop neighbor expansion)
  - Added 3-retry exponential backoff for STRING API
  - Added validation: node count < 500 triggers warning + file freshness check

Reference:
  - STRING database (PMID: 36573205, Szklarczyk NAR 2023)
  - species 9606 for human gene symbols

Input:
  - Seed pool gene list (human symbols)
  
Output:
  - string_ppi.tsv: PPI edge list
  - ppi_topology.json: Node topology features
  - stage5.log: Run log
"""

import os
import sys
import warnings
import logging
import json
import time
from datetime import datetime, timedelta

warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import networkx as nx

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import RESULTS_DIR, CUPROPTOSIS_GENES, BCP_TARGETS
from scripts.utils import setup_logger

STAGE_DIR = os.path.join(RESULTS_DIR, "stage5_string_ppi")
os.makedirs(STAGE_DIR, exist_ok=True)

logger = setup_logger("stage5", os.path.join(STAGE_DIR, "stage5.log"))

STRING_CONFIDENCE = 0.4  # medium confidence (PMID: 36573205)
STRING_SPECIES = "9606"  # FIX:[P0-3][Human species for human gene symbols]
MAX_RETRIES = 3
MIN_NODE_COUNT = 500  # FIX:[P0-3][minimum node validation threshold]
CACHE_MAX_AGE_DAYS = 7  # FIX:[P0-3][cache freshness threshold]


def build_seed_genes():
    """Build seed pool gene union (human symbols)
    
    FIX:[v7][优先读取Stage4过滤后的种子池，确保PPI紧邻过滤传导至Stage5]
    """
    seeds = set()
    
    # FIX:[v7][try Stage4 filtered seed pool first]
    stage4_seed_file = os.path.join(RESULTS_DIR, "stage4_seed_wgcna", "seed_pool_genes.txt")
    if os.path.exists(stage4_seed_file):
        logger.info(f"  Loading Stage4 filtered seed pool: {stage4_seed_file}")
        with open(stage4_seed_file, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    seeds.add(line.upper())
        logger.info(f"  Stage4 filtered seed pool: {len(seeds)} genes")
    else:
        logger.warning("  Stage4 seed pool not found, falling back to manual seed pool")
        
        for g in CUPROPTOSIS_GENES:
            seeds.add(g.upper())
        
        for g in BCP_TARGETS:
            seeds.add(g.upper())
        
        ciri_genes = [
            "IL6","TNF","IL1B","CXCL8","PTGS2","NFKB1","MAPK1","MAPK3",
            "AKT1","MTOR","HIF1A","VEGFA","BCL2","CASP3","TP53",
            "JAK2","STAT3","RELA","FOS","JUN","MYC","CDKN1A",
            "SOD1","CAT","GPX1","NRF1","NFE2L2","KEAP1"
        ]
        for g in ciri_genes:
            seeds.add(g.upper())
    
    logger.info(f"  Seed pool (final): {len(seeds)} genes")
    return seeds


def fetch_ppi_from_string(gene_list):
    """
    Fetch PPI network from STRING database with retry
    
    FIX:[P0-3][species=9606 human + 3-retry exponential backoff]
    FIX:[v8][batch requests to avoid HTTP 414 URI Too Long]
    Reference: https://string-db.org/help/api/
    """
    logger.info("Fetching PPI from STRING database (species=9606 human)...")
    
    batch_size = 200
    all_ppi_dfs = []
    n_batches = (len(gene_list) + batch_size - 1) // batch_size
    
    for batch_idx in range(n_batches):
        batch_genes = gene_list[batch_idx * batch_size:(batch_idx + 1) * batch_size]
        logger.info(f"  Batch {batch_idx + 1}/{n_batches}: {len(batch_genes)} genes")
        
        url_base = (
            "https://string-db.org/api/tsv/interactions"
            f"?identifiers={'%0d'.join(batch_genes)}"
            f"&species={STRING_SPECIES}"
            f"&limit=1000000"
        )
        
        batch_success = False
        for attempt in range(MAX_RETRIES):
            try:
                ppi_df = pd.read_csv(url_base, sep='\t')
                all_ppi_dfs.append(ppi_df)
                logger.info(f"    Batch {batch_idx + 1} returned: {len(ppi_df)} interactions")
                batch_success = True
                break
            except Exception as e:
                wait_time = 2 ** attempt * 5
                logger.warning(f"    Attempt {attempt + 1} failed: {e}")
                if attempt < MAX_RETRIES - 1:
                    time.sleep(wait_time)
        
        if not batch_success:
            logger.warning(f"  Batch {batch_idx + 1} failed after {MAX_RETRIES} attempts, continuing...")
        
        if batch_idx < n_batches - 1:
            time.sleep(1.0)  # rate limiting
    
    if not all_ppi_dfs:
        logger.error("  All batches failed")
        return None
    
    merged_df = pd.concat(all_ppi_dfs, ignore_index=True).drop_duplicates()
    logger.info(f"  Total STRING interactions (merged): {len(merged_df)}")
    return merged_df


def _file_is_fresh(filepath, max_age_days=7):
    """Check if cached file is within max_age_days"""
    if not os.path.exists(filepath):
        return False
    mtime = datetime.fromtimestamp(os.path.getmtime(filepath))
    return (datetime.now() - mtime) < timedelta(days=max_age_days)


def load_local_ppi():
    """Load locally cached PPI network with freshness validation"""
    # FIX:[v8][优先使用已存在的string_ppi.tsv (之前v7下载的)]
    existing_ppi = os.path.join(STAGE_DIR, "string_ppi.tsv")
    if os.path.exists(existing_ppi):
        logger.info(f"  Using existing PPI file: {existing_ppi}")
        ppi_df = pd.read_csv(existing_ppi, sep='\t')
        logger.info(f"    Loaded {len(ppi_df)} interactions")
        return ppi_df
    
    ppi_file = os.path.join(STAGE_DIR, "string_ppi_raw.tsv")
    if _file_is_fresh(ppi_file, CACHE_MAX_AGE_DAYS):
        logger.info(f"  Using fresh cached PPI ({os.path.getmtime(ppi_file)})")
        ppi_df = pd.read_csv(ppi_file, sep='\t')
        return ppi_df
    
    # Try old stage9 PPI data
    old_file = os.path.join(RESULTS_DIR, "stage9_ppi_gat", "string_ppi.tsv")
    if os.path.exists(old_file):
        logger.info("  Loading legacy PPI from stage9_ppi_gat...")
        ppi_df = pd.read_csv(old_file, sep='\t')
        if 'combined_score' in ppi_df.columns and 'score' not in ppi_df.columns:
            ppi_df['score'] = ppi_df['combined_score'] / 1000.0
        return ppi_df
    
    logger.warning("  No local PPI data available")
    return None


def filter_ppi_network(ppi_df, seed_genes):
    """
    Filter PPI network
    
    FIX:[P0-3][at least one end in seed pool (1-hop neighbor expansion)]
    Previous: both ends must be in seed pool (overly restrictive)
    """
    logger.info("Filtering PPI network (1-hop neighbor expansion)...")
    
    # Keep edges with score >= STRING_CONFIDENCE
    ppi_df = ppi_df[ppi_df['score'] >= STRING_CONFIDENCE].copy()
    
    # Remove self-loops
    ppi_df = ppi_df[ppi_df["preferredName_A"] != ppi_df["preferredName_B"]]
    
    # FIX:[P0-3][at least one end in seed pool]
    seed_genes_upper = {g.upper() for g in seed_genes}
    mask_a = ppi_df["preferredName_A"].isin(seed_genes_upper)
    mask_b = ppi_df["preferredName_B"].isin(seed_genes_upper)
    ppi_df = ppi_df[mask_a | mask_b].copy()
    
    logger.info(f"  After filtering: {len(ppi_df)} edges")
    return ppi_df


def compute_topology(ppi_df):
    """Compute PPI network topology features"""
    logger.info("Computing network topology...")
    
    G = nx.Graph()
    
    for _, row in ppi_df.iterrows():
        G.add_edge(row["preferredName_A"], row["preferredName_B"], weight=row["score"])
    
    logger.info(f"  Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")
    
    # FIX:[P0-3][validation: node count < MIN_NODE_COUNT]
    if G.number_of_nodes() < MIN_NODE_COUNT:
        logger.warning(
            f"  PPI network has only {G.number_of_nodes()} nodes (< {MIN_NODE_COUNT}). "
            f"This may indicate species mismatch or API failure."
        )
    
    # Topology metrics
    degree = dict(G.degree())
    betweenness = nx.betweenness_centrality(G, weight='weight')
    closeness = nx.closeness_centrality(G)
    eigenvector = nx.eigenvector_centrality(G, max_iter=1000, tol=1e-6, weight='weight')
    
    # Connectivity analysis
    n_components = nx.number_connected_components(G)
    largest_cc = max(nx.connected_components(G), key=len)
    logger.info(f"  Connected components: {n_components}, Largest: {len(largest_cc)} nodes")
    
    # Save topology data
    topology = {}
    all_genes = set(degree.keys())
    
    for gene in all_genes:
        topology[gene] = {
            'degree': degree.get(gene, 0),
            'betweenness': betweenness.get(gene, 0.0),
            'closeness': closeness.get(gene, 0.0),
            'eigenvector': eigenvector.get(gene, 0.0),
            'in_largest_cc': gene in largest_cc
        }
    
    # Statistics
    deg_values = list(degree.values())
    bc_values = [v for v in betweenness.values() if v > 0]
    
    logger.info(f"  Degree: mean={np.mean(deg_values):.1f}, max={max(deg_values)}")
    logger.info(f"  Betweenness > 0: {len(bc_values)}/{len(betweenness)} nodes")
    
    return G, topology, ppi_df


def save_results(ppi_df, topology):
    """Save results"""
    logger.info("Saving results...")
    
    # PPI edge list
    ppi_file = os.path.join(STAGE_DIR, "string_ppi.tsv")
    ppi_df.to_csv(ppi_file, sep='\t', index=False)
    logger.info(f"  ✓ PPI network: {ppi_file} ({len(ppi_df)} edges)")
    
    # Topology features
    topo_file = os.path.join(STAGE_DIR, "ppi_topology.json")
    with open(topo_file, 'w', encoding='utf-8') as f:
        json.dump(topology, f, indent=2, ensure_ascii=False)
    logger.info(f"  ✓ Topology: {topo_file} ({len(topology)} genes)")
    
    # Node degree ranking
    degree_ranking = sorted(
        [(g, d['degree']) for g, d in topology.items()],
        key=lambda x: x[1], reverse=True
    )
    rank_df = pd.DataFrame(degree_ranking, columns=['Gene', 'Degree'])
    rank_df['Rank'] = range(1, len(rank_df) + 1)
    rank_file = os.path.join(STAGE_DIR, "node_degree_ranking.csv")
    rank_df.to_csv(rank_file, index=False)
    logger.info(f"  ✓ Degree ranking: {rank_file}")


def main():
    logger.info("=" * 60)
    logger.info("Stage 5: STRING PPI Network Construction (v8)")
    logger.info("=" * 60)
    logger.info(f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Species: {STRING_SPECIES} (Homo sapiens)")
    logger.info(f"Confidence threshold: {STRING_CONFIDENCE}")
    
    # 1. Seed pool
    seed_genes = build_seed_genes()
    
    # 2. Load PPI data (v8: 优先使用本地缓存，避免网络请求)
    logger.info("Loading PPI data (v8: using local cache)...")
    ppi_df = load_local_ppi()
    
    if ppi_df is None:
        logger.info("  No local cache, fetching from STRING API...")
        ppi_df = fetch_ppi_from_string(list(seed_genes))
    
    if ppi_df is None:
        logger.error("Unable to obtain PPI data")
        return
    
    # 3. Filter (v8: 使用更大的种子池2029基因)
    ppi_df = filter_ppi_network(ppi_df, seed_genes)
    
    # 4. Topology analysis
    G, topology, ppi_df = compute_topology(ppi_df)
    
    # 5. Save
    save_results(ppi_df, topology)
    
    logger.info("\n" + "=" * 60)
    logger.info("Stage 5 completed!")
    logger.info(f"  Nodes: {len(topology)}, Edges: {len(ppi_df)}")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
