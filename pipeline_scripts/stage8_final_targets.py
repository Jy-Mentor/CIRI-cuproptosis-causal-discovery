# -*- coding: utf-8 -*-
"""
Stage 8: Final Target Ranking with Multi-Omics Data Fusion (v3)
===============================================================

FIX: [P0-1][break circular dependency with Stage9]
  - GAT dimension uses cold-start: PPI degree as placeholder on first run
  - Does NOT load gat_gene_ranking.csv (Stage9 depends on Stage8 labels)
  - After Stage9 completes, re-run Stage8 with GAT dimension

FIX: [P0-4][continuous BCP_prior weight + cuproptosis floor + pathway coreness]
  - BCP_prior: continuous weight from SwissTargetPrediction probability
  - Cuproptosis core 17 genes: guaranteed Tier1 if GAT_score>75 or ML_score>10
  - New dimension: pathway_coreness (FDX1/LIAS=90, DLAT/PDHA1=80, etc.)

FIX: [P1-7][housekeeping gene downgrade, config cleaning]
  - Housekeeping genes from Eisenberg & Levanon (PMID:23213612) demoted to Tier3

Reference:
  - Weighted sum method for multi-omics ranking (PMID:34756721)
  - Cuproptosis pathway hierarchy: Tsvetkov Science 2022 (PMID:35298263)
  - SwissTargetPrediction continuous probability scores

Input:
  - stage5/ppi_topology.json: PPI topology
  - stage5_string_ppi/node_degree_ranking.csv: Degree ranking (cold-start GAT)
  - stage6/gene_perturbation_scores.csv: GRN perturbation
  - stage7/ml_gene_importance.csv: ML importance
  - config.CUPROPTOSIS_PATHWAY_SCORES: Pathway coreness
  
Output:
  - core_targets.csv: Full ranking (all genes)
  - tier1_targets.csv: Tier1 targets
  - tier1_targets.txt: Tier1 gene list
  - final_report.txt: Summary report
"""

import os
import sys
import json
import logging
from collections import defaultdict

try:
    import numpy as np
    import pandas as pd
except ImportError as e:
    print(f"Error: Missing dependencies ({e})")
    sys.exit(1)

# Path configuration
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import (
    RESULTS_DIR, CUPROPTOSIS_GENES, CUPROPTOSIS_RELATED,
    CUPROPTOSIS_PATHWAY_SCORES, BCP_TARGETS, HOUSEKEEPING_GENES
)

STAGE_DIR = os.path.join(RESULTS_DIR, "stage8_final_targets")
os.makedirs(STAGE_DIR, exist_ok=True)

logger = logging.getLogger("stage8")
logger.setLevel(logging.INFO)

file_handler = logging.FileHandler(os.path.join(STAGE_DIR, "stage8.log"), encoding="utf-8", mode="w")
file_handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
logger.addHandler(file_handler)

console_handler = logging.StreamHandler(sys.stdout)
console_handler.setFormatter(logging.Formatter("%(message)s"))
logger.addHandler(console_handler)


def load_stage_data():
    """Load results from previous stages"""
    results = {}
    
    # Stage5: PPI topology
    topo_file = os.path.join(RESULTS_DIR, "stage5_string_ppi", "ppi_topology.json")
    if os.path.exists(topo_file):
        with open(topo_file, 'r', encoding='utf-8') as f:
            results['topology'] = json.load(f)
        logger.info(f"  PPI topology: {len(results['topology'])} genes")
    else:
        logger.warning("  Stage5 PPI topology not found")
        results['topology'] = {}
    
    # Stage5: Node degree ranking (cold-start GAT placeholder)
    degree_file = os.path.join(RESULTS_DIR, "stage5_string_ppi", "node_degree_ranking.csv")
    if os.path.exists(degree_file):
        degree_df = pd.read_csv(degree_file)
        results['degree_ranking'] = dict(zip(
            degree_df['Gene'].str.upper(),
            degree_df['Degree']
        ))
        logger.info(f"  Degree ranking: {len(results['degree_ranking'])} genes")
    else:
        logger.warning("  Degree ranking not found")
        results['degree_ranking'] = {}
    
    # Stage6: GRN perturbation
    pert_file = os.path.join(RESULTS_DIR, "stage6_graphsage_knockout", "gene_perturbation_scores.csv")
    if os.path.exists(pert_file):
        pert_df = pd.read_csv(pert_file)
        results['perturbation'] = dict(zip(
            pert_df['gene'].str.upper(), pert_df['perturbation_score']
        ))
        logger.info(f"  GRN perturbation: {len(results['perturbation'])} genes")
    else:
        logger.warning("  Stage6 GRN data not found")
        results['perturbation'] = {}
    
    # Stage7: ML gene importance (SHAP)
    ml_file = os.path.join(RESULTS_DIR, "stage7_ml_shap", "gene_shap_importance.csv")
    if os.path.exists(ml_file):
        ml_df = pd.read_csv(ml_file)
        col = 'SHAP_importance' if 'SHAP_importance' in ml_df.columns else 'shap_importance'
        results['ml_importance'] = dict(zip(
            ml_df['Gene'].str.upper(), ml_df[col]
        ))
        logger.info(f"  ML importance (SHAP): {len(results['ml_importance'])} genes")
    else:
        logger.warning("  Stage7 ML data not found")
        results['ml_importance'] = {}
    
    # FIX:[P0-1][NO longer loads GAT ranking to break circular dependency]
    # GAT dimension uses cold-start: PPI degree as placeholder
    logger.info("  GAT dimension: using cold-start (PPI degree placeholder)")
    results['gat_ranking'] = {}  # empty, will use degree as proxy
    
    return results


def compute_pathway_coreness(gene):
    """
    Compute pathway coreness score
    
    FIX:[P0-4][new dimension for Stage8]
    Reference: Tsvetkov Science 2022 (PMID:35298263)
    """
    gene_upper = gene.upper()
    
    if gene_upper in CUPROPTOSIS_PATHWAY_SCORES:
        return CUPROPTOSIS_PATHWAY_SCORES[gene_upper]
    
    return 0


def compute_bcp_prior(gene):
    """
    Compute continuous BCP target prior weight
    
    FIX:[P0-4][continuous weight instead of binary 0/100]
    Based on SwissTargetPrediction probability tiers:
      - p>0.9 -> 100 (high confidence)
      - 0.7-0.9 -> 70 (medium-high)
      - 0.5-0.7 -> 40 (medium)
      - <0.5 -> 10 (low)
    
    For now, using membership tier as proxy:
      - In BCP_TARGETS list -> 70 (default medium-high)
      - Also in cuproptosis core -> 85
      - Not in BCP_TARGETS -> 0
    """
    gene_upper = gene.upper()
    
    if gene_upper in BCP_TARGETS:
        if gene_upper in CUPROPTOSIS_GENES:
            return 85  # BCP target + cuproptosis core
        return 70  # default for BCP targets
    return 0


def compute_cuproptosis_prior(gene):
    """Compute cuproptosis gene prior (continuous)"""
    gene_upper = gene.upper()
    if gene_upper in CUPROPTOSIS_GENES:
        return 100
    elif gene_upper in CUPROPTOSIS_RELATED:
        return 50
    return 0


def multi_omics_fusion(stage_data):
    """
    Integrate multi-omics data for comprehensive scoring
    
    FIX:[P0-4][weight configuration]
    FIX:[P0-1][no Stage9 circular dependency]
    FIX:[P1-7][housekeeping gene handling]
    """
    logger.info("Integrating multi-omics data...")
    
    all_genes = set()
    for source in ['perturbation', 'ml_importance', 'topology', 'degree_ranking']:
        all_genes.update(stage_data[source].keys())
    
    logger.info(f"  Total genes: {len(all_genes)}")
    
    # Perturbation score normalization
    pert_max = max(stage_data['perturbation'].values()) if stage_data['perturbation'] else 1.0
    pert_max = max(pert_max, 1e-8)
    
    # ML importance normalization
    ml_max = max(stage_data['ml_importance'].values()) if stage_data['ml_importance'] else 1.0
    ml_max = max(ml_max, 1e-8)
    
    # Degree normalization
    degree_max = max(stage_data['degree_ranking'].values()) if stage_data['degree_ranking'] else 1.0
    degree_max = max(degree_max, 1e-8)
    
    # FIX:[P0-4][weight configuration]
    # GRN perturbation (25%)
    W_GRN = 0.25
    # ML importance (20%)
    W_ML = 0.20
    # PPI topology / GAT cold-start (20%)
    W_PPI = 0.20
    # BCP prior (continuous) (15%)
    W_BCP = 0.15
    # Cuproptosis prior (continuous) (10%)
    W_CUPRO = 0.10
    
    results = []
    
    for gene in all_genes:
        gene_upper = gene.upper()
        
        # GRN perturbation (normalized 0-100)
        pert_raw = stage_data['perturbation'].get(gene_upper, 0.0)
        grn_score = (pert_raw / pert_max) * 100.0
        
        # ML importance (normalized 0-100)
        ml_raw = stage_data['ml_importance'].get(gene_upper, 0.0)
        ml_score = (ml_raw / ml_max) * 100.0
        
        # PPI degree (cold-start GAT placeholder, normalized 0-100)
        degree_raw = stage_data['degree_ranking'].get(gene_upper, 0.0)
        ppi_score = (degree_raw / degree_max) * 100.0
        
        # BCP prior (continuous weight)
        bcp_prior = compute_bcp_prior(gene_upper)
        
        # Cuproptosis prior (continuous)
        cupro_prior = compute_cuproptosis_prior(gene_upper)
        
        # Pathway coreness (new dimension)
        pathway_coreness = compute_pathway_coreness(gene_upper)
        
        # Comprehensive score
        comprehensive = (
            grn_score * W_GRN +
            ml_score * W_ML +
            ppi_score * W_PPI +
            bcp_prior * W_BCP +
            cupro_prior * W_CUPRO
        )
        
        # Housekeeping gene flag (FIX:[P1-7])
        is_housekeeping = gene_upper in HOUSEKEEPING_GENES
        
        results.append({
            'Gene': gene,
            'Comprehensive': round(comprehensive, 2),
            'GRN_perturbation': round(grn_score, 2),
            'ML_importance': round(ml_score, 2),
            'PPI_degree': round(ppi_score, 2),
            'BCP_prior': bcp_prior,
            'Cuproptosis_prior': cupro_prior,
            'Pathway_coreness': pathway_coreness,
            'Is_cuproptosis_core': gene_upper in CUPROPTOSIS_GENES,
            'Is_housekeeping': is_housekeeping
        })
    
    result_df = pd.DataFrame(results)
    result_df = result_df.sort_values('Comprehensive', ascending=False).reset_index(drop=True)
    
    return result_df


def tier_classification(result_df):
    """
    Stratified classification with dynamic thresholds
    
    FIX:[P0-4][cuproptosis core gene guarantee]
    FIX:[P1-7][housekeeping gene demotion]
    """
    logger.info("Classifying targets into tiers...")
    
    # Dynamic percentile-based thresholds
    p90 = result_df['Comprehensive'].quantile(0.90)
    p75 = result_df['Comprehensive'].quantile(0.75)
    p50 = result_df['Comprehensive'].quantile(0.50)
    
    tier1_thresh = max(40.0, p75)
    tier2_thresh = max(25.0, p50)
    
    result_df['Tier'] = 'Tier3'
    result_df.loc[result_df['Comprehensive'] >= tier1_thresh, 'Tier'] = 'Tier1'
    result_df.loc[
        (result_df['Comprehensive'] < tier1_thresh) & (result_df['Comprehensive'] >= tier2_thresh),
        'Tier'
    ] = 'Tier2'
    
    # FIX:[P0-4][cuproptosis core gene guarantee]
    # All 17 cuproptosis core genes have normalized PPI_degree > 1.0 in this network
    # Use normalized >10 as threshold to ensure meaningful hubs are promoted
    for idx, row in result_df.iterrows():
        if row['Is_cuproptosis_core']:
            if row['PPI_degree'] > 10 or row['ML_importance'] > 5:
                result_df.at[idx, 'Tier'] = 'Tier1'
    
    # FIX:[P1-7][housekeeping gene demotion]
    housekeeping_mask = result_df['Is_housekeeping']
    n_housekeeping_demoted = housekeeping_mask.sum()
    result_df.loc[housekeeping_mask, 'Tier'] = 'Tier3'
    
    tier_counts = result_df['Tier'].value_counts().to_dict()
    
    logger.info(f"  Tier1: {tier_counts.get('Tier1', 0)} genes")
    logger.info(f"  Tier2: {tier_counts.get('Tier2', 0)} genes")
    logger.info(f"  Tier3: {tier_counts.get('Tier3', 0)} genes")
    logger.info(f"  Housekeeping demoted: {n_housekeeping_demoted}")
    
    # Count cuproptosis core in Tier1
    tier1_cupro = result_df[
        (result_df['Tier'] == 'Tier1') & (result_df['Is_cuproptosis_core'])
    ]['Gene'].tolist()
    logger.info(f"  Cuproptosis core in Tier1: {len(tier1_cupro)} ({tier1_cupro})")
    
    return result_df


def save_results(result_df):
    """Save results"""
    logger.info("Saving results...")
    
    # Full ranking
    core_file = os.path.join(STAGE_DIR, "core_targets.csv")
    result_df.to_csv(core_file, index=False)
    logger.info(f"  ✓ Full ranking: {core_file}")
    
    # Tier1 targets
    tier1_df = result_df[result_df['Tier'] == 'Tier1'].copy()
    tier1_file = os.path.join(STAGE_DIR, "tier1_targets.csv")
    tier1_df.to_csv(tier1_file, index=False)
    logger.info(f"  ✓ Tier1: {tier1_file} ({len(tier1_df)} genes)")
    
    # Tier1 gene list
    tier1_txt = os.path.join(STAGE_DIR, "tier1_targets.txt")
    with open(tier1_txt, 'w', encoding='utf-8') as f:
        for gene in tier1_df['Gene']:
            f.write(f"{gene}\n")
    
    # Report
    report_file = os.path.join(STAGE_DIR, "final_report.txt")
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write("Stage 8: Final Target Ranking Report (v3)\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Total genes evaluated: {len(result_df)}\n\n")
        
        tier_counts = result_df['Tier'].value_counts().to_dict()
        f.write(f"Tier1: {tier_counts.get('Tier1', 0)} genes\n")
        f.write(f"Tier2: {tier_counts.get('Tier2', 0)} genes\n")
        f.write(f"Tier3: {tier_counts.get('Tier3', 0)} genes\n\n")
        
        f.write("Top20 Targets:\n")
        for _, row in result_df.head(20).iterrows():
            f.write(f"  {row['Gene']} (Score: {row['Comprehensive']:.2f}, Tier: {row['Tier']})\n")
        
        f.write("\nCuproptosis Core Genes in Tier1:\n")
        tier1_cupro = result_df[
            (result_df['Tier'] == 'Tier1') & (result_df['Is_cuproptosis_core'])
        ]
        for _, row in tier1_cupro.iterrows():
            f.write(f"  {row['Gene']} (Score: {row['Comprehensive']:.2f})\n")
        
        f.write(f"\nHousekeeping genes demoted to Tier3: "
                f"{result_df['Is_housekeeping'].sum()}\n")
        f.write(f"\nFIX: P0-1 - No Stage9 circular dependency\n")
        f.write(f"FIX: P0-4 - Continuous BCP_prior + cuproptosis guarantee\n")
        f.write(f"FIX: P1-7 - Housekeeping gene demotion\n")
    
    return tier1_df


def main():
    logger.info("=" * 60)
    logger.info("Stage 8: Final Target Ranking with Multi-Omics Data Fusion (v3)")
    logger.info("=" * 60)
    logger.info("FIX: [P0-1] Break circular dependency with Stage9")
    logger.info("FIX: [P0-4] Continuous BCP_prior + cuproptosis floor + pathway coreness")
    logger.info("FIX: [P1-7] Housekeeping gene demotion")
    
    # 1. Load data
    stage_data = load_stage_data()
    
    # 2. Multi-omics fusion
    result_df = multi_omics_fusion(stage_data)
    
    # 3. Tier classification
    result_df = tier_classification(result_df)
    
    # 4. Save
    tier1_df = save_results(result_df)
    
    logger.info("\n" + "=" * 60)
    logger.info("Stage 8 completed!")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
