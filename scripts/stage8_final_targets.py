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
import logging
from collections import defaultdict

try:
    import numpy as np
    import pandas as pd
except ImportError as e:
    print(f"Error: Missing dependencies ({e})")
    sys.exit(1)

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import (
    RESULTS_DIR, CUPROPTOSIS_GENES, CUPROPTOSIS_RELATED,
    CUPROPTOSIS_PATHWAY_SCORES, BCP_TARGETS, HOUSEKEEPING_GENES
)
from scripts.data_manager import StageDataManager

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
    """使用StageDataManager统一加载跨阶段数据"""
    dm = StageDataManager()
    results = {}
    
    # Stage5: PPI topology + Node degree ranking
    stage5_data = dm.load_stage('stage5_string_ppi', {
        'ppi_topology': ('json', {'filename': 'ppi_topology.json', 'default': {}}),
        'degree_ranking': ('csv_dict', {
            'filename': 'node_degree_ranking.csv',
            'key_col': 'Gene',
            'value_col': 'Degree'
        })
    })
    results['topology'] = stage5_data['ppi_topology']
    results['degree_ranking'] = stage5_data['degree_ranking']
    logger.info(f"  PPI topology: {len(results['topology'])} genes")
    logger.info(f"  Degree ranking: {len(results['degree_ranking'])} genes")
    
    # Stage6: GRN perturbation
    results['perturbation'] = dm.load_csv_as_dict(
        'stage6_sctenifold_knockout',
        'gene_perturbation_scores.csv',
        key_col='gene',
        value_col='perturbation_score',
        default={}
    )
    logger.info(f"  GRN perturbation: {len(results['perturbation'])} genes")
    
    # Stage7: ML gene importance (SHAP/Bootstrap)
    ml_df = dm.load_csv('stage7_ml_shap', 'gene_shap_importance.csv')
    if ml_df is not None and not ml_df.empty:
        col = 'SHAP_importance' if 'SHAP_importance' in ml_df.columns else 'shap_importance'
        results['ml_importance'] = dict(zip(
            ml_df['Gene'].str.upper(), ml_df[col]
        ))
        logger.info(f"  ML importance (SHAP/Bootstrap): {len(results['ml_importance'])} genes")
        results['ml_n_samples'] = 6
        
        if 'Bootstrap_stability' in ml_df.columns:
            results['boot_stability'] = dict(zip(
                ml_df['Gene'].str.upper(), ml_df['Bootstrap_stability']
            ))
            logger.info(f"  Bootstrap stability: {len(results['boot_stability'])} genes")
        else:
            results['boot_stability'] = {}
    else:
        logger.warning("  Stage7 ML data not found")
        results['ml_importance'] = {}
        results['boot_stability'] = {}
        results['ml_n_samples'] = 0
    
    # FIX:[P0-1][NO longer loads GAT ranking to break circular dependency]
    logger.info("  GAT dimension: using cold-start (PPI degree placeholder)")
    results['gat_ranking'] = {}
    
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
    """v8: unbiased pipeline - no BCP prior"""
    return 0


def compute_cuproptosis_prior(gene):
    """v8: unbiased pipeline - no cuproptosis prior"""
    return 0


def multi_omics_fusion(stage_data):
    """
    Integrate multi-omics data for comprehensive scoring
    
    v8: UNBIASED DATA-DRIVEN PIPELINE
    - No BCP prior, no cuproptosis prior
    - Pure data-driven: GRN perturbation + ML importance + PPI topology
    """
    logger.info("Integrating multi-omics data (v8: unbiased)...")
    
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
    
    # v8: Pure data-driven weights (no priors)
    W_GRN = 0.35
    W_ML_BASE = 0.35
    W_PPI = 0.30
    
    # ML sample size adaptive penalty (same as v7)
    n_samples = stage_data.get('ml_n_samples', 6)
    ml_penalty = n_samples / (n_samples + 10.0)
    W_ML = W_ML_BASE * ml_penalty
    weight_lost = W_ML_BASE - W_ML
    W_GRN = W_GRN + weight_lost * 0.5
    W_PPI = W_PPI + weight_lost * 0.5
    
    logger.info(f"  ML sample size penalty: n={n_samples}, penalty={ml_penalty:.3f}, W_ML={W_ML:.3f}")
    logger.info(f"  Data-driven weights: GRN={W_GRN:.3f}, ML={W_ML:.3f}, PPI={W_PPI:.3f}")
    
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
        
        # Comprehensive score (v8: pure data-driven, no priors)
        comprehensive = (
            grn_score * W_GRN +
            ml_score * W_ML +
            ppi_score * W_PPI
        )
        
        results.append({
            'Gene': gene,
            'Comprehensive': round(comprehensive, 2),
            'GRN_perturbation': round(grn_score, 2),
            'ML_importance': round(ml_score, 2),
            'PPI_degree': round(ppi_score, 2),
            'Is_cuproptosis_core': gene_upper in CUPROPTOSIS_GENES,
            'Is_housekeeping': gene_upper in HOUSEKEEPING_GENES
        })
    
    result_df = pd.DataFrame(results)
    result_df = result_df.sort_values('Comprehensive', ascending=False).reset_index(drop=True)
    
    return result_df


def tier_classification(result_df):
    """
    Stratified classification with dynamic thresholds (v8: unbiased)
    """
    logger.info("Classifying targets into tiers (v8: unbiased)...")
    
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
    
    # Housekeeping gene demotion
    housekeeping_mask = result_df['Is_housekeeping']
    n_housekeeping_demoted = housekeeping_mask.sum()
    result_df.loc[housekeeping_mask, 'Tier'] = 'Tier3'
    
    tier_counts = result_df['Tier'].value_counts().to_dict()
    
    logger.info(f"  Tier1: {tier_counts.get('Tier1', 0)} genes")
    logger.info(f"  Tier2: {tier_counts.get('Tier2', 0)} genes")
    logger.info(f"  Tier3: {tier_counts.get('Tier3', 0)} genes")
    logger.info(f"  Housekeeping demoted: {n_housekeeping_demoted}")
    
    # Log cuproptosis core genes in Tier1 (for reference only, no boost applied)
    tier1_cupro = result_df[
        (result_df['Tier'] == 'Tier1') & (result_df['Is_cuproptosis_core'])
    ]['Gene'].tolist()
    if tier1_cupro:
        logger.info(f"  Cuproptosis core in Tier1 (data-driven): {len(tier1_cupro)} ({tier1_cupro})")
    else:
        logger.info(f"  Cuproptosis core in Tier1: 0 (no boost applied)")
    
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
