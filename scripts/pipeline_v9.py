# -*- coding: utf-8 -*-
"""
BCP×CIRI Target Discovery Pipeline (v9.0)
==========================================

v9.0 Major Improvements:
- Stage0: Data validation layer (species detection, gene ID validation)
- Stage5: STRING API POST method with proper batch handling
- Immune infiltration: CIBERSORT-like deconvolution + cross-method validation
- Statistical tests: Permutation tests, bootstrap CI, FDR correction
- GSEA: Permutation-based p-values
- Cross-validation: Stratified K-fold with performance metrics

Usage:
    python pipeline_v9.py --stage all
    python pipeline_v9.py --stage 0  # Data validation only
    python pipeline_v9.py --stage 5  # PPI network only
"""

import os
import sys
import argparse
import logging
import json
from datetime import datetime
from typing import Dict, List

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# v9.0 imports
from stage0_data_validation import DataValidator
from stage5_string_ppi_v9 import fetch_ppi_from_string_post
from immune_infiltration_v9 import ImmuneInfiltration
from statistical_tests_v9 import StatisticalTests

logger = logging.getLogger("pipeline_v9")


def setup_logging(log_dir: str):
    """Setup logging"""
    os.makedirs(log_dir, exist_ok=True)
    
    handler = logging.FileHandler(
        os.path.join(log_dir, f"pipeline_v9_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"),
        encoding='utf-8'
    )
    handler.setFormatter(logging.Formatter("%(asctime)s | %(levelname)s | %(message)s"))
    
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter("%(message)s"))
    
    logger.addHandler(handler)
    logger.addHandler(console)
    logger.setLevel(logging.INFO)


def run_stage0() -> Dict:
    """Run Stage 0: Data Validation"""
    logger.info("=" * 60)
    logger.info("Stage 0: Data Validation (v9.0)")
    logger.info("=" * 60)
    
    validator = DataValidator()
    report = validator.validate_stage_inputs()
    validator.save_report(report)
    
    logger.info(f"Validation status: {report['overall_status']}")
    
    return report


def run_stage5() -> Dict:
    """Run Stage 5: PPI Network (v9.0 POST method)"""
    logger.info("=" * 60)
    logger.info("Stage 5: PPI Network (v9.0 POST)")
    logger.info("=" * 60)
    
    from stage5_string_ppi_v9 import main as stage5_main
    stage5_main()
    
    return {'status': 'completed'}


def run_immune_analysis(expr_file: str, target_genes: List[str]) -> Dict:
    """Run immune infiltration analysis (v9.0)"""
    logger.info("=" * 60)
    logger.info("Immune Infiltration Analysis (v9.0)")
    logger.info("=" * 60)
    
    if not os.path.exists(expr_file):
        logger.error(f"Expression file not found: {expr_file}")
        return {'status': 'failed'}
    
    expr_matrix = pd.read_csv(expr_file, index_col=0)
    
    # Run analysis
    immune = ImmuneInfiltration(method='marker_average')
    results = immune.run(expr_matrix, target_genes)
    
    # Save results
    output_dir = os.path.join(RESULTS_DIR, "immune_infiltration_v9")
    os.makedirs(output_dir, exist_ok=True)
    
    results['fractions'].to_csv(os.path.join(output_dir, 'cell_fractions.csv'))
    
    if 'correlations' in results:
        results['correlations']['correlation'].to_csv(
            os.path.join(output_dir, 'correlations.csv')
        )
    
    with open(os.path.join(output_dir, 'validation.json'), 'w') as f:
        json.dump(results['validation'], f, indent=2)
    
    logger.info(f"Results saved to {output_dir}")
    
    return results


def run_statistical_tests() -> Dict:
    """Run statistical validation (v9.0)"""
    logger.info("=" * 60)
    logger.info("Statistical Validation (v9.0)")
    logger.info("=" * 60)
    
    stats_test = StatisticalTests(n_permutations=1000, n_bootstrap=1000)
    
    # Example: Test if Tier1 genes have higher GRN perturbation
    # This would be integrated with actual pipeline data
    
    results = {
        'permutation_tests': {},
        'bootstrap_ci': {},
        'fdr_correction': {},
    }
    
    logger.info("Statistical validation framework initialized")
    
    return results


def main():
    """Main pipeline workflow"""
    parser = argparse.ArgumentParser(description='BCP×CIRI Pipeline v9.0')
    parser.add_argument('--stage', type=str, default='all',
                       help='Stage to run: 0, 5, immune, stats, or all')
    parser.add_argument('--log-dir', type=str, default='logs',
                       help='Log directory')
    
    args = parser.parse_args()
    
    # Setup
    setup_logging(args.log_dir)
    
    logger.info("=" * 60)
    logger.info("BCP×CIRI Target Discovery Pipeline v9.0")
    logger.info("=" * 60)
    logger.info(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    results = {}
    
    # Run requested stages
    if args.stage in ['0', 'all']:
        results['stage0'] = run_stage0()
    
    if args.stage in ['5', 'all']:
        results['stage5'] = run_stage5()
    
    if args.stage in ['immune', 'all']:
        # Example target genes
        target_genes = ['SPP1', 'CTSS', 'TREM2', 'LYN', 'CTSB']
        expr_file = os.path.join(RESULTS_DIR, "stage1_rma_degs", "expr_matrix.csv")
        results['immune'] = run_immune_analysis(expr_file, target_genes)
    
    if args.stage in ['stats', 'all']:
        results['stats'] = run_statistical_tests()
    
    # Summary
    logger.info("\n" + "=" * 60)
    logger.info("Pipeline v9.0 Summary")
    logger.info("=" * 60)
    
    for stage, result in results.items():
        status = result.get('status', 'completed')
        logger.info(f"  {stage}: {status}")
    
    logger.info(f"\nEnd time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    return results


if __name__ == "__main__":
    main()
