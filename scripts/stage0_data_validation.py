# -*- coding: utf-8 -*-
"""
Stage 0: Data Validation Layer (v9.0)
=====================================

NEW in v9.0:
- Species detection from gene identifiers
- Gene ID format validation
- Expression matrix integrity checks
- Cross-species mapping verification

Input:
  - Raw expression matrices
  - Gene lists from various stages
  
Output:
  - validation_report.json
  - Data quality flags
  - Species assignments
"""

import os
import sys
import json
import logging
import re
from collections import Counter
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import RESULTS_DIR

STAGE_DIR = os.path.join(RESULTS_DIR, "stage0_validation")
os.makedirs(STAGE_DIR, exist_ok=True)

logger = logging.getLogger("stage0")
logger.setLevel(logging.INFO)

file_handler = logging.FileHandler(os.path.join(STAGE_DIR, "stage0.log"), encoding="utf-8", mode="w")
file_handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
logger.addHandler(file_handler)

console_handler = logging.StreamHandler(sys.stdout)
console_handler.setFormatter(logging.Formatter("%(message)s"))
logger.addHandler(console_handler)


# Species-specific gene patterns
SPECIES_PATTERNS = {
    'human': {
        'pattern': r'^[A-Z][A-Z0-9]+$',
        'examples': ['TP53', 'BRCA1', 'GAPDH', 'ACTB'],
        'case': 'upper',
    },
    'mouse': {
        'pattern': r'^[A-Z][a-z0-9]+$',
        'examples': ['Trp53', 'Brca1', 'Gapdh', 'Actb'],
        'case': 'title',
    },
    'rat': {
        'pattern': r'^[A-Z][a-z0-9]+$',
        'examples': ['Tp53', 'Brca1', 'Gapdh', 'Actb'],
        'case': 'title',
    },
}

# Common species-specific gene prefixes/suffixes
SPECIES_MARKERS = {
    'human': [],
    'mouse': [],
    'rat': [],
}


class DataValidator:
    """v9.0 Data validation layer"""
    
    def __init__(self):
        self.validation_results = {}
        
    def detect_species(self, gene_list: List[str], sample_size: int = 100) -> Tuple[str, float]:
        """
        Detect species from gene identifier patterns
        
        Returns:
            (species, confidence_score)
        """
        if not gene_list:
            return 'unknown', 0.0
            
        sample = gene_list[:min(sample_size, len(gene_list))]
        
        scores = {}
        for species, patterns in SPECIES_PATTERNS.items():
            matches = 0
            for gene in sample:
                gene_clean = gene.strip()
                if re.match(patterns['pattern'], gene_clean):
                    matches += 1
            scores[species] = matches / len(sample)
        
        best_species = max(scores, key=scores.get)
        confidence = scores[best_species]
        
        logger.info(f"Species detection: {best_species} (confidence: {confidence:.2%})")
        logger.info(f"  Scores: {scores}")
        
        return best_species, confidence
    
    def validate_gene_ids(self, gene_list: List[str], expected_species: str) -> Dict:
        """
        Validate gene IDs match expected species format
        
        Returns:
            Validation report dict
        """
        if expected_species not in SPECIES_PATTERNS:
            logger.warning(f"Unknown expected species: {expected_species}")
            return {'valid': False, 'error': 'unknown_species'}
        
        pattern = SPECIES_PATTERNS[expected_species]['pattern']
        valid_count = sum(1 for g in gene_list if re.match(pattern, g.strip()))
        
        validation_rate = valid_count / len(gene_list) if gene_list else 0
        
        result = {
            'valid': validation_rate >= 0.8,  # 80% threshold
            'validation_rate': validation_rate,
            'total_genes': len(gene_list),
            'valid_genes': valid_count,
            'expected_species': expected_species,
        }
        
        logger.info(f"Gene ID validation: {valid_count}/{len(gene_list)} ({validation_rate:.2%})")
        
        return result
    
    def check_expression_matrix(self, expr_file: str) -> Dict:
        """
        Check expression matrix integrity
        
        Returns:
            Integrity report
        """
        logger.info(f"Checking expression matrix: {expr_file}")
        
        report = {
            'file_exists': False,
            'readable': False,
            'n_genes': 0,
            'n_samples': 0,
            'missing_values': 0,
            'negative_values': False,
            'species_detected': 'unknown',
            'species_confidence': 0.0,
        }
        
        if not os.path.exists(expr_file):
            logger.error(f"File not found: {expr_file}")
            return report
        
        report['file_exists'] = True
        
        try:
            # Try to read as CSV first
            df = pd.read_csv(expr_file, index_col=0)
            report['readable'] = True
        except Exception as e:
            logger.error(f"Cannot read file: {e}")
            return report
        
        report['n_genes'] = len(df)
        report['n_samples'] = len(df.columns)
        report['missing_values'] = df.isnull().sum().sum()
        report['negative_values'] = (df < 0).any().any()
        
        # Detect species from gene IDs
        gene_list = df.index.tolist()
        species, confidence = self.detect_species(gene_list)
        report['species_detected'] = species
        report['species_confidence'] = confidence
        
        logger.info(f"  Genes: {report['n_genes']}, Samples: {report['n_samples']}")
        logger.info(f"  Missing values: {report['missing_values']}")
        logger.info(f"  Negative values: {report['negative_values']}")
        
        return report
    
    def verify_cross_species_mapping(self, original_genes: List[str], 
                                     mapped_genes: List[str],
                                     mapping_file: str) -> Dict:
        """
        Verify cross-species mapping quality
        
        Returns:
            Mapping statistics
        """
        logger.info("Verifying cross-species mapping...")
        
        original_set = set(original_genes)
        mapped_set = set(mapped_genes)
        
        report = {
            'original_count': len(original_set),
            'mapped_count': len(mapped_set),
            'mapping_rate': len(mapped_set) / len(original_set) if original_set else 0,
            'unmapped_genes': list(original_set - mapped_set),
        }
        
        logger.info(f"  Original: {report['original_count']}")
        logger.info(f"  Mapped: {report['mapped_count']}")
        logger.info(f"  Mapping rate: {report['mapping_rate']:.2%}")
        
        return report
    
    def validate_stage_inputs(self) -> Dict:
        """
        Validate inputs for all stages
        
        Returns:
            Complete validation report
        """
        logger.info("=" * 60)
        logger.info("Stage 0: Data Validation (v9.0)")
        logger.info("=" * 60)
        
        report = {
            'stage1_bulk': {},
            'stage2_single_cell': {},
            'stage3_homology': {},
            'overall_status': 'pending',
        }
        
        # Validate Stage1 output
        stage1_file = os.path.join(RESULTS_DIR, "stage1_rma_degs", "deg_results.csv")
        if os.path.exists(stage1_file):
            df = pd.read_csv(stage1_file)
            genes = df['Gene'].tolist() if 'Gene' in df.columns else []
            species, conf = self.detect_species(genes)
            report['stage1_bulk'] = {
                'species_detected': species,
                'confidence': conf,
                'n_degs': len(genes),
            }
        
        # Validate Stage3 output (should be human)
        stage3_file = os.path.join(RESULTS_DIR, "stage3_enrichment", "human_homolog_degs.csv")
        if os.path.exists(stage3_file):
            df = pd.read_csv(stage3_file)
            genes = df['Gene'].tolist() if 'Gene' in df.columns else df.iloc[:, 0].tolist()
            species, conf = self.detect_species(genes)
            report['stage3_homology'] = {
                'species_detected': species,
                'confidence': conf,
                'n_genes': len(genes),
                'expected_species': 'human',
                'validation': self.validate_gene_ids(genes, 'human'),
            }
        
        # Overall status
        if report['stage3_homology'].get('validation', {}).get('valid', False):
            report['overall_status'] = 'passed'
        else:
            report['overall_status'] = 'warning'
        
        return report
    
    def save_report(self, report: Dict):
        """Save validation report"""
        report_file = os.path.join(STAGE_DIR, "validation_report.json")
        with open(report_file, 'w', encoding='utf-8') as f:
            json.dump(report, f, indent=2, ensure_ascii=False)
        logger.info(f"Validation report saved: {report_file}")


def main():
    """Main validation workflow"""
    logger.info("Starting v9.0 data validation...")
    
    validator = DataValidator()
    report = validator.validate_stage_inputs()
    validator.save_report(report)
    
    # Summary
    logger.info("\n" + "=" * 60)
    logger.info("Validation Summary")
    logger.info("=" * 60)
    logger.info(f"Overall status: {report['overall_status'].upper()}")
    
    if report['overall_status'] == 'passed':
        logger.info("✓ All validations passed. Pipeline can proceed.")
    elif report['overall_status'] == 'warning':
        logger.warning("⚠ Some validations failed. Review warnings before proceeding.")
    else:
        logger.error("✗ Critical validation errors. Pipeline halted.")
    
    return report


if __name__ == "__main__":
    main()
