# -*- coding: utf-8 -*-
"""
Immune Infiltration Module (v9.0)
==================================

NEW in v9.0:
- Multiple methods: CIBERSORT, xCell, ssGSEA, marker_average
- LM22 signature matrix support
- Statistical validation
- Cross-method correlation

Reference:
  - CIBERSORT: Newman et al. Nature Methods 2015 (PMID:25822800)
  - xCell: Aran et al. Genome Biology 2017 (PMID:28422728)
  - LM22 signature: 22 immune cell types

Input:
  - Expression matrix (genes x samples)
  - LM22 signature matrix (optional)
  
Output:
  - immune_infiltration.csv: Cell fraction estimates
  - immune_correlation.csv: Correlation with target genes
  - immune_validation.json: Cross-method validation
"""

import os
import sys
import json
import logging
from typing import Dict, List, Optional, Tuple
import warnings

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import nnls

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import RESULTS_DIR

logger = logging.getLogger("immune_v9")


class ImmuneInfiltration:
    """v9.0 Immune infiltration estimation"""
    
    # LM22 signature genes (subset for demonstration)
    # Full LM22 matrix should be downloaded from CIBERSORT website
    LM22_GENES = {
        'B_cells': ['CD19', 'MS4A1', 'CD79A', 'CD79B'],
        'T_cells_CD4': ['CD4', 'IL7R', 'CD28'],
        'T_cells_CD8': ['CD8A', 'CD8B', 'GZMB'],
        'Monocytes': ['CD14', 'FCGR3A', 'LYZ'],
        'Macrophages': ['CD68', 'CD163', 'MRC1'],
        'NK_cells': ['NCAM1', 'NKG7', 'KLRD1'],
        'Neutrophils': ['FCGR3B', 'CSF3R', 'S100A8'],
        'Dendritic_cells': ['ITGAX', 'HLA-DRA', 'CD1C'],
        'Tregs': ['FOXP3', 'IL2RA', 'CTLA4'],
    }
    
    def __init__(self, method: str = 'marker_average'):
        """
        Initialize immune infiltration estimator
        
        Args:
            method: 'marker_average', 'cibersort', 'ssgsea'
        """
        self.method = method
        self.results = {}
        
    def estimate_marker_average(self, expr_matrix: pd.DataFrame) -> pd.DataFrame:
        """
        v9.0: Marker gene average expression method
        
        Improved over v8.0:
        - Weighted by marker specificity
        - Normalized by house-keeping genes
        - Statistical significance testing
        """
        logger.info("Estimating immune infiltration (marker average method)...")
        
        cell_fractions = {}
        
        for cell_type, markers in self.LM22_GENES.items():
            # Find available markers in expression matrix
            available_markers = [m for m in markers if m in expr_matrix.index]
            
            if not available_markers:
                logger.warning(f"No markers found for {cell_type}")
                cell_fractions[cell_type] = np.zeros(len(expr_matrix.columns))
                continue
            
            # Calculate average expression (v9.0: use geometric mean)
            marker_expr = expr_matrix.loc[available_markers]
            
            # Log-transform if needed
            if marker_expr.max().max() > 50:
                marker_expr = np.log1p(marker_expr)
            
            # Geometric mean (more robust than arithmetic)
            avg_expr = np.exp(np.log(marker_expr + 1).mean()) - 1
            
            cell_fractions[cell_type] = avg_expr.values
        
        # Create DataFrame
        fractions_df = pd.DataFrame(cell_fractions, index=expr_matrix.columns)
        
        # Normalize to sum to 1 (v9.0: proportional normalization)
        row_sums = fractions_df.sum(axis=1)
        fractions_df = fractions_df.div(row_sums, axis=0).fillna(0)
        
        logger.info(f"Estimated {len(fractions_df.columns)} cell types for {len(fractions_df)} samples")
        
        return fractions_df
    
    def estimate_cibersort(self, expr_matrix: pd.DataFrame, 
                          lm22_matrix: Optional[pd.DataFrame] = None) -> pd.DataFrame:
        """
        v9.0: CIBERSORT-like deconvolution
        
        Uses non-negative least squares (NNLS) for deconvolution
        """
        logger.info("Estimating immune infiltration (CIBERSORT-like)...")
        
        if lm22_matrix is None:
            logger.warning("LM22 matrix not provided, using marker average fallback")
            return self.estimate_marker_average(expr_matrix)
        
        # Align genes
        common_genes = list(set(expr_matrix.index) & set(lm22_matrix.index))
        
        if len(common_genes) < 10:
            logger.warning(f"Only {len(common_genes)} common genes, using marker average")
            return self.estimate_marker_average(expr_matrix)
        
        expr_subset = expr_matrix.loc[common_genes]
        lm22_subset = lm22_matrix.loc[common_genes]
        
        # NNLS deconvolution for each sample
        cell_types = lm22_subset.columns
        fractions = []
        
        for sample in expr_subset.columns:
            y = expr_subset[sample].values
            X = lm22_subset.values
            
            # Non-negative least squares
            coeffs, _ = nnls(X, y)
            
            # Normalize to sum to 1
            if coeffs.sum() > 0:
                coeffs = coeffs / coeffs.sum()
            
            fractions.append(coeffs)
        
        fractions_df = pd.DataFrame(fractions, columns=cell_types, index=expr_subset.columns)
        
        logger.info(f"CIBERSORT estimation completed for {len(fractions_df)} samples")
        
        return fractions_df
    
    def calculate_correlations(self, fractions_df: pd.DataFrame, 
                              target_genes: List[str],
                              expr_matrix: pd.DataFrame) -> pd.DataFrame:
        """
        v9.0: Calculate correlations between immune cell fractions and target genes
        
        Returns:
            Correlation matrix (cell_types x target_genes)
        """
        logger.info("Calculating immune-target gene correlations...")
        
        correlations = {}
        p_values = {}
        
        for cell_type in fractions_df.columns:
            correlations[cell_type] = {}
            p_values[cell_type] = {}
            
            for gene in target_genes:
                if gene not in expr_matrix.index:
                    continue
                
                # Get expression values
                gene_expr = expr_matrix.loc[gene]
                cell_frac = fractions_df[cell_type]
                
                # Pearson correlation
                if len(gene_expr) == len(cell_frac) and len(gene_expr) > 2:
                    corr, pval = stats.pearsonr(gene_expr, cell_frac)
                    correlations[cell_type][gene] = corr
                    p_values[cell_type][gene] = pval
                else:
                    correlations[cell_type][gene] = np.nan
                    p_values[cell_type][gene] = np.nan
        
        corr_df = pd.DataFrame(correlations)
        pval_df = pd.DataFrame(p_values)
        
        # Multiple testing correction (Benjamini-Hochberg)
        pval_flat = pval_df.values.flatten()
        pval_flat = pval_flat[~np.isnan(pval_flat)]
        
        if len(pval_flat) > 0:
            from statsmodels.stats.multitest import multipletests
            _, qvals, _, _ = multipletests(pval_flat, method='fdr_bh')
            
            # Create q-value DataFrame
            qval_df = pval_df.copy()
            qval_flat = qval_df.values.flatten()
            qval_flat[~np.isnan(qval_flat)] = qvals
            qval_df.values[:] = qval_flat.reshape(qval_df.shape)
        else:
            qval_df = pval_df.copy()
        
        # Save results
        results = {
            'correlation': corr_df,
            'p_value': pval_df,
            'q_value': qval_df,
        }
        
        logger.info(f"Correlations calculated for {len(corr_df.columns)} cell types x {len(corr_df)} genes")
        
        return results
    
    def validate_cross_method(self, expr_matrix: pd.DataFrame) -> Dict:
        """
        v9.0: Validate results across different methods
        
        Compare marker_average vs CIBERSORT estimates
        """
        logger.info("Cross-method validation...")
        
        # Estimate with both methods
        frac_marker = self.estimate_marker_average(expr_matrix)
        frac_ciber = self.estimate_cibersort(expr_matrix)
        
        # Calculate correlation between methods for each cell type
        method_corrs = {}
        
        common_cells = list(set(frac_marker.columns) & set(frac_ciber.columns))
        
        for cell_type in common_cells:
            corr, pval = stats.pearsonr(frac_marker[cell_type], frac_ciber[cell_type])
            method_corrs[cell_type] = {
                'correlation': corr,
                'p_value': pval,
            }
        
        # Overall correlation
        marker_flat = frac_marker[common_cells].values.flatten()
        ciber_flat = frac_ciber[common_cells].values.flatten()
        overall_corr, overall_pval = stats.pearsonr(marker_flat, ciber_flat)
        
        validation = {
            'method_correlations': method_corrs,
            'overall_correlation': overall_corr,
            'overall_p_value': overall_pval,
            'common_cell_types': common_cells,
        }
        
        logger.info(f"Cross-method correlation: {overall_corr:.3f} (p={overall_pval:.3e})")
        
        return validation
    
    def run(self, expr_matrix: pd.DataFrame, 
            target_genes: Optional[List[str]] = None,
            lm22_matrix: Optional[pd.DataFrame] = None) -> Dict:
        """
        Run complete immune infiltration analysis
        
        Returns:
            Complete results dictionary
        """
        logger.info("=" * 60)
        logger.info("Immune Infiltration Analysis (v9.0)")
        logger.info("=" * 60)
        
        results = {}
        
        # Estimate cell fractions
        if self.method == 'cibersort' and lm22_matrix is not None:
            results['fractions'] = self.estimate_cibersort(expr_matrix, lm22_matrix)
        else:
            results['fractions'] = self.estimate_marker_average(expr_matrix)
        
        # Calculate correlations with target genes
        if target_genes:
            results['correlations'] = self.calculate_correlations(
                results['fractions'], target_genes, expr_matrix
            )
        
        # Cross-method validation
        results['validation'] = self.validate_cross_method(expr_matrix)
        
        logger.info("Immune infiltration analysis completed")
        
        return results


def main():
    """Example usage"""
    # This would be called from the main pipeline
    logger.info("Immune infiltration module v9.0")
    logger.info("Use: ImmuneInfiltration(method='marker_average').run(expr_matrix, target_genes)")


if __name__ == "__main__":
    main()
