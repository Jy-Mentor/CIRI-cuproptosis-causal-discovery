# -*- coding: utf-8 -*-
"""
Statistical Tests Module (v9.0)
================================

NEW in v9.0:
- Permutation tests for significance
- Bootstrap confidence intervals
- Benjamini-Hochberg FDR correction
- Cross-validation performance metrics

Reference:
  - Efron & Tibshirani, An Introduction to the Bootstrap, 1994
  - Benjamini & Hochberg, J R Stat Soc B, 1995

Functions:
  - permutation_test: Calculate empirical p-values
  - bootstrap_ci: Bootstrap confidence intervals
  - fdr_correction: Multiple testing correction
  - cv_performance: Cross-validation performance metrics
"""

import os
import sys
import logging
from typing import Dict, List, Tuple, Optional, Callable
import warnings

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.metrics import roc_auc_score, precision_recall_curve, average_precision_score

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

logger = logging.getLogger("statistical_v9")


class StatisticalTests:
    """v9.0 Statistical testing framework"""
    
    def __init__(self, n_permutations: int = 1000, n_bootstrap: int = 1000, 
                 random_state: int = 42):
        """
        Initialize statistical tests
        
        Args:
            n_permutations: Number of permutations for permutation test
            n_bootstrap: Number of bootstrap iterations
            random_state: Random seed for reproducibility
        """
        self.n_permutations = n_permutations
        self.n_bootstrap = n_bootstrap
        self.random_state = random_state
        np.random.seed(random_state)
        
    def permutation_test(self, data1: np.ndarray, data2: np.ndarray,
                        statistic_func: Callable = np.mean,
                        alternative: str = 'two-sided') -> Dict:
        """
        v9.0: Permutation test for comparing two groups
        
        Args:
            data1: First group data
            data2: Second group data
            statistic_func: Function to compute test statistic
            alternative: 'two-sided', 'greater', 'less'
            
        Returns:
            Dictionary with observed_stat, p_value, null_distribution
        """
        logger.info(f"Running permutation test (n={self.n_permutations})...")
        
        # Observed statistic
        observed_stat = statistic_func(data1) - statistic_func(data2)
        
        # Pooled data
        pooled = np.concatenate([data1, data2])
        n1 = len(data1)
        n = len(pooled)
        
        # Permutation
        null_distribution = np.zeros(self.n_permutations)
        
        for i in range(self.n_permutations):
            # Shuffle and split
            np.random.shuffle(pooled)
            perm1 = pooled[:n1]
            perm2 = pooled[n1:]
            
            null_distribution[i] = statistic_func(perm1) - statistic_func(perm2)
        
        # Calculate p-value
        if alternative == 'two-sided':
            p_value = np.mean(np.abs(null_distribution) >= np.abs(observed_stat))
        elif alternative == 'greater':
            p_value = np.mean(null_distribution >= observed_stat)
        elif alternative == 'less':
            p_value = np.mean(null_distribution <= observed_stat)
        else:
            raise ValueError(f"Unknown alternative: {alternative}")
        
        result = {
            'observed_statistic': float(observed_stat),
            'p_value': float(p_value),
            'null_distribution': null_distribution,
            'n_permutations': self.n_permutations,
            'alternative': alternative,
        }
        
        logger.info(f"  Observed: {observed_stat:.4f}, p-value: {p_value:.4e}")
        
        return result
    
    def bootstrap_ci(self, data: np.ndarray, 
                    statistic_func: Callable = np.mean,
                    confidence_level: float = 0.95) -> Dict:
        """
        v9.0: Bootstrap confidence interval
        
        Args:
            data: Input data
            statistic_func: Function to compute statistic
            confidence_level: Confidence level (e.g., 0.95)
            
        Returns:
            Dictionary with statistic, ci_lower, ci_upper, bootstrap_distribution
        """
        logger.info(f"Computing bootstrap CI (n={self.n_bootstrap})...")
        
        n = len(data)
        observed_stat = statistic_func(data)
        
        # Bootstrap
        bootstrap_stats = np.zeros(self.n_bootstrap)
        
        for i in range(self.n_bootstrap):
            # Resample with replacement
            sample = np.random.choice(data, size=n, replace=True)
            bootstrap_stats[i] = statistic_func(sample)
        
        # Percentile CI
        alpha = 1 - confidence_level
        ci_lower = np.percentile(bootstrap_stats, alpha/2 * 100)
        ci_upper = np.percentile(bootstrap_stats, (1 - alpha/2) * 100)
        
        # Bias-corrected CI (optional)
        z0 = stats.norm.ppf(np.mean(bootstrap_stats < observed_stat))
        
        result = {
            'observed_statistic': float(observed_stat),
            'ci_lower': float(ci_lower),
            'ci_upper': float(ci_upper),
            'confidence_level': confidence_level,
            'bootstrap_distribution': bootstrap_stats,
            'bias_correction': float(z0),
        }
        
        logger.info(f"  Statistic: {observed_stat:.4f}")
        logger.info(f"  CI: [{ci_lower:.4f}, {ci_upper:.4f}]")
        
        return result
    
    def fdr_correction(self, p_values: np.ndarray, 
                      method: str = 'benjamini_hochberg') -> Dict:
        """
        v9.0: Multiple testing correction
        
        Args:
            p_values: Array of p-values
            method: 'benjamini_hochberg' or 'bonferroni'
            
        Returns:
            Dictionary with corrected_p_values, rejected
        """
        logger.info(f"FDR correction ({method})...")
        
        p_values = np.array(p_values)
        n = len(p_values)
        
        if method == 'benjamini_hochberg':
            # Sort p-values
            sorted_indices = np.argsort(p_values)
            sorted_pvals = p_values[sorted_indices]
            
            # BH procedure
            corrected = np.zeros(n)
            for i in range(n):
                corrected[sorted_indices[i]] = min(
                    sorted_pvals[i] * n / (i + 1),
                    1.0
                )
            
            # Ensure monotonicity
            for i in range(n-2, -1, -1):
                corrected[sorted_indices[i]] = min(
                    corrected[sorted_indices[i]],
                    corrected[sorted_indices[i+1]]
                )
            
        elif method == 'bonferroni':
            corrected = np.minimum(p_values * n, 1.0)
        else:
            raise ValueError(f"Unknown method: {method}")
        
        rejected = corrected < 0.05
        
        result = {
            'original_p_values': p_values,
            'corrected_p_values': corrected,
            'rejected': rejected,
            'n_tests': n,
            'n_rejected': np.sum(rejected),
            'method': method,
        }
        
        logger.info(f"  Tests: {n}, Rejected: {np.sum(rejected)}")
        
        return result
    
    def cv_performance(self, X: np.ndarray, y: np.ndarray, 
                      model, cv_folds: int = 5) -> Dict:
        """
        v9.0: Cross-validation performance metrics
        
        Args:
            X: Feature matrix
            y: Target vector
            model: Scikit-learn model
            cv_folds: Number of CV folds
            
        Returns:
            Dictionary with performance metrics
        """
        logger.info(f"Cross-validation ({cv_folds}-fold)...")
        
        # Stratified K-Fold
        skf = StratifiedKFold(n_splits=cv_folds, shuffle=True, 
                             random_state=self.random_state)
        
        # Metrics
        auc_scores = []
        ap_scores = []
        fold_results = []
        
        for fold, (train_idx, val_idx) in enumerate(skf.split(X, y)):
            X_train, X_val = X[train_idx], X[val_idx]
            y_train, y_val = y[train_idx], y[val_idx]
            
            # Train
            model.fit(X_train, y_train)
            
            # Predict
            y_pred_proba = model.predict_proba(X_val)[:, 1]
            
            # Metrics
            auc = roc_auc_score(y_val, y_pred_proba)
            ap = average_precision_score(y_val, y_pred_proba)
            
            auc_scores.append(auc)
            ap_scores.append(ap)
            
            fold_results.append({
                'fold': fold + 1,
                'auc': auc,
                'average_precision': ap,
                'n_train': len(train_idx),
                'n_val': len(val_idx),
            })
            
            logger.info(f"  Fold {fold+1}: AUC={auc:.3f}, AP={ap:.3f}")
        
        # Summary statistics
        result = {
            'auc_mean': float(np.mean(auc_scores)),
            'auc_std': float(np.std(auc_scores)),
            'auc_ci': (
                float(np.percentile(auc_scores, 2.5)),
                float(np.percentile(auc_scores, 97.5))
            ),
            'ap_mean': float(np.mean(ap_scores)),
            'ap_std': float(np.std(ap_scores)),
            'fold_results': fold_results,
            'cv_folds': cv_folds,
        }
        
        logger.info(f"  Mean AUC: {result['auc_mean']:.3f} ± {result['auc_std']:.3f}")
        
        return result
    
    def gsea_permutation(self, ranked_genes: List[str], 
                        gene_set: List[str],
                        n_permutations: int = 1000) -> Dict:
        """
        v9.0: GSEA with permutation test
        
        Args:
            ranked_genes: List of genes ranked by statistic
            gene_set: Gene set to test
            n_permutations: Number of permutations
            
        Returns:
            Dictionary with ES, NES, p-value
        """
        logger.info(f"GSEA permutation test (n={n_permutations})...")
        
        # Calculate observed enrichment score
        gene_set = set(gene_set)
        n = len(ranked_genes)
        n_set = len(gene_set)
        
        # Ranking positions
        positions = [i for i, g in enumerate(ranked_genes) if g in gene_set]
        
        if not positions:
            return {'es': 0, 'nes': 0, 'p_value': 1.0}
        
        # Kolmogorov-Smirnov-like statistic
        def calc_es(positions, n):
            if not positions:
                return 0
            
            # ECDF
            ecdf = np.zeros(n)
            for pos in positions:
                ecdf[pos:] += 1
            ecdf = ecdf / len(positions)
            
            # Null ECDF (uniform)
            null_ecdf = np.arange(1, n+1) / n
            
            # Max deviation
            es = np.max(np.abs(ecdf - null_ecdf))
            return es
        
        observed_es = calc_es(positions, n)
        
        # Permutation
        null_es = np.zeros(n_permutations)
        all_genes = list(ranked_genes)
        
        for i in range(n_permutations):
            # Random gene set of same size
            random_set = np.random.choice(all_genes, size=n_set, replace=False)
            random_positions = [j for j, g in enumerate(ranked_genes) if g in random_set]
            null_es[i] = calc_es(random_positions, n)
        
        # Normalize ES
        mean_null = np.mean(null_es)
        std_null = np.std(null_es)
        
        if std_null > 0:
            nes = (observed_es - mean_null) / std_null
        else:
            nes = observed_es
        
        # P-value
        p_value = np.mean(null_es >= observed_es)
        
        result = {
            'es': float(observed_es),
            'nes': float(nes),
            'p_value': float(p_value),
            'n_permutations': n_permutations,
            'n_genes': n,
            'n_gene_set': n_set,
        }
        
        logger.info(f"  ES: {observed_es:.3f}, NES: {nes:.3f}, p: {p_value:.4e}")
        
        return result


def main():
    """Example usage"""
    logger.info("Statistical Tests Module v9.0")
    
    # Example: Permutation test
    np.random.seed(42)
    data1 = np.random.normal(0, 1, 50)
    data2 = np.random.normal(0.5, 1, 50)
    
    stats_test = StatisticalTests(n_permutations=1000)
    result = stats_test.permutation_test(data1, data2)
    
    print(f"Permutation test p-value: {result['p_value']:.4f}")


if __name__ == "__main__":
    main()
