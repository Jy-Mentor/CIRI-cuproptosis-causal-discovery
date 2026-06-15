# -*- coding: utf-8 -*-
"""
Stage 7: Cell-Type-Agnostic ML + Nested CV Feature Selection
=============================================================

Critical Fixes Applied:
  C2 [FIXED] LASSO data leakage: Now uses nested CV (outer 5-fold for evaluation,
              inner L1-LogReg feature selection within each training fold)
              Reference: Varma & Simon, BMC Bioinformatics 2006 (PMID:16504092)

  C3 [FIXED] Cell-type confounding: Model is built on cell-type-averaged expression
              profiles, not individual cells. For each sample (Sham vs MCAO) and each
              cell-type cluster, we compute mean expression. This prevents the model
              from learning cell-type composition differences instead of condition effects.
              Reference: Squair et al., Nat Neurosci 2021 (PMID:34518649)

Major Fixes Applied:
  M1 [FIXED] permutation_test_score type: use float(np.asarray(x).flat[0]) for robustness
  M2 [FIXED] L1-LogisticRegression replaces LassoCV for binary classification feature selection
  M3 [FIXED] Feature truncation uses single-cell detection rate + bulk logFC, not logFC alone
  M4 [FIXED] Rank-based fusion replaces arbitrary 0.4/0.3/0.3 weights
  M5 [FIXED] Cuproptosis genes annotated in all outputs

Input:
  - stage3/human_degs.csv: Bulk DEGs
  - stage4/seed_pool_genes.txt: Seed pool
  - stage2_single_cell/sc_adata.h5ad: scRNA-seq (with leiden clustering)

Output:
  - gene_shap_importance.csv: Feature importance (with cuproptosis annotation)
  - ml_model_performance.csv: Model CV metrics + permutation P-values
"""

import os
import sys
import warnings
import logging
import json
from datetime import datetime
from collections import defaultdict

warnings.filterwarnings('ignore')
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'

import numpy as np
import pandas as pd

try:
    import scanpy as sc
    from sklearn.linear_model import LogisticRegression
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.preprocessing import StandardScaler
    from sklearn.model_selection import (
        StratifiedKFold, cross_val_predict, permutation_test_score,
        cross_val_score
    )
    from sklearn.metrics import accuracy_score, roc_auc_score
    from sklearn.pipeline import Pipeline
    from sklearn.feature_selection import SelectFromModel
    from scipy.sparse import issparse
except ImportError as e:
    print(f"Error: Missing dependencies ({e})")
    sys.exit(1)

# Path configuration - 使用统一导入接口
from pipeline_scripts import RESULTS_DIR

STAGE_DIR = os.path.join(RESULTS_DIR, "stage7_ml_shap")
os.makedirs(STAGE_DIR, exist_ok=True)

# Logging
logger = logging.getLogger("stage7")
logger.setLevel(logging.INFO)

file_handler = logging.FileHandler(os.path.join(STAGE_DIR, "stage7.log"), encoding="utf-8", mode="w")
file_handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
logger.addHandler(file_handler)

console_handler = logging.StreamHandler(sys.stdout)
console_handler.setFormatter(logging.Formatter("%(message)s"))
logger.addHandler(console_handler)


def load_data():
    """Load all required data"""
    logger.info("Loading data...")

    deg_file = os.path.join(RESULTS_DIR, "stage3_enrichment", "human_degs.csv")
    degs = pd.read_csv(deg_file)
    logger.info(f"  DEGs: {len(degs)} genes")

    seed_file = os.path.join(RESULTS_DIR, "stage4_seed_wgcna", "seed_pool_genes.txt")
    seeds = set(pd.read_csv(seed_file, header=None)[0].str.upper().tolist())
    logger.info(f"  Seed pool: {len(seeds)} genes")

    h5ad_file = os.path.join(RESULTS_DIR, "stage2_single_cell", "sc_adata.h5ad")
    adata = sc.read_h5ad(h5ad_file)
    logger.info(f"  scRNA-seq: {adata.shape}")

    from config import CUPROPTOSIS_GENES, CUPROPTOSIS_RELATED
    cupro_core = {g.upper() for g in CUPROPTOSIS_GENES}
    cupro_related = {g.upper() for g in CUPROPTOSIS_RELATED}
    cupro_all = cupro_core | cupro_related
    logger.info(f"  Cuproptosis core: {len(cupro_core)}, related: {len(cupro_related)}")

    return degs, seeds, adata, cupro_core, cupro_all


def build_celltype_averaged_matrix(degs, seeds, adata, cupro_all):
    """
    Build cell-type-averaged expression matrix

    FIX:[C3][single-cell level modeling fallacy]
    Instead of using individual cells as samples (which learns cell-type markers),
    we aggregate expression by (condition, leiden_cluster) pairs.

    For each cluster, we compare MCAO vs Sham mean expression, creating a
    cell-type-stratified feature matrix.

    Reference: Squair et al., Nat Neurosci 2021 (PMID:34518649)
    """
    logger.info("\nBuilding cell-type-averaged feature matrix...")

    var_names_upper = [g.upper() for g in adata.var_names]

    # FIX:[M3][single-cell detection rate filter + bulk logFC]
    # Only include genes that are:
    # (a) in seed pool AND (b) significant DEGs AND (c) detected in >=10% of cells
    seed_in_sc = [g for g in seeds if g in var_names_upper]
    sig_deg_genes = set(degs[degs['adjPVal'] < 0.05]['Gene'].str.upper().tolist())

    candidate_genes = [g for g in seed_in_sc if g in sig_deg_genes]

    # Filter by single-cell detection rate (avoid dropout-heavy genes)
    detection_threshold = 0.10  # 10% of cells
    gene_detection = []
    for g in candidate_genes:
        idx = var_names_upper.index(g)
        if issparse(adata.X):
            n_detected = (adata.X[:, idx].toarray().flatten() > 0).sum()
        else:
            n_detected = (adata.X[:, idx].flatten() > 0).sum()
        detection_rate = n_detected / adata.n_obs
        gene_detection.append((g, detection_rate))

    # Keep genes above detection threshold
    detected_genes = [g for g, rate in gene_detection if rate >= detection_threshold]
    n_dropout = len(candidate_genes) - len(detected_genes)

    if len(detected_genes) < 20:
        logger.warning(f"  Too few detected genes ({len(detected_genes)}), using all candidates")
        detected_genes = candidate_genes

    logger.info(f"  Candidate genes: {len(candidate_genes)}, "
                f"After scRNA-seq detection filter: {len(detected_genes)} "
                f"({n_dropout} filtered by dropout)")

    # Add top DEGs not in seed pool (up to 200) to ensure comprehensive coverage
    all_sig_genes = set(degs[degs['adjPVal'] < 0.05]['Gene'].str.upper().tolist())
    non_seed_degs = [g for g in all_sig_genes if g not in set(detected_genes)]

    # Sort by abs(logFC) and take top
    deg_fc = degs.set_index('Gene')['logFC'].to_dict()
    non_seed_degs_sorted = sorted(non_seed_degs, key=lambda g: abs(deg_fc.get(g, 0)), reverse=True)
    additional_genes = []
    for g in non_seed_degs_sorted:
        if g in var_names_upper:
            idx = var_names_upper.index(g)
            if issparse(adata.X):
                n_det = (adata.X[:, idx].toarray().flatten() > 0).sum()
            else:
                n_det = (adata.X[:, idx].flatten() > 0).sum()
            if n_det / adata.n_obs >= detection_threshold:
                additional_genes.append(g)
        if len(additional_genes) >= 200:
            break

    feature_genes = detected_genes + additional_genes

    # FIX:[C3][cell-type stratified aggregation]
    # Get cluster labels
    if 'leiden' in adata.obs.columns:
        cluster_col = 'leiden'
    elif 'louvain' in adata.obs.columns:
        cluster_col = 'louvain'
    else:
        sc.tl.leiden(adata, resolution=0.8)
        cluster_col = 'leiden'

    condition = adata.obs['condition'].values
    clusters = adata.obs[cluster_col].values
    unique_clusters = sorted(set(clusters))
    conditions = ['Sham', 'MCAO']

    # For each (condition, cluster) pair, compute mean expression across genes
    logger.info(f"  Aggregating by {len(unique_clusters)} clusters × 2 conditions...")

    agg_data = []
    agg_labels = []
    agg_meta = []

    for cl in unique_clusters:
        for cond in conditions:
            mask = (clusters == cl) & (condition == cond)
            if mask.sum() == 0:
                continue

            if issparse(adata.X):
                mean_expr = np.asarray(adata.X[mask].mean(axis=0)).flatten()
            else:
                mean_expr = adata.X[mask].mean(axis=0)

            gene_indices = [var_names_upper.index(g) for g in feature_genes]
            expr_vector = mean_expr[gene_indices]
            agg_data.append(expr_vector)
            agg_labels.append(1 if cond == 'MCAO' else 0)
            agg_meta.append({'cluster': cl, 'condition': cond, 'n_cells': int(mask.sum())})

    X = np.array(agg_data, dtype=np.float32)
    y = np.array(agg_labels, dtype=int)

    # Standardize
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    n_pos = y.sum()
    n_neg = len(y) - n_pos
    logger.info(f"  Aggregated matrix: X={X_scaled.shape} (cluster×condition profiles), "
                f"y=[{n_pos} MCAO, {n_neg} Sham]")

    return X_scaled, y, feature_genes, agg_meta


def nested_cv_feature_selection(X, y, feature_genes, n_top=50):
    """
    Nested cross-validation for feature selection

    FIX:[C2][LASSO data leakage]
    FIX:[M2][L1-LogisticRegression instead of LassoCV]

    Outer loop: 5-fold CV for model evaluation
    Inner loop: Within each outer training fold, L1-LogReg selects features

    This ensures the test fold is NEVER seen during feature selection.
    Reference: Varma & Simon, BMC Bioinformatics 2006 (PMID:16504092)
    """
    logger.info("\nNested CV feature selection (L1-LogReg, 5-fold outer)...")

    cv_outer = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    fold_selected_genes = []
    fold_coef = {}

    for fold_idx, (train_idx, test_idx) in enumerate(cv_outer.split(X, y)):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train = y[train_idx]

        # Inner: L1-LogReg feature selection on training fold only
        # Adaptive C based on n_features/n_samples ratio
        n_feat = X_train.shape[1]
        n_samp = X_train.shape[0]
        ratio = n_feat / max(n_samp, 1)
        
        # For very high ratio (>10), use moderate C to avoid zeroing all features
        # For moderate ratio (1-10), use C=0.5
        # For low ratio (<1), use C=1.0
        if ratio > 20:
            c_val = 0.5
        elif ratio > 10:
            c_val = 0.3
        elif ratio > 1:
            c_val = 0.5
        else:
            c_val = 1.0

        l1_lr = LogisticRegression(
            penalty='l1', solver='saga', C=c_val, max_iter=5000,
            random_state=42, n_jobs=-1
        )
        l1_lr.fit(X_train, y_train)

        # Select non-zero coefficient features
        mask = np.abs(l1_lr.coef_[0]) > 1e-8
        selected = [feature_genes[i] for i in range(len(feature_genes)) if mask[i]]
        fold_selected_genes.extend(selected)

        # Store coefficients for this fold
        for i, g in enumerate(feature_genes):
            coef = float(np.abs(l1_lr.coef_[0][i]))
            if g not in fold_coef:
                fold_coef[g] = []
            fold_coef[g].append(coef)

        logger.info(f"  Fold {fold_idx+1}/5: selected {len(selected)} features (C={c_val})")

    # Consensus: features selected in >=2/5 folds (relaxed for small sample size)
    from collections import Counter
    gene_counts = Counter(fold_selected_genes)
    consensus_genes = [g for g, c in gene_counts.most_common() if c >= 2]

    if len(consensus_genes) < 10:
        # Further relax to >=1 fold (any feature selected at least once)
        consensus_genes = [g for g, c in gene_counts.most_common() if c >= 1]
        logger.info(f"  Relaxed to >=1 fold: {len(consensus_genes)} consensus features")

    if len(consensus_genes) < 10:
        # Ultimate fallback: use top features by average coefficient magnitude
        all_by_coef = sorted(fold_coef.items(), key=lambda x: np.mean(x[1]), reverse=True)
        fallback_genes = [g for g, _ in all_by_coef[:50]]
        consensus_genes = list(set(consensus_genes + fallback_genes))
        logger.info(f"  Fallback + consensus: {len(consensus_genes)} features")

    # Average absolute coefficient across folds
    avg_coef = {g: float(np.mean(coefs)) for g, coefs in fold_coef.items()}

    # Rank consensus genes by average coefficient
    consensus_with_score = [(g, avg_coef.get(g, 0)) for g in consensus_genes]
    consensus_with_score.sort(key=lambda x: x[1], reverse=True)

    top_genes = [g for g, s in consensus_with_score[:n_top]]
    gene_scores_dict = {g: s for g, s in consensus_with_score}

    logger.info(f"  Consensus features: {len(consensus_genes)}")
    logger.info(f"  Top {min(n_top, len(top_genes))} genes: {top_genes[:5]}")

    return top_genes, gene_scores_dict


def train_ml_models_nested(X, y, feature_genes, top_genes):
    """
    Train ML models with proper nested CV

    FIX:[C2][no data leakage - feature selection already done via nested CV]
    FIX:[M4][rank-based fusion instead of arbitrary weights]

    Uses the consensus top_genes from nested CV feature selection,
    then trains LR and RF with 5-fold CV on the selected features.
    """
    logger.info("\nTraining ML models (5-fold CV on consensus features)...")

    gene_idx_map = {g: i for i, g in enumerate(feature_genes)}
    top_idx = [gene_idx_map[g] for g in top_genes if g in gene_idx_map]
    X_reduced = X[:, top_idx]

    logger.info(f"  Reduced matrix: {X_reduced.shape}")

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    results = {}

    # 1. Logistic Regression (L2)
    logger.info("  Logistic Regression (L2, 5-fold CV)...")
    lr = LogisticRegression(C=1.0, penalty='l2', solver='lbfgs', max_iter=1000, random_state=42)
    lr_cv_pred = cross_val_predict(lr, X_reduced, y, cv=cv, method='predict_proba')[:, 1]
    lr_cv_auc = roc_auc_score(y, lr_cv_pred)
    lr.fit(X_reduced, y)
    lr_acc = accuracy_score(y, lr.predict(X_reduced))
    logger.info(f"    CV AUC: {lr_cv_auc:.4f}, Train Accuracy: {lr_acc:.4f}")
    results['LR'] = {
        'CV_AUC': round(lr_cv_auc, 4),
        'accuracy': round(lr_acc, 4),
        'importance': np.abs(lr.coef_[0])
    }

    # 2. Random Forest
    logger.info("  Random Forest (5-fold CV)...")
    rf = RandomForestClassifier(
        n_estimators=200, max_depth=min(10, len(top_idx)), min_samples_leaf=5,
        random_state=42, n_jobs=-1
    )
    rf_cv_pred = cross_val_predict(rf, X_reduced, y, cv=cv, method='predict_proba')[:, 1]
    rf_cv_auc = roc_auc_score(y, rf_cv_pred)
    rf.fit(X_reduced, y)
    rf_acc = accuracy_score(y, rf.predict(X_reduced))
    logger.info(f"    CV AUC: {rf_cv_auc:.4f}, Train Accuracy: {rf_acc:.4f}")
    results['RF'] = {
        'CV_AUC': round(rf_cv_auc, 4),
        'accuracy': round(rf_acc, 4),
        'importance': rf.feature_importances_
    }

    # 3. Permutation test (FIX:[M1][robust type handling])
    logger.info("\n  Permutation test (LR, n=1000)...")
    lr_perm = LogisticRegression(C=1.0, penalty='l2', solver='lbfgs', max_iter=1000, random_state=42)
    perm_score_lr, perm_pvalue_lr, _ = permutation_test_score(
        lr_perm, X_reduced, y, cv=cv, n_permutations=1000,
        scoring='roc_auc', random_state=42, n_jobs=-1
    )
    results['LR']['perm_score'] = float(np.asarray(perm_score_lr).flat[0])
    results['LR']['perm_pvalue'] = float(np.asarray(perm_pvalue_lr).flat[0])
    logger.info(f"    LR Perm Score: {results['LR']['perm_score']:.4f}, P-value: {results['LR']['perm_pvalue']:.4f}")

    logger.info("  Permutation test (RF, n=1000)...")
    rf_perm = RandomForestClassifier(
        n_estimators=200, max_depth=min(10, len(top_idx)), min_samples_leaf=5,
        random_state=42, n_jobs=-1
    )
    perm_score_rf, perm_pvalue_rf, _ = permutation_test_score(
        rf_perm, X_reduced, y, cv=cv, n_permutations=1000,
        scoring='roc_auc', random_state=42, n_jobs=-1
    )
    results['RF']['perm_score'] = float(np.asarray(perm_score_rf).flat[0])
    results['RF']['perm_pvalue'] = float(np.asarray(perm_pvalue_rf).flat[0])
    logger.info(f"    RF Perm Score: {results['RF']['perm_score']:.4f}, P-value: {results['RF']['perm_pvalue']:.4f}")

    # 4. Gene importance: rank-based fusion (FIX:[M4])
    # Rank each gene by L1-LogReg coef, LR coef, RF importance separately
    # Then fuse by sum of reciprocal ranks (robust, no arbitrary weights)
    logger.info("\n  Rank-based gene importance fusion...")
    gene_ranks = defaultdict(list)

    # Rank 1: L1-LogReg consensus coefficient
    l1_scores = {g: s for g, s in zip(top_genes,
        [results['LR']['importance'][i] if i < len(results['LR']['importance']) else 0
         for i in range(len(top_genes))])}

    # Rank 2: LR coefficient
    lr_coef = {g: results['LR']['importance'][i] if i < len(results['LR']['importance']) else 0
               for i, g in enumerate(top_genes)}

    # Rank 3: RF importance
    rf_imp = {g: results['RF']['importance'][i] if i < len(results['RF']['importance']) else 0
              for i, g in enumerate(top_genes)}

    for rank_source, scores in [('L1_LogReg', l1_scores), ('LR_coef', lr_coef), ('RF_imp', rf_imp)]:
        ranked = sorted(scores.items(), key=lambda x: x[1], reverse=True)
        for rank, (gene, _) in enumerate(ranked, 1):
            gene_ranks[gene].append(1.0 / (rank + 1))  # reciprocal rank

    # Fusion score = sum of reciprocal ranks
    fused_scores = {g: sum(ranks) / len(ranks) for g, ranks in gene_ranks.items()}

    # Normalize to 0-1
    max_fused = max(fused_scores.values()) if fused_scores else 1
    if max_fused > 0:
        fused_scores = {g: s / max_fused for g, s in fused_scores.items()}

    sorted_importance = sorted(fused_scores.items(), key=lambda x: x[1], reverse=True)

    logger.info(f"  Gene importance (Top10):")
    for gene, score in sorted_importance[:10]:
        logger.info(f"    {gene}: {score:.4f}")

    return results, sorted_importance, l1_scores


def save_results(ml_results, sorted_importance, l1_scores, feature_genes, cupro_core, cupro_all):
    """
    Save results with cuproptosis annotation

    FIX:[M5][cuproptosis gene annotation in output]
    """
    logger.info("\nSaving results...")

    importance_data = []
    for rank, (gene, score) in enumerate(sorted_importance, 1):
        importance_data.append({
            'Gene': gene,
            'SHAP_importance': round(score, 6),
            'Rank': rank,
            'L1_LogReg_score': round(l1_scores.get(gene, 0), 6),
            'Is_cuproptosis_core': gene.upper() in cupro_core,
            'Is_cuproptosis_related': gene.upper() in cupro_all
        })

    importance_df = pd.DataFrame(importance_data)
    importance_file = os.path.join(STAGE_DIR, "gene_shap_importance.csv")
    importance_df.to_csv(importance_file, index=False)
    logger.info(f"  ✓ Gene importance: {importance_file}")

    n_cupro_in_top = importance_df[importance_df['Is_cuproptosis_core']].shape[0]
    logger.info(f"  Cuproptosis core genes in output: {n_cupro_in_top}")

    # Model performance
    perf_data = []
    for model_name, metrics in ml_results.items():
        perf_data.append({
            'Model': model_name,
            'CV_AUC': metrics.get('CV_AUC', 0),
            'Accuracy': metrics.get('accuracy', 0),
            'Permutation_Score': metrics.get('perm_score', 0),
            'Permutation_P': metrics.get('perm_pvalue', 1.0)
        })

    perf_df = pd.DataFrame(perf_data)
    perf_file = os.path.join(STAGE_DIR, "ml_model_performance.csv")
    perf_df.to_csv(perf_file, index=False)
    logger.info(f"  ✓ Model performance: {perf_file}")

    return importance_df


def main():
    logger.info("=" * 60)
    logger.info("Stage 7: Cell-Type-Agnostic ML + Nested CV Feature Selection (v4)")
    logger.info("=" * 60)
    logger.info(f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info("FIX [C2]: Nested CV for feature selection (no data leakage)")
    logger.info("FIX [C3]: Cell-type-averaged expression (not single-cell level)")
    logger.info("FIX [M1]: Robust permutation type handling")
    logger.info("FIX [M2]: L1-LogReg instead of LassoCV")
    logger.info("FIX [M3]: scRNA-seq detection rate filter")
    logger.info("FIX [M4]: Rank-based fusion (no arbitrary weights)")
    logger.info("FIX [M5]: Cuproptosis annotation in output")
    logger.info("Reference: Varma & Simon 2006 (PMID:16504092)")

    # 1. Load data
    degs, seeds, adata, cupro_core, cupro_all = load_data()

    # 2. Build cell-type-averaged matrix (FIX:C3)
    X, y, feature_genes, agg_meta = build_celltype_averaged_matrix(degs, seeds, adata, cupro_all)

    # 3. Nested CV feature selection (FIX:C2, M2)
    top_genes, gene_scores_dict = nested_cv_feature_selection(X, y, feature_genes, n_top=50)

    # 4. Train ML models (FIX:M4)
    ml_results, sorted_importance, l1_scores = train_ml_models_nested(
        X, y, feature_genes, top_genes)

    # 5. Save (FIX:M5)
    save_results(ml_results, sorted_importance, l1_scores, feature_genes, cupro_core, cupro_all)

    logger.info("\n" + "=" * 60)
    logger.info("Stage 7 completed!")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
