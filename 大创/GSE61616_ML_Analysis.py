#!/usr/bin/env python3
"""
GSE61616 Expression Data Processing
Extract expression values for 8 Hub genes and perform ML analysis
"""

import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore')

def smart_install(package):
    import subprocess
    import sys
    try:
        __import__(package)
    except ImportError:
        subprocess.check_call([sys.executable, "-m", "pip", "install", package, "-q"])

for pkg in ["scikit-learn", "matplotlib", "openpyxl"]:
    smart_install(pkg)

from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LassoCV, LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFECV
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.metrics import roc_auc_score, roc_curve, roc_curve
import matplotlib.pyplot as plt

HUB_GENES_8 = ['Stat5a', 'Hspa5', 'Hmox1', 'Il6', 'Mapk1', 'Rela', 'Stat3', 'Mdm2']

RAT_HUB_GENES_MAPPING = {
    'Stat5a': ['1368231_at'],
    'Hspa5': ['1370283_at'],
    'Hmox1': ['1370080_at'],
    'Il6': ['1369191_at'],
    'Mapk1': ['1367697_at'],
    'Rela': ['1372853_at'],
    'Stat3': ['1370224_at'],
    'Mdm2': ['1383288_at']
}

SAMPLE_GROUPS = {
    'GSM1509422': ('Sham', 0), 'GSM1509423': ('Sham', 1), 'GSM1509424': ('Sham', 2),
    'GSM1509425': ('Sham', 3), 'GSM1509426': ('Sham', 4),
    'GSM1509427': ('Model', 0), 'GSM1509428': ('Model', 1), 'GSM1509429': ('Model', 2),
    'GSM1509430': ('Model', 3), 'GSM1509431': ('Model', 4),
    'GSM1509432': ('XST', 0), 'GSM1509433': ('XST', 1), 'GSM1509434': ('XST', 2),
    'GSM1509435': ('XST', 3), 'GSM1509436': ('XST', 4)
}

def load_series_matrix(file_path):
    """Load GSE61616 series matrix file"""
    print(f"Loading: {file_path}")

    with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()

    header_line = None
    data_start = None
    data_end = None

    for i, line in enumerate(lines):
        if line.startswith('!series_matrix_table_begin'):
            data_start = i + 1
        elif line.startswith('!series_matrix_table_end'):
            data_end = i
            break
        elif '"ID_REF"' in line and header_line is None:
            header_line = i

    print(f"Header at line: {header_line}")
    print(f"Data from line {data_start} to {data_end}")

    headers = lines[header_line].strip().split('\t')
    headers = [h.strip('"') for h in headers]

    data_lines = []
    for i in range(data_start, data_end):
        parts = lines[i].strip().split('\t')
        parts = [p.strip('"') for p in parts]
        data_lines.append(parts)

    df = pd.DataFrame(data_lines, columns=headers)
    df = df.set_index('ID_REF')

    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    print(f"Loaded {len(df)} probes x {len(df.columns)} samples")
    return df

def extract_hub_gene_expression(expr_df, hub_mapping):
    """Extract expression values for Hub genes"""
    print("\nExtracting Hub gene expression...")

    hub_expr = {}

    for gene, probes in hub_mapping.items():
        probe_values = []
        for probe in probes:
            if probe in expr_df.index:
                values = expr_df.loc[probe].values
                probe_values.append(values)

        if probe_values:
            hub_expr[gene] = np.mean(probe_values, axis=0)
            print(f"  {gene}: {len(probe_values)} probe(s) found")
        else:
            print(f"  {gene}: NO probes found!")

    result_df = pd.DataFrame(hub_expr, index=expr_df.columns)
    return result_df

def prepare_ml_data(hub_expr_df, groups):
    """Prepare data for ML analysis (Sham vs Model)"""
    print("\nPreparing ML data (Sham vs Model)...")

    X_data = []
    y_data = []
    sample_ids = []

    for sample_id, (group, _) in groups.items():
        if sample_id in hub_expr_df.index and group in ['Sham', 'Model']:
            X_data.append(hub_expr_df.loc[sample_id].values)
            y_data.append(1 if group == 'Model' else 0)
            sample_ids.append(f"{group}_{sample_id[-4:]}")

    X = pd.DataFrame(X_data, columns=hub_expr_df.columns, index=sample_ids)
    y = np.array(y_data)

    print(f"  X shape: {X.shape}")
    print(f"  y shape: {y.shape}")
    print(f"  Model (case): {sum(y)}, Sham (ctrl): {len(y) - sum(y)}")

    return X, y

def perform_lasso_selection(X, y):
    """LASSO feature selection"""
    print("\n" + "="*60)
    print("Step 1: LASSO Regression Feature Selection")
    print("="*60)

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    lassocv = LassoCV(cv=5, random_state=42, max_iter=10000)
    lassocv.fit(X_scaled, y)

    coefficients = pd.DataFrame({
        'Gene': X.columns,
        'Coefficient': lassocv.coef_,
        'Abs_Coef': np.abs(lassocv.coef_)
    }).sort_values('Abs_Coef', ascending=False)

    selected = coefficients[coefficients['Coefficient'] != 0]['Gene'].tolist()

    print(f"\nLASSO selected {len(selected)} genes:")
    for _, row in coefficients.iterrows():
        status = "✓" if row['Coefficient'] != 0 else "✗"
        print(f"  {status} {row['Gene']}: {row['Coefficient']:.4f}")

    return selected, coefficients

def perform_svm_rfe(X, y):
    """SVM-RFE feature selection"""
    print("\n" + "="*60)
    print("Step 2: SVM-RFE Feature Selection")
    print("="*60)

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    svm_rfe = RFECV(
        estimator=SVC(kernel='linear', random_state=42),
        step=1,
        cv=StratifiedKFold(5, shuffle=True, random_state=42),
        scoring='roc_auc',
        min_features_to_select=1
    )
    svm_rfe.fit(X_scaled, y)

    ranking = pd.DataFrame({
        'Gene': X.columns,
        'Ranking': svm_rfe.ranking_,
        'Selected': svm_rfe.support_
    }).sort_values('Ranking')

    selected = ranking[ranking['Selected']]['Gene'].tolist()

    print(f"\nSVM-RFE selected {len(selected)} genes:")
    for _, row in ranking.iterrows():
        status = "✓" if row['Selected'] else "✗"
        print(f"  {status} {row['Gene']}: Rank {row['Ranking']}")

    print(f"\nOptimal number of features: {svm_rfe.n_features_}")

    return selected, ranking

def perform_rf_selection(X, y):
    """Random Forest feature selection"""
    print("\n" + "="*60)
    print("Step 3: Random Forest Feature Selection")
    print("="*60)

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    rf = RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=-1)
    rf.fit(X_scaled, y)

    importance = pd.DataFrame({
        'Gene': X.columns,
        'Importance': rf.feature_importances_
    }).sort_values('Importance', ascending=False)

    threshold = importance['Importance'].mean()
    selected = importance[importance['Importance'] >= threshold]['Gene'].tolist()

    print(f"\nRF selected {len(selected)} genes (>= mean {threshold:.4f}):")
    for _, row in importance.iterrows():
        status = "✓" if row['Importance'] >= threshold else "✗"
        print(f"  {status} {row['Gene']}: {row['Importance']:.4f}")

    return selected, importance

def find_intersection(lasso_genes, svm_genes, rf_genes):
    """Find intersection of three algorithms"""
    print("\n" + "="*60)
    print("Step 4: Finding Intersection")
    print("="*60)

    intersection = set(lasso_genes) & set(svm_genes) & set(rf_genes)

    print(f"\nLASSO: {sorted(lasso_genes)}")
    print(f"SVM-RFE: {sorted(svm_genes)}")
    print(f"RF: {sorted(rf_genes)}")
    print(f"\nIntersection: {sorted(intersection)}")

    return sorted(intersection)

def perform_roc_validation(X, y, selected_genes):
    """ROC validation"""
    print("\n" + "="*60)
    print("Step 5: ROC Validation")
    print("="*60)

    if len(selected_genes) == 0:
        print("No genes selected!")
        return None, None

    gene_indices = [list(X.columns).index(g) for g in selected_genes if g in X.columns]
    X_selected = X.iloc[:, gene_indices]

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_selected)

    model = LogisticRegression(random_state=42, max_iter=1000)
    cv = StratifiedKFold(5, shuffle=True, random_state=42)
    cv_scores = cross_val_score(model, X_scaled, y, cv=cv, scoring='roc_auc')

    model.fit(X_scaled, y)
    y_pred_proba = model.predict_proba(X_scaled)[:, 1]

    auc = roc_auc_score(y, y_pred_proba)

    print(f"\nROC Results:")
    print(f"  Selected genes: {selected_genes}")
    print(f"  AUC: {auc:.4f}")
    print(f"  CV AUC: {cv_scores.mean():.4f} ± {cv_scores.std():.4f}")
    print(f"  AUC > 0.7: {'✓ PASS' if auc > 0.7 else '✗ FAIL'}")

    fpr, tpr, _ = roc_curve(y, y_pred_proba)

    plt.figure(figsize=(8, 6))
    plt.plot(fpr, tpr, 'b-', linewidth=2, label=f'ROC (AUC = {auc:.3f})')
    plt.plot([0, 1], [0, 1], 'k--', linewidth=1, label='Random')
    plt.xlabel('False Positive Rate', fontsize=12)
    plt.ylabel('True Positive Rate', fontsize=12)
    plt.title('ROC Curve - Intersection Genes\n(Sham vs Model)', fontsize=14, fontweight='bold')
    plt.legend(loc='lower right', fontsize=11)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('ROC_GSE61616.pdf', format='pdf', dpi=300)
    plt.savefig('ROC_GSE61616.png', format='png', dpi=300)
    plt.close()

    return auc, cv_scores

def save_results(X, y, lasso_df, svm_df, rf_df, intersection, auc, cv_scores):
    """Save results"""
    output_file = 'GSE61616_ML_Results.xlsx'

    X_with_group = X.copy()
    X_with_group['Group'] = ['Model' if yy == 1 else 'Sham' for yy in y]
    X_with_group.to_excel(output_file, sheet_name='Expression_Data')

    lasso_df.to_excel(output_file, sheet_name='LASSO_Coefficients', index=False)
    svm_df.to_excel(output_file, sheet_name='SVM_RFE_Ranking', index=False)
    rf_df.to_excel(output_file, sheet_name='RF_Importance', index=False)

    pd.DataFrame({'Intersection_Genes': intersection}).to_excel(
        output_file, sheet_name='Intersection_Genes', index=False)

    if auc is not None:
        pd.DataFrame({
            'Metric': ['AUC', 'CV_AUC_Mean', 'CV_AUC_Std', 'Threshold'],
            'Value': [auc, cv_scores.mean(), cv_scores.std(), 0.7]
        }).to_excel(output_file, sheet_name='ROC_Results', index=False)

    print(f"\nResults saved to: {output_file}")

def main():
    print("="*70)
    print("GSE61616 Expression Data Analysis")
    print("Machine Learning: LASSO + SVM-RFE + RF")
    print("="*70)

    data_file = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\大创\GSE61616_series_matrix.txt"

    expr_df = load_series_matrix(data_file)

    hub_expr = extract_hub_gene_expression(expr_df, RAT_HUB_GENES_MAPPING)

    X, y = prepare_ml_data(hub_expr, SAMPLE_GROUPS)

    print("\nHub Gene Expression Summary:")
    print(hub_expr.describe())

    lasso_genes, lasso_df = perform_lasso_selection(X, y)
    svm_genes, svm_df = perform_svm_rfe(X, y)
    rf_genes, rf_df = perform_rf_selection(X, y)

    intersection = find_intersection(lasso_genes, svm_genes, rf_genes)

    auc, cv_scores = perform_roc_validation(X, y, intersection)

    save_results(X, y, lasso_df, svm_df, rf_df, intersection, auc, cv_scores)

    print("\n" + "="*70)
    print("Analysis Complete!")
    print("="*70)

if __name__ == "__main__":
    main()
