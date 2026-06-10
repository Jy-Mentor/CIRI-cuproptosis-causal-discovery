# -*- coding: utf-8 -*-
"""
多算法融合桥梁靶点预测系统 (v1.0)
=====================================
网络药理学第一阶段：药物-靶点-疾病桥梁基因发现
使用 13 种 ML 算法独立预测，5 折分层 CV，RRF 集成排名

输入数据：
  - 增强基因特征矩阵 (1,072 维)
  - 药物靶点列表 (BCP)
  - 疾病基因列表 (CIRI)
  - (可选) GAT 桥梁基因用于对比

输出：
  - algorithm_metrics.csv      各算法 DT/DG 任务 AUROC/AUPRC
  - ml_bridge_genes_all.csv    所有未知基因的预测详情
  - top20_bridge_genes_ml.csv  Top-20 桥梁基因
  - feature_importance.csv     特征重要性
"""

import os
import sys
import time
import warnings
import traceback
from functools import partial

import numpy as np
import pandas as pd

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.decomposition import PCA
from sklearn.feature_selection import SequentialFeatureSelector
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.utils.class_weight import compute_class_weight
from scipy.stats import spearmanr

warnings.filterwarnings("ignore")

# ───────────────────────────── 全局配置 ─────────────────────────────
SEED = 42
N_FOLDS = 5
RRF_K = 60  # RRF 平滑常数

# 输入数据路径
DATA_DIR = r"D:/反向网络药理学/GAT拓展维度/cache"
FEATURE_PATH = os.path.join(DATA_DIR, "enhanced_gene_features.csv")
DRUG_TARGETS_PATH = r"C:/Users/Jy-Mentor-7/Desktop/GAT/drug_targets.txt"
DISEASE_GENES_PATH = r"C:/Users/Jy-Mentor-7/Desktop/GAT/disease_genes.txt"
SUBGRAPH_GENES_PATH = r"C:/Users/Jy-Mentor-7/Desktop/GAT/subgraph_genes.txt"
GAT_BRIDGE_PATH = r"C:/Users/Jy-Mentor-7/Desktop/GAT/top20_bridge_genes.csv"

# 输出目录
OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "ml_output")
os.makedirs(OUTPUT_DIR, exist_ok=True)

np.random.seed(SEED)

# ───────────────────────────── 日志工具 ─────────────────────────────
def log(msg: str):
    timestamp = time.strftime("%H:%M:%S")
    print(f"[{timestamp}] {msg}", flush=True)

# ───────────────────────────── 1. 数据加载与标签构建 ─────────────────────────────
def load_data():
    """加载特征矩阵、药物靶点、疾病基因，构建标签列"""
    log("=" * 60)
    log("STEP 1: 数据加载与标签构建")
    log("=" * 60)

    # 1a. 加载增强特征矩阵
    log(f"加载特征矩阵: {FEATURE_PATH}")
    feat_df = pd.read_csv(FEATURE_PATH)
    log(f"  特征矩阵: {feat_df.shape[0]} genes × {feat_df.shape[1] - 1} features")

    # 统一基因名为大写，设为索引
    feat_df['gene_symbol'] = feat_df['gene_symbol'].str.upper()
    feat_df = feat_df.drop_duplicates(subset='gene_symbol', keep='first')
    feat_df = feat_df.set_index('gene_symbol')

    # 填充缺失值
    n_missing = feat_df.isnull().sum().sum()
    if n_missing > 0:
        log(f"  发现 {n_missing} 个缺失值，使用列均值填充")
        feat_df = feat_df.fillna(feat_df.mean())
        feat_df = feat_df.fillna(0.0)  # 整列为 NaN 的极端情况

    all_genes = set(feat_df.index)
    feature_cols = list(feat_df.columns)
    log(f"  有效基因数: {len(all_genes)}, 特征维度: {len(feature_cols)}")

    # 1b. 加载药物靶点
    log(f"加载药物靶点: {DRUG_TARGETS_PATH}")
    with open(DRUG_TARGETS_PATH, 'r') as f:
        drug_targets_raw = set(line.strip().upper() for line in f if line.strip())
    drug_targets = drug_targets_raw & all_genes
    log(f"  药物靶点 (BCP): {len(drug_targets_raw)} 原始 → {len(drug_targets)} 命中特征矩阵")

    # 1c. 加载疾病基因
    log(f"加载疾病基因: {DISEASE_GENES_PATH}")
    with open(DISEASE_GENES_PATH, 'r') as f:
        disease_genes_raw = set(line.strip().upper() for line in f if line.strip())
    disease_genes = disease_genes_raw & all_genes
    log(f"  疾病基因 (CIRI): {len(disease_genes_raw)} 原始 → {len(disease_genes)} 命中特征矩阵")

    # 1d. (可选) 限制到 subgraph_genes
    if os.path.exists(SUBGRAPH_GENES_PATH):
        log(f"加载子图基因列表: {SUBGRAPH_GENES_PATH}")
        with open(SUBGRAPH_GENES_PATH, 'r') as f:
            subgraph_genes = set(line.strip().upper() for line in f if line.strip())
            # 跳过标题行
            if 'GENE' in subgraph_genes:
                subgraph_genes.discard('GENE')
        kept_genes = all_genes & subgraph_genes
        log(f"  子图基因过滤: {len(all_genes)} → {len(kept_genes)}")
        feat_df = feat_df.loc[list(kept_genes)]
        all_genes = set(feat_df.index)
        drug_targets = drug_targets & all_genes
        disease_genes = disease_genes & all_genes

    # 1e. 构建标签
    feat_df['is_drug_target'] = feat_df.index.isin(drug_targets).astype(int)
    feat_df['is_disease_gene'] = feat_df.index.isin(disease_genes).astype(int)

    n_dt = feat_df['is_drug_target'].sum()
    n_dg = feat_df['is_disease_gene'].sum()
    n_both = ((feat_df['is_drug_target'] == 1) & (feat_df['is_disease_gene'] == 1)).sum()
    n_unknown = len(feat_df) - n_dt - n_dg + n_both

    log(f"  标签统计:")
    log(f"    is_drug_target=1: {n_dt} ({100*n_dt/len(feat_df):.2f}%)")
    log(f"    is_disease_gene=1: {n_dg} ({100*n_dg/len(feat_df):.2f}%)")
    log(f"    同时为 DT 和 DG: {n_both}")
    log(f"    完全未知基因 (候选池): {n_unknown}")

    return feat_df, feature_cols, all_genes

# ───────────────────────────── 2. 定义算法库 ─────────────────────────────
def build_algorithms():
    """构建 11+ 种分类器，返回 (name, classifier) 列表"""
    algos = []

    # 1. L1 逻辑回归 (LASSO)
    algos.append(("L1_LR_LASSO", LogisticRegression(
        penalty='l1', solver='saga', C=0.1,
        class_weight='balanced', max_iter=5000, random_state=SEED)))

    # 2. L2 逻辑回归 (Ridge)
    algos.append(("L2_LR_Ridge", LogisticRegression(
        penalty='l2', C=1.0,
        class_weight='balanced', max_iter=5000, random_state=SEED)))

    # 3. 弹性网络逻辑回归
    algos.append(("ElasticNet_LR", LogisticRegression(
        penalty='elasticnet', solver='saga', l1_ratio=0.5,
        class_weight='balanced', max_iter=5000, random_state=SEED)))

    # 4. 随机森林
    algos.append(("RandomForest", RandomForestClassifier(
        n_estimators=200, class_weight='balanced',
        random_state=SEED, n_jobs=-1)))

    # 5. 梯度提升机
    algos.append(("GradientBoosting", GradientBoostingClassifier(
        n_estimators=200, learning_rate=0.1, random_state=SEED)))

    # 6. XGBoost
    try:
        import xgboost as xgb
        algos.append(("XGBoost", xgb.XGBClassifier(
            n_estimators=200, learning_rate=0.1,
            eval_metric='logloss', random_state=SEED, verbosity=0)))
    except ImportError:
        log("  [WARN] xgboost 未安装，跳过")

    # 7. LightGBM
    try:
        import lightgbm as lgb
        algos.append(("LightGBM", lgb.LGBMClassifier(
            n_estimators=200, learning_rate=0.1,
            class_weight='balanced', random_state=SEED, verbose=-1)))
    except ImportError:
        log("  [WARN] lightgbm 未安装，跳过")

    # 8. SVM (RBF 核)
    algos.append(("SVM_RBF", SVC(
        kernel='rbf', probability=True,
        class_weight='balanced', random_state=SEED)))

    # 9. glmBoost (GB with deviance loss)
    algos.append(("glmBoost", GradientBoostingClassifier(
        loss='log_loss', n_estimators=200,
        learning_rate=0.1, random_state=SEED)))

    # 10. plsRglm (PCA + LR)
    algos.append(("PCA_LR", "plsRglm_proxy"))  # will be handled separately

    # 11. Stepglm (SFS + LR)
    algos.append(("Stepglm", "stepglm_proxy"))  # will be handled separately

    log(f"构建了 {len(algos)} 种算法 (含 2 种代理算法)")
    return algos

# ───────────────────────────── 3. 严格 5 折分层 CV（防泄漏） ─────────────────────────────
def run_cv_for_task(X_base, y, gene_index, algo_specs, task_name, unknown_mask):
    """
    对某一任务 (DT/DG) 运行所有算法的 5 折分层 CV。
    
    关键防泄漏：
      - StandardScaler 在每折训练集上 fit，然后 transform 验证集和未知基因
      - 未知基因在每折均预测，最后取 5 折平均值
    
    返回:
      - all_probas: dict[algo_name] = np.array(n_genes)  （每个基因的 CV 概率）
      - metrics: dict[algo_name] = {'auroc': ..., 'auprc': ...}
    """
    n_genes = X_base.shape[0]
    n_labeled = (~unknown_mask).sum()

    all_probas = {}
    metrics = {}

    skf = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=SEED)
    labeled_idx = np.where(~unknown_mask)[0]
    y_labeled = y[~unknown_mask]

    for algo_name, clf_raw in algo_specs:
        log(f"  [{task_name}] {algo_name} ...")

        probas = np.full(n_genes, np.nan, dtype=np.float64)
        fold_aurocs = []
        fold_auprcs = []

        # plsRglm: PCA(15) + LR
        is_pls = (clf_raw == "plsRglm_proxy")
        # Stepglm: SFS(20) + LR
        is_sfs = (clf_raw == "stepglm_proxy")

        for fold, (train_idx, val_idx) in enumerate(skf.split(labeled_idx, y_labeled)):
            train_pos = labeled_idx[train_idx]
            val_pos = labeled_idx[val_idx]

            X_train = X_base[train_pos].copy()
            y_train = y[train_pos]
            X_val = X_base[val_pos].copy()
            y_val = y[val_pos]

            # 标准缩放（仅在训练集上 fit）
            scaler = StandardScaler()
            X_train_scaled = scaler.fit_transform(X_train)
            X_val_scaled = scaler.transform(X_val)

            # 对未知基因也进行缩放
            X_unknown_scaled = scaler.transform(X_base[unknown_mask])

            if is_pls:
                # PCA(15) + LR
                pca = PCA(n_components=min(15, X_train_scaled.shape[1]), random_state=SEED)
                X_train_tr = pca.fit_transform(X_train_scaled)
                X_val_tr = pca.transform(X_val_scaled)
                X_unk_tr = pca.transform(X_unknown_scaled)
                clf = LogisticRegression(
                    class_weight='balanced', max_iter=5000, random_state=SEED)
            elif is_sfs:
                # SFS(20) + LR
                n_feats = min(20, X_train_scaled.shape[1])
                sfs = SequentialFeatureSelector(
                    LogisticRegression(class_weight='balanced', max_iter=2000, random_state=SEED),
                    n_features_to_select=n_feats, direction='forward',
                    scoring='roc_auc', cv=2, n_jobs=-1)
                try:
                    sfs.fit(X_train_scaled, y_train)
                    X_train_tr = sfs.transform(X_train_scaled)
                    X_val_tr = sfs.transform(X_val_scaled)
                    X_unk_tr = sfs.transform(X_unknown_scaled)
                except Exception:
                    # Fallback: use PCA
                    pca = PCA(n_components=20, random_state=SEED)
                    X_train_tr = pca.fit_transform(X_train_scaled)
                    X_val_tr = pca.transform(X_val_scaled)
                    X_unk_tr = pca.transform(X_unknown_scaled)
                clf = LogisticRegression(
                    class_weight='balanced', max_iter=5000, random_state=SEED)
            else:
                X_train_tr = X_train_scaled
                X_val_tr = X_val_scaled
                X_unk_tr = X_unknown_scaled
                clf = clf_raw

            try:
                clf.fit(X_train_tr, y_train)
            except Exception as e:
                log(f"    [WARN] Fold {fold+1} fit 失败: {e}")
                continue

            # 验证集预测
            y_val_prob = clf.predict_proba(X_val_tr)[:, 1]
            probas[val_pos] = y_val_prob

            # 未知基因预测
            y_unk_prob = clf.predict_proba(X_unk_tr)[:, 1]
            unknown_pos = np.where(unknown_mask)[0]
            for j, upos in enumerate(unknown_pos):
                if np.isnan(probas[upos]):
                    probas[upos] = y_unk_prob[j] / N_FOLDS
                else:
                    probas[upos] += y_unk_prob[j] / N_FOLDS

            # 评估指标
            try:
                fold_aurocs.append(roc_auc_score(y_val, y_val_prob))
                fold_auprcs.append(average_precision_score(y_val, y_val_prob))
            except Exception:
                pass

        # 检查是否有未知基因未被覆盖（理论上不会，但防御性检查）
        if np.isnan(probas[unknown_mask]).any():
            n_nan = np.isnan(probas[unknown_mask]).sum()
            log(f"    [WARN] {n_nan} 个未知基因未被覆盖，使用全局均值填充")
            global_mean = np.nanmean(probas[~unknown_mask])
            unknown_pos = np.where(unknown_mask)[0]
            for upos in unknown_pos:
                if np.isnan(probas[upos]):
                    probas[upos] = global_mean

        all_probas[algo_name] = probas

        mean_auroc = np.mean(fold_aurocs) if fold_aurocs else 0
        mean_auprc = np.mean(fold_auprcs) if fold_auprcs else 0
        metrics[algo_name] = {'auroc': mean_auroc, 'auprc': mean_auprc}
        log(f"    AUROC={mean_auroc:.4f}, AUPRC={mean_auprc:.4f}")

    return all_probas, metrics

# ───────────────────────────── 4. 计算桥梁得分 ─────────────────────────────
def compute_bridge_scores(drug_probas, disease_probas, algo_names):
    """bridge_score = P_drug × P_disease"""
    bridge = {}
    for algo in algo_names:
        dp = drug_probas.get(algo)
        disp = disease_probas.get(algo)
        if dp is not None and disp is not None:
            bridge[algo] = dp * disp
    return bridge

# ───────────────────────────── 5. RRF 集成排名 ─────────────────────────────
def rrf_integrate(bridge_scores, gene_list, unknown_mask):
    """
    对每种算法的桥梁得分排名，使用 RRF 融合。
    RRF_score(g) = Σ 1/(k + rank_i(g))
    """
    n_genes = len(gene_list)
    rrf_sum = np.zeros(n_genes, dtype=np.float64)

    for algo, scores in bridge_scores.items():
        # 仅对未知基因排名
        unknown_scores = scores[unknown_mask]
        # 降序排名 (最高分 rank=1)
        order = np.argsort(-unknown_scores)  # 降序
        ranks = np.empty(len(unknown_scores), dtype=np.float64)
        ranks[order] = np.arange(1, len(unknown_scores) + 1)
        # 加权倒数
        rrf_sum[unknown_mask] += 1.0 / (RRF_K + ranks)

    return rrf_sum

# ───────────────────────────── 6. 特征重要性 ─────────────────────────────
def extract_feature_importance(algo_name, clf, feature_cols, task_name):
    """提取树模型和线性模型的特征重要性"""
    importances = None
    try:
        if hasattr(clf, 'feature_importances_'):
            importances = clf.feature_importances_
        elif hasattr(clf, 'coef_'):
            importances = np.abs(clf.coef_).flatten()
    except Exception:
        pass

    if importances is None:
        return None

    # 确保长度匹配
    if len(importances) != len(feature_cols):
        return None

    df = pd.DataFrame({
        'algorithm': algo_name,
        'task': task_name,
        'feature': feature_cols,
        'importance': importances,
    })
    return df.sort_values('importance', ascending=False)

# ───────────────────────────── 主流程 ─────────────────────────────
def main():
    start_time = time.time()
    log("多算法融合桥梁靶点预测系统启动")
    log(f"随机种子: {SEED}, 折数: {N_FOLDS}, RRF_k={RRF_K}")

    # ── 1. 加载数据 ──
    feat_df, feature_cols, all_genes = load_data()
    X = feat_df[feature_cols].values.astype(np.float64)
    y_dt = feat_df['is_drug_target'].values.astype(int)
    y_dg = feat_df['is_disease_gene'].values.astype(int)
    gene_index = feat_df.index.tolist()

    unknown_mask = (y_dt == 0) & (y_dg == 0)
    n_total = len(gene_index)
    n_unknown = unknown_mask.sum()
    n_dt = y_dt.sum()
    n_dg = y_dg.sum()
    log(f"总基因: {n_total}, DT+: {n_dt}, DG+: {n_dg}, 未知: {n_unknown}")

    # ── 2. 构建算法库 ──
    log("=" * 60)
    log("STEP 2: 构建算法库")
    log("=" * 60)
    algo_specs = build_algorithms()
    algo_names = [a[0] for a in algo_specs]

    # ── 3. 药物靶点任务 CV ──
    log("=" * 60)
    log("STEP 3: 药物可靶向性预测 (DT Task)")
    log("=" * 60)
    drug_probas, dt_metrics = run_cv_for_task(
        X, y_dt, gene_index, algo_specs, "DT", unknown_mask)

    # ── 4. 疾病基因任务 CV ──
    log("=" * 60)
    log("STEP 4: 疾病相关性预测 (DG Task)")
    log("=" * 60)
    disease_probas, dg_metrics = run_cv_for_task(
        X, y_dg, gene_index, algo_specs, "DG", unknown_mask)

    # ── 5. 计算桥梁得分 ──
    log("=" * 60)
    log("STEP 5: 计算桥梁得分")
    log("=" * 60)
    bridge_scores = compute_bridge_scores(drug_probas, disease_probas, algo_names)
    log(f"  共 {len(bridge_scores)} 种算法产生桥梁得分")

    # ── 6. RRF 集成 ──
    log("=" * 60)
    log("STEP 6: RRF 集成排名")
    log("=" * 60)
    rrf = rrf_integrate(bridge_scores, gene_index, unknown_mask)
    # 对未知基因排序
    unknown_idx = np.where(unknown_mask)[0]
    rrf_sort_order = unknown_idx[np.argsort(-rrf[unknown_mask])]
    rrf_ranks = np.full(n_total, np.nan)
    for rank, idx in enumerate(rrf_sort_order):
        rrf_ranks[idx] = rank + 1

    # ── 7. 构建输出 DataFrame ──
    log("=" * 60)
    log("STEP 7: 构建输出文件")
    log("=" * 60)

    # 7a. 所有未知基因详情
    df_all = pd.DataFrame({'gene_symbol': gene_index})
    for algo in algo_names:
        if algo in drug_probas:
            df_all[f'{algo}_drug_prob'] = drug_probas[algo]
        if algo in disease_probas:
            df_all[f'{algo}_disease_prob'] = disease_probas[algo]
        if algo in bridge_scores:
            df_all[f'{algo}_bridge_score'] = bridge_scores[algo]
    df_all['RRF_score'] = rrf
    df_all['final_rank'] = rrf_ranks
    # 保留仅未知基因
    df_unknown = df_all[unknown_mask].sort_values('RRF_score', ascending=False)
    unknown_path = os.path.join(OUTPUT_DIR, "ml_bridge_genes_all.csv")
    df_unknown.to_csv(unknown_path, index=False, encoding='utf-8-sig')
    log(f"  所有未知基因详情: {unknown_path} ({len(df_unknown)} genes)")

    # 7b. Top-20
    top20 = df_unknown.head(20)[['gene_symbol'] + 
        [c for c in df_unknown.columns if c.startswith('RRF') or c.startswith('final')]]
    # 加入每个算法的桥梁得分
    bridge_cols = [f'{a}_bridge_score' for a in algo_names if f'{a}_bridge_score' in df_unknown.columns]
    top20 = df_unknown.head(20)[['gene_symbol', 'RRF_score', 'final_rank'] + bridge_cols]
    top20_path = os.path.join(OUTPUT_DIR, "top20_bridge_genes_ml.csv")
    top20.to_csv(top20_path, index=False, encoding='utf-8-sig')
    log(f"  Top-20 桥梁基因: {top20_path}")

    # 打印 Top-10
    log("\n" + "=" * 60)
    log("TOP-10 桥梁基因 (RRF 融合)")
    log("=" * 60)
    for i, (_, row) in enumerate(top20.head(10).iterrows()):
        log(f"  {int(row['final_rank']):>4d}. {row['gene_symbol']:<12s}  RRF={row['RRF_score']:.6f}")

    # 7c. 算法评估指标
    metrics_rows = []
    for algo in algo_names:
        if algo in dt_metrics:
            metrics_rows.append({
                'algorithm': algo,
                'task': 'DT',
                'auroc': dt_metrics[algo]['auroc'],
                'auprc': dt_metrics[algo]['auprc'],
            })
        if algo in dg_metrics:
            metrics_rows.append({
                'algorithm': algo,
                'task': 'DG',
                'auroc': dg_metrics[algo]['auroc'],
                'auprc': dg_metrics[algo]['auprc'],
            })
    df_metrics = pd.DataFrame(metrics_rows)
    metrics_path = os.path.join(OUTPUT_DIR, "algorithm_metrics.csv")
    df_metrics.to_csv(metrics_path, index=False, encoding='utf-8-sig')
    log(f"  算法指标: {metrics_path}")

    # ── 8. GAT 对比 ──
    log("=" * 60)
    log("STEP 8: GAT 对比分析")
    log("=" * 60)
    if os.path.exists(GAT_BRIDGE_PATH):
        try:
            gat_df = pd.read_csv(GAT_BRIDGE_PATH)
            gat_genes = list(gat_df['gene_symbol'].str.upper())
            log(f"  GAT Top-20 基因: {gat_genes}")

            # 获取这些基因在我们的未知基因中的 RRF 排名
            overlap = []
            for g in gat_genes:
                if g in df_unknown['gene_symbol'].values:
                    rank_info = df_unknown[df_unknown['gene_symbol'] == g]
                    rrf_rank = rank_info['final_rank'].values[0]
                    overlap.append((g, int(rrf_rank)))
                else:
                    log(f"    GAT基因 {g} 不在未知候选池中")

            n_overlap = len(overlap)
            log(f"  Top-20 交集: {n_overlap}")
            if overlap:
                for g, r in overlap[:5]:
                    log(f"    {g}: GAT Top-20, ML RRF rank={r}")

            # 计算斯皮尔曼等级相关
            # 对于两个列表中共有的基因，比较它们的排名
            ml_ranks = {}
            for i, (_, row) in enumerate(df_unknown.iterrows()):
                ml_ranks[row['gene_symbol']] = i + 1  # 1-based rank

            gat_ranks_list = []
            ml_ranks_list = []
            for rank_gat, g in enumerate(gat_genes):
                if g in ml_ranks:
                    gat_ranks_list.append(rank_gat + 1)
                    ml_ranks_list.append(ml_ranks[g])
            
            if len(gat_ranks_list) >= 5:
                corr, pval = spearmanr(gat_ranks_list, ml_ranks_list)
                log(f"  Spearman ρ = {corr:.4f} (p = {pval:.4f})")
            else:
                log(f"  共同基因不足 5 个，跳过 Spearman 计算")
        except Exception as e:
            log(f"  [WARN] GAT 对比失败: {e}")
    else:
        log(f"  GAT 桥梁基因文件未找到，跳过对比")

    # ── 9. 特征重要性 ──
    log("=" * 60)
    log("STEP 9: 特征重要性提取")
    log("=" * 60)
    # 对每种算法的最终全局模型提取重要性
    all_importances = []
    for algo_name, clf_raw in algo_specs:
        if clf_raw in ("plsRglm_proxy", "stepglm_proxy"):
            continue
        # 在全量有标签数据上训练获取重要性
        for task_name, y_task in [("DT", y_dt), ("DG", y_dg)]:
            try:
                scaler = StandardScaler()
                X_scaled = scaler.fit_transform(X)
                clf = clf_raw
                clf.fit(X_scaled, y_task)
                imp_df = extract_feature_importance(algo_name, clf, feature_cols, task_name)
                if imp_df is not None:
                    all_importances.append(imp_df.head(50))  # Top-50 per model
            except Exception:
                pass

    if all_importances:
        df_imp = pd.concat(all_importances, ignore_index=True)
        imp_path = os.path.join(OUTPUT_DIR, "feature_importance.csv")
        df_imp.to_csv(imp_path, index=False, encoding='utf-8-sig')
        log(f"  特征重要性: {imp_path} ({len(df_imp)} entries)")
    else:
        log("  无特征重要性可提取")

    # ── 完成 ──
    elapsed = time.time() - start_time
    log("\n" + "=" * 60)
    log(f"完成! 总耗时: {elapsed/60:.1f} 分钟")
    log(f"输出目录: {OUTPUT_DIR}")
    log("=" * 60)

if __name__ == "__main__":
    main()