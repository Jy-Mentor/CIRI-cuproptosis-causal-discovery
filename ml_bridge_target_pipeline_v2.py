# -*- coding: utf-8 -*-
"""
多算法融合桥梁靶点预测系统 (v2.0 — 特征工程×基础算法 笛卡尔积)
=================================================================
网络药理学第一阶段：药物-靶点-疾病桥梁基因发现

核心设计（论文范式）：
  特征工程策略 × 基础分类器 = 笛卡尔积组合模型
  → 5-fold 分层 CV + 防泄漏 → 桥梁得分 (P_drug × P_disease) → RRF 集成排名

特征工程策略 (5种)：
  1. raw         — 原始特征 (1,072维)
  2. pca_10      — PCA 降维至 10维
  3. pca_50      — PCA 降维至 50维  
  4. lasso_sel   — LassoCV 特征选择 (~80维)
  5. pls_10      — PLS 降维至 10维 (偏最小二乘)

基础分类器 (7种)：
  L1_LR (liblinear), L2_LR (lbfgs), ElasticNet_SGD, RF, XGBoost, LightGBM, GB

组合总数: 5 × 7 = 35 种模型

v2.0 修复 (vs v1.0):
  - [FIX] 致命bug: 未知基因概率累加错误地除以 N_FOLDS 两次
  - [NEW] 笛卡尔积架构: 特征工程 × 分类器 = 论文级模型组合
  - [NEW] PLS 特征提取策略
  - [FIX] probas 初始化为 0.0 而非 np.nan
  - [FIX] 有标签基因 NaN 处理
  - [FIX] SFS 替换为 SelectFromModel(LassoCV)
  - [IMPROVE] 特征工程对象在每折训练集上严格 fit/transform

输入数据：
  - 增强基因特征矩阵 (1,072 维)
  - 药物靶点列表 (BCP) / 疾病基因列表 (CIRI)

输出：
  - algorithm_metrics.csv         各组合模型 DT/DG 任务的 AUROC/AUPRC
  - ml_bridge_genes_all.csv       所有未知基因的预测详情
  - top20_bridge_genes_ml.csv     Top-20 桥梁基因
  - feature_importance.csv        特征重要性 (树模型/LASSO)
"""

import os
import sys
import time
import warnings
import traceback
from collections import defaultdict

import numpy as np
import pandas as pd

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression, LassoCV, SGDClassifier
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.decomposition import PCA
from sklearn.cross_decomposition import PLSRegression
from sklearn.feature_selection import SelectFromModel
from sklearn.metrics import roc_auc_score, average_precision_score
from scipy.stats import spearmanr

warnings.filterwarnings("ignore")

# ═══════════════════════════════════════════════════════════════════════════════
# 全局配置
# ═══════════════════════════════════════════════════════════════════════════════
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
OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "ml_output_v2")
os.makedirs(OUTPUT_DIR, exist_ok=True)

np.random.seed(SEED)


# ═══════════════════════════════════════════════════════════════════════════════
# 日志工具
# ═══════════════════════════════════════════════════════════════════════════════
def log(msg: str):
    timestamp = time.strftime("%H:%M:%S")
    print(f"[{timestamp}] {msg}", flush=True)


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 1: 数据加载与标签构建
# ═══════════════════════════════════════════════════════════════════════════════
def load_data():
    log("=" * 70)
    log("STEP 1: 数据加载与标签构建")
    log("=" * 70)

    # 1a. 增强特征矩阵
    log(f"加载特征矩阵: {FEATURE_PATH}")
    feat_df = pd.read_csv(FEATURE_PATH)
    feat_df['gene_symbol'] = feat_df['gene_symbol'].str.upper()
    feat_df = feat_df.drop_duplicates(subset='gene_symbol', keep='first')
    feat_df = feat_df.set_index('gene_symbol')

    n_missing = feat_df.isnull().sum().sum()
    if n_missing > 0:
        log(f"  填充 {n_missing} 个缺失值 (列均值)")
        feat_df = feat_df.fillna(feat_df.mean())
        feat_df = feat_df.fillna(0.0)

    all_genes = set(feat_df.index)
    feature_cols = list(feat_df.columns)
    log(f"  特征矩阵: {len(all_genes)} genes × {len(feature_cols)} features")

    # 1b. 药物靶点
    log(f"加载药物靶点: {DRUG_TARGETS_PATH}")
    with open(DRUG_TARGETS_PATH, 'r') as f:
        drug_targets_raw = set(line.strip().upper() for line in f if line.strip())
    drug_targets = drug_targets_raw & all_genes
    log(f"  BCP 药物靶点: {len(drug_targets_raw)} → {len(drug_targets)} 命中")

    # 1c. 疾病基因
    log(f"加载疾病基因: {DISEASE_GENES_PATH}")
    with open(DISEASE_GENES_PATH, 'r') as f:
        disease_genes_raw = set(line.strip().upper() for line in f if line.strip())
    disease_genes = disease_genes_raw & all_genes
    log(f"  CIRI 疾病基因: {len(disease_genes_raw)} → {len(disease_genes)} 命中")

    # 1d. 可选：限制到子图基因
    if os.path.exists(SUBGRAPH_GENES_PATH):
        log(f"加载子图基因: {SUBGRAPH_GENES_PATH}")
        with open(SUBGRAPH_GENES_PATH, 'r') as f:
            sg = set(line.strip().upper() for line in f if line.strip())
            sg.discard('GENE')
        kept = all_genes & sg
        if kept:
            log(f"  子图过滤: {len(all_genes)} → {len(kept)}")
            feat_df = feat_df.loc[list(kept)]
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

    log(f"  标签统计: DT+={n_dt} ({100*n_dt/len(feat_df):.1f}%), "
        f"DG+={n_dg} ({100*n_dg/len(feat_df):.1f}%), "
        f"Both={n_both}, Unknown={n_unknown}")

    return feat_df, feature_cols, all_genes


# ═══════════════════════════════════════════════════════════════════════════════
# 模块级 Builder 函数 (供 build_model_combinations 和 main 的 STEP 7e 共用)
# ═══════════════════════════════════════════════════════════════════════════════

# ── 特征工程策略 ──
def fe_raw():
    """原始特征 — 不做任何变换"""
    return None

def fe_pca_10():
    return PCA(n_components=10, random_state=SEED)

def fe_pca_50():
    return PCA(n_components=50, random_state=SEED)

def fe_lasso_sel():
    """LassoCV 特征选择 — 自动选择非零系数特征"""
    return SelectFromModel(
        LassoCV(cv=2, n_alphas=20, random_state=SEED, max_iter=2000, n_jobs=-1),
        max_features=50
    )

def fe_pls_10():
    """PLS 偏最小二乘 — 同时利用 X 和 y 的协方差进行降维"""
    return PLSRegression(n_components=10, scale=False)

FE_MAP = {
    "raw":       fe_raw,
    "pca_10":    fe_pca_10,
    "pca_50":    fe_pca_50,
    "lasso_sel": fe_lasso_sel,
    "pls_10":    fe_pls_10,
}

# 需要 y 的有监督特征工程 (PCA 不需要)
SUPERVISED_FE = {PLSRegression, SelectFromModel}

# ── 基础分类器 ──
def clf_l1_lr():
    return LogisticRegression(
        penalty='l1', solver='liblinear', C=0.1,
        class_weight='balanced', max_iter=5000, random_state=SEED)

def clf_l2_lr():
    return LogisticRegression(
        penalty='l2', C=1.0,
        class_weight='balanced', max_iter=5000, random_state=SEED)

def clf_elasticnet_lr():
    return SGDClassifier(
        loss='log_loss', penalty='elasticnet', alpha=0.001, l1_ratio=0.5,
        class_weight='balanced', max_iter=2000, random_state=SEED)

def clf_rf():
    return RandomForestClassifier(
        n_estimators=200, class_weight='balanced',
        random_state=SEED, n_jobs=-1)

def clf_gb():
    """GradientBoostingClassifier — class_weight='balanced' 处理样本不平衡"""
    return GradientBoostingClassifier(
        n_estimators=100, learning_rate=0.1,
        class_weight='balanced', random_state=SEED)

def clf_xgb():
    import xgboost as xgb
    return xgb.XGBClassifier(
        n_estimators=200, learning_rate=0.1,
        eval_metric='logloss', random_state=SEED, verbosity=0)

def clf_lgb():
    import lightgbm as lgb
    return lgb.LGBMClassifier(
        n_estimators=200, learning_rate=0.1,
        class_weight='balanced', random_state=SEED, verbose=-1)

CLF_MAP = {
    "L1_LR":        clf_l1_lr,
    "L2_LR":        clf_l2_lr,
    "ElasticNet_LR": clf_elasticnet_lr,
    "RF":           clf_rf,
    "GB":           clf_gb,
}

# XGBoost (可选)
try:
    import xgboost
    CLF_MAP["XGBoost"] = clf_xgb
except ImportError:
    pass

# LightGBM (可选)
try:
    import lightgbm
    CLF_MAP["LightGBM"] = clf_lgb
except ImportError:
    pass


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 2: 构建特征工程 × 基础分类器的笛卡尔积模型组合
# ═══════════════════════════════════════════════════════════════════════════════
def build_model_combinations():
    """
    返回 (combinations, fe_map, clf_map):
      - combinations: list of (model_name, fe_builder_fn, clf_builder_fn)
      - fe_map: {fe_name: fe_builder_fn}
      - clf_map: {clf_name: clf_builder_fn}
    """
    feature_engineerings = list(FE_MAP.items())
    classifiers = list(CLF_MAP.items())

    combinations = []
    for fe_name, fe_builder in feature_engineerings:
        for clf_name, clf_builder in classifiers:
            model_name = f"{fe_name}__{clf_name}"
            combinations.append((model_name, fe_builder, clf_builder))

    log(f"  特征工程策略: {len(feature_engineerings)}")
    log(f"  基础分类器: {len(classifiers)}")
    log(f"  笛卡尔积组合: {len(combinations)} 种模型")
    return combinations, FE_MAP, CLF_MAP


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 3: 严格5折分层CV + 特征工程防泄漏
# ═══════════════════════════════════════════════════════════════════════════════
def run_cv_for_task(X_base, y, gene_index, model_combos, task_name, unknown_mask):
    """
    对某任务 (DT/DG) 运行所有模型组合的 5 折分层 CV。

    防泄漏关键设计:
      1. StandardScaler 在每折训练集上 fit → transform 验证集 + 未知基因
      2. 特征工程 (PCA/PLS/SelectFromModel) 同样在训练集上 fit → transform
      3. 未知基因在每折均预测, 概率直接累加 → 循环结束后统一除以 N_FOLDS

    【v2.0 修复】probas 初始化为 0.0, 累加原始概率, 最后统一除以 N_FOLDS
    """
    n_genes = X_base.shape[0]

    all_probas = {}      # {model_name: np.array(shape=n_genes)}
    metrics = {}         # {model_name: {'auroc': ..., 'auprc': ...}}

    skf = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=SEED)
    labeled_idx = np.where(~unknown_mask)[0]
    y_labeled = y[~unknown_mask]
    unknown_pos = np.where(unknown_mask)[0]
    n_unknown = len(unknown_pos)

    total_models = len(model_combos)

    for m_idx, (model_name, fe_builder, clf_builder) in enumerate(model_combos):
        log(f"  [{task_name}] [{m_idx+1}/{total_models}] {model_name} ...")

        # 初始化: 0.0 (而非 np.nan), 用于累加
        probas = np.zeros(n_genes, dtype=np.float64)
        # 标记有标签基因是否已被赋值 (每个有标签基因只在一折中出现)
        labeled_assigned = np.zeros(n_genes, dtype=bool)

        fold_aurocs = []
        fold_auprcs = []
        fold_success = 0

        for fold, (train_fold_idx, val_fold_idx) in enumerate(skf.split(labeled_idx, y_labeled)):
            train_pos = labeled_idx[train_fold_idx]
            val_pos = labeled_idx[val_fold_idx]

            X_train = X_base[train_pos].copy()
            y_train = y[train_pos]
            X_val = X_base[val_pos].copy()
            y_val = y[val_pos]

            # --- 1. StandardScaler (仅在训练集 fit) ---
            scaler = StandardScaler()
            X_train_s = scaler.fit_transform(X_train)
            X_val_s = scaler.transform(X_val)
            X_unk_s = scaler.transform(X_base[unknown_pos])

            # --- 2. 特征工程 (仅在训练集 fit) ---
            try:
                fe_obj = fe_builder() if fe_builder is not None else None
            except Exception as e:
                log(f"    [SKIP] Fold {fold+1}: 特征工程创建失败: {e}")
                continue

            if fe_obj is not None:
                try:
                    # 有监督方法 (PLS/LassoCV) 需要 y, 无监督方法 (PCA) 不需要
                    if isinstance(fe_obj, tuple(SUPERVISED_FE)):
                        fe_obj.fit(X_train_s, y_train)
                    else:
                        fe_obj.fit(X_train_s)
                    X_train_fe = fe_obj.transform(X_train_s)
                    X_val_fe = fe_obj.transform(X_val_s)
                    X_unk_fe = fe_obj.transform(X_unk_s)
                except Exception as e:
                    log(f"    [WARN] Fold {fold+1}: 特征工程 fit/transform 失败: {e}, 回退到原始特征")
                    X_train_fe = X_train_s
                    X_val_fe = X_val_s
                    X_unk_fe = X_unk_s

                # 处理 PLS 输出的特殊情况 (transform 返回 X_scores 元组)
                if isinstance(fe_obj, PLSRegression):
                    if isinstance(X_train_fe, tuple):
                        X_train_fe = X_train_fe[0]
                    if isinstance(X_val_fe, tuple):
                        X_val_fe = X_val_fe[0]
                    if isinstance(X_unk_fe, tuple):
                        X_unk_fe = X_unk_fe[0]
            else:
                X_train_fe = X_train_s
                X_val_fe = X_val_s
                X_unk_fe = X_unk_s

            # --- 3. 训练分类器 ---
            try:
                clf = clf_builder()
                clf.fit(X_train_fe, y_train)
            except Exception as e:
                log(f"    [SKIP] Fold {fold+1}: 分类器 fit 失败: {e}")
                continue

            # --- 4. 预测 ---
            try:
                # 验证集 (有标签基因) — 直接赋值
                y_val_prob = clf.predict_proba(X_val_fe)[:, 1]
                probas[val_pos] = y_val_prob
                labeled_assigned[val_pos] = True

                # 未知基因 — 累加原始概率 (不除以 N_FOLDS)
                y_unk_prob = clf.predict_proba(X_unk_fe)[:, 1]
                probas[unknown_pos] += y_unk_prob

                fold_success += 1
            except Exception as e:
                log(f"    [WARN] Fold {fold+1}: predict 失败: {e}")
                continue

            # --- 5. 评估 ---
            try:
                fold_aurocs.append(roc_auc_score(y_val, y_val_prob))
                fold_auprcs.append(average_precision_score(y_val, y_val_prob))
            except Exception:
                log(f"    [WARN] Fold {fold+1}: 评估指标计算失败 (可能是单类)")

        # --- 折循环结束: 对未知基因取平均 ---
        if fold_success > 0:
            probas[unknown_pos] /= fold_success
        else:
            log(f"    [ERR] 所有折均失败, 该模型无有效预测")
            all_probas[model_name] = probas
            metrics[model_name] = {'auroc': 0.0, 'auprc': 0.0}
            continue

        # 检查是否有未赋值的标签基因 → 用均值填充
        if not labeled_assigned[~unknown_mask].all():
            n_unassigned = (~labeled_assigned[~unknown_mask]).sum()
            log(f"    [WARN] {n_unassigned} 个标签基因未被任何折覆盖, 用全局均值填充")
            global_mean = probas[~unknown_mask][labeled_assigned[~unknown_mask]].mean()
            unassigned_idx = np.where(~unknown_mask & ~labeled_assigned)[0]
            probas[unassigned_idx] = global_mean

        all_probas[model_name] = probas

        mean_auroc = np.mean(fold_aurocs) if fold_aurocs else 0.0
        mean_auprc = np.mean(fold_auprcs) if fold_auprcs else 0.0
        metrics[model_name] = {'auroc': mean_auroc, 'auprc': mean_auprc}
        log(f"    AUROC={mean_auroc:.4f}, AUPRC={mean_auprc:.4f}  ({fold_success}/{N_FOLDS} folds)")

    return all_probas, metrics


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 4: 桥梁得分 = P_drug × P_disease
# ═══════════════════════════════════════════════════════════════════════════════
def compute_bridge_scores(drug_probas, disease_probas):
    bridge = {}
    all_model_names = set(drug_probas.keys()) & set(disease_probas.keys())
    for model in all_model_names:
        bridge[model] = drug_probas[model] * disease_probas[model]
    return bridge


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 5: RRF 集成排名
# ═══════════════════════════════════════════════════════════════════════════════
def rrf_integrate(bridge_scores, n_genes, unknown_mask):
    """RRF: RRF_score(g) = Σ 1/(k + rank_i(g))"""
    rrf_sum = np.zeros(n_genes, dtype=np.float64)

    for model, scores in bridge_scores.items():
        uk_scores = scores[unknown_mask]
        order = np.argsort(-uk_scores)  # 降序
        ranks = np.empty(len(uk_scores), dtype=np.float64)
        ranks[order] = np.arange(1, len(uk_scores) + 1)
        rrf_sum[unknown_mask] += 1.0 / (RRF_K + ranks)

    return rrf_sum


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 6: 特征重要性
# ═══════════════════════════════════════════════════════════════════════════════
def extract_importance(model_name, clf, feature_cols):
    """提取树模型和线性模型的特征重要性 (在变换后的特征空间)"""
    try:
        if hasattr(clf, 'feature_importances_'):
            imp = clf.feature_importances_
        elif hasattr(clf, 'coef_'):
            imp = np.abs(clf.coef_).flatten()
        else:
            return None

        if len(imp) != len(feature_cols):
            return None

        df = pd.DataFrame({
            'model': model_name,
            'feature': feature_cols,
            'importance': imp,
        }).sort_values('importance', ascending=False).head(30)
        return df
    except Exception:
        return None


# ═══════════════════════════════════════════════════════════════════════════════
# 主流程
# ═══════════════════════════════════════════════════════════════════════════════
def main():
    start_time = time.time()
    log("═" * 70)
    log("  多算法融合桥梁靶点预测系统 v2.0")
    log("  特征工程 × 基础分类器 = 论文级笛卡尔积组合")
    log(f"  随机种子: {SEED}, {N_FOLDS}折CV, RRF_k={RRF_K}")
    log("═" * 70)

    # ── STEP 1: 加载数据 ──
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
    log(f"  总基因: {n_total}, DT+={n_dt}, DG+={n_dg}, Unknown={n_unknown}")

    # ── STEP 2: 构建模型组合 ──
    log("\n" + "=" * 70)
    log("STEP 2: 构建 特征工程 × 分类器 笛卡尔积组合")
    log("=" * 70)
    model_combos, fe_map, clf_map = build_model_combinations()
    model_names = [m[0] for m in model_combos]
    log(f"  共 {len(model_combos)} 种模型组合")

    # ── STEP 3: DT Task CV ──
    log("\n" + "=" * 70)
    log("STEP 3: 药物可靶向性预测 (DT Task)")
    log("=" * 70)
    drug_probas, dt_metrics = run_cv_for_task(
        X, y_dt, gene_index, model_combos, "DT", unknown_mask)

    # ── STEP 4: DG Task CV ──
    log("\n" + "=" * 70)
    log("STEP 4: 疾病相关性预测 (DG Task)")
    log("=" * 70)
    disease_probas, dg_metrics = run_cv_for_task(
        X, y_dg, gene_index, model_combos, "DG", unknown_mask)

    # ── STEP 5: 桥梁得分 ──
    log("\n" + "=" * 70)
    log("STEP 5: 计算桥梁得分 (P_drug × P_disease)")
    log("=" * 70)
    bridge_scores = compute_bridge_scores(drug_probas, disease_probas)
    log(f"  有效组合: {len(bridge_scores)}/{len(model_combos)}")

    # ── STEP 6: RRF 集成 ──
    log("\n" + "=" * 70)
    log("STEP 6: RRF 倒数排名融合")
    log("=" * 70)
    rrf = rrf_integrate(bridge_scores, n_total, unknown_mask)
    unknown_idx = np.where(unknown_mask)[0]
    rrf_sort_order = unknown_idx[np.argsort(-rrf[unknown_mask])]
    rrf_ranks = np.full(n_total, np.nan)
    for rank, idx in enumerate(rrf_sort_order):
        rrf_ranks[idx] = rank + 1

    # ── STEP 7: 输出 ──
    log("\n" + "=" * 70)
    log("STEP 7: 生成输出文件")
    log("=" * 70)

    # 7a. 所有未知基因详情
    df_all = pd.DataFrame({'gene_symbol': gene_index})
    for model in model_names:
        if model in drug_probas:
            df_all[f'{model}_drug_prob'] = drug_probas[model]
        if model in disease_probas:
            df_all[f'{model}_disease_prob'] = disease_probas[model]
        if model in bridge_scores:
            df_all[f'{model}_bridge_score'] = bridge_scores[model]
    df_all['RRF_score'] = rrf
    df_all['final_rank'] = rrf_ranks

    df_unknown = df_all[unknown_mask].sort_values('RRF_score', ascending=False)
    unknown_path = os.path.join(OUTPUT_DIR, "ml_bridge_genes_all.csv")
    df_unknown.to_csv(unknown_path, index=False, encoding='utf-8-sig')
    log(f"  所有未知基因: {unknown_path} ({len(df_unknown)} genes)")

    # 7b. Top-20
    top20 = df_unknown.head(20)
    top20_path = os.path.join(OUTPUT_DIR, "top20_bridge_genes_ml.csv")
    top20.to_csv(top20_path, index=False, encoding='utf-8-sig')
    log(f"  Top-20: {top20_path}")

    log("\n" + "=" * 70)
    log("TOP-10 桥梁基因 (RRF 融合)")
    log("=" * 70)
    for i, (_, row) in enumerate(top20.head(10).iterrows()):
        log(f"  {int(row['final_rank']):>4d}. {row['gene_symbol']:<12s}  RRF={row['RRF_score']:.6f}")

    # 7c. 算法指标
    metrics_rows = []
    for model in model_names:
        if model in dt_metrics:
            metrics_rows.append({
                'model': model, 'task': 'DT',
                'auroc': dt_metrics[model]['auroc'],
                'auprc': dt_metrics[model]['auprc'],
            })
        if model in dg_metrics:
            metrics_rows.append({
                'model': model, 'task': 'DG',
                'auroc': dg_metrics[model]['auroc'],
                'auprc': dg_metrics[model]['auprc'],
            })
    df_metrics = pd.DataFrame(metrics_rows)
    metrics_path = os.path.join(OUTPUT_DIR, "algorithm_metrics.csv")
    df_metrics.to_csv(metrics_path, index=False, encoding='utf-8-sig')
    log(f"  算法指标: {metrics_path}")

    # 7d. 打印最佳模型
    log("\n--- 各特征工程策略的 DT 最佳模型 ---")
    for fe in ["raw", "pca_10", "pca_50", "lasso_sel", "pls_10"]:
        fe_models = [r for r in metrics_rows if r['model'].startswith(fe + "__") and r['task'] == 'DT']
        if fe_models:
            best = max(fe_models, key=lambda x: x['auroc'])
            log(f"  {fe}: {best['model']} AUROC={best['auroc']:.4f}")

    # ── STEP 7e: 特征重要性提取 ──
    log("\n--- 特征重要性 (raw 特征空间, 最佳模型) ---")
    importance_dfs = []
    # DT 最佳模型
    if metrics_rows:
        dt_best_row = max([r for r in metrics_rows if r['task'] == 'DT'], key=lambda x: x['auroc'])
        dg_best_row = max([r for r in metrics_rows if r['task'] == 'DG'], key=lambda x: x['auroc'])
        best_models = {'DT': dt_best_row['model'], 'DG': dg_best_row['model']}
        log(f"  DT最佳: {best_models['DT']} | DG最佳: {best_models['DG']}")

        for task_name, best_model_name in best_models.items():
            try:
                fe_name, clf_name = best_model_name.split("__", 1)
                y_task = y_dt if task_name == 'DT' else y_dg
                labeled = y_task != 0
                X_lbl, y_lbl = X[labeled], y_task[labeled]

                if len(np.unique(y_lbl)) < 2:
                    continue

                scaler = StandardScaler()
                X_s = scaler.fit_transform(X_lbl)

                # 使用 build_model_combinations 返回的 fe_map / clf_map
                fe_builder = fe_map.get(fe_name)
                clf_builder = clf_map.get(clf_name)
                if fe_builder is None or clf_builder is None:
                    log(f"    [WARN] {best_model_name}: 未找到 builder 函数")
                    continue

                fe_obj = fe_builder()
                if fe_obj is not None:
                    if isinstance(fe_obj, tuple(SUPERVISED_FE)):
                        fe_obj.fit(X_s, y_lbl)
                    else:
                        fe_obj.fit(X_s)
                    X_fe = fe_obj.transform(X_s)
                    if isinstance(X_fe, tuple):
                        X_fe = X_fe[0]
                else:
                    X_fe = X_s

                clf = clf_builder()
                clf.fit(X_fe, y_lbl)

                # For raw features, extract importance on original feature names
                if fe_name == "raw":
                    imp_df = extract_importance(best_model_name, clf, feature_cols)
                    if imp_df is not None:
                        imp_df['task'] = task_name
                        importance_dfs.append(imp_df)
                        log(f"    {best_model_name}: Top-5 = {list(imp_df['feature'].head(5))}")
                else:
                    log(f"    {best_model_name}: 特征变换空间({fe_name}), 重要性跳过")
            except Exception as e:
                log(f"    [WARN] {best_model_name} 重要性提取失败: {e}")

    if importance_dfs:
        imp_all = pd.concat(importance_dfs, ignore_index=True)
        imp_path = os.path.join(OUTPUT_DIR, "feature_importance.csv")
        imp_all.to_csv(imp_path, index=False, encoding='utf-8-sig')
        log(f"  特征重要性: {imp_path}")

    # ── STEP 8: GAT 对比 ──
    log("\n" + "=" * 70)
    log("STEP 8: GAT 桥梁基因对比")
    log("=" * 70)
    if os.path.exists(GAT_BRIDGE_PATH):
        try:
            gat_df = pd.read_csv(GAT_BRIDGE_PATH)
            gat_genes = list(gat_df['gene_symbol'].str.upper())
            log(f"  GAT Top-20: {gat_genes}")

            overlap = []
            for g in gat_genes:
                if g in df_unknown['gene_symbol'].values:
                    rank_info = df_unknown[df_unknown['gene_symbol'] == g]
                    rrf_rank = rank_info['final_rank'].values[0]
                    overlap.append((g, int(rrf_rank)))

            n_overlap = len(overlap)
            log(f"  GAT 基因在 ML 候选池中: {n_overlap}/20")
            for g, r in overlap[:5]:
                log(f"    {g}: GAT Top-20, ML RRF rank={r}")

            ml_ranks = {}
            for i, (_, row) in enumerate(df_unknown.iterrows()):
                ml_ranks[row['gene_symbol']] = i + 1

            gat_ranks_list, ml_ranks_list = [], []
            for rank_gat, g in enumerate(gat_genes):
                if g in ml_ranks:
                    gat_ranks_list.append(rank_gat + 1)
                    ml_ranks_list.append(ml_ranks[g])

            if len(gat_ranks_list) >= 5:
                corr, pval = spearmanr(gat_ranks_list, ml_ranks_list)
                log(f"  Spearman ρ = {corr:.4f} (p = {pval:.4f})")
            else:
                log(f"  共同基因不足5个，跳过 Spearman")
        except Exception as e:
            log(f"  [WARN] GAT 对比失败: {e}")
    else:
        log(f"  GAT 桥梁基因文件未找到，跳过")

    # ── 完成 ──
    elapsed = time.time() - start_time
    log("\n" + "=" * 70)
    log(f"[OK] 完成! 总耗时: {elapsed/60:.1f} 分钟")
    log(f"  输出目录: {OUTPUT_DIR}")
    log(f"  模型组合数: {len(model_combos)}")
    log(f"  每组合 5折CV × 2任务 = {len(model_combos)*5*2} 次训练")
    log("=" * 70)


if __name__ == "__main__":
    main()