# -*- coding: utf-8 -*-
"""
多算法融合桥梁靶点预测系统 (v4.0 — 四算法共识 + SHAP分层重要性 + 环境暴露 + 外部验证)
============================================================================
网络药理学第一阶段：药物-靶点-疾病桥梁基因发现

v4.0 核心升级 (vs v3.0):
  1. 四算法交集筛选: LASSO/SVM/Tree/Linear 四族各选Top30, 取交集 → 核心共识靶点
     [参考 CdCl2→AS 论文 PMID:42113191]
  2. RF特征重要性三层分层: 原始重要性 + Permutation + SHAP
     [参考暴露-疾病网络毒理框架]
  3. 环境暴露特征: TaRGET II Toxi BPA暴露 FPKM → log2FC 特征拼接
  4. 外部交叉验证: GSE61616 大鼠CIRI DEGs 跨数据集一致性
     [参考 MLGANN 鲁棒性测试]

保持自 v3.0 的: 16 FE × 11 CLF 笛卡尔积 | N_FOLDS=5 | RRF | joblib并行 | GPU检测
"""

import os
import sys
import json
import time
import warnings
import traceback
import glob
from collections import defaultdict

import numpy as np
import pandas as pd

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression, LassoCV, SGDClassifier, PassiveAggressiveClassifier
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier, ExtraTreesClassifier
from sklearn.decomposition import PCA
from sklearn.cross_decomposition import PLSRegression
from sklearn.feature_selection import SelectFromModel, SelectKBest, f_classif
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC
from sklearn.inspection import permutation_importance

from sklearn.metrics import roc_auc_score, average_precision_score
from scipy.stats import spearmanr

# SHAP (可选)
try:
    import shap
    _SHAP_AVAILABLE = True
except ImportError:
    _SHAP_AVAILABLE = False

# mygene (可选, 用于鼠标↔人同源基因映射)
try:
    import mygene
    _MYGENE_AVAILABLE = True
except ImportError:
    _MYGENE_AVAILABLE = False

# 并行化
from joblib import Parallel, delayed

warnings.filterwarnings("ignore")
from sklearn.exceptions import ConvergenceWarning
warnings.filterwarnings("ignore", category=ConvergenceWarning)

# ═══════════════════════════════════════════════════════════════════════════════
# 全局配置
# ═══════════════════════════════════════════════════════════════════════════════
SEED = 42
N_FOLDS = 5
RRF_K = 60
N_JOBS = 4  # 限制并行数，避免 GPU 显存溢出
GPU_ENABLED = True  # v4 GPU 加速 (XGBoost/LightGBM 可在 GPU 上运行)

# ---------------------------------------------------------------------------
# GPU 检测 (同 v3.0, 略去详细注释)
# ---------------------------------------------------------------------------
_XGB_GPU_PARAMS = {}
_XGB_GPU_AVAILABLE = False
_XGB_GPU_INFO = "未检测"

def detect_gpu_xgb(enabled=True):
    global _XGB_GPU_PARAMS, _XGB_GPU_AVAILABLE, _XGB_GPU_INFO
    if not enabled:
        _XGB_GPU_AVAILABLE = False
        _XGB_GPU_PARAMS = {'tree_method': 'hist', 'predictor': 'cpu_predictor'}
        _XGB_GPU_INFO = "显存安全模式: 回退 CPU"
        log(f"  [GPU] XGBoost: {_XGB_GPU_INFO}")
        return
    try:
        import xgboost as xgb
        params = {'tree_method': 'gpu_hist', 'predictor': 'gpu_predictor',
                  'n_estimators': 1, 'verbosity': 0, 'n_jobs': 1}
        X_dummy = np.random.randn(10, 5).astype(np.float32)
        y_dummy = np.random.randint(0, 2, 10)
        clf = xgb.XGBClassifier(**params)
        clf.fit(X_dummy, y_dummy)
        _XGB_GPU_AVAILABLE = True
        _XGB_GPU_PARAMS = {'tree_method': 'gpu_hist', 'predictor': 'gpu_predictor'}
        _XGB_GPU_INFO = "GPU 加速启用"
        log(f"  [GPU] XGBoost GPU 加速: 启用 [OK]")
    except Exception as e:
        _XGB_GPU_AVAILABLE = False
        _XGB_GPU_PARAMS = {'tree_method': 'hist', 'predictor': 'cpu_predictor'}
        _XGB_GPU_INFO = f"回退CPU ({str(e).split(chr(10))[0]})"
        log(f"  [GPU] XGBoost GPU: 不可用 → {_XGB_GPU_INFO}")

_LGB_GPU_PARAMS = {}
_LGB_GPU_AVAILABLE = False
_LGB_GPU_INFO = "未检测"

def detect_gpu_lgb(enabled=True):
    global _LGB_GPU_PARAMS, _LGB_GPU_AVAILABLE, _LGB_GPU_INFO
    if not enabled:
        _LGB_GPU_AVAILABLE = False
        _LGB_GPU_PARAMS = {}
        _LGB_GPU_INFO = "显存安全模式: 回退 CPU"
        log(f"  [GPU] LightGBM: {_LGB_GPU_INFO}")
        return
    try:
        import lightgbm as lgb
        params = {'device': 'gpu', 'gpu_platform_id': 0, 'gpu_device_id': 0,
                  'n_estimators': 1, 'verbose': -1}
        X_dummy = np.random.randn(10, 5).astype(np.float32)
        y_dummy = np.random.randint(0, 2, 10)
        clf = lgb.LGBMClassifier(**params)
        clf.fit(X_dummy, y_dummy)
        _LGB_GPU_AVAILABLE = True
        _LGB_GPU_PARAMS = {'device': 'gpu', 'gpu_platform_id': 0, 'gpu_device_id': 0}
        _LGB_GPU_INFO = "GPU 加速启用"
        log(f"  [GPU] LightGBM GPU 加速: 启用 [OK]")
    except Exception as e:
        _LGB_GPU_AVAILABLE = False
        _LGB_GPU_PARAMS = {}
        _LGB_GPU_INFO = f"回退CPU ({str(e).split(chr(10))[0]})"
        log(f"  [GPU] LightGBM GPU: 不可用 → {_LGB_GPU_INFO}")

def detect_gpu_environment():
    log("-" * 50)
    log("GPU 环境检测")
    log("-" * 50)
    effective_jobs = os.cpu_count() if N_JOBS == -1 else N_JOBS
    log(f"  并行核心: {N_JOBS} → 有效 {effective_jobs} 进程")
    try:
        import torch
        cuda_avail = torch.cuda.is_available()
        if cuda_avail:
            gpu_name = torch.cuda.get_device_name(0)
            total_vram = torch.cuda.get_device_properties(0).total_memory / 1e9
            log(f"  PyTorch CUDA: [OK] ({gpu_name}, {total_vram:.1f}GB VRAM)")
        else:
            log(f"  PyTorch CUDA: 不可用")
    except ImportError:
        log(f"  PyTorch CUDA: 未安装")
    if not GPU_ENABLED:
        log("  [GPU] GPU_ENABLED=False, 所有模型使用 CPU")
        return
    gpu_parallel_threshold = 4
    if effective_jobs > gpu_parallel_threshold:
        log(f"  [WARN] 显存安全: N_JOBS={effective_jobs} > {gpu_parallel_threshold}, GPU回退CPU")
        detect_gpu_xgb(enabled=False)
        detect_gpu_lgb(enabled=False)
    else:
        detect_gpu_xgb(enabled=True)
        detect_gpu_lgb(enabled=True)

# 数据路径
DATA_DIR = r"D:/反向网络药理学/GAT拓展维度/cache"
FEATURE_PATH = os.path.join(DATA_DIR, "enhanced_gene_features.csv")
DRUG_TARGETS_PATH = r"C:/Users/Jy-Mentor-7/Desktop/GAT/drug_targets.txt"
DISEASE_GENES_PATH = r"C:/Users/Jy-Mentor-7/Desktop/GAT/disease_genes.txt"
SUBGRAPH_GENES_PATH = r"C:/Users/Jy-Mentor-7/Desktop/GAT/subgraph_genes.txt"
GAT_BRIDGE_PATH = r"C:/Users/Jy-Mentor-7/Desktop/GAT/top20_bridge_genes.csv"
TOXI_FPKM_DIR = r"D:/反向网络药理学/GAT拓展维度/Toxi/rna_fpkm"
GSE61616_DEG_PATH = r"c:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\大创\GSE61616_human_homologs_DEGs.tsv"

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "ml_output_v4")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# 本地同源映射缓存路径 (避免 mygene/Ensembl API 失败)
ORTHO_CACHE_PATH = os.path.join(OUTPUT_DIR, "mouse_to_human_orthologs.csv")

np.random.seed(SEED)


# ═══════════════════════════════════════════════════════════════════════════════
# 日志
# ═══════════════════════════════════════════════════════════════════════════════
def log(msg: str):
    timestamp = time.strftime("%H:%M:%S")
    print(f"[{timestamp}] {msg}", flush=True)


def _save_ortholog_cache(ensmusg_to_symbol, mouse_to_human):
    """保存 ENSMUSG → 小鼠基因符号 → 人类基因符号 映射到本地 CSV"""
    try:
        rows = []
        for eid, msym in ensmusg_to_symbol.items():
            hsym = mouse_to_human.get(msym, '')
            rows.append({'ensmusg': eid, 'mouse_symbol': msym, 'human_symbol': hsym})
        pd.DataFrame(rows).to_csv(ORTHO_CACHE_PATH, index=False, encoding='utf-8-sig')
        log(f"  [Cache] 已保存 {len(rows)} 条同源映射到 {ORTHO_CACHE_PATH}")
    except Exception as e:
        log(f"  [Cache] 保存失败: {e}")


detect_gpu_environment()


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 1a: 主特征数据加载
# ═══════════════════════════════════════════════════════════════════════════════
def load_data():
    log("=" * 70)
    log("STEP 1a: 主特征数据加载与标签构建")
    log("=" * 70)

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

    log(f"加载药物靶点: {DRUG_TARGETS_PATH}")
    with open(DRUG_TARGETS_PATH, 'r') as f:
        drug_targets_raw = set(line.strip().upper() for line in f if line.strip())
    drug_targets = drug_targets_raw & all_genes
    log(f"  BCP 药物靶点: {len(drug_targets_raw)} → {len(drug_targets)} 命中")

    log(f"加载疾病基因: {DISEASE_GENES_PATH}")
    with open(DISEASE_GENES_PATH, 'r') as f:
        disease_genes_raw = set(line.strip().upper() for line in f if line.strip())
    disease_genes = disease_genes_raw & all_genes
    log(f"  CIRI 疾病基因: {len(disease_genes_raw)} → {len(disease_genes)} 命中")

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
# STEP 1b: Toxi 环境暴露特征提取 [v4.0 NEW]
# ═══════════════════════════════════════════════════════════════════════════════
def extract_toxi_features(drug_targets_human, all_genes_set):
    """
    从 toxirna_enhanced_features.csv 加载 TaRGET II 多源特征
    - toxirna_feature_extractor.py 预计算: FPKM (16 PCA) + ATAC (16 PCA) = 32 维
    - 小鼠基因符号 → 人同源映射 (本地 ORTHO_CACHE_PATH 优先)
    返回: (toxi_feat_df, toxi_cols) 或 (None, []) 如果失败。
    """
    log("\n" + "=" * 70)
    log("STEP 1b: Toxi 环境暴露特征 (TaRGET II FPKM+ATAC, 32维 PCA)")
    log("=" * 70)

    toxi_csv = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "toxirna_enhanced_features.csv")
    if not os.path.exists(toxi_csv):
        log(f"  [SKIP] 未找到预计算特征: {toxi_csv}")
        return None, []

    # ── 加载预计算特征 ──
    feat_raw = pd.read_csv(toxi_csv, index_col=0)  # index = mouse gene symbol
    feat_raw.index = feat_raw.index.str.upper()
    mouse_symbols = set(feat_raw.index)
    toxi_cols = [c for c in feat_raw.columns if c.startswith("toxirna_")]
    log(f"  预计算特征: {len(feat_raw)} 小鼠基因 × {len(toxi_cols)} 维")

    # ── 加载鼠→人同源映射 ──
    mouse_to_human = {}
    if os.path.exists(ORTHO_CACHE_PATH):
        cache_df = pd.read_csv(ORTHO_CACHE_PATH)
        for _, row in cache_df.iterrows():
            msym = str(row.get('mouse_symbol', '')).strip().upper()
            hsym = str(row.get('human_symbol', '')).strip().upper()
            if msym and hsym:
                mouse_to_human[msym] = hsym
        log(f"  同源映射: {len(mouse_to_human)} 对 (本地缓存)")

    # ── 对齐到 all_genes_set ──
    aligned = pd.DataFrame(0.0, index=[g.upper() for g in all_genes_set],
                           columns=toxi_cols, dtype=np.float32)
    matched = 0
    for msym, hsym in mouse_to_human.items():
        if msym in mouse_symbols and hsym.upper() in aligned.index:
            aligned.loc[hsym.upper()] = feat_raw.loc[msym, toxi_cols].values
            matched += 1

    log(f"  映射匹配: {matched}/{len(aligned)} 基因")
    log(f"  特征维度: {len(toxi_cols)} 列 (FPKM PCA×16 + ATAC PCA×16)")
    if matched > 0:
        top5 = aligned.abs().sum(axis=1).sort_values(ascending=False).head(5)
        log(f"  Top-5 非零基因: {list(zip(top5.index, top5.values.round(4)))}")

    return aligned, toxi_cols


# ═══════════════════════════════════════════════════════════════════════════════
# 特征工程 Builder (同 v3.0)
# ═══════════════════════════════════════════════════════════════════════════════
def fe_raw():
    return None

def fe_pca_10():
    return PCA(n_components=10, random_state=SEED)
def fe_pca_20():
    return PCA(n_components=20, random_state=SEED)
def fe_pca_30():
    return PCA(n_components=30, random_state=SEED)
def fe_pca_50():
    return PCA(n_components=50, random_state=SEED)
def fe_pca_80():
    return PCA(n_components=80, random_state=SEED)
def fe_pca_100():
    return PCA(n_components=100, random_state=SEED)

def fe_lasso_sel(max_features=50):
    return SelectFromModel(
        LassoCV(cv=2, n_alphas=20, random_state=SEED, max_iter=2000, n_jobs=-1),
        max_features=max_features)
def fe_lasso_20():
    return fe_lasso_sel(20)
def fe_lasso_30():
    return fe_lasso_sel(30)
def fe_lasso_50():
    return fe_lasso_sel(50)
def fe_lasso_80():
    return fe_lasso_sel(80)

def fe_pls_5():
    return PLSRegression(n_components=5, scale=False)
def fe_pls_10():
    return PLSRegression(n_components=10, scale=False)
def fe_pls_15():
    return PLSRegression(n_components=15, scale=False)

def fe_kbest_50():
    return SelectKBest(f_classif, k=50)
def fe_kbest_100():
    return SelectKBest(f_classif, k=100)

FE_MAP = {
    "raw": fe_raw,
    "pca_10": fe_pca_10, "pca_20": fe_pca_20, "pca_30": fe_pca_30,
    "pca_50": fe_pca_50, "pca_80": fe_pca_80, "pca_100": fe_pca_100,
    "lasso_20": fe_lasso_20, "lasso_30": fe_lasso_30,
    "lasso_50": fe_lasso_50, "lasso_80": fe_lasso_80,
    "pls_5": fe_pls_5, "pls_10": fe_pls_10, "pls_15": fe_pls_15,
    "kbest_50": fe_kbest_50, "kbest_100": fe_kbest_100,
}

SUPERVISED_FE = {PLSRegression, SelectFromModel, SelectKBest}


# ═══════════════════════════════════════════════════════════════════════════════
# 分类器 Builder (同 v3.0)
# ═══════════════════════════════════════════════════════════════════════════════
def clf_l1_lr():
    return LogisticRegression(penalty='l1', solver='liblinear', C=0.1,
        class_weight='balanced', max_iter=5000, random_state=SEED)
def clf_l2_lr():
    return LogisticRegression(penalty='l2', C=1.0,
        class_weight='balanced', max_iter=5000, random_state=SEED)
def clf_elasticnet_lr():
    return SGDClassifier(loss='log_loss', penalty='elasticnet', alpha=0.001,
        l1_ratio=0.5, class_weight='balanced', max_iter=2000, random_state=SEED)
def clf_rf():
    return RandomForestClassifier(n_estimators=200, class_weight='balanced',
        random_state=SEED, n_jobs=1)
def clf_gb():
    return GradientBoostingClassifier(n_estimators=100, learning_rate=0.1,
        random_state=SEED)
def clf_xgb():
    import xgboost as xgb
    params = {'n_estimators': 200, 'learning_rate': 0.1,
        'eval_metric': 'logloss', 'random_state': SEED, 'verbosity': 0, 'n_jobs': 1}
    params.update(_XGB_GPU_PARAMS)
    return xgb.XGBClassifier(**params)
def clf_lgb():
    import lightgbm as lgb
    params = {'n_estimators': 200, 'learning_rate': 0.1,
        'class_weight': 'balanced', 'random_state': SEED, 'verbose': -1}
    params.update(_LGB_GPU_PARAMS)
    return lgb.LGBMClassifier(**params)
def clf_nb():
    return GaussianNB()
def clf_svc():
    from sklearn.calibration import CalibratedClassifierCV
    base_svc = SVC(C=1.0, kernel='rbf', class_weight='balanced',
                   random_state=SEED, max_iter=5000, probability=False)
    return CalibratedClassifierCV(base_svc, cv=3, method='sigmoid')
def clf_extratrees():
    return ExtraTreesClassifier(n_estimators=200, class_weight='balanced',
        random_state=SEED, n_jobs=1)
def clf_pac():
    from sklearn.calibration import CalibratedClassifierCV
    base_pac = PassiveAggressiveClassifier(
        C=0.1, class_weight='balanced', max_iter=2000, random_state=SEED)
    return CalibratedClassifierCV(base_pac, cv=3, method='sigmoid')

CLF_MAP = {
    "L1_LR": clf_l1_lr, "L2_LR": clf_l2_lr, "ElasticNet_LR": clf_elasticnet_lr,
    "RF": clf_rf, "GB": clf_gb, "NB": clf_nb, "SVC": clf_svc,
    "ExtraTrees": clf_extratrees, "PAC": clf_pac,
}
try:
    import xgboost
    CLF_MAP["XGBoost"] = clf_xgb
except ImportError:
    pass
try:
    import lightgbm
    CLF_MAP["LightGBM"] = clf_lgb
except ImportError:
    pass

NO_CLASS_WEIGHT_CLFS = {GradientBoostingClassifier, GaussianNB}


def compute_sample_weight(y):
    n_pos = y.sum()
    n_neg = len(y) - n_pos
    return np.where(y == 1, len(y) / (2 * max(n_pos, 1)), len(y) / (2 * max(n_neg, 1)))


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 2: 构建模型组合 (同 v3.0)
# ═══════════════════════════════════════════════════════════════════════════════
def build_model_combinations():
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
# STEP 3-4: 并行 CV (同 v3.0)
# ═══════════════════════════════════════════════════════════════════════════════
def train_single_model(model_name, fe_builder, clf_builder, X_base, y,
                       unknown_mask, model_idx):
    warnings.filterwarnings("ignore", category=ConvergenceWarning)
    warnings.filterwarnings("ignore", category=FutureWarning)
    warnings.filterwarnings("ignore", message="X does not have valid feature names")

    seed = SEED + model_idx
    n_genes = X_base.shape[0]
    probas = np.zeros(n_genes, dtype=np.float64)
    labeled_assigned = np.zeros(n_genes, dtype=bool)

    skf = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=seed)
    labeled_idx = np.where(~unknown_mask)[0]
    y_labeled = y[~unknown_mask]
    unknown_pos = np.where(unknown_mask)[0]

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

        scaler = StandardScaler()
        X_train_s = scaler.fit_transform(X_train)
        X_val_s = scaler.transform(X_val)
        X_unk_s = scaler.transform(X_base[unknown_pos])

        try:
            fe_obj = fe_builder() if fe_builder is not None else None
        except Exception:
            fe_obj = None

        if fe_obj is not None:
            try:
                if isinstance(fe_obj, tuple(SUPERVISED_FE)):
                    fe_obj.fit(X_train_s, y_train)
                else:
                    fe_obj.fit(X_train_s)
                X_train_fe = fe_obj.transform(X_train_s)
                X_val_fe = fe_obj.transform(X_val_s)
                X_unk_fe = fe_obj.transform(X_unk_s)
            except Exception:
                X_train_fe = X_train_s
                X_val_fe = X_val_s
                X_unk_fe = X_unk_s

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

        try:
            clf = clf_builder()
            if isinstance(clf, tuple(NO_CLASS_WEIGHT_CLFS)):
                clf.fit(X_train_fe, y_train, sample_weight=compute_sample_weight(y_train))
            else:
                clf.fit(X_train_fe, y_train)
        except Exception:
            continue

        try:
            y_val_prob = clf.predict_proba(X_val_fe)[:, 1]
            probas[val_pos] = y_val_prob
            labeled_assigned[val_pos] = True
            y_unk_prob = clf.predict_proba(X_unk_fe)[:, 1]
            probas[unknown_pos] += y_unk_prob
            fold_success += 1
        except Exception:
            continue

        try:
            fold_aurocs.append(roc_auc_score(y_val, y_val_prob))
            fold_auprcs.append(average_precision_score(y_val, y_val_prob))
        except Exception:
            pass

    if fold_success > 0:
        probas[unknown_pos] /= fold_success
    else:
        return model_name, probas, {'auroc': 0.0, 'auprc': 0.0}

    if not labeled_assigned[~unknown_mask].all():
        global_mean = probas[~unknown_mask][labeled_assigned[~unknown_mask]].mean()
        unassigned_idx = np.where(~unknown_mask & ~labeled_assigned)[0]
        probas[unassigned_idx] = global_mean

    mean_auroc = np.mean(fold_aurocs) if fold_aurocs else 0.0
    mean_auprc = np.mean(fold_auprcs) if fold_auprcs else 0.0
    return model_name, probas, {'auroc': mean_auroc, 'auprc': mean_auprc}


def run_cv_parallel(X_base, y, model_combos, task_name, unknown_mask):
    n_models = len(model_combos)
    log(f"  [{task_name}] 启动 {n_models} 个模型的并行训练 ({N_JOBS} 核)...")
    results = Parallel(n_jobs=N_JOBS, verbose=5)(
        delayed(train_single_model)(name, fe, clf, X_base, y, unknown_mask, idx)
        for idx, (name, fe, clf) in enumerate(model_combos)
    )
    all_probas = {}
    all_metrics = {}
    for model_name, probas, metric in results:
        all_probas[model_name] = probas
        all_metrics[model_name] = metric
        log(f"  [{task_name}] {model_name}: AUROC={metric['auroc']:.4f}, AUPRC={metric['auprc']:.4f}")
    return all_probas, all_metrics


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 5-6: 桥梁得分 & RRF (同 v3.0)
# ═══════════════════════════════════════════════════════════════════════════════
def compute_bridge_scores(drug_probas, disease_probas):
    bridge = {}
    all_model_names = set(drug_probas.keys()) & set(disease_probas.keys())
    for model in all_model_names:
        bridge[model] = drug_probas[model] * disease_probas[model]
    return bridge

def rrf_integrate(bridge_scores, n_genes, unknown_mask):
    rrf_sum = np.zeros(n_genes, dtype=np.float64)
    for model, scores in bridge_scores.items():
        uk_scores = scores[unknown_mask]
        order = np.argsort(-uk_scores)
        ranks = np.empty(len(uk_scores), dtype=np.float64)
        ranks[order] = np.arange(1, len(uk_scores) + 1)
        rrf_sum[unknown_mask] += 1.0 / (RRF_K + ranks)
    return rrf_sum


# ═══════════════════════════════════════════════════════════════════════════════
# [v4.0 NEW] STEP 7: 四算法家族交集共识筛选
#   参考: CdCl2→AS 论文 (PMID:42113191) 四算法交集策略
# ═══════════════════════════════════════════════════════════════════════════════
ALGORITHM_FAMILIES = {
    "LASSO":  ["L1_LR", "ElasticNet_LR"],
    "SVM":    ["SVC"],
    "Tree":   ["RF", "ExtraTrees", "GB", "XGBoost", "LightGBM"],
    "Linear": ["L2_LR", "PAC", "NB"],
}

def four_algorithm_consensus(bridge_scores, model_names, gene_index, unknown_mask,
                              drug_probas, disease_probas, top_k=30):
    """
    四算法家族共识筛选:
      1) 将 176 模型按 4 算法族分组
      2) 每族内取所有模型的平均桥梁得分
      3) 每族独立选出 top_k 个未知基因
      4) 取 4 族交集 → "核心跨算法共识靶点"
    返回: (consensus_df, family_results)
    """
    log("\n" + "=" * 70)
    log("STEP 7: 四算法家族交集共识筛选 (LASSO/SVM/Tree/Linear)")
    log("=" * 70)

    unknown_pos = np.where(unknown_mask)[0]
    family_top_genes = {}  # {family: set of top_k gene symbols}

    log("--- 各算法族 Top-30 ---")
    for fname, clf_list in ALGORITHM_FAMILIES.items():
        # 筛选属于该家族的模型
        family_models = []
        for mname in model_names:
            for target_clf in clf_list:
                if mname.endswith(f"__{target_clf}"):
                    if mname in bridge_scores:
                        family_models.append(mname)
                    break

        if not family_models:
            log(f"  [{fname}] 无可用模型, 跳过")
            family_top_genes[fname] = set()
            continue

        # 计算家族内所有模型的平均桥梁得分 (仅未知基因)
        family_avg = np.zeros(len(unknown_pos), dtype=np.float64)
        for mname in family_models:
            family_avg += bridge_scores[mname][unknown_pos]
        family_avg /= len(family_models)

        # 选 top_k
        top_idx = np.argsort(-family_avg)[:top_k]
        top_genes = {gene_index[unknown_pos[i]] for i in top_idx}
        family_top_genes[fname] = top_genes

        # 输出 top-5 示例
        top5 = [gene_index[unknown_pos[i]] for i in top_idx[:5]]
        log(f"  [{fname}] {len(family_models)} 个模型 → Top-5: {top5}")

    # 四族交集
    all_sets = [s for s in family_top_genes.values() if s]
    if len(all_sets) < 2:
        log("  [WARN] 不足2族有效, 无法计算交集")
        return None, family_top_genes

    consensus = set.intersection(*all_sets)
    log(f"\n--- 四算法交集 ---")
    log(f"  交集基因数: {len(consensus)}")
    log(f"  各家族: LASSO={len(family_top_genes.get('LASSO',set()))}, "
        f"SVM={len(family_top_genes.get('SVM',set()))}, "
        f"Tree={len(family_top_genes.get('Tree',set()))}, "
        f"Linear={len(family_top_genes.get('Linear',set()))}")

    if len(consensus) >= 5:
        consensus_list = sorted(consensus)[:10]
        log(f"  示例: {consensus_list}")

    # 返回 DataFrame 格式 (与 df_unknown 格式一致)
    if len(consensus) == 0:
        log("  [WARN] 四族交集为空!")

    return consensus, family_top_genes


# ═══════════════════════════════════════════════════════════════════════════════
# [v4.0 NEW] STEP 8: RF 特征重要性三层分层
#   参考: 暴露-疾病网络毒理框架
# ═══════════════════════════════════════════════════════════════════════════════
def stratified_rf_importance(X, y_dt, y_dg, feature_cols, fe_map, clf_map,
                              metrics_rows, top_k=20):
    """
    RF 特征重要性三层分层:
      L1) 内置特征重要性 (feature_importances_)
      L2) Permutation importance (特征打乱后性能下降)
      L3) SHAP 值 (正 SHAP > 0 的样本比例作为稳定性指标)
    三维度取并集 top_k 特征。

    返回: (stratified_df, union_top_features) 或 (None, set())
    """
    log("\n" + "=" * 70)
    log("STEP 8: RF 特征重要性三层分层 (原始 + Permutation + SHAP)")
    log("=" * 70)

    if not _SHAP_AVAILABLE:
        log("  [SKIP] SHAP 未安装 (pip install shap)")
        # 仍然尝试 L1+L2
    else:
        log(f"  SHAP v{shap.__version__} [OK]")

    # 选择最佳 RF 模型 — 优先 raw__RF (特征空间可解释)
    dt_rf_models = [r for r in metrics_rows
                    if r['task'] == 'DT' and '__RF' in r['model']]
    if not dt_rf_models:
        log("  [SKIP] 无 RF 模型")
        return None, set()

    # 优先使用 raw__RF, 否则用 AUROC 最高的 RF
    raw_rf = [r for r in dt_rf_models if r['model'] == 'raw__RF']
    if raw_rf:
        best_rf = raw_rf[0]
        log(f"  使用 raw__RF (AUROC={best_rf['auroc']:.4f}), 特征维度可解释")
    else:
        best_rf = max(dt_rf_models, key=lambda x: x['auroc'])
        log(f"  最佳 RF 模型: {best_rf['model']} (AUROC={best_rf['auroc']:.4f})")
    model_name = best_rf['model']
    fe_name, clf_name = model_name.split("__", 1)

    # 在全部有标签数据上训练
    y_task = y_dt
    scaler = StandardScaler()
    X_s = scaler.fit_transform(X)

    fe_builder = fe_map.get(fe_name)
    if fe_builder is not None:
        fe_obj = fe_builder()
        if fe_obj is not None:
            if isinstance(fe_obj, tuple(SUPERVISED_FE)):
                fe_obj.fit(X_s, y_task)
            else:
                fe_obj.fit(X_s)
            X_fe = fe_obj.transform(X_s)
            if isinstance(X_fe, tuple):
                X_fe = X_fe[0]
        else:
            X_fe = X_s
    else:
        X_fe = X_s

    clf = RandomForestClassifier(n_estimators=200, class_weight='balanced',
                                  random_state=SEED, n_jobs=1)
    clf.fit(X_fe, y_task)

    # 如果特征维度 ≠ 原始特征数, 跳过 (变换后无解释性)
    if X_fe.shape[1] != len(feature_cols) or fe_name != "raw":
        log(f"  [SKIP] 非 raw 特征空间 ({fe_name}), 维度={X_fe.shape[1]} ≠ {len(feature_cols)}")
        return None, set()

    n_features = len(feature_cols)
    top_sets = {}

    # ── L1: 内置特征重要性 ──
    l1_imp = clf.feature_importances_
    l1_rank = np.argsort(-l1_imp)[:top_k]
    top_sets['L1'] = set(l1_rank)
    log(f"  L1 (内置重要性): Top-5 = {[feature_cols[i] for i in l1_rank[:5]]}")

    # ── L2: Permutation Importance ──
    try:
        perm = permutation_importance(clf, X_fe, y_task, n_repeats=5,
                                       random_state=SEED, n_jobs=1,
                                       scoring='roc_auc')
        l2_imp = perm.importances_mean
        l2_rank = np.argsort(-l2_imp)[:top_k]
        top_sets['L2'] = set(l2_rank)
        log(f"  L2 (Permutation): Top-5 = {[feature_cols[i] for i in l2_rank[:5]]}")
    except Exception as e:
        log(f"  L2 (Permutation): 失败 ({e}), 跳过")
        top_sets['L2'] = set()

    # ── L3: SHAP 值 ──
    if _SHAP_AVAILABLE:
        try:
            # 采样避免内存爆炸
            n_samples = min(500, X_fe.shape[0])
            idx_sample = np.random.choice(X_fe.shape[0], n_samples, replace=False)
            X_sample = X_fe[idx_sample]

            explainer = shap.TreeExplainer(clf, check_additivity=False)
            shap_raw = explainer.shap_values(X_sample)

            # 兼容 shap 多版本返回值:
            #   < 0.42.0: list of arrays (二分类: [neg, pos]) 或 2D array
            #   >= 0.42.0: Explanation object (values 属性)
            if hasattr(shap_raw, 'values'):
                shap_values = shap_raw.values
            else:
                shap_values = shap_raw

            # shap_values shape: (n_samples, n_features) 或 (n_samples, n_features, 2)
            if isinstance(shap_values, list) and len(shap_values) == 2:
                sv_pos = shap_values[1]  # 正类 SHAP
            elif len(shap_values.shape) == 3:
                sv_pos = shap_values[:, :, 1]
            else:
                sv_pos = shap_values

            # 稳定性指标: SHAP > 0 的样本比例
            l3_stability = np.mean(sv_pos > 0, axis=0)
            l3_rank = np.argsort(-l3_stability)[:top_k]
            top_sets['L3'] = set(l3_rank)
            log(f"  L3 (SHAP稳定性): Top-5 = {[feature_cols[i] for i in l3_rank[:5]]}")
        except Exception as e:
            log(f"  L3 (SHAP): 失败 ({e}), 跳过")
            top_sets['L3'] = set()
    else:
        top_sets['L3'] = set()

    # 三维度并集
    all_valid = [s for s in top_sets.values() if s]
    if not all_valid:
        return None, set()

    union_idx = sorted(set.union(*all_valid))
    log(f"\n  三层并集: {len(union_idx)} 个特征")
    log(f"  L1={len(top_sets['L1'])}, L2={len(top_sets.get('L2', set()))}, "
        f"L3={len(top_sets.get('L3', set()))}, Union={len(union_idx)}")

    # 构建详细报告
    rows = []
    for idx in union_idx:
        row = {
            'feature': feature_cols[idx],
            'L1_importance': l1_imp[idx],
            'L2_permutation': l2_imp[idx] if idx in top_sets.get('L2', set()) else np.nan,
        }
        if 'L3' in top_sets and idx in top_sets['L3']:
            row['L3_SHAP_stability'] = l3_stability[idx]
        else:
            row['L3_SHAP_stability'] = np.nan
        row['in_L1'] = idx in top_sets.get('L1', set())
        row['in_L2'] = idx in top_sets.get('L2', set())
        row['in_L3'] = idx in top_sets.get('L3', set())
        rows.append(row)

    df = pd.DataFrame(rows).sort_values('L1_importance', ascending=False)
    return df, union_idx


# ═══════════════════════════════════════════════════════════════════════════════
# [v4.0 NEW] STEP 9: 外部交叉验证 (GSE61616 大鼠 CIRI)
#   参考: MLGANN 鲁棒性测试
# ═══════════════════════════════════════════════════════════════════════════════
def external_validation_gse61616(df_unknown, gene_index):
    """
    使用 GSE61616 (大鼠 CIRI 模型) 作为外部验证集。
    计算 ML RRF 排名与 GSE61616 |logFC| 的 Spearman 相关性。
    返回: (spearman_rho, spearman_pvalue, n_overlap)
    """
    log("\n" + "=" * 70)
    log("STEP 9: 外部交叉验证 (GSE61616 大鼠 CIRI 模型)")
    log("=" * 70)

    if not os.path.exists(GSE61616_DEG_PATH):
        log(f"  [SKIP] GSE61616 文件未找到: {GSE61616_DEG_PATH}")
        return None, None, 0

    try:
        gse = pd.read_csv(GSE61616_DEG_PATH, sep='\t')
    except Exception:
        try:
            gse = pd.read_csv(GSE61616_DEG_PATH)
        except Exception as e:
            log(f"  [SKIP] 无法读取 GSE61616: {e}")
            return None, None, 0

    # 标准化列名
    col_map = {}
    for c in gse.columns:
        cl = c.lower()
        if 'human' in cl and 'symbol' in cl:
            col_map['gene'] = c
        elif c.lower() == 'human_symbol':
            col_map['gene'] = c
        elif 'logfc' in cl and 'adj' not in cl:
            col_map['logfc'] = c
        elif 'adj.p.val' in cl or 'padj' in cl or 'fdr' in cl:
            col_map['padj'] = c

    if 'gene' not in col_map:
        # 尝试第一列
        col_map['gene'] = gse.columns[0]
    if 'logfc' not in col_map:
        # 找 logFC 列
        for c in gse.columns:
            if c.lower() in ['logfc', 'log2fc', 'log2foldchange']:
                col_map['logfc'] = c
                break

    if 'gene' not in col_map or 'logfc' not in col_map:
        log(f"  [SKIP] 找不到 gene/logFC 列: {gse.columns.tolist()}")
        return None, None, 0

    gse['gene'] = gse[col_map['gene']].str.upper()
    gse['abs_logFC'] = gse[col_map['logfc']].abs()

    # 去重 (保留最大 |logFC|)
    gse = gse.sort_values('abs_logFC', ascending=False).drop_duplicates(subset='gene')

    # 交集: GSE61616 基因 ∩ ML 未知基因
    ml_genes_set = set(df_unknown['gene_symbol'])
    gse_genes_set = set(gse['gene'])
    common = ml_genes_set & gse_genes_set

    log(f"  GSE61616 总基因: {len(gse)}, 与 ML 交集: {len(common)}")

    if len(common) < 10:
        log("  [SKIP] 交集基因不足10个")
        return None, None, len(common)

    # Spearman: ML RRF score vs GSE61616 |logFC|
    df_u = df_unknown.set_index('gene_symbol')
    gse_ranked = gse[gse['gene'].isin(common)].set_index('gene')

    rrf_vals = []
    logfc_vals = []
    for g in common:
        rrf_vals.append(df_u.loc[g, 'RRF_score'])
        logfc_vals.append(gse_ranked.loc[g, 'abs_logFC'])

    rho, pval = spearmanr(rrf_vals, logfc_vals)
    log(f"  Spearman ρ = {rho:.4f} (p = {pval:.4e}), n = {len(common)}")

    # 可选: 筛选跨数据集稳健的基因 (ML高分 + GSE61616 |logFC| > 中位数)
    ml_thr = np.median(rrf_vals)
    gse_thr = np.median(logfc_vals)
    robust = {g for g, r, l in zip(common, rrf_vals, logfc_vals)
              if r > ml_thr and l > gse_thr}
    log(f"  跨数据集稳健基因: {len(robust)} (ML RRF > {ml_thr:.4f} & GSE |logFC| > {gse_thr:.2f})")
    if robust:
        log(f"  示例: {sorted(robust)[:10]}")

    return rho, pval, len(common)


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 10: 特征重要性 (基础版, 同 v3.0)
# ═══════════════════════════════════════════════════════════════════════════════
def extract_importance(model_name, clf, feature_cols):
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
            'model': model_name, 'feature': feature_cols, 'importance': imp,
        }).sort_values('importance', ascending=False).head(30)
        return df
    except Exception:
        return None


# ═══════════════════════════════════════════════════════════════════════════════
# 主流程 (v4.0)
# ═══════════════════════════════════════════════════════════════════════════════
def main():
    start_time = time.time()
    log("═" * 70)
    log("  多算法融合桥梁靶点预测系统 v4.0")
    log("  - 四算法家族交集共识  - SHAP 重要性分层")
    log("  - 环境暴露特征(Toxi) - 外部验证(GSE61616)")
    log(f"  随机种子: {SEED}, {N_FOLDS}折CV, RRF_k={RRF_K}")
    log(f"  并行核心: {N_JOBS}, GPU加速: {GPU_ENABLED}")
    log("═" * 70)

    # ── STEP 1a: 加载主特征 ──
    feat_df, feature_cols, all_genes = load_data()
    X_base = feat_df[feature_cols].values.astype(np.float64)
    y_dt = feat_df['is_drug_target'].values.astype(int)
    y_dg = feat_df['is_disease_gene'].values.astype(int)
    gene_index = feat_df.index.tolist()
    drug_targets = {g for g in gene_index if feat_df.loc[g, 'is_drug_target'] == 1}

    # ── STEP 1b: Toxi 环境暴露特征 [v4.0 NEW] ──
    toxi_df, toxi_cols = extract_toxi_features(drug_targets, all_genes)
    if toxi_df is not None and toxi_cols:
        # 拼接 toxi 特征到 X_base
        toxi_matrix = toxi_df.loc[gene_index, toxi_cols].values.astype(np.float64)
        X_base = np.hstack([X_base, toxi_matrix])
        feature_cols = feature_cols + toxi_cols
        log(f"  [Toxi] 特征拼接完成: {X_base.shape[1]} 维 (新增 {len(toxi_cols)} 列)")
    else:
        log("  [Toxi] 跳过环境暴露特征")

    # ── STEP 2: 构建模型组合 ──
    log("\n" + "=" * 70)
    log("STEP 2: 构建 特征工程 × 分类器 笛卡尔积组合")
    log("=" * 70)
    model_combos, fe_map, clf_map = build_model_combinations()
    model_names = [m[0] for m in model_combos]
    n_models = len(model_combos)
    log(f"  共 {n_models} 种模型组合")

    unknown_mask = (y_dt == 0) & (y_dg == 0)
    n_total = len(gene_index)
    n_unknown = unknown_mask.sum()
    log(f"  基因: {n_total} (DT+={y_dt.sum()}, DG+={y_dg.sum()}, Unknown={n_unknown})")

    # ── STEP 3: DT Task ──
    log("\n" + "=" * 70)
    log("STEP 3: 药物可靶向性预测 (DT Task) — 并行")
    log("=" * 70)
    drug_probas, dt_metrics = run_cv_parallel(X_base, y_dt, model_combos, "DT", unknown_mask)

    # ── STEP 4: DG Task ──
    log("\n" + "=" * 70)
    log("STEP 4: 疾病相关性预测 (DG Task) — 并行")
    log("=" * 70)
    disease_probas, dg_metrics = run_cv_parallel(X_base, y_dg, model_combos, "DG", unknown_mask)

    # ── STEP 5: 桥梁得分 ──
    log("\n" + "=" * 70)
    log("STEP 5: 计算桥梁得分 (P_drug × P_disease)")
    log("=" * 70)
    bridge_scores = compute_bridge_scores(drug_probas, disease_probas)
    log(f"  有效组合: {len(bridge_scores)}/{n_models}")

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

    # ── STEP 7: 四算法家族共识 [v4.0 NEW] ──
    consensus_genes, family_results = four_algorithm_consensus(
        bridge_scores, model_names, gene_index, unknown_mask,
        drug_probas, disease_probas)

    # ── STEP 8: RF 重要性分层 [v4.0 NEW] ──
    # 收集指标行
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

    stratified_df, union_features = stratified_rf_importance(
        X_base, y_dt, y_dg, feature_cols, fe_map, clf_map, metrics_rows)

    # ── STEP 9: 输出文件 ──
    log("\n" + "=" * 70)
    log("STEP 9: 生成输出文件")
    log("=" * 70)

    # 9a. 所有未知基因
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

    # 标记共识基因
    if consensus_genes is not None and len(consensus_genes) > 0:
        df_all['four_algo_consensus'] = df_all['gene_symbol'].isin(consensus_genes).astype(int)
        log(f"  四算法共识标记: {df_all['four_algo_consensus'].sum()} 个基因")

    df_unknown = df_all[unknown_mask].sort_values('RRF_score', ascending=False)
    unknown_path = os.path.join(OUTPUT_DIR, "ml_bridge_genes_all.csv")
    df_unknown.to_csv(unknown_path, index=False, encoding='utf-8-sig')
    log(f"  所有未知基因: {unknown_path} ({len(df_unknown)} genes)")

    # 9b. Top-20
    top20 = df_unknown.head(20)
    top20_path = os.path.join(OUTPUT_DIR, "top20_bridge_genes_ml.csv")
    top20.to_csv(top20_path, index=False, encoding='utf-8-sig')
    log(f"  Top-20: {top20_path}")

    log("\n" + "=" * 70)
    log("TOP-10 桥梁基因 (RRF 融合)")
    log("=" * 70)
    for i, (_, row) in enumerate(top20.head(10).iterrows()):
        consensus_tag = " *" if (consensus_genes and row['gene_symbol'] in consensus_genes) else ""
        log(f"  {int(row['final_rank']):>4d}. {row['gene_symbol']:<12s}  RRF={row['RRF_score']:.6f}{consensus_tag}")

    # 9c. 算法指标
    df_metrics = pd.DataFrame(metrics_rows)
    metrics_path = os.path.join(OUTPUT_DIR, "algorithm_metrics.csv")
    df_metrics.to_csv(metrics_path, index=False, encoding='utf-8-sig')
    log(f"  算法指标: {metrics_path}")

    # 9d. 各策略最佳模型
    log("\n--- 各特征工程策略的 DT 最佳模型 ---")
    for fe in sorted(FE_MAP.keys()):
        fe_models = [r for r in metrics_rows if r['model'].startswith(fe + "__") and r['task'] == 'DT']
        if fe_models:
            best = max(fe_models, key=lambda x: x['auroc'])
            log(f"  {fe}: {best['model']} AUROC={best['auroc']:.4f}")

    # 9e. 分层特征重要性
    if stratified_df is not None:
        imp_path = os.path.join(OUTPUT_DIR, "rf_stratified_importance.csv")
        stratified_df.to_csv(imp_path, index=False, encoding='utf-8-sig')
        log(f"\n  分层特征重要性: {imp_path} ({len(stratified_df)} features)")

    # 9f. 四算法共识基因
    if consensus_genes is not None and len(consensus_genes) > 0:
        consensus_path = os.path.join(OUTPUT_DIR, "four_algo_consensus_genes.csv")
        consensus_list = sorted(consensus_genes)
        # 附加 RRF 信息
        consensus_rows = []
        for g in consensus_list:
            if g in df_unknown['gene_symbol'].values:
                row_data = df_unknown[df_unknown['gene_symbol'] == g].iloc[0]
                consensus_rows.append({
                    'gene_symbol': g,
                    'RRF_score': row_data['RRF_score'],
                    'final_rank': row_data['final_rank'],
                })
        if consensus_rows:
            pd.DataFrame(consensus_rows).to_csv(consensus_path, index=False, encoding='utf-8-sig')
            log(f"  四算法共识基因: {consensus_path} ({len(consensus_rows)} genes)")

    # ── STEP 10: 外部交叉验证 [v4.0 NEW] ──
    rho_gse, pval_gse, n_overlap_gse = external_validation_gse61616(df_unknown, gene_index)

    # ── STEP 11: GAT 对比 (同 v3.0) ──
    log("\n" + "=" * 70)
    log("STEP 11: GAT 桥梁基因对比")
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
        log("  GAT 桥梁基因文件未找到，跳过")

    # ── 完成 ──
    elapsed = time.time() - start_time
    log("\n" + "=" * 70)
    log(f"[OK] 完成! 总耗时: {elapsed/60:.1f} 分钟")
    log(f"  输出目录: {OUTPUT_DIR}")
    log(f"  模型组合数: {n_models}")
    log(f"  每组合 5折CV × 2任务 = {n_models*5*2} 次训练")
    n_cons = len(consensus_genes) if consensus_genes else 0
    log(f"  四算法共识基因数: {n_cons}")
    if rho_gse is not None:
        log(f"  GSE61616 Spearman ρ = {rho_gse:.4f} (p = {pval_gse:.4e})")
    log("=" * 70)


if __name__ == "__main__":
    main()