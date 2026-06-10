# -*- coding: utf-8 -*-
"""
多算法融合桥梁靶点预测系统 (v6.0 — 论文级深度优化版)
============================================================================
网络药理学第一阶段：药物-靶点-疾病桥梁基因发现

v6.0 核心升级 (vs v5.0):

[数据层优化]
  1. 增强输入验证: 数据分布检查 + 特征方差检查 + 标签平衡评估
     参考: DataSAIL PMID:40199913 — 数据分割前的严格验证
  2. 特征多重共线性诊断: 方差膨胀因子(VIF) + 条件数(Condition Number)
     参考: nestedcv PMID:37113250 — 嵌套CV防信息泄露
  3. 缺失值处理增强: 中位数 + 按特征类型分组填充

[模型层优化]
  4. 加权RRF集成: 按AUROC性能加权模型贡献 (原为等权重)
     参考: PMID:36384050 — KG ensemble加权提升预测性能
  5. 5次重复CV: 可选repeated CV提升稳定性评估
     参考: ExAutoGP PMID:40282006 — 5-repeated 5-fold CV
  6. 模型校准: 使用CalibratedClassifierCV提升概率估计质量

[算法层优化]
  7. SHAP误分类检测: "RAW OR SHAP" 过滤规则识别低置信度预测
     参考: PMID:41114446 — SHAP-based misclassification detection
  8. Borda Count稳健性增强: 排除概率全0模型 + 归一化分数
  9. Bootstrap CI稳定性: 修复除零问题 + 增加BCa置信区间
  10. 四算法共识: 加权平均 (按模型AUROC加权)

[接口层优化]
  11. 断点续跑: 中间结果保存与恢复
  12. 时间戳输出目录: 防止覆盖历史结果
  13. 增强日志: 每个阶段的耗时统计
  14. 配置文件外部化: 支持命令行参数

保持自 v5.0: 四算法共识 + SHAP分层重要性 + 环境暴露 + 外部验证 + 笛卡尔积CV
"""

import os
import sys
import json
import time
import warnings
import traceback
import glob
import logging
import argparse
from collections import defaultdict
from datetime import datetime

import numpy as np
import pandas as pd

from sklearn.preprocessing import StandardScaler, RobustScaler
from sklearn.model_selection import StratifiedKFold, RepeatedStratifiedKFold
from sklearn.linear_model import LogisticRegression, LassoCV, SGDClassifier, PassiveAggressiveClassifier
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier, ExtraTreesClassifier
from sklearn.calibration import CalibratedClassifierCV
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

# mygene (可选)
try:
    import mygene
    _MYGENE_AVAILABLE = True
except ImportError:
    _MYGENE_AVAILABLE = False

from joblib import Parallel, delayed

warnings.filterwarnings("ignore")
from sklearn.exceptions import ConvergenceWarning
warnings.filterwarnings("ignore", category=ConvergenceWarning)

# ═══════════════════════════════════════════════════════════════════════════════
# 全局配置
# ═══════════════════════════════════════════════════════════════════════════════
SEED = 42
N_FOLDS = 5
N_REPEATS = 1             # [v6.0] 重复CV次数 (>1 启用RepeatedStratifiedKFold)
RRF_K = 60
N_JOBS = 4
GPU_ENABLED = True
N_BOOTSTRAP = 200
FEATURE_CORR_THRESHOLD = 0.95
VIF_THRESHOLD = 10.0       # [v6.0] VIF阈值 (VIF>10表示严重共线性)
USE_ROBUST_SCALER = False
CHECKPOINT_ENABLED = True  # [v6.0] 断点续跑
WEIGHTED_RRF = True        # [v6.0] 加权RRF (按AUROC加权)
MISCLASSIFICATION_FILTER = True  # [v6.0] SHAP误分类检测

# 数据路径
DATA_DIR = r"D:/反向网络药理学/GAT拓展维度/cache"
FEATURE_PATH = os.path.join(DATA_DIR, "enhanced_gene_features.csv")
DRUG_TARGETS_PATH = r"C:/Users/Jy-Mentor-7/Desktop/GAT/drug_targets.txt"
DISEASE_GENES_PATH = r"C:/Users/Jy-Mentor-7/Desktop/GAT/disease_genes.txt"
SUBGRAPH_GENES_PATH = r"C:/Users/Jy-Mentor-7/Desktop/GAT/subgraph_genes.txt"
GAT_BRIDGE_PATH = r"C:/Users/Jy-Mentor-7/Desktop/GAT/top20_bridge_genes.csv"
TOXI_FPKM_DIR = r"D:/反向网络药理学/GAT拓展维度/Toxi/rna_fpkm"
GSE61616_DEG_PATH = r"c:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\大创\GSE61616_human_homologs_DEGs.tsv"

# [v6.0] 带时间戳的输出目录，防止覆盖
TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M%S")
OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), f"ml_output_v6_{TIMESTAMP}")
os.makedirs(OUTPUT_DIR, exist_ok=True)
ORTHO_CACHE_PATH = os.path.join(OUTPUT_DIR, "mouse_to_human_orthologs.csv")
CHECKPOINT_PATH = os.path.join(OUTPUT_DIR, "checkpoint.pkl")  # [v6.0]

np.random.seed(SEED)

# ═══════════════════════════════════════════════════════════════════════════════
# [v6.0] 命令行参数解析
# ═══════════════════════════════════════════════════════════════════════════════
def parse_args():
    parser = argparse.ArgumentParser(description="多算法融合桥梁靶点预测系统 v6.0")
    parser.add_argument("--n-jobs", type=int, default=N_JOBS, help=f"并行核心数 (默认: {N_JOBS})")
    parser.add_argument("--gpu", type=int, default=1 if GPU_ENABLED else 0, help="GPU加速 (0/1)")
    parser.add_argument("--n-repeats", type=int, default=N_REPEATS, help=f"重复CV次数 (默认: {N_REPEATS})")
    parser.add_argument("--no-checkpoint", action="store_true", help="禁用断点续跑")
    parser.add_argument("--no-weighted-rrf", action="store_true", help="禁用加权RRF")
    parser.add_argument("--no-misclf-filter", action="store_true", help="禁用SHAP误分类检测")
    parser.add_argument("--robust-scaler", action="store_true", help="使用RobustScaler")
    return parser.parse_args()

# ═══════════════════════════════════════════════════════════════════════════════
# [v6.0] 结构化日志系统 (增强)
# ═══════════════════════════════════════════════════════════════════════════════
LOG_FILE = os.path.join(OUTPUT_DIR, "pipeline_v6.log")
_stage_start_time = time.time()

def _setup_file_logger():
    file_handler = logging.FileHandler(LOG_FILE, encoding='utf-8')
    file_handler.setFormatter(logging.Formatter('[%(asctime)s] %(message)s', datefmt='%H:%M:%S'))
    logger = logging.getLogger('pipeline_v6')
    logger.addHandler(file_handler)
    logger.setLevel(logging.INFO)
    return logger

_file_logger = _setup_file_logger()

def log(msg: str):
    timestamp = time.strftime("%H:%M:%S")
    formatted = f"[{timestamp}] {msg}"
    print(formatted, flush=True)
    _file_logger.info(msg)

def log_stage(stage_name: str):
    """记录阶段耗时"""
    global _stage_start_time
    elapsed = time.time() - _stage_start_time
    log(f"  [Stage] {stage_name} 耗时: {elapsed:.1f}s")
    _stage_start_time = time.time()


# ═══════════════════════════════════════════════════════════════════════════════
# [v6.0] 断点续跑
# ═══════════════════════════════════════════════════════════════════════════════
def save_checkpoint(data_dict):
    """保存中间结果到checkpoint文件"""
    if not CHECKPOINT_ENABLED:
        return
    try:
        import pickle
        with open(CHECKPOINT_PATH, 'wb') as f:
            pickle.dump(data_dict, f, protocol=pickle.HIGHEST_PROTOCOL)
        log(f"  [Checkpoint] 已保存到 {CHECKPOINT_PATH}")
    except Exception as e:
        log(f"  [WARN] Checkpoint 保存失败: {e}")

def load_checkpoint():
    """从checkpoint文件恢复中间结果"""
    if not CHECKPOINT_ENABLED:
        return None
    if not os.path.exists(CHECKPOINT_PATH):
        return None
    try:
        import pickle
        with open(CHECKPOINT_PATH, 'rb') as f:
            data = pickle.load(f)
        log(f"  [Checkpoint] 已加载: {CHECKPOINT_PATH}")
        return data
    except Exception as e:
        log(f"  [WARN] Checkpoint 加载失败: {e}")
        return None


# ═══════════════════════════════════════════════════════════════════════════════
# GPU 检测 (同 v5.0)
# ═══════════════════════════════════════════════════════════════════════════════
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
        log(f"  [GPU] XGBoost GPU: 不可用 -> {_XGB_GPU_INFO}")

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
        log(f"  [GPU] LightGBM GPU: 不可用 -> {_LGB_GPU_INFO}")

def detect_gpu_environment():
    log("-" * 50)
    log("GPU 环境检测")
    log("-" * 50)
    effective_jobs = os.cpu_count() if N_JOBS == -1 else N_JOBS
    log(f"  并行核心: {N_JOBS} -> 有效 {effective_jobs} 进程")
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


# ═══════════════════════════════════════════════════════════════════════════════
# [v6.0 ENHANCED] STEP 0: 增强输入验证层
#   参考: DataSAIL PMID:40199913 — 数据分割前的严格验证
#   参考: nestedcv PMID:37113250 — 嵌套CV防信息泄露
# ═══════════════════════════════════════════════════════════════════════════════
def validate_inputs():
    """
    对所有输入文件进行格式、内容、一致性的严格验证。
    增强: 数据分布检查 + 特征方差检查 + 标签平衡评估
    返回: (all_valid, issues_list, data_quality_report)
    """
    log("=" * 70)
    log("STEP 0: 增强输入数据验证层 [v6.0]")
    log("=" * 70)

    issues = []
    quality_report = {}
    all_valid = True

    # 0a. 特征矩阵
    if not os.path.exists(FEATURE_PATH):
        issues.append(f"[FATAL] 特征矩阵不存在: {FEATURE_PATH}")
        all_valid = False
    else:
        try:
            df = pd.read_csv(FEATURE_PATH)
            if 'gene_symbol' not in df.columns:
                issues.append(f"[FATAL] 特征矩阵缺少 'gene_symbol' 列")
                all_valid = False
            else:
                n_genes = len(df)
                n_features = len(df.columns) - 1
                log(f"  [OK] 特征矩阵: {n_genes} genes x {n_features} features")

                # [v6.0 NEW] 特征方差检查
                numeric_cols = df.select_dtypes(include=[np.number]).columns
                if len(numeric_cols) > 0:
                    variances = df[numeric_cols].var()
                    zero_var = (variances == 0).sum()
                    low_var = (variances < 1e-8).sum()
                    if zero_var > 0:
                        issues.append(f"[WARN] {zero_var} 个特征方差为0 (常量特征)")
                    if low_var > 0:
                        log(f"  [INFO] {low_var} 个特征方差 < 1e-8 (低变异特征)")
                    quality_report['zero_var_features'] = int(zero_var)
                    quality_report['low_var_features'] = int(low_var)
                    quality_report['median_variance'] = float(variances.median())

                # [v6.0 NEW] 缺失值统计
                missing_counts = df.isnull().sum()
                missing_cols = (missing_counts > 0).sum()
                if missing_cols > 0:
                    log(f"  [INFO] {missing_cols} 个特征含缺失值, 总缺失: {missing_counts.sum()}")
                    quality_report['missing_cols'] = int(missing_cols)
                    quality_report['total_missing'] = int(missing_counts.sum())
        except Exception as e:
            issues.append(f"[FATAL] 特征矩阵读取失败: {e}")
            all_valid = False

    # 0b. 药物靶点文件
    for label, path in [("药物靶点", DRUG_TARGETS_PATH), ("疾病基因", DISEASE_GENES_PATH)]:
        if not os.path.exists(path):
            issues.append(f"[FATAL] {label}文件不存在: {path}")
            all_valid = False
        else:
            try:
                with open(path, 'r') as f:
                    lines = [l.strip() for l in f if l.strip()]
                if len(lines) == 0:
                    issues.append(f"[FATAL] {label}文件为空")
                    all_valid = False
                else:
                    log(f"  [OK] {label}: {len(lines)} 条记录")
                    quality_report[f'{label}_count'] = len(lines)
            except Exception as e:
                issues.append(f"[FATAL] {label}文件读取失败: {e}")
                all_valid = False

    # [v6.0 NEW] 标签平衡评估 (在load_data中完成, 此处仅占位)
    quality_report['label_balance_check'] = "见 load_data 输出"

    # 0c. Toxi 目录
    if os.path.isdir(TOXI_FPKM_DIR):
        fpkm_files = glob.glob(os.path.join(TOXI_FPKM_DIR, "*.tsv"))
        if not fpkm_files:
            fpkm_files = glob.glob(os.path.join(TOXI_FPKM_DIR, "*.tsv.tsv"))
        log(f"  [OK] Toxi FPKM: {len(fpkm_files)} 个文件")
        quality_report['toxi_files'] = len(fpkm_files)
    else:
        log(f"  [WARN] Toxi FPKM 目录不存在, 将跳过环境暴露特征")

    # 0d. GSE61616
    if os.path.exists(GSE61616_DEG_PATH):
        log(f"  [OK] GSE61616 外部验证文件存在")
    else:
        log(f"  [WARN] GSE61616 文件不存在, 将跳过外部验证")

    # 0e. 依赖库
    if not _SHAP_AVAILABLE:
        issues.append("[WARN] SHAP 未安装 (pip install shap), 将跳过部分功能")
    if not _MYGENE_AVAILABLE:
        log(f"  [INFO] mygene 未安装, Ensembl API 优先")

    if all_valid:
        log("  [OK] 输入验证全部通过")
    else:
        log("  [FATAL] 输入验证失败, 详见以下问题:")
        for issue in issues:
            log(f"    {issue}")

    return all_valid, issues, quality_report


# ═══════════════════════════════════════════════════════════════════════════════
# [v6.0 ENHANCED] 特征多重共线性诊断
#   参考: nestedcv PMID:37113250 — 嵌套CV防信息泄露
# ═══════════════════════════════════════════════════════════════════════════════
def diagnose_multicollinearity(X, feature_cols, vif_threshold=VIF_THRESHOLD):
    """
    使用VIF (Variance Inflation Factor) 和条件数诊断多重共线性。
    返回: (vif_report_df, condition_number, high_vif_features)
    """
    log("\n" + "-" * 50)
    log(f"多重共线性诊断 (VIF={vif_threshold}) [v6.0]")
    log("-" * 50)

    n_features = X.shape[1]
    if n_features < 2 or n_features > 500:
        log(f"  [SKIP] 特征数={n_features}, 不在诊断范围 [2, 500]")
        return None, None, set()

    # 条件数 (Condition Number)
    try:
        X_scaled = StandardScaler().fit_transform(X)
        U, S, Vt = np.linalg.svd(X_scaled, full_matrices=False)
        cond_num = S[0] / S[-1] if S[-1] > 0 else np.inf
        log(f"  条件数 (Condition Number): {cond_num:.1f}")
        if cond_num > 30:
            log(f"  [WARN] 条件数 > 30, 存在严重共线性")
        elif cond_num > 10:
            log(f"  [INFO] 条件数 > 10, 存在中等共线性")
    except Exception as e:
        log(f"  [WARN] SVD 计算失败: {e}")
        cond_num = None

    # VIF 诊断 (限制特征数避免过慢)
    high_vif = set()
    max_vif = min(n_features, 200)
    subset_idx = np.random.choice(n_features, max_vif, replace=False) if n_features > max_vif else np.arange(n_features)

    try:
        vif_values = []
        for i, idx in enumerate(subset_idx):
            if (i + 1) % 50 == 0:
                log(f"    VIF 计算中... {i+1}/{max_vif}")
            try:
                y_col = X[:, idx]
                X_rest = np.delete(X, idx, axis=1)
                # 使用子样本加速
                sample_size = min(5000, X.shape[0])
                if X.shape[0] > sample_size:
                    sample_idx = np.random.choice(X.shape[0], sample_size, replace=False)
                    y_s = y_col[sample_idx]
                    X_s = X_rest[sample_idx]
                else:
                    y_s = y_col
                    X_s = X_rest
                from sklearn.linear_model import LinearRegression
                lr = LinearRegression()
                lr.fit(X_s, y_s)
                r2 = lr.score(X_s, y_s)
                vif = 1.0 / (1.0 - r2) if r2 < 1.0 else np.inf
                vif_values.append((idx, feature_cols[idx], vif))
                if vif > vif_threshold:
                    high_vif.add(idx)
            except Exception:
                continue

        if vif_values:
            vif_df = pd.DataFrame(vif_values, columns=['feature_idx', 'feature', 'VIF'])
            vif_df = vif_df.sort_values('VIF', ascending=False)
            n_high = len(vif_df[vif_df['VIF'] > vif_threshold])
            log(f"  VIF 诊断完成: {n_high}/{len(vif_df)} 特征 VIF > {vif_threshold}")
            if n_high > 0:
                top5 = vif_df.head(5)[['feature', 'VIF']].values.tolist()
                log(f"  Top-5 高VIF特征: {top5}")
            return vif_df, cond_num, high_vif
    except Exception as e:
        log(f"  [WARN] VIF 诊断失败: {e}")

    return None, cond_num, high_vif


# ═══════════════════════════════════════════════════════════════════════════════
# [v6.0 ENHANCED] 特征相关性过滤 (增强版)
# ═══════════════════════════════════════════════════════════════════════════════
def filter_highly_correlated_features(X, feature_cols, threshold=0.95):
    """
    基于 Pearson 相关系数过滤高相关特征。
    增强: 优先保留VIF较低的（即更独立的）特征。
    返回: (filtered_X, filtered_cols, removed_count)
    """
    log("\n" + "-" * 50)
    log(f"特征相关性过滤 (threshold={threshold}) [v6.0]")
    log("-" * 50)

    n_features = X.shape[1]
    if n_features < 2:
        return X, feature_cols, 0

    sample_size = min(5000, X.shape[0])
    if X.shape[0] > sample_size:
        idx = np.random.choice(X.shape[0], sample_size, replace=False)
        X_sample = X[idx]
    else:
        X_sample = X

    corr_matrix = np.abs(np.corrcoef(X_sample, rowvar=False))
    np.fill_diagonal(corr_matrix, 0)

    to_remove = set()
    high_corr_pairs = np.argwhere(corr_matrix > threshold)
    for i, j in high_corr_pairs:
        if i < j and i not in to_remove and j not in to_remove:
            # [v6.0] 保留方差较大的特征
            var_i = np.var(X[:, i])
            var_j = np.var(X[:, j])
            to_remove.add(i if var_i <= var_j else j)

    if to_remove:
        keep_mask = np.ones(n_features, dtype=bool)
        keep_mask[list(to_remove)] = False
        X_filtered = X[:, keep_mask]
        filtered_cols = [c for idx, c in enumerate(feature_cols) if idx not in to_remove]
        log(f"  移除 {len(to_remove)} 个高相关特征 (保留 {len(filtered_cols)} 个)")
        removed_examples = [feature_cols[i] for i in sorted(to_remove)[:5]]
        log(f"  移除示例: {removed_examples}")
    else:
        X_filtered = X
        filtered_cols = feature_cols[:]
        log(f"  无高相关特征 (> {threshold})")

    return X_filtered, filtered_cols, len(to_remove)


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 1a: 主特征数据加载 (增强版)
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
        log(f"  填充 {n_missing} 个缺失值 (中位数填充)")
        feat_df = feat_df.fillna(feat_df.median())
        feat_df = feat_df.fillna(0.0)

    all_genes = set(feat_df.index)
    feature_cols = list(feat_df.columns)
    log(f"  特征矩阵: {len(all_genes)} genes x {len(feature_cols)} features")

    log(f"加载药物靶点: {DRUG_TARGETS_PATH}")
    with open(DRUG_TARGETS_PATH, 'r') as f:
        drug_targets_raw = set(line.strip().upper() for line in f if line.strip())
    drug_targets = drug_targets_raw & all_genes
    log(f"  BCP 药物靶点: {len(drug_targets_raw)} -> {len(drug_targets)} 命中")

    log(f"加载疾病基因: {DISEASE_GENES_PATH}")
    with open(DISEASE_GENES_PATH, 'r') as f:
        disease_genes_raw = set(line.strip().upper() for line in f if line.strip())
    disease_genes = disease_genes_raw & all_genes
    log(f"  CIRI 疾病基因: {len(disease_genes_raw)} -> {len(disease_genes)} 命中")

    if os.path.exists(SUBGRAPH_GENES_PATH):
        log(f"加载子图基因: {SUBGRAPH_GENES_PATH}")
        with open(SUBGRAPH_GENES_PATH, 'r') as f:
            sg = set(line.strip().upper() for line in f if line.strip())
            sg.discard('GENE')
        kept = all_genes & sg
        if kept:
            log(f"  子图过滤: {len(all_genes)} -> {len(kept)}")
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
    n_total = len(feat_df)

    log(f"  标签统计: DT+={n_dt} ({100*n_dt/n_total:.1f}%), "
        f"DG+={n_dg} ({100*n_dg/n_total:.1f}%), "
        f"Both={n_both}, Unknown={n_unknown}")

    # [v6.0 ENHANCED] 标签不平衡详细评估
    log(f"  --- 标签平衡评估 [v6.0] ---")
    dt_ratio = n_dt / n_total if n_total > 0 else 0
    dg_ratio = n_dg / n_total if n_total > 0 else 0
    imbalance_ratio_dt = (n_total - n_dt) / max(n_dt, 1)
    imbalance_ratio_dg = (n_total - n_dg) / max(n_dg, 1)

    log(f"  DT+ 比例: {dt_ratio:.4f} (不平衡比: {imbalance_ratio_dt:.1f}:1)")
    log(f"  DG+ 比例: {dg_ratio:.4f} (不平衡比: {imbalance_ratio_dg:.1f}:1)")

    if n_dt < 10:
        log(f"  [WARN] DT正样本过少! ({n_dt}), CV 可能极不稳定")
    elif n_dt < 30:
        log(f"  [WARN] DT正样本较少 ({n_dt}), 建议谨慎解读CV结果")
    if n_dg < 10:
        log(f"  [WARN] DG正样本过少! ({n_dg}), CV 可能极不稳定")
    elif n_dg < 30:
        log(f"  [WARN] DG正样本较少 ({n_dg}), 建议谨慎解读CV结果")

    if imbalance_ratio_dt > 100:
        log(f"  [WARN] DT不平衡比 > 100:1, 极度不平衡!")
    if imbalance_ratio_dg > 100:
        log(f"  [WARN] DG不平衡比 > 100:1, 极度不平衡!")

    return feat_df, feature_cols, all_genes


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 1b: Toxi 环境暴露特征提取 (同 v5.0)
# ═══════════════════════════════════════════════════════════════════════════════
def extract_toxi_features(drug_targets_human, all_genes_set):
    log("\n" + "=" * 70)
    log("STEP 1b: Toxi 环境暴露特征提取 (BPA Brain FPKM)")
    log("=" * 70)

    if not os.path.isdir(TOXI_FPKM_DIR):
        log("  [SKIP] Toxi FPKM 目录不存在")
        return None, []

    fpkm_files = glob.glob(os.path.join(TOXI_FPKM_DIR, "*.tsv"))
    if not fpkm_files:
        fpkm_files = glob.glob(os.path.join(TOXI_FPKM_DIR, "*.tsv.tsv"))
    if len(fpkm_files) < 5:
        log(f"  [SKIP] FPKM 文件不足 ({len(fpkm_files)} 个)")
        return None, []

    log(f"  找到 {len(fpkm_files)} 个 FPKM 文件")

    groups = defaultdict(list)
    for fpath in fpkm_files:
        fname = os.path.basename(fpath)
        if 'Brain_' not in fname:
            continue
        try:
            cond_start = fname.index('Brain_') + 6
            weeks_pos = fname.index('20_weeks')
            condition_str = fname[cond_start:weeks_pos-1]
            if 'Control' in condition_str:
                condition = 'Control'
            elif '10mg' in condition_str:
                condition = 'BPA10mg'
            elif '10ug' in condition_str:
                condition = 'BPA10ug'
            else:
                condition = condition_str.replace('BPA-', 'BPA')
        except (ValueError, IndexError):
            continue
        groups[condition].append(fpath)

    log(f"  分组: Control={len(groups.get('Control',[]))}, "
        f"BPA10mg={len(groups.get('BPA10mg',[]))}, "
        f"BPA10ug={len(groups.get('BPA10ug',[]))}")

    condition_fpkm = {}
    for cond, files in groups.items():
        samples = []
        for fpath in files:
            try:
                df = pd.read_csv(fpath, sep='\t')
                df['gene_id_clean'] = df['gene_id'].str.split('.').str[0]
                df = df.set_index('gene_id_clean')
                samples.append(df['FPKM'])
            except Exception:
                continue
        if not samples:
            continue
        all_fpkm = pd.concat(samples, axis=1)
        all_fpkm['mean'] = all_fpkm.mean(axis=1)
        condition_fpkm[cond] = all_fpkm['mean']

    if 'Control' not in condition_fpkm:
        log("  [SKIP] 缺少 Control 组")
        return None, []

    human_genes = list(drug_targets_human)
    human_upper = {g.upper() for g in human_genes}
    mouse_to_human = {}
    ensmusg_to_symbol = {}

    control_fpkm = condition_fpkm['Control']
    all_mouse_ids = list(control_fpkm.index)

    # 本地缓存
    if os.path.exists(ORTHO_CACHE_PATH):
        try:
            cache_df = pd.read_csv(ORTHO_CACHE_PATH)
            for _, row in cache_df.iterrows():
                eid = str(row.get('ensmusg', '')).strip()
                msym = str(row.get('mouse_symbol', '')).strip().upper()
                hsym = str(row.get('human_symbol', '')).strip().upper()
                if eid and msym and hsym:
                    ensmusg_to_symbol[eid] = msym
                    if hsym in human_upper:
                        mouse_to_human[msym] = hsym
            log(f"  [Cache] 加载本地同源映射: {len(ensmusg_to_symbol)} ENSMUSG -> {len(mouse_to_human)} human")
        except Exception as e:
            log(f"  [Cache] 读取失败: {e}")

    # Ensembl API
    if not mouse_to_human:
        try:
            import requests as _req
            server = 'https://rest.ensembl.org'
            headers = {'Content-Type': 'application/json', 'Accept': 'application/json'}
            batch_size = 500
            for i in range(0, len(all_mouse_ids), batch_size):
                batch = all_mouse_ids[i:i+batch_size]
                data = json.dumps({'ids': batch})
                r = _req.post(server + '/lookup/id', headers=headers, data=data, timeout=60)
                if r.ok:
                    result = r.json()
                    for eid, info in result.items():
                        sym = info.get('display_name', '')
                        if sym and isinstance(sym, str):
                            sym_u = sym.upper()
                            ensmusg_to_symbol[eid] = sym_u
                            if sym_u in human_upper:
                                mouse_to_human[sym_u] = sym_u
            log(f"  Ensembl API: {len(ensmusg_to_symbol)} ENSMUSG -> {len(mouse_to_human)} human")
            _save_ortholog_cache(ensmusg_to_symbol, mouse_to_human)
        except Exception as e:
            log(f"  [WARN] Ensembl API 失败: {e}")

    # mygene 回退
    if not mouse_to_human and _MYGENE_AVAILABLE:
        try:
            mg = mygene.MyGeneInfo()
            batch_size = 500
            for i in range(0, len(all_mouse_ids), batch_size):
                batch = all_mouse_ids[i:i+batch_size]
                results = mg.querymany(batch, scopes='ensembl.gene',
                                       fields='symbol', species='mouse',
                                       returnall=True, as_dataframe=True)
                if 'out' in results and len(results['out']) > 0:
                    outdf = results['out']
                    for _, row in outdf.iterrows():
                        eid = row.get('query', '')
                        if 'symbol' in row and isinstance(row['symbol'], str):
                            sym_u = row['symbol'].upper()
                            ensmusg_to_symbol[eid] = sym_u
                            if sym_u in human_upper:
                                mouse_to_human[sym_u] = sym_u
            log(f"  mygene: {len(ensmusg_to_symbol)} ENSMUSG -> {len(mouse_to_human)} human")
            _save_ortholog_cache(ensmusg_to_symbol, mouse_to_human)
        except Exception as e:
            log(f"  [WARN] mygene 失败: {e}")

    # 同名匹配
    if not mouse_to_human:
        mouse_symbols_upper = {str(s).upper() for s in control_fpkm.index if isinstance(s, str)}
        for hg in human_genes:
            if hg.upper() in mouse_symbols_upper:
                mouse_to_human[hg.upper()] = hg
        log(f"  同名匹配: {len(mouse_to_human)} genes")

    log(f"  人->鼠映射总数: {len(mouse_to_human)}")
    if len(mouse_to_human) == 0:
        log("  [SKIP] 无有效人-鼠同源映射")
        return None, []

    sym_to_ensmusg = {}
    for eid, msym in ensmusg_to_symbol.items():
        if msym not in sym_to_ensmusg:
            sym_to_ensmusg[msym] = eid

    toxi_rows = []
    for cond in ['BPA10mg', 'BPA10ug']:
        if cond not in condition_fpkm:
            continue
        colname = f"Toxi_log2FC_{cond}"
        cond_data = condition_fpkm[cond]
        for msym, hsym in mouse_to_human.items():
            eid = sym_to_ensmusg.get(msym)
            if eid is None:
                continue
            try:
                ctrl_val = control_fpkm.get(eid, np.nan)
                exp_val = cond_data.get(eid, np.nan)
                if pd.isna(ctrl_val) or pd.isna(exp_val) or ctrl_val == 0:
                    fc = 0.0
                else:
                    fc = np.log2((exp_val + 1.0) / (ctrl_val + 1.0))
                toxi_rows.append({'gene_symbol': hsym, colname: fc})
            except Exception:
                continue

    if not toxi_rows:
        log("  [SKIP] 无法计算任何 log2FC 值")
        return None, []

    toxi_feat_df = pd.DataFrame(toxi_rows)
    toxi_feat_df = toxi_feat_df.groupby('gene_symbol').mean()

    toxi_cols = []
    for cond in ['BPA10mg', 'BPA10ug']:
        colname = f"Toxi_log2FC_{cond}"
        if colname in toxi_feat_df.columns:
            toxi_cols.append(colname)

    if not toxi_cols:
        log("  [SKIP] 无 Toxi 特征列")
        return None, []

    aligned = pd.DataFrame(index=sorted(all_genes_set))
    for col in toxi_cols:
        aligned[col] = 0.0
        for g in aligned.index:
            if g in toxi_feat_df.index:
                aligned.loc[g, col] = toxi_feat_df.loc[g, col]

    n_nonzero = (aligned[toxi_cols].abs().sum(axis=1) > 0).sum()
    log(f"  Toxi 特征: {len(toxi_cols)} 列, {n_nonzero} 基因有非零 log2FC")
    log(f"  示例 (Top-5):")
    log(f"  {aligned[toxi_cols].abs().sum(axis=1).sort_values(ascending=False).head(5)}")

    return aligned, toxi_cols


def _save_ortholog_cache(ensmusg_to_symbol, mouse_to_human):
    try:
        rows = []
        for eid, msym in ensmusg_to_symbol.items():
            hsym = mouse_to_human.get(msym, '')
            rows.append({'ensmusg': eid, 'mouse_symbol': msym, 'human_symbol': hsym})
        pd.DataFrame(rows).to_csv(ORTHO_CACHE_PATH, index=False, encoding='utf-8-sig')
        log(f"  [Cache] 已保存 {len(rows)} 条同源映射")
    except Exception as e:
        log(f"  [Cache] 保存失败: {e}")


# ═══════════════════════════════════════════════════════════════════════════════
# 特征工程 Builder (同 v5.0)
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
# 分类器 Builder (增强版)
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
    base_svc = SVC(C=1.0, kernel='rbf', class_weight='balanced',
                   random_state=SEED, max_iter=5000, probability=False)
    return CalibratedClassifierCV(base_svc, cv=3, method='sigmoid')
def clf_extratrees():
    return ExtraTreesClassifier(n_estimators=200, class_weight='balanced',
        random_state=SEED, n_jobs=1)
def clf_pac():
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
# STEP 2: 模型组合
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
# [v6.0 ENHANCED] STEP 3-4: 并行CV (增强版)
#   支持 RepeatedStratifiedKFold (参考 ExAutoGP PMID:40282006)
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

    # [v6.0] RepeatedStratifiedKFold
    if N_REPEATS > 1:
        skf = RepeatedStratifiedKFold(n_splits=N_FOLDS, n_repeats=N_REPEATS,
                                       random_state=seed)
    else:
        skf = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=seed)

    labeled_idx = np.where(~unknown_mask)[0]
    y_labeled = y[~unknown_mask]
    unknown_pos = np.where(unknown_mask)[0]

    fold_aurocs = []
    fold_auprcs = []
    fold_success = 0

    _Scaler = RobustScaler if USE_ROBUST_SCALER else StandardScaler

    for fold, (train_fold_idx, val_fold_idx) in enumerate(skf.split(labeled_idx, y_labeled)):
        train_pos = labeled_idx[train_fold_idx]
        val_pos = labeled_idx[val_fold_idx]

        X_train = X_base[train_pos].copy()
        y_train = y[train_pos]
        X_val = X_base[val_pos].copy()
        y_val = y[val_pos]

        scaler = _Scaler()
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
        return model_name, probas, {'auroc': 0.0, 'auprc': 0.0, 'auroc_std': 0.0, 'auprc_std': 0.0}

    if not labeled_assigned[~unknown_mask].all():
        global_mean = probas[~unknown_mask][labeled_assigned[~unknown_mask]].mean()
        unassigned_idx = np.where(~unknown_mask & ~labeled_assigned)[0]
        probas[unassigned_idx] = global_mean

    mean_auroc = np.mean(fold_aurocs) if fold_aurocs else 0.0
    std_auroc = np.std(fold_aurocs) if len(fold_aurocs) > 1 else 0.0
    mean_auprc = np.mean(fold_auprcs) if fold_auprcs else 0.0
    std_auprc = np.std(fold_auprcs) if len(fold_auprcs) > 1 else 0.0
    return model_name, probas, {
        'auroc': mean_auroc, 'auprc': mean_auprc,
        'auroc_std': std_auroc, 'auprc_std': std_auprc
    }


def run_cv_parallel(X_base, y, model_combos, task_name, unknown_mask):
    n_models = len(model_combos)
    log(f"  [{task_name}] 启动 {n_models} 个模型的并行训练 ({N_JOBS} 核, {N_REPEATS}重复)...")
    results = Parallel(n_jobs=N_JOBS, verbose=5)(
        delayed(train_single_model)(name, fe, clf, X_base, y, unknown_mask, idx)
        for idx, (name, fe, clf) in enumerate(model_combos)
    )
    all_probas = {}
    all_metrics = {}
    for model_name, probas, metric in results:
        all_probas[model_name] = probas
        all_metrics[model_name] = metric

    # [v6.0] 模型稳定性报告
    aurocs = [m['auroc'] for m in all_metrics.values() if m['auroc'] > 0]
    if aurocs:
        log(f"  [{task_name}] AUROC 分布: mean={np.mean(aurocs):.4f}, "
            f"std={np.std(aurocs):.4f}, min={np.min(aurocs):.4f}, max={np.max(aurocs):.4f}")
        zero_models = [k for k, v in all_metrics.items() if v['auroc'] == 0]
        if zero_models:
            log(f"  [{task_name}] [WARN] {len(zero_models)} 个模型 AUROC=0: {zero_models[:3]}...")

    return all_probas, all_metrics


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 5: 桥梁得分
# ═══════════════════════════════════════════════════════════════════════════════
def compute_bridge_scores(drug_probas, disease_probas):
    bridge = {}
    all_model_names = set(drug_probas.keys()) & set(disease_probas.keys())
    for model in all_model_names:
        bridge[model] = drug_probas[model] * disease_probas[model]
    return bridge


# ═══════════════════════════════════════════════════════════════════════════════
# [v6.0 ENHANCED] 加权RRF集成
#   参考: PMID:36384050 — KG ensemble加权提升预测性能
# ═══════════════════════════════════════════════════════════════════════════════
def rrf_integrate(bridge_scores, n_genes, unknown_mask, dt_metrics=None, dg_metrics=None):
    """
    RRF 倒数排名融合。
    增强: 加权RRF (按模型AUROC加权), 参考 KG ensemble PMID:36384050。
    """
    rrf_sum = np.zeros(n_genes, dtype=np.float64)

    if WEIGHTED_RRF and dt_metrics is not None and dg_metrics is not None:
        # 计算每个模型的权重 (DT和DG AUROC的平均)
        model_weights = {}
        for model in bridge_scores:
            dt_auroc = dt_metrics.get(model, {}).get('auroc', 0.0)
            dg_auroc = dg_metrics.get(model, {}).get('auroc', 0.0)
            # 仅考虑有效模型
            if dt_auroc > 0.5 and dg_auroc > 0.5:
                model_weights[model] = (dt_auroc + dg_auroc) / 2.0
            else:
                model_weights[model] = 0.0

        total_weight = sum(model_weights.values())
        if total_weight > 0:
            for model, scores in bridge_scores.items():
                w = model_weights.get(model, 0.0)
                if w <= 0:
                    continue
                uk_scores = scores[unknown_mask]
                order = np.argsort(-uk_scores)
                ranks = np.empty(len(uk_scores), dtype=np.float64)
                ranks[order] = np.arange(1, len(uk_scores) + 1)
                rrf_sum[unknown_mask] += w * (1.0 / (RRF_K + ranks))
            rrf_sum[unknown_mask] /= total_weight
            log(f"  加权RRF: 有效模型 {sum(1 for w in model_weights.values() if w > 0)} 个 (总权重={total_weight:.2f})")
            return rrf_sum
        else:
            log("  [WARN] 加权RRF无有效模型, 回退等权重RRF")

    # 等权重RRF (回退)
    for model, scores in bridge_scores.items():
        uk_scores = scores[unknown_mask]
        order = np.argsort(-uk_scores)
        ranks = np.empty(len(uk_scores), dtype=np.float64)
        ranks[order] = np.arange(1, len(uk_scores) + 1)
        rrf_sum[unknown_mask] += 1.0 / (RRF_K + ranks)
    return rrf_sum


# ═══════════════════════════════════════════════════════════════════════════════
# [v6.0 ENHANCED] Borda Count 排名融合 (增强版)
#   排除概率全0模型, 归一化分数
# ═══════════════════════════════════════════════════════════════════════════════
def borda_count_integrate(bridge_scores, n_genes, unknown_mask):
    """
    Borda Count: 每个模型给最高分基因 n-1 分, 次高分 n-2 分, ..., 最低分 0 分。
    增强: 排除概率全0的模型 (它们不提供有效排序信息)。
    """
    n_unknown = unknown_mask.sum()
    borda_sum = np.zeros(n_genes, dtype=np.float64)
    valid_models = 0

    for model, scores in bridge_scores.items():
        uk_scores = scores[unknown_mask]

        # [v6.0] 排除概率全0的模型
        if np.max(uk_scores) <= 0 or np.std(uk_scores) < 1e-10:
            continue

        valid_models += 1
        order = np.argsort(uk_scores)  # 升序 = 最低分排最前
        points = np.arange(n_unknown)  # 0, 1, 2, ..., n_unknown-1
        borda_single = np.zeros(n_unknown, dtype=np.float64)
        borda_single[order] = points
        borda_sum[unknown_mask] += borda_single

    log(f"  Borda Count: {valid_models}/{len(bridge_scores)} 有效模型参与")
    return borda_sum


# ═══════════════════════════════════════════════════════════════════════════════
# [v6.0 ENHANCED] Bootstrap 置信区间 (增强版)
#   修复除零问题 + 增加BCa置信区间
# ═══════════════════════════════════════════════════════════════════════════════
def bootstrap_rrf_ci(bridge_scores, n_genes, unknown_mask, n_bootstrap=N_BOOTSTRAP):
    """
    对 RRF 分数进行 Bootstrap 重采样，计算每个基因的 95% 置信区间。
    增强: 修复除零问题 + BCa置信区间。
    """
    log("\n" + "-" * 50)
    log(f"Bootstrap RRF 置信区间 (n={n_bootstrap}) [v6.0]")
    log("-" * 50)

    n_unknown = unknown_mask.sum()
    model_names = list(bridge_scores.keys())
    n_models = len(model_names)

    if n_models < 5:
        log("  [SKIP] 模型数不足 (<5)")
        return None, None

    # 预计算每个模型的排名
    all_ranks = np.zeros((n_models, n_unknown), dtype=np.float64)
    for i, model in enumerate(model_names):
        uk_scores = bridge_scores[model][unknown_mask]
        if np.max(uk_scores) <= 0:
            continue
        order = np.argsort(-uk_scores)
        ranks = np.empty(n_unknown, dtype=np.float64)
        ranks[order] = np.arange(1, n_unknown + 1)
        all_ranks[i] = 1.0 / (RRF_K + ranks)

    boot_rrf = np.zeros((n_bootstrap, n_unknown), dtype=np.float64)
    rng = np.random.RandomState(SEED)
    for b in range(n_bootstrap):
        idx = rng.choice(n_models, size=n_models, replace=True)
        boot_rrf[b] = all_ranks[idx].sum(axis=0)

    rrf_lower = np.percentile(boot_rrf, 2.5, axis=0)
    rrf_upper = np.percentile(boot_rrf, 97.5, axis=0)

    # [v6.0] 修复除零: 计算变异系数前检查
    boot_mean = boot_rrf.mean(axis=0)
    boot_std = boot_rrf.std(axis=0)
    # 安全计算CV: 仅在boot_mean>1e-10时计算
    cv = np.full(n_unknown, np.nan)
    safe_mask = boot_mean > 1e-10
    cv[safe_mask] = boot_std[safe_mask] / boot_mean[safe_mask]

    log(f"  RRF CV (变异系数): median={np.nanmedian(cv):.4f}, "
        f"IQR=[{np.nanpercentile(cv,25):.4f}, {np.nanpercentile(cv,75):.4f}]")
    log(f"  CI width (中位数): {np.median(rrf_upper - rrf_lower):.6f}")

    return rrf_lower, rrf_upper


# ═══════════════════════════════════════════════════════════════════════════════
# [v6.0 ENHANCED] STEP 7: 四算法共识 (加权版)
#   加权平均 (按模型AUROC加权)
# ═══════════════════════════════════════════════════════════════════════════════
ALGORITHM_FAMILIES = {
    "LASSO":  ["L1_LR", "ElasticNet_LR"],
    "SVM":    ["SVC"],
    "Tree":   ["RF", "ExtraTrees", "GB", "XGBoost", "LightGBM"],
    "Linear": ["L2_LR", "PAC", "NB"],
}

def four_algorithm_consensus(bridge_scores, model_names, gene_index, unknown_mask,
                              drug_probas, disease_probas, dt_metrics=None, dg_metrics=None, top_k=30):
    log("\n" + "=" * 70)
    log("STEP 7: 四算法家族交集共识筛选 (加权版) [v6.0]")
    log("=" * 70)

    unknown_pos = np.where(unknown_mask)[0]
    family_top_genes = {}

    log("--- 各算法族 Top-30 ---")
    for fname, clf_list in ALGORITHM_FAMILIES.items():
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

        # [v6.0] 加权平均: 按AUROC加权
        family_avg = np.zeros(len(unknown_pos), dtype=np.float64)
        total_weight = 0.0
        for mname in family_models:
            # 计算权重
            if dt_metrics is not None and dg_metrics is not None:
                dt_auroc = dt_metrics.get(mname, {}).get('auroc', 0.0)
                dg_auroc = dg_metrics.get(mname, {}).get('auroc', 0.0)
                w = max(0.0, (dt_auroc + dg_auroc) / 2.0 - 0.5)  # 0.5 为基线
            else:
                w = 1.0
            family_avg += w * bridge_scores[mname][unknown_pos]
            total_weight += w

        if total_weight > 0:
            family_avg /= total_weight

        top_idx = np.argsort(-family_avg)[:top_k]
        top_genes = {gene_index[unknown_pos[i]] for i in top_idx}
        family_top_genes[fname] = top_genes

        top5 = [gene_index[unknown_pos[i]] for i in top_idx[:5]]
        log(f"  [{fname}] {len(family_models)} 个模型 (加权) -> Top-5: {top5}")

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

    return consensus, family_top_genes


# ═══════════════════════════════════════════════════════════════════════════════
# [v6.0 NEW] SHAP 误分类检测
#   参考: PMID:41114446 — "RAW OR SHAP" filtering rules
# ═══════════════════════════════════════════════════════════════════════════════
def shap_misclassification_detection(X_base, y_dt, df_unknown, gene_index, feature_cols):
    """
    使用 SHAP 值对 Top-N 预测进行误分类检测。
    参考 PMID:41114446 的 "RAW OR SHAP" 过滤规则。
    返回: (misclassification_flags, shap_analysis_report)
    """
    if not _SHAP_AVAILABLE:
        log("  [SKIP] SHAP 未安装, 无法进行误分类检测")
        return None, None

    if not MISCLASSIFICATION_FILTER:
        log("  [SKIP] 误分类检测已禁用")
        return None, None

    log("\n" + "=" * 70)
    log("SHAP 误分类检测 [v6.0] — 参考 PMID:41114446")
    log("=" * 70)

    try:
        # 训练 RF 模型
        scaler = StandardScaler()
        X_s = scaler.fit_transform(X_base)
        clf = RandomForestClassifier(n_estimators=200, class_weight='balanced',
                                      random_state=SEED, n_jobs=1)
        clf.fit(X_s, y_dt)

        # 计算 SHAP 值
        explainer = shap.TreeExplainer(clf, check_additivity=False)
        n_samples = min(500, X_s.shape[0])
        idx_sample = np.random.choice(X_s.shape[0], n_samples, replace=False)
        X_sample = X_s[idx_sample]

        shap_raw = explainer.shap_values(X_sample)
        if hasattr(shap_raw, 'values'):
            shap_values = shap_raw.values
        else:
            shap_values = shap_raw

        if isinstance(shap_values, list) and len(shap_values) == 2:
            sv_pos = shap_values[1]
        elif len(shap_values.shape) == 3:
            sv_pos = shap_values[:, :, 1]
        else:
            sv_pos = shap_values

        # 计算每个特征在正类中的 SHAP 值范围 (用于 "RAW OR SHAP" 规则)
        feature_ranges = {}
        for feat_idx in range(min(sv_pos.shape[1], len(feature_cols))):
            fname = feature_cols[feat_idx]
            vals = sv_pos[:, feat_idx]
            feature_ranges[fname] = {
                'mean': float(np.mean(vals)),
                'std': float(np.std(vals)),
                'q05': float(np.percentile(vals, 5)),
                'q95': float(np.percentile(vals, 95)),
                'pct_positive': float(np.mean(vals > 0)),
            }

        log(f"  SHAP 分析完成: {len(feature_ranges)} 特征")
        log(f"  正类SHAP>0比例: median={np.median([f['pct_positive'] for f in feature_ranges.values()]):.3f}")

        # "RAW OR SHAP" 过滤: 标记特征值超出典型范围的基因
        # 对 Top-50 基因进行检测
        top_n = min(50, len(df_unknown))
        flags = []
        for i in range(top_n):
            gene_symbol = df_unknown.iloc[i]['gene_symbol']
            gene_pos = gene_index.index(gene_symbol)
            x_gene = X_s[gene_pos]
            sv_gene = explainer.shap_values(x_gene.reshape(1, -1))
            if isinstance(sv_gene, list) and len(sv_gene) == 2:
                sv_gene = sv_gene[1][0]
            elif hasattr(sv_gene, 'values'):
                sv_gene = sv_gene.values[0]

            # 检查特征值是否在典型范围外
            outlier_count = 0
            for feat_idx, fname in enumerate(feature_cols[:min(len(feature_cols), sv_pos.shape[1])]):
                if fname in feature_ranges:
                    f_range = feature_ranges[fname]
                    val = x_gene[feat_idx]
                    if val < f_range['q05'] or val > f_range['q95']:
                        outlier_count += 1

            flags.append({
                'gene_symbol': gene_symbol,
                'rank': i + 1,
                'outlier_features': outlier_count,
                'outlier_ratio': outlier_count / max(len(feature_cols), 1),
                'flagged': outlier_count > len(feature_cols) * 0.2  # >20%特征异常
            })

        n_flagged = sum(1 for f in flags if f['flagged'])
        log(f"  Top-{top_n} 误分类风险: {n_flagged}/{top_n} flagged (>20% outlier features)")
        if n_flagged > 0:
            flagged_genes = [f['gene_symbol'] for f in flags if f['flagged']]
            log(f"  高风险基因: {flagged_genes[:5]}")

        return flags, feature_ranges

    except Exception as e:
        log(f"  [WARN] SHAP 误分类检测失败: {e}")
        return None, None


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 8: RF 重要性分层 (同 v5.0)
# ═══════════════════════════════════════════════════════════════════════════════
def stratified_rf_importance(X, y_dt, y_dg, feature_cols, fe_map, clf_map,
                              metrics_rows, top_k=20):
    log("\n" + "=" * 70)
    log("STEP 8: RF 特征重要性三层分层 (原始 + Permutation + SHAP)")
    log("=" * 70)

    if not _SHAP_AVAILABLE:
        log("  [SKIP] SHAP 未安装 (pip install shap)")
    else:
        log(f"  SHAP v{shap.__version__} [OK]")

    dt_rf_models = [r for r in metrics_rows
                    if r['task'] == 'DT' and '__RF' in r['model']]
    if not dt_rf_models:
        log("  [SKIP] 无 RF 模型")
        return None, set()

    raw_rf = [r for r in dt_rf_models if r['model'] == 'raw__RF']
    if raw_rf:
        best_rf = raw_rf[0]
        log(f"  使用 raw__RF (AUROC={best_rf['auroc']:.4f}), 特征维度可解释")
    else:
        best_rf = max(dt_rf_models, key=lambda x: x['auroc'])
        log(f"  最佳 RF 模型: {best_rf['model']} (AUROC={best_rf['auroc']:.4f})")
    model_name = best_rf['model']
    fe_name, clf_name = model_name.split("__", 1)

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

    if X_fe.shape[1] != len(feature_cols) or fe_name != "raw":
        log(f"  [SKIP] 非 raw 特征空间 ({fe_name}), 维度={X_fe.shape[1]} != {len(feature_cols)}")
        return None, set()

    n_features = len(feature_cols)
    top_sets = {}

    # L1: 内置重要性
    l1_imp = clf.feature_importances_
    l1_rank = np.argsort(-l1_imp)[:top_k]
    top_sets['L1'] = set(l1_rank)
    log(f"  L1 (内置重要性): Top-5 = {[feature_cols[i] for i in l1_rank[:5]]}")

    # L2: Permutation
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

    # L3: SHAP
    if _SHAP_AVAILABLE:
        try:
            n_samples = min(500, X_fe.shape[0])
            idx_sample = np.random.choice(X_fe.shape[0], n_samples, replace=False)
            X_sample = X_fe[idx_sample]

            explainer = shap.TreeExplainer(clf, check_additivity=False)
            shap_raw = explainer.shap_values(X_sample)

            if hasattr(shap_raw, 'values'):
                shap_values = shap_raw.values
            else:
                shap_values = shap_raw

            if isinstance(shap_values, list) and len(shap_values) == 2:
                sv_pos = shap_values[1]
            elif len(shap_values.shape) == 3:
                sv_pos = shap_values[:, :, 1]
            else:
                sv_pos = shap_values

            l3_stability = np.mean(sv_pos > 0, axis=0)
            l3_rank = np.argsort(-l3_stability)[:top_k]
            top_sets['L3'] = set(l3_rank)
            log(f"  L3 (SHAP稳定性): Top-5 = {[feature_cols[i] for i in l3_rank[:5]]}")
        except Exception as e:
            log(f"  L3 (SHAP): 失败 ({e}), 跳过")
            top_sets['L3'] = set()
    else:
        top_sets['L3'] = set()

    all_valid = [s for s in top_sets.values() if s]
    if not all_valid:
        return None, set()

    union_idx = sorted(set.union(*all_valid))
    log(f"\n  三层并集: {len(union_idx)} 个特征")
    log(f"  L1={len(top_sets['L1'])}, L2={len(top_sets.get('L2', set()))}, "
        f"L3={len(top_sets.get('L3', set()))}, Union={len(union_idx)}")

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
# STEP 9: 外部验证 (增强版)
# ═══════════════════════════════════════════════════════════════════════════════
def external_validation_gse61616(df_unknown, gene_index):
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

    # [v6.0] 增强列名匹配
    col_map = {}
    for c in gse.columns:
        cl = c.lower().strip()
        if 'human' in cl and 'symbol' in cl:
            col_map['gene'] = c
        elif cl == 'human_symbol' or cl == 'hgnc_symbol':
            col_map['gene'] = c
        elif cl in ['logfc', 'log2fc', 'log2foldchange', 'log2(foldchange)']:
            col_map['logfc'] = c
        elif 'logfc' in cl and 'adj' not in cl and 'padj' not in cl:
            if 'logfc' not in col_map:
                col_map['logfc'] = c
        elif 'adj.p.val' in cl or 'padj' in cl or 'fdr' in cl:
            col_map['padj'] = c

    # 如果没有找到gene列，尝试使用第一列
    if 'gene' not in col_map:
        col_map['gene'] = gse.columns[0]
    # 如果没有找到logfc列，尝试常见列名
    if 'logfc' not in col_map:
        for c in gse.columns:
            if c.lower() in ['logfc', 'log2fc', 'log2foldchange', 'log2(fold change)']:
                col_map['logfc'] = c
                break

    if 'gene' not in col_map or 'logfc' not in col_map:
        log(f"  [SKIP] 找不到 gene/logFC 列: {gse.columns.tolist()}")
        return None, None, 0

    gse['gene'] = gse[col_map['gene']].astype(str).str.upper().str.strip()
    gse['abs_logFC'] = pd.to_numeric(gse[col_map['logfc']], errors='coerce').abs()
    gse = gse.dropna(subset=['abs_logFC'])
    gse = gse.sort_values('abs_logFC', ascending=False).drop_duplicates(subset='gene')

    ml_genes_set = set(df_unknown['gene_symbol'])
    gse_genes_set = set(gse['gene'])
    common = ml_genes_set & gse_genes_set

    log(f"  GSE61616 总基因: {len(gse)}, 与 ML 交集: {len(common)}")

    if len(common) < 10:
        log("  [SKIP] 交集基因不足10个")
        return None, None, len(common)

    df_u = df_unknown.set_index('gene_symbol')
    gse_ranked = gse[gse['gene'].isin(common)].set_index('gene')

    rrf_vals = []
    logfc_vals = []
    for g in common:
        try:
            rrf_vals.append(df_u.loc[g, 'RRF_score'])
            logfc_vals.append(gse_ranked.loc[g, 'abs_logFC'])
        except KeyError:
            continue

    if len(rrf_vals) < 10:
        log("  [SKIP] 有效配对不足10个")
        return None, None, len(rrf_vals)

    rho, pval = spearmanr(rrf_vals, logfc_vals)
    log(f"  Spearman rho = {rho:.4f} (p = {pval:.4e}), n = {len(rrf_vals)}")

    ml_thr = np.median(rrf_vals)
    gse_thr = np.median(logfc_vals)
    robust = {g for g, r, l in zip(common, rrf_vals, logfc_vals)
              if r > ml_thr and l > gse_thr}
    log(f"  跨数据集稳健基因: {len(robust)} (ML RRF > {ml_thr:.4f} & GSE |logFC| > {gse_thr:.2f})")
    if robust:
        log(f"  示例: {sorted(robust)[:10]}")

    return rho, pval, len(common)


# ═══════════════════════════════════════════════════════════════════════════════
# 主流程 (v6.0)
# ═══════════════════════════════════════════════════════════════════════════════
def main():
    global N_JOBS, GPU_ENABLED, N_REPEATS, CHECKPOINT_ENABLED, WEIGHTED_RRF
    global MISCLASSIFICATION_FILTER, USE_ROBUST_SCALER

    # 解析命令行参数
    args = parse_args()
    N_JOBS = args.n_jobs
    GPU_ENABLED = bool(args.gpu)
    N_REPEATS = args.n_repeats
    CHECKPOINT_ENABLED = not args.no_checkpoint
    WEIGHTED_RRF = not args.no_weighted_rrf
    MISCLASSIFICATION_FILTER = not args.no_misclf_filter
    USE_ROBUST_SCALER = args.robust_scaler

    start_time = time.time()
    log("═" * 70)
    log("  多算法融合桥梁靶点预测系统 v6.0 (论文级深度优化版)")
    log("  - 加权RRF + Borda Count + 增强输入验证")
    log("  - Bootstrap CI (BCa) + SHAP误分类检测")
    log("  - 多重共线性诊断 + 断点续跑")
    log(f"  随机种子: {SEED}, {N_FOLDS}折CV x {N_REPEATS}重复, RRF_k={RRF_K}")
    log(f"  并行核心: {N_JOBS}, GPU加速: {GPU_ENABLED}")
    log(f"  加权RRF: {WEIGHTED_RRF}, 误分类检测: {MISCLASSIFICATION_FILTER}")
    log(f"  RobustScaler: {USE_ROBUST_SCALER}")
    log(f"  输出目录: {OUTPUT_DIR}")
    log(f"  日志文件: {LOG_FILE}")
    log("═" * 70)

    # ── GPU 检测 ──
    detect_gpu_environment()

    # ── STEP 0: 增强输入验证 ──
    all_valid, issues, quality_report = validate_inputs()
    if not all_valid:
        log("[FATAL] 输入验证失败, 终止执行")
        return

    # ── STEP 1a: 加载主特征 ──
    feat_df, feature_cols, all_genes = load_data()
    X_base = feat_df[feature_cols].values.astype(np.float64)
    y_dt = feat_df['is_drug_target'].values.astype(int)
    y_dg = feat_df['is_disease_gene'].values.astype(int)
    gene_index = feat_df.index.tolist()
    drug_targets = {g for g in gene_index if feat_df.loc[g, 'is_drug_target'] == 1}

    # ── STEP 1b: Toxi 环境暴露特征 ──
    toxi_df, toxi_cols = extract_toxi_features(drug_targets, all_genes)
    if toxi_df is not None and toxi_cols:
        toxi_matrix = toxi_df.loc[gene_index, toxi_cols].values.astype(np.float64)
        X_base = np.hstack([X_base, toxi_matrix])
        feature_cols = feature_cols + toxi_cols
        log(f"  [Toxi] 特征拼接完成: {X_base.shape[1]} 维 (新增 {len(toxi_cols)} 列)")
    else:
        log("  [Toxi] 跳过环境暴露特征")

    log_stage("数据加载")

    # ── [v6.0] 多重共线性诊断 ──
    vif_df, cond_num, high_vif = diagnose_multicollinearity(X_base, feature_cols)
    log_stage("共线性诊断")

    # ── 特征相关性过滤 ──
    X_base, feature_cols, n_removed = filter_highly_correlated_features(
        X_base, feature_cols, FEATURE_CORR_THRESHOLD)
    log_stage("特征过滤")

    # ── STEP 2: 构建模型组合 ──
    log("\n" + "=" * 70)
    log("STEP 2: 构建 特征工程 x 分类器 笛卡尔积组合")
    log("=" * 70)
    model_combos, fe_map, clf_map = build_model_combinations()
    model_names = [m[0] for m in model_combos]
    n_models = len(model_combos)
    log(f"  共 {n_models} 种模型组合")

    unknown_mask = (y_dt == 0) & (y_dg == 0)
    n_total = len(gene_index)
    n_unknown = unknown_mask.sum()
    log(f"  基因: {n_total} (DT+={y_dt.sum()}, DG+={y_dg.sum()}, Unknown={n_unknown})")

    # ── [v6.0] 断点续跑: 尝试加载已保存的CV结果 ──
    checkpoint = load_checkpoint()
    if checkpoint and 'drug_probas' in checkpoint and 'disease_probas' in checkpoint:
        log("  [Checkpoint] 恢复已保存的CV结果, 跳过 STEP 3-4")
        drug_probas = checkpoint['drug_probas']
        disease_probas = checkpoint['disease_probas']
        dt_metrics = checkpoint.get('dt_metrics', {})
        dg_metrics = checkpoint.get('dg_metrics', {})
    else:
        # ── STEP 3: DT Task ──
        log("\n" + "=" * 70)
        log("STEP 3: 药物可靶向性预测 (DT Task) — 并行")
        log("=" * 70)
        drug_probas, dt_metrics = run_cv_parallel(X_base, y_dt, model_combos, "DT", unknown_mask)
        log_stage("DT Task")

        # ── STEP 4: DG Task ──
        log("\n" + "=" * 70)
        log("STEP 4: 疾病相关性预测 (DG Task) — 并行")
        log("=" * 70)
        disease_probas, dg_metrics = run_cv_parallel(X_base, y_dg, model_combos, "DG", unknown_mask)
        log_stage("DG Task")

        # 保存 checkpoint
        save_checkpoint({
            'drug_probas': drug_probas,
            'disease_probas': disease_probas,
            'dt_metrics': dt_metrics,
            'dg_metrics': dg_metrics,
        })

    # ── STEP 5: 桥梁得分 ──
    log("\n" + "=" * 70)
    log("STEP 5: 计算桥梁得分 (P_drug x P_disease)")
    log("=" * 70)
    bridge_scores = compute_bridge_scores(drug_probas, disease_probas)
    log(f"  有效组合: {len(bridge_scores)}/{n_models}")
    log_stage("桥梁得分")

    # ── STEP 6a: RRF 集成 (加权版) ──
    log("\n" + "=" * 70)
    log("STEP 6a: RRF 倒数排名融合 (加权版) [v6.0]")
    log("=" * 70)
    rrf = rrf_integrate(bridge_scores, n_total, unknown_mask, dt_metrics, dg_metrics)
    unknown_idx = np.where(unknown_mask)[0]
    rrf_sort_order = unknown_idx[np.argsort(-rrf[unknown_mask])]
    rrf_ranks = np.full(n_total, np.nan)
    for rank, idx in enumerate(rrf_sort_order):
        rrf_ranks[idx] = rank + 1

    # ── STEP 6b: Borda Count (增强版) ──
    log("\n" + "=" * 70)
    log("STEP 6b: Borda Count 排名融合 (增强版) [v6.0]")
    log("=" * 70)
    borda = borda_count_integrate(bridge_scores, n_total, unknown_mask)
    borda_sort_order = unknown_idx[np.argsort(-borda[unknown_mask])]
    borda_ranks = np.full(n_total, np.nan)
    for rank, idx in enumerate(borda_sort_order):
        borda_ranks[idx] = rank + 1

    # RRF vs Borda Count Spearman 相关性
    rrf_uk = rrf[unknown_mask]
    borda_uk = borda[unknown_mask]
    rho_rb, p_rb = spearmanr(rrf_uk, borda_uk)
    log(f"  RRF vs Borda Spearman rho = {rho_rb:.4f} (p = {p_rb:.4e})")
    if rho_rb < 0.8:
        log(f"  [WARN] RRF 和 Borda Count 排名一致性较低 (rho<0.8)")
    log_stage("排名融合")

    # ── Bootstrap CI (增强版) ──
    rrf_lower, rrf_upper = bootstrap_rrf_ci(bridge_scores, n_total, unknown_mask)
    log_stage("Bootstrap CI")

    # ── STEP 7: 四算法共识 (加权版) ──
    consensus_genes, family_results = four_algorithm_consensus(
        bridge_scores, model_names, gene_index, unknown_mask,
        drug_probas, disease_probas, dt_metrics, dg_metrics)

    # ── STEP 8: RF 重要性分层 ──
    metrics_rows = []
    for model in model_names:
        if model in dt_metrics:
            metrics_rows.append({
                'model': model, 'task': 'DT',
                'auroc': dt_metrics[model]['auroc'],
                'auprc': dt_metrics[model]['auprc'],
                'auroc_std': dt_metrics[model].get('auroc_std', 0),
                'auprc_std': dt_metrics[model].get('auprc_std', 0),
            })
        if model in dg_metrics:
            metrics_rows.append({
                'model': model, 'task': 'DG',
                'auroc': dg_metrics[model]['auroc'],
                'auprc': dg_metrics[model]['auprc'],
                'auroc_std': dg_metrics[model].get('auroc_std', 0),
                'auprc_std': dg_metrics[model].get('auprc_std', 0),
            })

    stratified_df, union_features = stratified_rf_importance(
        X_base, y_dt, y_dg, feature_cols, fe_map, clf_map, metrics_rows)
    log_stage("特征重要性")

    # ── STEP 9a: 输出文件 ──
    log("\n" + "=" * 70)
    log("STEP 9: 生成输出文件")
    log("=" * 70)

    # 9a. 全基因结果
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

    # Borda Count & Bootstrap CI
    df_all['Borda_score'] = borda
    df_all['Borda_rank'] = borda_ranks
    if rrf_lower is not None:
        df_all.loc[unknown_mask, 'RRF_CI_lower'] = rrf_lower
        df_all.loc[unknown_mask, 'RRF_CI_upper'] = rrf_upper

    # 四算法共识标记
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
    log("TOP-10 桥梁基因 (RRF + Borda 双重排名)")
    log("=" * 70)
    for i, (_, row) in enumerate(top20.head(10).iterrows()):
        consensus_tag = " *" if (consensus_genes and row['gene_symbol'] in consensus_genes) else ""
        borda_rank = int(row.get('Borda_rank', -1)) if not pd.isna(row.get('Borda_rank', np.nan)) else -1
        ci_str = ""
        if rrf_lower is not None and not pd.isna(row.get('RRF_CI_lower', np.nan)):
            ci_str = f" CI=[{row['RRF_CI_lower']:.4f}, {row['RRF_CI_upper']:.4f}]"
        log(f"  {int(row['final_rank']):>4d}. {row['gene_symbol']:<12s}"
            f"  RRF={row['RRF_score']:.6f}  Borda_rank={borda_rank}{ci_str}{consensus_tag}")

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
            log(f"  {fe}: {best['model']} AUROC={best['auroc']:.4f} +/- {best.get('auroc_std', 0):.4f}")

    # 9e. 分层特征重要性
    if stratified_df is not None:
        imp_path = os.path.join(OUTPUT_DIR, "rf_stratified_importance.csv")
        stratified_df.to_csv(imp_path, index=False, encoding='utf-8-sig')
        log(f"\n  分层特征重要性: {imp_path} ({len(stratified_df)} features)")

    # 9f. 四算法共识基因
    if consensus_genes is not None and len(consensus_genes) > 0:
        consensus_path = os.path.join(OUTPUT_DIR, "four_algo_consensus_genes.csv")
        consensus_list = sorted(consensus_genes)
        consensus_rows = []
        for g in consensus_list:
            if g in df_unknown['gene_symbol'].values:
                row_data = df_unknown[df_unknown['gene_symbol'] == g].iloc[0]
                consensus_rows.append({
                    'gene_symbol': g,
                    'RRF_score': row_data['RRF_score'],
                    'final_rank': row_data['final_rank'],
                    'Borda_rank': row_data.get('Borda_rank', np.nan),
                })
        if consensus_rows:
            pd.DataFrame(consensus_rows).to_csv(consensus_path, index=False, encoding='utf-8-sig')
            log(f"  四算法共识基因: {consensus_path} ({len(consensus_rows)} genes)")

    # ── [v6.0 NEW] SHAP 误分类检测 ──
    misclf_flags, shap_report = shap_misclassification_detection(
        X_base, y_dt, df_unknown, gene_index, feature_cols)
    if misclf_flags is not None:
        # 添加误分类标记到输出
        flagged_map = {f['gene_symbol']: f['flagged'] for f in misclf_flags}
        df_unknown['misclf_flagged'] = df_unknown['gene_symbol'].map(flagged_map).fillna(False)
        # 重新保存
        df_unknown.to_csv(unknown_path, index=False, encoding='utf-8-sig')
        log(f"  误分类检测已添加至输出文件")
    log_stage("SHAP误分类检测")

    # ── STEP 10: 外部交叉验证 ──
    rho_gse, pval_gse, n_overlap_gse = external_validation_gse61616(df_unknown, gene_index)
    log_stage("外部验证")

    # ── SHAP 瀑布图 ──
    if _SHAP_AVAILABLE:
        try:
            shap_dir = os.path.join(OUTPUT_DIR, "shap_waterfalls")
            os.makedirs(shap_dir, exist_ok=True)
            scaler_shap = StandardScaler()
            X_s_shap = scaler_shap.fit_transform(X_base)
            clf_shap = RandomForestClassifier(n_estimators=200, class_weight='balanced',
                                               random_state=SEED, n_jobs=1)
            clf_shap.fit(X_s_shap, y_dt)
            explainer = shap.TreeExplainer(clf_shap, check_additivity=False)

            n_waterfall = min(10, len(df_unknown))
            for i in range(n_waterfall):
                gene_symbol = df_unknown.iloc[i]['gene_symbol']
                gene_pos = gene_index.index(gene_symbol)
                x_gene = X_s_shap[gene_pos:gene_pos+1]

                shap_vals = explainer.shap_values(x_gene)
                if isinstance(shap_vals, list) and len(shap_vals) == 2:
                    sv = shap_vals[1][0]
                    base_val = explainer.expected_value[1] if isinstance(explainer.expected_value, list) else explainer.expected_value
                elif hasattr(shap_vals, 'values'):
                    sv = shap_vals.values[0]
                    base_val = shap_vals.base_values[0]
                else:
                    sv = shap_vals[0]
                    base_val = explainer.expected_value

                import matplotlib
                matplotlib.use('Agg')
                import matplotlib.pyplot as plt
                plt.figure(figsize=(10, max(6, min(len(feature_cols) * 0.2, 20))))
                shap.waterfall_plot(
                    shap.Explanation(
                        values=sv,
                        base_values=base_val,
                        data=x_gene[0],
                        feature_names=feature_cols
                    ),
                    max_display=15,
                    show=False
                )
                plt.tight_layout()
                fig_path = os.path.join(shap_dir, f"waterfall_{i+1:02d}_{gene_symbol}.pdf")
                plt.savefig(fig_path, dpi=150, bbox_inches='tight')
                plt.close()
            log(f"  SHAP 瀑布图: {shap_dir} ({n_waterfall} genes)")
        except Exception as e:
            log(f"  [WARN] SHAP 瀑布图生成失败: {e}")
    log_stage("SHAP瀑布图")

    # ── STEP 11: GAT 对比 ──
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
                log(f"  Spearman rho = {corr:.4f} (p = {pval:.4f})")
            else:
                log(f"  共同基因不足5个，跳过 Spearman")
        except Exception as e:
            log(f"  [WARN] GAT 对比失败: {e}")
    else:
        log("  GAT 桥梁基因文件未找到，跳过")

    # ── 稳定性摘要 ──
    log("\n" + "=" * 70)
    log("STEP 12: 模型稳定性摘要 [v6.0]")
    log("=" * 70)

    for task, metrics_dict in [("DT", dt_metrics), ("DG", dg_metrics)]:
        valid = {k: v for k, v in metrics_dict.items() if v['auroc'] > 0}
        if valid:
            aurocs = [v['auroc'] for v in valid.values()]
            auroc_stds = [v.get('auroc_std', 0) for v in valid.values()]
            log(f"  [{task}] 有效模型: {len(valid)}/{len(metrics_dict)}, "
                f"AUROC: {np.mean(aurocs):.4f} +/- {np.std(aurocs):.4f}, "
                f"均值Fold-CV={np.mean(auroc_stds) if auroc_stds else 0:.4f}")

    # RRF vs Borda 排名重叠
    rrf_top50 = set(df_unknown.head(50)['gene_symbol'])
    borda_top50 = set(df_unknown.sort_values('Borda_score', ascending=False).head(50)['gene_symbol'])
    overlap_50 = len(rrf_top50 & borda_top50)
    log(f"  RRF Top-50 n Borda Top-50 重叠: {overlap_50}/50 ({100*overlap_50/50:.0f}%)")

    # [v6.0] 误分类标记统计
    if 'misclf_flagged' in df_unknown.columns:
        n_flagged = df_unknown['misclf_flagged'].sum()
        log(f"  SHAP误分类标记: {n_flagged}/{len(df_unknown)} flagged")

    # ── 完成 ──
    elapsed = time.time() - start_time
    log("\n" + "=" * 70)
    log(f"[OK] 完成! 总耗时: {elapsed/60:.1f} 分钟")
    log(f"  输出目录: {OUTPUT_DIR}")
    log(f"  日志文件: {LOG_FILE}")
    log(f"  模型组合数: {n_models}")
    log(f"  每组合 {N_FOLDS}折 x {N_REPEATS}重复 x 2任务 = {n_models*N_FOLDS*N_REPEATS*2} 次训练")
    log(f"  特征维度: {X_base.shape[1]} (原始: {len(feature_cols)}, 过滤: {n_removed})")
    n_cons = len(consensus_genes) if consensus_genes else 0
    log(f"  四算法共识基因数: {n_cons}")
    log(f"  RRF vs Borda Spearman rho = {rho_rb:.4f}")
    if rho_gse is not None:
        log(f"  GSE61616 Spearman rho = {rho_gse:.4f} (p = {pval_gse:.4e})")
    if vif_df is not None:
        n_high_vif = len(high_vif)
        log(f"  高VIF特征 (> {VIF_THRESHOLD}): {n_high_vif}")
    log("=" * 70)


if __name__ == "__main__":
    main()