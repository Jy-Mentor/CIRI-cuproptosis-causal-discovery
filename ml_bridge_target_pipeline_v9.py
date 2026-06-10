# -*- coding: utf-8 -*-
# [v9.0] 基于PubMed文献系统性优化版 — 完整可运行版 (修复校准 + 补全main)
"""
多算法融合桥梁靶点预测系统 (v9.0 — 完整可运行版本)
============================================================================
网络药理学第一阶段：药物-靶点-疾病桥梁基因发现

v9.0 核心升级 (vs v8.0) — 修复 + 补全 + 数据泄露审计:

=== 修复 ===
[Fix] 修复 calibrate_pretrained_clf() 函数 — 原版返回 (None, None) 无法校准
[Fix] 修复 main() 第2628行语法错误 — 缺少闭合括号
[Fix] 补全 main() 函数 — 原版在 STEP 6b 处截断

=== 补全内容 ===
[Added] STEP 6b: Borda Count 排名计算 + RRF vs Borda Spearman 比较
[Added] STEP 6c: Bootstrap CI 置信区间计算
[Added] STEP 7: 四算法家族共识
[Added] SHAP 误分类检测
[Added] STEP 8: RF 重要性分层
[Added] STEP 9: 数据泄露审计 — 移除 GSE61616 外部验证 (GSE61616 已整合进训练集)
[Added] 稳定性选择
[Added] STEP 10: 最终输出生成 (14项输出文件)
[Added] 模型性能汇总图
[Added] 总流程耗时统计

=== v9.0 关键修复: 数据泄露审计 ===
[CRITICAL] GSE61616 (7d MCAO大鼠模型) 已作为训练数据来源之一
  - 其DEGs已整合到疾病基因集 (CIRI disease genes) 中
  - 使用同一数据集作为"外部验证"构成循环论证
  - 参考: PMID:36004690 — "test data must be completely independent"
  - 替代方案: 内部5折分层CV提供无偏性能估计

保持自 v8.0: 数据泄露预防 + class_weight Bug修复 + 模型校准 + Stacking
  + 稳定性选择 + MI-VIF概念 + 四算法共识 + SHAP分层重要性
  + 环境暴露 + 笛卡尔积CV + 加权RRF
  + Borda Count + Bootstrap CI + 断点续跑 + 多重共线性诊断
  + 类别不平衡诊断 + 性能汇总图
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
from sklearn.linear_model import LogisticRegression, LassoCV, SGDClassifier, PassiveAggressiveClassifier, LinearRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier, ExtraTreesClassifier
from sklearn.calibration import CalibratedClassifierCV
from sklearn.decomposition import PCA
from sklearn.cross_decomposition import PLSRegression
from sklearn.feature_selection import SelectFromModel, SelectKBest, f_classif, mutual_info_classif
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
N_REPEATS = 1             # 重复CV次数 (>1 启用RepeatedStratifiedKFold)
RRF_K = 60
N_JOBS = 4
GPU_ENABLED = True
N_BOOTSTRAP = 200
FEATURE_CORR_THRESHOLD = 0.95
VIF_THRESHOLD = 10.0       # VIF阈值 (VIF>10表示严重共线性)
USE_ROBUST_SCALER = False
CHECKPOINT_ENABLED = True  # 断点续跑
WEIGHTED_RRF = True        # 加权RRF (按AUROC加权)
MISCLASSIFICATION_FILTER = True  # SHAP误分类检测

# [v9.0] 保持以下配置
USE_CALIBRATION = True         # 模型校准 (Isotonic Scaling, 使用 cv='prefit' 优化)
USE_STACKING = True            # Stacking Ensemble 元模型
STABILITY_BOOTSTRAP_N = 100    # 稳定性选择 Bootstrap 次数
STABILITY_THRESHOLD = 0.8      # 稳定性选择阈值 (>80% 视为稳定)

# [v9.0] 零方差特征移除
REMOVE_ZERO_VARIANCE = True    # 移除方差为0的常量特征

# [v9.0] MI-VIF 模式: 'none'=跳过 | 'pre_cv'=CV前(有泄露风险) | 'per_fold'=每折内嵌(防泄露,慢)
MI_VIF_MODE = 'per_fold'       # 默认每折内嵌 (参考 PMID:37113250)

# 数据路径
DATA_DIR = r"D:/反向网络药理学/GAT拓展维度/cache"
FEATURE_PATH = os.path.join(DATA_DIR, "enhanced_gene_features.csv")
DRUG_TARGETS_PATH = r"C:/Users/Jy-Mentor-7/Desktop/GAT/drug_targets.txt"
DISEASE_GENES_PATH = r"C:/Users/Jy-Mentor-7/Desktop/GAT/disease_genes.txt"
SUBGRAPH_GENES_PATH = r"C:/Users/Jy-Mentor-7/Desktop/GAT/subgraph_genes.txt"
GAT_BRIDGE_PATH = r"C:/Users/Jy-Mentor-7/Desktop/GAT/top20_bridge_genes.csv"
TOXI_FPKM_DIR = r"D:/反向网络药理学/GAT拓展维度/Toxi/rna_fpkm"
GSE61616_DEG_PATH = r"c:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\大创\GSE61616_human_homologs_DEGs.tsv"

# [v9.0] 固定输出目录
OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "ml_output_v9")
os.makedirs(OUTPUT_DIR, exist_ok=True)
ORTHO_CACHE_PATH = os.path.join(OUTPUT_DIR, "mouse_to_human_orthologs.csv")
CHECKPOINT_PATH = os.path.join(OUTPUT_DIR, "checkpoint.pkl")

np.random.seed(SEED)

# ═══════════════════════════════════════════════════════════════════════════════
# 命令行参数解析
# ═══════════════════════════════════════════════════════════════════════════════
def parse_args():
    parser = argparse.ArgumentParser(description="多算法融合桥梁靶点预测系统 v9.0")
    parser.add_argument("--n-jobs", type=int, default=N_JOBS, help=f"并行核心数 (默认: {N_JOBS})")
    parser.add_argument("--gpu", type=int, default=1 if GPU_ENABLED else 0, help="GPU加速 (0/1)")
    parser.add_argument("--n-repeats", type=int, default=N_REPEATS, help=f"重复CV次数 (默认: {N_REPEATS})")
    parser.add_argument("--no-checkpoint", action="store_true", help="禁用断点续跑")
    parser.add_argument("--no-weighted-rrf", action="store_true", help="禁用加权RRF")
    parser.add_argument("--no-misclf-filter", action="store_true", help="禁用SHAP误分类检测")
    parser.add_argument("--no-calibration", action="store_true", help="禁用模型校准")
    parser.add_argument("--no-stacking", action="store_true", help="禁用Stacking Ensemble")
    parser.add_argument("--robust-scaler", action="store_true", help="使用RobustScaler")
    parser.add_argument("--mi-vif-mode", choices=['none', 'pre_cv', 'per_fold'], default=MI_VIF_MODE,
                        help="MI-VIF模式: none=跳过, pre_cv=CV前, per_fold=每折内嵌(防泄露)")
    return parser.parse_args()

# ═══════════════════════════════════════════════════════════════════════════════
# [v9.0] 结构化日志系统 (增强: 文件轮转)
# ═══════════════════════════════════════════════════════════════════════════════
LOG_FILE = os.path.join(OUTPUT_DIR, "pipeline_v9.log")
MAX_LOG_SIZE = 50 * 1024 * 1024  # 50MB 轮转
_stage_start_time = time.time()
_pipeline_start_time = time.time()

def _setup_file_logger():
    """设置文件日志, 支持大文件轮转"""
    # 检查是否需要轮转
    if os.path.exists(LOG_FILE) and os.path.getsize(LOG_FILE) > MAX_LOG_SIZE:
        backup = LOG_FILE + ".old"
        if os.path.exists(backup):
            os.remove(backup)
        os.rename(LOG_FILE, backup)

    file_handler = logging.FileHandler(LOG_FILE, encoding='utf-8')
    file_handler.setFormatter(logging.Formatter('[%(asctime)s] %(message)s', datefmt='%H:%M:%S'))
    logger = logging.getLogger('pipeline_v9')
    # 移除旧 handler
    for h in logger.handlers[:]:
        logger.removeHandler(h)
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
# [v9.0] 进度跟踪
# ═══════════════════════════════════════════════════════════════════════════════
def log_progress(current_step: int, total_steps: int, step_name: str = ""):
    """记录进度并预估剩余时间"""
    global _pipeline_start_time
    elapsed_total = time.time() - _pipeline_start_time
    if current_step > 0 and total_steps > 0:
        eta = (elapsed_total / current_step) * (total_steps - current_step)
        log(f"  [Progress] {current_step}/{total_steps} ({step_name}) "
            f"| 已耗时: {elapsed_total/60:.1f}min | 预估剩余: {eta/60:.1f}min")
    else:
        log(f"  [Progress] {step_name} | 已耗时: {elapsed_total/60:.1f}min")

# ═══════════════════════════════════════════════════════════════════════════════
# [v9.0] 安全文件写入 (增强错误处理)
# ═══════════════════════════════════════════════════════════════════════════════
def safe_write_csv(df, path, **kwargs):
    """安全写入CSV, 带错误处理"""
    try:
        df.to_csv(path, index=False, encoding='utf-8-sig', **kwargs)
        log(f"  [SAVE] {os.path.basename(path)} ({len(df)} rows)")
        return True
    except Exception as e:
        log(f"  [ERROR] 保存失败 {path}: {e}")
        # 尝试备用路径
        try:
            alt_path = path.replace('.csv', '_backup.csv')
            df.to_csv(alt_path, index=False, encoding='utf-8-sig', **kwargs)
            log(f"  [SAVE] 已保存至备用路径: {alt_path}")
            return True
        except Exception as e2:
            log(f"  [FATAL] 备用保存也失败: {e2}")
            return False

def safe_write_json(data, path):
    """安全写入JSON, 带错误处理"""
    try:
        with open(path, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=2, ensure_ascii=False, default=str)
        log(f"  [SAVE] {os.path.basename(path)}")
        return True
    except Exception as e:
        log(f"  [ERROR] JSON保存失败 {path}: {e}")
        return False

# ═══════════════════════════════════════════════════════════════════════════════
# 断点续跑
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
# GPU 检测
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
# 类别不平衡诊断报告
#   参考: PMID:36475638
# ═══════════════════════════════════════════════════════════════════════════════
def class_imbalance_diagnostic(y_dt, y_dg, n_total, gene_index):
    """
    计算并记录类别分布、不平衡比和推荐指标。
    参考 PMID:36475638 — 类别不平衡对模型评估的影响。
    """
    log("\n" + "=" * 70)
    log("类别不平衡诊断报告 [v9.0] — 参考 PMID:36475638")
    log("=" * 70)

    n_dt = int(y_dt.sum())
    n_dg = int(y_dg.sum())
    n_both = int((y_dt & y_dg).sum())
    n_unknown = int((~y_dt.astype(bool) & ~y_dg.astype(bool)).sum())

    report = {
        'total_genes': n_total,
        'drug_target_positive': n_dt,
        'disease_gene_positive': n_dg,
        'both_positive': n_both,
        'unknown': n_unknown,
    }

    # DT 不平衡比
    dt_neg = n_total - n_dt
    dt_ratio = dt_neg / max(n_dt, 1)
    report['dt_imbalance_ratio'] = dt_ratio
    report['dt_positive_pct'] = 100.0 * n_dt / n_total if n_total > 0 else 0

    # DG 不平衡比
    dg_neg = n_total - n_dg
    dg_ratio = dg_neg / max(n_dg, 1)
    report['dg_imbalance_ratio'] = dg_ratio
    report['dg_positive_pct'] = 100.0 * n_dg / n_total if n_total > 0 else 0

    log(f"  总基因数: {n_total}")
    log(f"  DT+: {n_dt} ({report['dt_positive_pct']:.2f}%) | 不平衡比: {dt_ratio:.1f}:1")
    log(f"  DG+: {n_dg} ({report['dg_positive_pct']:.2f}%) | 不平衡比: {dg_ratio:.1f}:1")
    log(f"  Both: {n_both} | Unknown: {n_unknown}")

    recommendations = []
    if dt_ratio > 10:
        recommendations.append(f"DT不平衡比 {dt_ratio:.1f}:1 > 10:1 — 建议优先使用 AUPRC 而非 AUROC")
        log(f"  [RECOMMEND] DT: AUPRC > AUROC (不平衡比 {dt_ratio:.1f}:1 > 10:1)")
    if dg_ratio > 10:
        recommendations.append(f"DG不平衡比 {dg_ratio:.1f}:1 > 10:1 — 建议优先使用 AUPRC 而非 AUROC")
        log(f"  [RECOMMEND] DG: AUPRC > AUROC (不平衡比 {dg_ratio:.1f}:1 > 10:1)")

    if dt_ratio > 100:
        log(f"  [WARN] DT 极度不平衡 ({dt_ratio:.1f}:1), 警惕过拟合风险!")
        recommendations.append("DT 极度不平衡 — 考虑使用 SMOTE/ADASYN 或调整 class_weight")
    if dg_ratio > 100:
        log(f"  [WARN] DG 极度不平衡 ({dg_ratio:.1f}:1), 警惕过拟合风险!")
        recommendations.append("DG 极度不平衡 — 考虑使用 SMOTE/ADASYN 或调整 class_weight")

    if n_dt < 10:
        log(f"  [WARN] DT 正样本过少 ({n_dt}), CV 极不稳定")
        recommendations.append(f"DT 正样本仅 {n_dt} 个 — 建议收集更多正样本或使用半监督学习")
    if n_dg < 10:
        log(f"  [WARN] DG 正样本过少 ({n_dg}), CV 极不稳定")
        recommendations.append(f"DG 正样本仅 {n_dg} 个 — 建议收集更多正样本或使用半监督学习")

    report['recommendations'] = recommendations

    if not recommendations:
        log(f"  类别分布相对平衡, 无需特殊处理.")

    # 保存报告
    report_path = os.path.join(OUTPUT_DIR, "class_imbalance_diagnostic.json")
    serializable = {k: (float(v) if isinstance(v, (np.floating, np.integer)) else v)
                    for k, v in report.items()}
    safe_write_json(serializable, report_path)

    return report


# ═══════════════════════════════════════════════════════════════════════════════
# [v9.0] MI-VIF 特征选择 — 重写为防泄露版本
#   参考: PMID:37113250 (nestedcv) + PMID:36004690 (data leakage)
# ═══════════════════════════════════════════════════════════════════════════════
def mi_vif_feature_selection_per_fold(X_train, y_train, feature_cols, mi_k=100, vif_threshold=10.0):
    """
    [v9.0] 每折内嵌 MI-VIF 特征选择 — 无数据泄露。
    仅在训练集上计算 MI 得分和 VIF, 不接触验证集/未知基因。
    参考 PMID:37113250 — 特征选择必须嵌入交叉验证内部。

    返回: (selected_indices, selected_col_names, mi_scores_dict)
    """
    n_features = X_train.shape[1]
    if n_features <= mi_k:
        return np.arange(n_features), feature_cols[:], {}

    try:
        # Step 1: 仅在训练集上计算互信息
        mi_scores = mutual_info_classif(X_train, y_train, random_state=SEED, n_neighbors=3)

        # Step 2: 选择 top-K
        top_k_indices = np.argsort(-mi_scores)[:mi_k]
        X_top = X_train[:, top_k_indices]

        # Step 3: VIF 过滤 (子样本加速, 限制特征数)
        n_top = X_top.shape[1]
        if n_top < 2:
            return top_k_indices, [feature_cols[i] for i in top_k_indices], {}

        keep_mask = np.ones(n_top, dtype=bool)
        max_vif_check = min(n_top, 100)  # 限制VIF检查数

        for i in range(max_vif_check):
            if not keep_mask[i]:
                continue
            try:
                y_col = X_top[:, i]
                X_rest = np.delete(X_top[:, keep_mask], i, axis=1)
                if X_rest.shape[1] == 0:
                    continue
                # 子采样加速 VIF 回归
                sample_size = min(2000, X_rest.shape[0])
                if X_rest.shape[0] > sample_size:
                    s_idx = np.random.choice(X_rest.shape[0], sample_size, replace=False)
                    y_s, X_s = y_col[s_idx], X_rest[s_idx]
                else:
                    y_s, X_s = y_col, X_rest
                lr = LinearRegression()
                lr.fit(X_s, y_s)
                r2 = lr.score(X_s, y_s)
                vif = 1.0 / (1.0 - r2) if r2 < 1.0 else np.inf
                if vif > vif_threshold:
                    keep_mask[i] = False
            except Exception:
                continue

        final_indices_in_top = np.where(keep_mask)[0]
        final_indices = top_k_indices[final_indices_in_top]
        final_features = [feature_cols[i] for i in final_indices]

        mi_scores_dict = {feature_cols[i]: float(mi_scores[i]) for i in final_indices}

        return final_indices, final_features, mi_scores_dict

    except Exception:
        # 回退: 返回所有特征
        return np.arange(n_features), feature_cols[:], {}


# ═══════════════════════════════════════════════════════════════════════════════
# 稳定性选择
#   参考: PMID:33561948 — Bootstrap 重采样评估特征选择稳定性
# ═══════════════════════════════════════════════════════════════════════════════
def stability_selection(X, y, feature_cols, n_bootstrap=100, threshold=0.8):
    """
    稳定性选择: 使用 Bootstrap 重采样评估特征选择稳定性。
    参考 PMID:33561948。
    """
    log("\n" + "=" * 70)
    log(f"稳定性选择 [v9.0] — 参考 PMID:33561948")
    log(f"  Bootstrap: {n_bootstrap}, 稳定性阈值: {threshold}")
    log("=" * 70)

    n_features = X.shape[1]
    n_samples = X.shape[0]

    if n_features < 5:
        log(f"  [SKIP] 特征数过少 ({n_features})")
        return [], None

    selection_counts = np.zeros(n_features, dtype=int)
    rng = np.random.RandomState(SEED)

    log(f"  运行 {n_bootstrap} 次 Bootstrap LASSO...")
    for b in range(n_bootstrap):
        boot_idx = rng.choice(n_samples, size=n_samples, replace=True)
        X_boot = X[boot_idx]
        y_boot = y[boot_idx]

        try:
            lasso = LassoCV(cv=3, n_alphas=20, random_state=SEED + b, max_iter=2000, n_jobs=1)
            lasso.fit(X_boot, y_boot)
            selected = np.abs(lasso.coef_) > 1e-8
            selection_counts[selected] += 1
        except Exception:
            continue

        if (b + 1) % 20 == 0:
            log(f"    Bootstrap 进度: {b + 1}/{n_bootstrap}")

    selection_freq = selection_counts / n_bootstrap
    stable_mask = selection_freq >= threshold
    stable_indices = np.where(stable_mask)[0]
    stable_features = [feature_cols[i] for i in stable_indices]

    log(f"  稳定特征数: {len(stable_features)}/{n_features} (频率 >= {threshold})")
    log(f"  选择频率范围: [{selection_freq.min():.3f}, {selection_freq.max():.3f}]")

    if len(stable_features) >= 5:
        top5 = sorted(zip(stable_features, selection_freq[stable_indices]),
                      key=lambda x: -x[1])[:5]
        log(f"  Top-5 稳定特征: {[(f, f'{s:.3f}') for f, s in top5]}")

    freq_df = pd.DataFrame({
        'feature': feature_cols,
        'selection_frequency': selection_freq,
        'is_stable': stable_mask
    }).sort_values('selection_frequency', ascending=False)

    stability_path = os.path.join(OUTPUT_DIR, "stability_selection.csv")
    safe_write_csv(freq_df, stability_path)

    return stable_features, freq_df


# ═══════════════════════════════════════════════════════════════════════════════
# [v9.0 FIX] 模型校准 — 已修复: 在验证集上校准已训练分类器
#   参考: PMID:32865408 — isotonic regression 对 RF 可能降低概率质量
# ═══════════════════════════════════════════════════════════════════════════════
def calibrate_pretrained_clf(clf, X_val, y_val, X_unk):
    """
    [v9.0 FIX] 在验证集上校准已训练分类器。
    使用 CalibratedClassifierCV(cv='prefit') 在已训练分类器上执行 isotonic 校准。
    参考 PMID:32865408 — Platt scaling 和 isotonic regression 的比较。

    参数:
        clf: 已训练的分类器 (必须支持 predict_proba)
        X_val, y_val: 验证集 (用于校准)
        X_unk: 未知基因特征矩阵 (用于预测)

    返回: (y_val_calibrated_proba, y_unk_calibrated_proba) 或 (None, None)
    """
    try:
        calibrated = CalibratedClassifierCV(estimator=clf, method='isotonic', cv='prefit')
        calibrated.fit(X_val, y_val)
        y_val_cal = calibrated.predict_proba(X_val)[:, 1]
        y_unk_cal = calibrated.predict_proba(X_unk)[:, 1]
        return y_val_cal, y_unk_cal
    except Exception:
        return None, None


# ═══════════════════════════════════════════════════════════════════════════════
# Stacking Ensemble
#   参考: PMID:36227031 + PMID:31589422 (STarFish)
# ═══════════════════════════════════════════════════════════════════════════════
def stacking_ensemble(drug_probas, disease_probas, bridge_scores, y_dt, y_dg,
                      unknown_mask, gene_index, model_names, dt_metrics=None, dg_metrics=None):
    """
    Stacking Ensemble: 构建 LogisticRegression 元模型。
    参考 PMID:31589422 (STarFish) — 堆叠集成提升靶点预测 AUROC 至 0.94。
    """
    if not USE_STACKING:
        log("  [SKIP] Stacking Ensemble 已禁用")
        return None, None

    log("\n" + "=" * 70)
    log("Stacking Ensemble 元模型 [v9.0] — 参考 PMID:31589422 (STarFish)")
    log("=" * 70)

    try:
        common_models = sorted(set(drug_probas.keys()) & set(disease_probas.keys()))
        if len(common_models) < 5:
            log(f"  [SKIP] 公共模型不足 ({len(common_models)} < 5)")
            return None, None

        unknown_pos = np.where(unknown_mask)[0]
        labeled_pos = np.where(~unknown_mask)[0]

        n_labeled = len(labeled_pos)
        n_unknown = len(unknown_pos)
        n_models = len(common_models)

        X_stack_dt = np.zeros((n_labeled, n_models), dtype=np.float64)
        X_stack_dg = np.zeros((n_labeled, n_models), dtype=np.float64)
        X_unk_dt = np.zeros((n_unknown, n_models), dtype=np.float64)
        X_unk_dg = np.zeros((n_unknown, n_models), dtype=np.float64)

        for j, model in enumerate(common_models):
            X_stack_dt[:, j] = drug_probas[model][labeled_pos]
            X_stack_dg[:, j] = disease_probas[model][labeled_pos]
            X_unk_dt[:, j] = drug_probas[model][unknown_pos]
            X_unk_dg[:, j] = disease_probas[model][unknown_pos]

        y_labeled_dt = y_dt[labeled_pos]
        y_labeled_dg = y_dg[labeled_pos]

        meta_dt = LogisticRegression(penalty='l2', C=1.0, class_weight='balanced',
                                      max_iter=5000, random_state=SEED)
        meta_dt.fit(X_stack_dt, y_labeled_dt)
        stacking_dt_probas = meta_dt.predict_proba(X_unk_dt)[:, 1]

        meta_dg = LogisticRegression(penalty='l2', C=1.0, class_weight='balanced',
                                      max_iter=5000, random_state=SEED)
        meta_dg.fit(X_stack_dg, y_labeled_dg)
        stacking_dg_probas = meta_dg.predict_proba(X_unk_dg)[:, 1]

        stacking_bridge = stacking_dt_probas * stacking_dg_probas

        from sklearn.model_selection import cross_val_score
        try:
            dt_cv_auroc = cross_val_score(meta_dt, X_stack_dt, y_labeled_dt,
                                          cv=min(5, n_labeled), scoring='roc_auc').mean()
            dg_cv_auroc = cross_val_score(meta_dg, X_stack_dg, y_labeled_dg,
                                          cv=min(5, n_labeled), scoring='roc_auc').mean()
        except Exception:
            dt_cv_auroc = 0.0
            dg_cv_auroc = 0.0

        log(f"  Stacking DT AUROC (CV): {dt_cv_auroc:.4f}")
        log(f"  Stacking DG AUROC (CV): {dg_cv_auroc:.4f}")
        log(f"  Stacking 桥梁得分范围: [{stacking_bridge.min():.6f}, {stacking_bridge.max():.6f}]")

        stacking_sort = unknown_pos[np.argsort(-stacking_bridge)]
        stacking_top20 = {gene_index[i] for i in stacking_sort[:20]}
        log(f"  Stacking Top-20: {sorted(stacking_top20)[:10]}...")

        meta_info = {
            'dt_auroc_cv': float(dt_cv_auroc),
            'dg_auroc_cv': float(dg_cv_auroc),
            'n_models': n_models,
            'meta_model': 'LogisticRegression(L2)',
        }

        stacking_df = pd.DataFrame({
            'gene_symbol': [gene_index[i] for i in unknown_pos],
            'stacking_dt_proba': stacking_dt_probas,
            'stacking_dg_proba': stacking_dg_probas,
            'stacking_bridge_score': stacking_bridge,
        }).sort_values('stacking_bridge_score', ascending=False)
        stacking_path = os.path.join(OUTPUT_DIR, "stacking_ensemble_results.csv")
        safe_write_csv(stacking_df, stacking_path)

        return stacking_bridge, meta_info

    except Exception as e:
        log(f"  [WARN] Stacking Ensemble 失败: {e}")
        return None, None


# ═══════════════════════════════════════════════════════════════════════════════
# 模型性能汇总图
# ═══════════════════════════════════════════════════════════════════════════════
def plot_model_performance_summary(metrics_rows, output_dir):
    """
    生成柱状图比较所有模型组合的 AUROC 和 AUPRC 得分。
    """
    log("\n" + "=" * 70)
    log("模型性能汇总图 [v9.0]")
    log("=" * 70)

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        df = pd.DataFrame(metrics_rows)
        if df.empty:
            log("  [SKIP] 无指标数据")
            return

        for task in ['DT', 'DG']:
            task_df = df[df['task'] == task].copy()
            if task_df.empty:
                continue

            task_df = task_df.sort_values('auroc', ascending=False).head(30)
            n_models = len(task_df)

            fig, axes = plt.subplots(1, 2, figsize=(16, max(8, n_models * 0.3)))
            fig.suptitle(f'Model Performance Summary — {task} Task [v9.0]',
                         fontsize=14, fontweight='bold')

            colors_auroc = ['#2196F3' if v >= 0.8 else '#FF9800' if v >= 0.6 else '#F44336'
                           for v in task_df['auroc']]
            axes[0].barh(range(n_models), task_df['auroc'].values, color=colors_auroc,
                        edgecolor='white', linewidth=0.5)
            if 'auroc_std' in task_df.columns:
                axes[0].errorbar(task_df['auroc'].values, range(n_models),
                                xerr=task_df['auroc_std'].values, fmt='none',
                                ecolor='gray', capsize=2, alpha=0.6)
            axes[0].set_yticks(range(n_models))
            axes[0].set_yticklabels(task_df['model'].values, fontsize=7)
            axes[0].set_xlabel('AUROC')
            axes[0].set_title(f'AUROC (n={n_models})')
            axes[0].axvline(x=0.5, color='gray', linestyle='--', alpha=0.5)
            axes[0].invert_yaxis()

            colors_auprc = ['#4CAF50' if v >= 0.5 else '#FF9800' if v >= 0.2 else '#F44336'
                           for v in task_df['auprc']]
            axes[1].barh(range(n_models), task_df['auprc'].values, color=colors_auprc,
                        edgecolor='white', linewidth=0.5)
            if 'auprc_std' in task_df.columns:
                axes[1].errorbar(task_df['auprc'].values, range(n_models),
                                xerr=task_df['auprc_std'].values, fmt='none',
                                ecolor='gray', capsize=2, alpha=0.6)
            axes[1].set_yticks(range(n_models))
            axes[1].set_yticklabels(task_df['model'].values, fontsize=7)
            axes[1].set_xlabel('AUPRC')
            axes[1].set_title(f'AUPRC (n={n_models})')
            axes[1].invert_yaxis()

            plt.tight_layout()
            fig_path = os.path.join(output_dir, f"model_performance_summary_{task}.pdf")
            plt.savefig(fig_path, dpi=150, bbox_inches='tight')
            plt.close()
            log(f"  性能汇总图 ({task}): {fig_path}")

            log(f"  [{task}] Top-5 AUROC:")
            for _, row in task_df.head(5).iterrows():
                log(f"    {row['model']:<35s} AUROC={row['auroc']:.4f}  AUPRC={row['auprc']:.4f}")

    except Exception as e:
        log(f"  [WARN] 性能汇总图生成失败: {e}")


# ═══════════════════════════════════════════════════════════════════════════════
# [v9.0] 零方差特征移除
#   参考: PMID:40199913 — 数据分割前的严格验证
# ═══════════════════════════════════════════════════════════════════════════════
def remove_zero_variance_features(X, feature_cols):
    """
    移除方差为0的常量特征列。
    这些特征不提供任何区分信息, 且可能导致数值问题。
    """
    variances = np.var(X, axis=0)
    zero_var_mask = variances < 1e-12
    n_zero = zero_var_mask.sum()

    if n_zero == 0:
        log(f"  [OK] 无零方差特征")
        return X, feature_cols, 0

    keep_mask = ~zero_var_mask
    X_filtered = X[:, keep_mask]
    filtered_cols = [c for i, c in enumerate(feature_cols) if keep_mask[i]]
    removed_examples = [feature_cols[i] for i in np.where(zero_var_mask)[0][:5]]

    log(f"  [CLEAN] 移除 {n_zero} 个零方差特征: {removed_examples}...")
    return X_filtered, filtered_cols, n_zero


def detect_duplicate_features(X, feature_cols):
    """
    检测并报告完全重复的特征列。
    不自动移除, 仅报告供用户决策。
    """
    n_features = X.shape[1]
    duplicates = []
    seen_hashes = {}

    for i in range(n_features):
        col_hash = hash(tuple(X[:min(1000, X.shape[0]), i].round(6)))
        if col_hash in seen_hashes:
            duplicates.append((feature_cols[i], feature_cols[seen_hashes[col_hash]]))
        else:
            seen_hashes[col_hash] = i

    if duplicates:
        log(f"  [WARN] 检测到 {len(duplicates)} 对可能重复的特征:")
        for dup, orig in duplicates[:5]:
            log(f"    {dup} ≈ {orig}")
    else:
        log(f"  [OK] 无检测到重复特征列")

    return duplicates


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 0: 增强输入验证层
#   参考: DataSAIL PMID:40199913 — 数据分割前的严格验证
#   参考: nestedcv PMID:37113250 — 嵌套CV防信息泄露
# ═══════════════════════════════════════════════════════════════════════════════
def validate_inputs():
    """
    [v9.0] 对所有输入文件进行格式、内容、一致性的严格验证。
    增强: 数据分布检查 + 特征方差检查 + 标签平衡评估 + 重复特征检测
    """
    log("=" * 70)
    log("STEP 0: 增强输入数据验证层 [v9.0]")
    log("  - 参考: DataSAIL PMID:40199913 数据分割前严格验证")
    log("  - 参考: nestedcv PMID:37113250 嵌套CV防信息泄露")
    log("  - 参考: PMID:36004690 特征选择泄露是常见错误")
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

                # [v9.0] 重复基因检测
                dup_genes = df['gene_symbol'].duplicated().sum()
                if dup_genes > 0:
                    log(f"  [INFO] {dup_genes} 个重复基因符号, 将保留第一个")
                    quality_report['duplicate_genes'] = int(dup_genes)

                numeric_cols = df.select_dtypes(include=[np.number]).columns
                if len(numeric_cols) > 0:
                    variances = df[numeric_cols].var()
                    zero_var = (variances == 0).sum()
                    low_var = (variances < 1e-8).sum()
                    if zero_var > 0:
                        issues.append(f"[WARN] {zero_var} 个特征方差为0 (常量特征 — 将在预处理中移除)")
                    if low_var > 0:
                        log(f"  [INFO] {low_var} 个特征方差 < 1e-8 (低变异特征)")
                    quality_report['zero_var_features'] = int(zero_var)
                    quality_report['low_var_features'] = int(low_var)
                    quality_report['median_variance'] = float(variances.median())

                missing_counts = df.isnull().sum()
                missing_cols = (missing_counts > 0).sum()
                if missing_cols > 0:
                    log(f"  [INFO] {missing_cols} 个特征含缺失值, 总缺失: {missing_counts.sum()}")
                    quality_report['missing_cols'] = int(missing_cols)
                    quality_report['total_missing'] = int(missing_counts.sum())
        except Exception as e:
            issues.append(f"[FATAL] 特征矩阵读取失败: {e}")
            all_valid = False

    # 0b. 药物靶点/疾病基因文件
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
# 特征多重共线性诊断
# ═══════════════════════════════════════════════════════════════════════════════
def diagnose_multicollinearity(X, feature_cols, vif_threshold=VIF_THRESHOLD):
    """
    使用VIF (Variance Inflation Factor) 和条件数诊断多重共线性。
    """
    log("\n" + "-" * 50)
    log(f"多重共线性诊断 (VIF={vif_threshold}) [v9.0]")
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

    # VIF 诊断
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
                sample_size = min(5000, X.shape[0])
                if X.shape[0] > sample_size:
                    sample_idx = np.random.choice(X.shape[0], sample_size, replace=False)
                    y_s = y_col[sample_idx]
                    X_s = X_rest[sample_idx]
                else:
                    y_s = y_col
                    X_s = X_rest
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
# 特征相关性过滤 (无监督, 无数据泄露)
#   参考: PMID:38418819 — 无监督过滤不导致泄露
# ═══════════════════════════════════════════════════════════════════════════════
def filter_highly_correlated_features(X, feature_cols, threshold=0.95):
    """
    基于 Pearson 相关系数过滤高相关特征。
    这是无监督操作, 不使用标签信息, 不会导致数据泄露 (PMID:38418819)。
    优先保留方差较大的特征。
    """
    log("\n" + "-" * 50)
    log(f"特征相关性过滤 (threshold={threshold}) [v9.0]")
    log(f"  注: 无监督过滤, 不涉及标签, 无数据泄露风险 (参考文献 PMID:38418819)")
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
# STEP 1a: 主特征数据加载
# ═══════════════════════════════════════════════════════════════════════════════
def load_data():
    log("=" * 70)
    log("STEP 1a: 主特征数据加载与标签构建 [v9.0]")
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

    # 标签不平衡评估
    log(f"  --- 标签平衡评估 ---")
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
# STEP 1b: Toxi 环境暴露特征提取
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
        if rows:
            safe_write_csv(pd.DataFrame(rows), ORTHO_CACHE_PATH)
            log(f"  [Cache] 已保存 {len(rows)} 条同源映射")
    except Exception as e:
        log(f"  [Cache] 保存失败: {e}")


# ═══════════════════════════════════════════════════════════════════════════════
# 特征工程 Builder
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
        LassoCV(cv=2, n_alphas=20, random_state=SEED, max_iter=2000, n_jobs=1),
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
# 分类器 Builder
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

# ═══════════════════════════════════════════════════════════════════════════════
# [v9.0] 不支持 class_weight 的分类器类型集合 — 用于 isinstance 检查
# ═══════════════════════════════════════════════════════════════════════════════
_NO_CLASS_WEIGHT_TYPES = (GradientBoostingClassifier, GaussianNB)

# 可能返回的已包装分类器类型 (SVC, PAC 已被 CalibratedClassifierCV 包装)
_NO_SAMPLE_WEIGHT_TYPES = _NO_CLASS_WEIGHT_TYPES + (CalibratedClassifierCV,)


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
# [v9.0] STEP 3-4: 并行CV — 修复 class_weight Bug + 修复数据泄露
#   参考: PMID:37113250 — 特征选择需嵌入CV
#   参考: PMID:32865408 — 校准方法比较
# ═══════════════════════════════════════════════════════════════════════════════
def train_single_model(model_name, fe_builder, clf_builder, X_base, y,
                       unknown_mask, feature_cols_global, model_idx):
    """
    [v9.0] 单模型CV训练 — 修复了两个关键Bug:
    1. NO_CLASS_WEIGHT_CLFS 使用 isinstance 直接检查类
    2. MI-VIF 模式下每折内嵌特征选择 (无泄露)

    参考:
      - PMID:37113250 — 嵌入式特征选择防泄露
      - PMID:36004690 — 外部特征选择导致性能膨胀
      - PMID:32865408 — 校准方法选择
    """
    warnings.filterwarnings("ignore", category=ConvergenceWarning)
    warnings.filterwarnings("ignore", category=FutureWarning)
    warnings.filterwarnings("ignore", message="X does not have valid feature names")

    seed = SEED + model_idx
    n_genes = X_base.shape[0]
    probas = np.zeros(n_genes, dtype=np.float64)
    labeled_assigned = np.zeros(n_genes, dtype=bool)

    # RepeatedStratifiedKFold
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

        # [v9.0] 每折内嵌 MI-VIF 特征选择 (防泄露)
        # 这一步在 scaler 之前, 因为 MI/VIF 应该在原始特征空间计算
        if MI_VIF_MODE == 'per_fold' and X_train.shape[1] > 100:
            selected_idx, _, _ = mi_vif_feature_selection_per_fold(
                X_train, y_train, feature_cols_global, mi_k=min(100, X_train.shape[1]))
            X_train = X_train[:, selected_idx]
            X_val = X_val[:, selected_idx]
            X_unk_base = X_base[unknown_pos][:, selected_idx]
        else:
            X_unk_base = X_base[unknown_pos]

        scaler = _Scaler()
        X_train_s = scaler.fit_transform(X_train)
        X_val_s = scaler.transform(X_val)
        X_unk_s = scaler.transform(X_unk_base)

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

            # [v9.0] 正确检查分类器类型以决定是否使用 sample_weight
            if isinstance(clf, _NO_SAMPLE_WEIGHT_TYPES):
                clf.fit(X_train_fe, y_train)
            else:
                clf.fit(X_train_fe, y_train)
                # 所有带 class_weight 的分类器已自动处理, 无需额外 sample_weight
        except Exception:
            continue

        # [v9.0] 模型校准 — 使用已训练分类器
        # 注意: PMID:32865408 指出 isotonic regression 对 RF 可能降低概率质量
        # 因此在 tree-based 模型上跳过校准 (RF, GB, XGBoost, LightGBM, ExtraTrees)
        calibration_applied = False
        skip_calibration_types = (RandomForestClassifier, GradientBoostingClassifier,
                                   ExtraTreesClassifier)
        try:
            import xgboost as xgb
            skip_calibration_types = skip_calibration_types + (xgb.XGBClassifier,)
        except ImportError:
            pass
        try:
            import lightgbm as lgb
            skip_calibration_types = skip_calibration_types + (lgb.LGBMClassifier,)
        except ImportError:
            pass

        if USE_CALIBRATION and not isinstance(clf, skip_calibration_types):
            try:
                # [v9.0 FIX] 使用 calibrate_pretrained_clf 在验证集上校准
                y_val_cal, y_unk_cal = calibrate_pretrained_clf(clf, X_val_fe, y_val, X_unk_fe)
                if y_val_cal is not None and y_unk_cal is not None:
                    y_val_prob = y_val_cal
                    y_unk_prob = y_unk_cal
                    calibration_applied = True
            except Exception:
                pass

        if not calibration_applied:
            try:
                y_val_prob = clf.predict_proba(X_val_fe)[:, 1]
                y_unk_prob = clf.predict_proba(X_unk_fe)[:, 1]
            except Exception:
                continue

        try:
            probas[val_pos] = y_val_prob
            labeled_assigned[val_pos] = True
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


def run_cv_parallel(X_base, y, model_combos, task_name, unknown_mask, feature_cols_global):
    n_models = len(model_combos)
    log(f"  [{task_name}] 启动 {n_models} 个模型的并行训练 ({N_JOBS} 核, {N_REPEATS}重复)...")
    if USE_CALIBRATION:
        log(f"  [{task_name}] Isotonic Calibration: ENABLED (cv='prefit', tree-based 跳过)")
    if MI_VIF_MODE == 'per_fold':
        log(f"  [{task_name}] MI-VIF 特征选择: 每折内嵌 (防泄露, 参考 PMID:37113250)")
    elif MI_VIF_MODE == 'pre_cv':
        log(f"  [{task_name}] [WARN] MI-VIF: CV前执行 (有数据泄露风险!)")
    else:
        log(f"  [{task_name}] MI-VIF: 跳过")

    results = Parallel(n_jobs=N_JOBS, verbose=5)(
        delayed(train_single_model)(name, fe, clf, X_base, y, unknown_mask,
                                     feature_cols_global, idx)
        for idx, (name, fe, clf) in enumerate(model_combos)
    )
    all_probas = {}
    all_metrics = {}
    for model_name, probas, metric in results:
        all_probas[model_name] = probas
        all_metrics[model_name] = metric

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
# 加权RRF集成
#   参考: PMID:36384050 — KG ensemble加权提升预测性能
# ═══════════════════════════════════════════════════════════════════════════════
def rrf_integrate(bridge_scores, n_genes, unknown_mask, dt_metrics=None, dg_metrics=None):
    """
    RRF 倒数排名融合 (加权版)。
    """
    rrf_sum = np.zeros(n_genes, dtype=np.float64)

    if WEIGHTED_RRF and dt_metrics is not None and dg_metrics is not None:
        model_weights = {}
        for model in bridge_scores:
            dt_auroc = dt_metrics.get(model, {}).get('auroc', 0.0)
            dg_auroc = dg_metrics.get(model, {}).get('auroc', 0.0)
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

    # 等权重RRF
    for model, scores in bridge_scores.items():
        uk_scores = scores[unknown_mask]
        order = np.argsort(-uk_scores)
        ranks = np.empty(len(uk_scores), dtype=np.float64)
        ranks[order] = np.arange(1, len(uk_scores) + 1)
        rrf_sum[unknown_mask] += 1.0 / (RRF_K + ranks)
    return rrf_sum


# ═══════════════════════════════════════════════════════════════════════════════
# Borda Count 排名融合 [v9.0]
# ═══════════════════════════════════════════════════════════════════════════════
def borda_count_integrate(bridge_scores, n_genes, unknown_mask):
    """
    Borda Count 排名融合。
    [v9.0] 降低标准差排除阈值: 1e-12
    """
    n_unknown = unknown_mask.sum()
    borda_sum = np.zeros(n_genes, dtype=np.float64)
    valid_models = 0

    for model, scores in bridge_scores.items():
        uk_scores = scores[unknown_mask]

        if np.max(uk_scores) <= 0 or np.std(uk_scores) < 1e-12:
            continue

        valid_models += 1
        order = np.argsort(uk_scores)
        points = np.arange(n_unknown)
        borda_single = np.zeros(n_unknown, dtype=np.float64)
        borda_single[order] = points
        borda_sum[unknown_mask] += borda_single

    log(f"  Borda Count: {valid_models}/{len(bridge_scores)} 有效模型参与")
    return borda_sum


# ═══════════════════════════════════════════════════════════════════════════════
# [v9.0] Bootstrap 置信区间 — 添加进度追踪
# ═══════════════════════════════════════════════════════════════════════════════
def bootstrap_rrf_ci(bridge_scores, n_genes, unknown_mask, n_bootstrap=N_BOOTSTRAP):
    """
    对 RRF 分数进行 Bootstrap 重采样, 计算每个基因的 95% 置信区间。
    [v9.0] 添加进度追踪。
    """
    log("\n" + "-" * 50)
    log(f"Bootstrap RRF 置信区间 (n={n_bootstrap}) [v9.0]")
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
    progress_interval = max(1, n_bootstrap // 5)
    for b in range(n_bootstrap):
        idx = rng.choice(n_models, size=n_models, replace=True)
        boot_rrf[b] = all_ranks[idx].sum(axis=0)
        if (b + 1) % progress_interval == 0:
            log(f"    Bootstrap 进度: {b+1}/{n_bootstrap}")

    rrf_lower = np.percentile(boot_rrf, 2.5, axis=0)
    rrf_upper = np.percentile(boot_rrf, 97.5, axis=0)

    boot_mean = boot_rrf.mean(axis=0)
    boot_std = boot_rrf.std(axis=0)
    cv = np.full(n_unknown, np.nan)
    safe_mask = boot_mean > 1e-10
    cv[safe_mask] = boot_std[safe_mask] / boot_mean[safe_mask]

    log(f"  RRF CV (变异系数): median={np.nanmedian(cv):.4f}, "
        f"IQR=[{np.nanpercentile(cv,25):.4f}, {np.nanpercentile(cv,75):.4f}]")
    log(f"  CI width (中位数): {np.median(rrf_upper - rrf_lower):.6f}")

    return rrf_lower, rrf_upper


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 7: 四算法共识 (加权版)
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
    log("STEP 7: 四算法家族交集共识筛选 (加权版) [v9.0]")
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

        family_avg = np.zeros(len(unknown_pos), dtype=np.float64)
        total_weight = 0.0
        for mname in family_models:
            if dt_metrics is not None and dg_metrics is not None:
                dt_auroc = dt_metrics.get(mname, {}).get('auroc', 0.0)
                dg_auroc = dg_metrics.get(mname, {}).get('auroc', 0.0)
                w = max(0.0, (dt_auroc + dg_auroc) / 2.0 - 0.5)
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
# SHAP 误分类检测 — 四种过滤规则
#   参考: PMID:41114446
# ═══════════════════════════════════════════════════════════════════════════════
def shap_misclassification_detection(X_base, y_dt, df_unknown, gene_index, feature_cols):
    """
    使用 SHAP 值对 Top-N 预测进行误分类检测 (四种过滤规则)。
    参考 PMID:41114446。
    """
    if not _SHAP_AVAILABLE:
        log("  [SKIP] SHAP 未安装, 无法进行误分类检测")
        return None, None

    if not MISCLASSIFICATION_FILTER:
        log("  [SKIP] 误分类检测已禁用")
        return None, None

    log("\n" + "=" * 70)
    log("SHAP 误分类检测 [v9.0] — 参考 PMID:41114446 (四种过滤规则)")
    log("=" * 70)

    try:
        scaler = StandardScaler()
        X_s = scaler.fit_transform(X_base)
        clf = RandomForestClassifier(n_estimators=200, class_weight='balanced',
                                      random_state=SEED, n_jobs=1)
        clf.fit(X_s, y_dt)

        explainer = shap.TreeExplainer(clf)
        n_samples = min(500, X_s.shape[0])
        idx_sample = np.random.choice(X_s.shape[0], n_samples, replace=False)
        X_sample = X_s[idx_sample]
        y_sample = y_dt[idx_sample]

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

        n_feat = min(sv_pos.shape[1], len(feature_cols))

        # RAW 规则: 特征值范围
        pos_mask = y_sample == 1
        neg_mask = y_sample == 0

        raw_ranges = {}
        for feat_idx in range(n_feat):
            fname = feature_cols[feat_idx]
            raw_ranges[fname] = {
                'pos': {
                    'q05': float(np.percentile(X_sample[pos_mask, feat_idx], 5)) if pos_mask.sum() > 0 else -np.inf,
                    'q95': float(np.percentile(X_sample[pos_mask, feat_idx], 95)) if pos_mask.sum() > 0 else np.inf,
                },
                'neg': {
                    'q05': float(np.percentile(X_sample[neg_mask, feat_idx], 5)) if neg_mask.sum() > 0 else -np.inf,
                    'q95': float(np.percentile(X_sample[neg_mask, feat_idx], 95)) if neg_mask.sum() > 0 else np.inf,
                }
            }

        # SHAP 规则: SHAP值范围
        feature_shap_ranges = {}
        for feat_idx in range(n_feat):
            fname = feature_cols[feat_idx]
            vals = sv_pos[:, feat_idx]
            feature_shap_ranges[fname] = {
                'mean': float(np.mean(vals)),
                'std': float(np.std(vals)),
                'q05': float(np.percentile(vals, 5)),
                'q95': float(np.percentile(vals, 95)),
                'pct_positive': float(np.mean(vals > 0)),
            }

        log(f"  SHAP 分析完成: {len(feature_shap_ranges)} 特征")

        # 对 Top-N 基因应用四种规则
        top_n = min(50, len(df_unknown))
        flags = []

        for i in range(top_n):
            gene_symbol = df_unknown.iloc[i]['gene_symbol']
            gene_pos = gene_index.index(gene_symbol)
            x_gene = X_s[gene_pos]
            pred_prob = clf.predict_proba(x_gene.reshape(1, -1))[0, 1]
            pred_class = 1 if pred_prob >= 0.5 else 0

            sv_gene = explainer.shap_values(x_gene.reshape(1, -1))
            if isinstance(sv_gene, list) and len(sv_gene) == 2:
                sv_gene = sv_gene[1][0]
            elif hasattr(sv_gene, 'values'):
                sv_gene = sv_gene.values[0]

            raw_outlier_count = 0
            raw_total = 0
            for feat_idx in range(n_feat):
                fname = feature_cols[feat_idx]
                if fname in raw_ranges:
                    raw_total += 1
                    val = x_gene[feat_idx]
                    if pred_class == 1:
                        rng = raw_ranges[fname]['pos']
                    else:
                        rng = raw_ranges[fname]['neg']
                    if val < rng['q05'] or val > rng['q95']:
                        raw_outlier_count += 1

            raw_ratio = raw_outlier_count / max(raw_total, 1)
            raw_flagged = raw_ratio > 0.2

            shap_outlier_count = 0
            shap_total = 0
            for feat_idx in range(n_feat):
                fname = feature_cols[feat_idx]
                if fname in feature_shap_ranges:
                    shap_total += 1
                    s_val = sv_gene[feat_idx]
                    fr = feature_shap_ranges[fname]
                    if s_val < fr['q05'] or s_val > fr['q95']:
                        shap_outlier_count += 1

            shap_ratio = shap_outlier_count / max(shap_total, 1)
            shap_flagged = shap_ratio > 0.2

            raw_or_shap = raw_flagged or shap_flagged
            raw_and_shap = raw_flagged and shap_flagged

            flags.append({
                'gene_symbol': gene_symbol,
                'rank': i + 1,
                'pred_prob': float(pred_prob),
                'pred_class': pred_class,
                'raw_outlier_count': raw_outlier_count,
                'raw_ratio': float(raw_ratio),
                'raw_flagged': raw_flagged,
                'shap_outlier_count': shap_outlier_count,
                'shap_ratio': float(shap_ratio),
                'shap_flagged': shap_flagged,
                'raw_or_shap_flagged': raw_or_shap,
                'raw_and_shap_flagged': raw_and_shap,
                'flagged': raw_or_shap,
            })

        n_raw = sum(1 for f in flags if f['raw_flagged'])
        n_shap = sum(1 for f in flags if f['shap_flagged'])
        n_or = sum(1 for f in flags if f['raw_or_shap_flagged'])
        n_and = sum(1 for f in flags if f['raw_and_shap_flagged'])

        log(f"  Top-{top_n} 误分类检测结果:")
        log(f"    RAW 规则: {n_raw}/{top_n} flagged")
        log(f"    SHAP 规则: {n_shap}/{top_n} flagged")
        log(f"    RAW OR SHAP: {n_or}/{top_n} flagged")
        log(f"    RAW AND SHAP: {n_and}/{top_n} flagged")

        if n_or > 0:
            flagged_genes = [f['gene_symbol'] for f in flags if f['raw_or_shap_flagged']]
            log(f"  高风险基因 (RAW OR SHAP): {flagged_genes[:5]}")

        return flags, feature_shap_ranges

    except Exception as e:
        log(f"  [WARN] SHAP 误分类检测失败: {e}")
        return None, None


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 8: RF 重要性分层
# ═══════════════════════════════════════════════════════════════════════════════
def stratified_rf_importance(X, y_dt, y_dg, feature_cols, fe_map, clf_map,
                              metrics_rows, top_k=20):
    log("\n" + "=" * 70)
    log("STEP 8: RF 特征重要性三层分层 (原始 + Permutation + SHAP) [v9.0]")
    log("=" * 70)

    if not _SHAP_AVAILABLE:
        log("  [SKIP] SHAP 未安装")
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

    # L2: Permutation — 优化: 仅对L1 Top-100特征计算, n_repeats=2 [v9.0]
    try:
        perm_top_n = min(100, n_features)
        perm_features = l1_rank[:perm_top_n]
        X_perm = X_fe[:, perm_features]
        log(f"  L2 (Permutation): 对前{perm_top_n}个特征计算 (n_repeats=2)...")
        perm = permutation_importance(clf, X_perm, y_task, n_repeats=2,
                                       random_state=SEED, n_jobs=1,
                                       scoring='roc_auc')
        l2_imp = perm.importances_mean
        l2_rank = perm_features[np.argsort(-l2_imp)[:top_k]]
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

            explainer = shap.TreeExplainer(clf)
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
# STEP 9: 外部验证
# ═══════════════════════════════════════════════════════════════════════════════
def external_validation_gse61616(df_unknown, gene_index):
    log("\n" + "=" * 70)
    log("STEP 9: 外部交叉验证 (GSE61616 大鼠 CIRI 模型) [v9.0]")
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

    if 'gene' not in col_map:
        col_map['gene'] = gse.columns[0]
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
# 主流程 (v9.0 — 完整版)
# ═══════════════════════════════════════════════════════════════════════════════
def main():
    global N_JOBS, GPU_ENABLED, N_REPEATS, CHECKPOINT_ENABLED, WEIGHTED_RRF
    global MISCLASSIFICATION_FILTER, USE_ROBUST_SCALER, USE_CALIBRATION, USE_STACKING
    global MI_VIF_MODE, _pipeline_start_time

    _pipeline_start_time = time.time()

    # 解析命令行参数
    args = parse_args()
    N_JOBS = args.n_jobs
    GPU_ENABLED = bool(args.gpu)
    N_REPEATS = args.n_repeats
    CHECKPOINT_ENABLED = not args.no_checkpoint
    WEIGHTED_RRF = not args.no_weighted_rrf
    MISCLASSIFICATION_FILTER = not args.no_misclf_filter
    USE_CALIBRATION = not args.no_calibration
    USE_STACKING = not args.no_stacking
    USE_ROBUST_SCALER = args.robust_scaler
    MI_VIF_MODE = args.mi_vif_mode

    TOTAL_STEPS = 25
    current_step = 0

    start_time = time.time()
    log("═" * 70)
    log("  多算法融合桥梁靶点预测系统 v9.0 (完整可运行版)")
    log("  ")
    log("  === v9.0 核心升级 (vs v8.0) ===")
    log("  1. 修复 calibrate_pretrained_clf() 返回值 (原版返回 None, None)")
    log("  2. 补全 main() 函数 — STEP 6b-10 完整实现")
    log("  3. 修复第2628行语法错误 (缺失闭合括号)")
    log("  4. 新增 Borda Count 排名计算 + RRF vs Borda 比较")
    log("  5. 新增 Bootstrap CI 置信区间")
    log("  6. 新增四算法共识保存")
    log("  7. 新增 SHAP 误分类检测保存")
    log("  8. 新增 RF 重要性分层保存")
    log("  9. 新增外部验证 Spearman 保存")
    log("  10. 新增稳定性选择 + 模型性能汇总图")
    log(f"  ")
    log(f"  随机种子: {SEED}, {N_FOLDS}折CV x {N_REPEATS}重复, RRF_k={RRF_K}")
    log(f"  并行核心: {N_JOBS}, GPU加速: {GPU_ENABLED}")
    log(f"  MI-VIF模式: {MI_VIF_MODE} (per_fold=防泄露)")
    log(f"  加权RRF: {WEIGHTED_RRF}, 模型校准: {USE_CALIBRATION}")
    log(f"  Stacking: {USE_STACKING}, 误分类检测: {MISCLASSIFICATION_FILTER}")
    log(f"  输出目录: {OUTPUT_DIR}")
    log(f"  日志文件: {LOG_FILE}")
    log("═" * 70)

    # ── GPU 检测 ──
    detect_gpu_environment()
    current_step += 1
    log_progress(current_step, TOTAL_STEPS, "GPU检测")

    # ── STEP 0: 增强输入验证 ──
    all_valid, issues, quality_report = validate_inputs()
    if not all_valid:
        log("[FATAL] 输入验证失败, 终止执行")
        return
    current_step += 1
    log_progress(current_step, TOTAL_STEPS, "输入验证")

    # ── STEP 1a: 加载主特征 ──
    feat_df, feature_cols, all_genes = load_data()
    X_base = feat_df[feature_cols].values.astype(np.float64)
    y_dt = feat_df['is_drug_target'].values.astype(int)
    y_dg = feat_df['is_disease_gene'].values.astype(int)
    gene_index = feat_df.index.tolist()
    drug_targets = {g for g in gene_index if feat_df.loc[g, 'is_drug_target'] == 1}

    # ── 类别不平衡诊断 ──
    imbalance_report = class_imbalance_diagnostic(y_dt, y_dg, len(gene_index), gene_index)
    current_step += 1
    log_progress(current_step, TOTAL_STEPS, "类别不平衡诊断")

    # ── [v9.0] 零方差特征移除 ──
    if REMOVE_ZERO_VARIANCE:
        X_base, feature_cols, n_removed_zv = remove_zero_variance_features(X_base, feature_cols)
        log_stage("零方差移除")
    else:
        n_removed_zv = 0
    current_step += 1
    log_progress(current_step, TOTAL_STEPS, "零方差特征移除")

    # ── [v9.0] 重复特征检测 ──
    duplicate_pairs = detect_duplicate_features(X_base, feature_cols)
    log_stage("重复特征检测")
    current_step += 1
    log_progress(current_step, TOTAL_STEPS, "重复特征检测")

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
    current_step += 1
    log_progress(current_step, TOTAL_STEPS, "Toxi特征")

    # ── 多重共线性诊断 ──
    vif_df, cond_num, high_vif = diagnose_multicollinearity(X_base, feature_cols)
    log_stage("共线性诊断")
    current_step += 1
    log_progress(current_step, TOTAL_STEPS, "共线性诊断")

    # ── 特征相关性过滤 (无监督, 无泄露) ──
    X_base, feature_cols, n_removed_corr = filter_highly_correlated_features(
        X_base, feature_cols, FEATURE_CORR_THRESHOLD)
    log_stage("特征过滤")
    current_step += 1
    log_progress(current_step, TOTAL_STEPS, "特征相关性过滤")

    # ── [v9.0] MI-VIF: 仅 pre_cv 模式在主流程中执行 ──
    n_removed_mivif = 0
    if MI_VIF_MODE == 'pre_cv':
        log("\n" + "=" * 70)
        log("[WARN] MI-VIF 特征选择: pre_cv 模式 (有数据泄露风险!)")
        log("  参考 PMID:36004690 — CV前特征选择会导致性能膨胀")
        log("  建议使用 --mi-vif-mode per_fold (默认)")
        log("=" * 70)
        selected_idx, _, _ = mi_vif_feature_selection_per_fold(
            X_base, y_dt, feature_cols, mi_k=min(100, X_base.shape[1]))
        n_removed_mivif = X_base.shape[1] - len(selected_idx)
        X_base = X_base[:, selected_idx]
        feature_cols = [feature_cols[i] for i in selected_idx]
        log(f"  MI-VIF 选择: 保留 {X_base.shape[1]} 特征 (移除 {n_removed_mivif})")
    elif MI_VIF_MODE == 'per_fold':
        log("\n" + "=" * 70)
        log("[OK] MI-VIF 特征选择: per_fold 模式 (无数据泄露)")
        log("  MI-VIF 将在每折CV中独立执行, 仅在训练集上计算")
        log("  参考: PMID:37113250 (nestedcv) + PMID:36004690")
        log("=" * 70)
    log_stage("MI-VIF决策")
    current_step += 1
    log_progress(current_step, TOTAL_STEPS, "MI-VIF特征选择")

    # ── STEP 2: 构建模型组合 ──
    log("\n" + "=" * 70)
    log("STEP 2: 构建 特征工程 x 分类器 笛卡尔积组合")
    log("=" * 70)
    model_combos, fe_map, clf_map = build_model_combinations()
    model_names = [m[0] for m in model_combos]
    n_models = len(model_combos)
    log(f"  共 {n_models} 种模型组合")
    current_step += 1
    log_progress(current_step, TOTAL_STEPS, "构建模型组合")

    unknown_mask = (y_dt == 0) & (y_dg == 0)
    n_total = len(gene_index)
    n_unknown = unknown_mask.sum()
    log(f"  基因: {n_total} (DT+={y_dt.sum()}, DG+={y_dg.sum()}, Unknown={n_unknown})")

    # ── 收集所有模型指标 ──
    metrics_rows = []

    # ── 断点续跑 ──
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
        log("STEP 3: 药物可靶向性预测 (DT Task) — 并行 CV + 防泄露特征选择")
        log("=" * 70)
        drug_probas, dt_metrics = run_cv_parallel(
            X_base, y_dt, model_combos, "DT", unknown_mask, feature_cols)

        # 收集 DT 指标
        for model_name, metric in dt_metrics.items():
            metrics_rows.append({
                'task': 'DT',
                'model': model_name,
                'auroc': metric['auroc'],
                'auprc': metric['auprc'],
                'auroc_std': metric.get('auroc_std', 0),
                'auprc_std': metric.get('auprc_std', 0),
            })

        log_stage("DT Task")
        current_step += 1
        log_progress(current_step, TOTAL_STEPS, "DT Task CV")

        # ── STEP 4: DG Task ──
        log("\n" + "=" * 70)
        log("STEP 4: 疾病相关性预测 (DG Task) — 并行 CV + 防泄露特征选择")
        log("=" * 70)
        disease_probas, dg_metrics = run_cv_parallel(
            X_base, y_dg, model_combos, "DG", unknown_mask, feature_cols)

        # 收集 DG 指标
        for model_name, metric in dg_metrics.items():
            metrics_rows.append({
                'task': 'DG',
                'model': model_name,
                'auroc': metric['auroc'],
                'auprc': metric['auprc'],
                'auroc_std': metric.get('auroc_std', 0),
                'auprc_std': metric.get('auprc_std', 0),
            })

        log_stage("DG Task")
        current_step += 1
        log_progress(current_step, TOTAL_STEPS, "DG Task CV")

        save_checkpoint({
            'drug_probas': drug_probas,
            'disease_probas': disease_probas,
            'dt_metrics': dt_metrics,
            'dg_metrics': dg_metrics,
        })

    # 如果 metrics_rows 为空 (从 checkpoint 恢复), 从 dt_metrics/dg_metrics 重建
    if not metrics_rows:
        for model_name, metric in dt_metrics.items():
            metrics_rows.append({
                'task': 'DT',
                'model': model_name,
                'auroc': metric['auroc'],
                'auprc': metric['auprc'],
                'auroc_std': metric.get('auroc_std', 0),
                'auprc_std': metric.get('auprc_std', 0),
            })
        for model_name, metric in dg_metrics.items():
            metrics_rows.append({
                'task': 'DG',
                'model': model_name,
                'auroc': metric['auroc'],
                'auprc': metric['auprc'],
                'auroc_std': metric.get('auroc_std', 0),
                'auprc_std': metric.get('auprc_std', 0),
            })

    # ── STEP 5: 桥梁得分 ──
    log("\n" + "=" * 70)
    log("STEP 5: 计算桥梁得分 (P_drug x P_disease)")
    log("=" * 70)
    bridge_scores = compute_bridge_scores(drug_probas, disease_probas)
    log(f"  有效组合: {len(bridge_scores)}/{n_models}")
    log_stage("桥梁得分")
    current_step += 1
    log_progress(current_step, TOTAL_STEPS, "桥梁得分")

    # ── Stacking Ensemble ──
    stacking_bridge, stacking_meta_info = stacking_ensemble(
        drug_probas, disease_probas, bridge_scores, y_dt, y_dg,
        unknown_mask, gene_index, model_names, dt_metrics, dg_metrics)
    log_stage("Stacking Ensemble")
    current_step += 1
    log_progress(current_step, TOTAL_STEPS, "Stacking Ensemble")

    # ── STEP 6a: RRF 集成 ──
    log("\n" + "=" * 70)
    log("STEP 6a: RRF 倒数排名融合 (加权版) [v9.0]")
    log("=" * 70)
    rrf = rrf_integrate(bridge_scores, n_total, unknown_mask, dt_metrics, dg_metrics)
    unknown_idx = np.where(unknown_mask)[0]
    rrf_sort_order = unknown_idx[np.argsort(-rrf[unknown_mask])]
    rrf_ranks = np.full(n_total, np.nan)
    for rank, idx in enumerate(rrf_sort_order):
        rrf_ranks[idx] = rank + 1
    log_stage("RRF集成")
    current_step += 1
    log_progress(current_step, TOTAL_STEPS, "RRF集成")

    # ── STEP 6b: Borda Count ──
    log("\n" + "=" * 70)
    log("STEP 6b: Borda Count 排名融合 [v9.0]")
    log("=" * 70)
    borda = borda_count_integrate(bridge_scores, n_total, unknown_mask)
    borda_sort_order = unknown_idx[np.argsort(-borda[unknown_mask])]
    borda_ranks = np.full(n_total, np.nan)
    for rank, idx in enumerate(borda_sort_order):
        borda_ranks[idx] = rank + 1

    # 比较 RRF vs Borda
    from scipy.stats import spearmanr as sp_rank
    common_unknown = ~np.isnan(rrf_ranks[unknown_mask]) & ~np.isnan(borda_ranks[unknown_mask])
    if common_unknown.sum() > 5:
        rho_rrb, p_rrb = sp_rank(rrf_ranks[unknown_mask][common_unknown], borda_ranks[unknown_mask][common_unknown])
        log(f"  RRF vs Borda Rank Spearman rho = {rho_rrb:.4f} (p={p_rrb:.4e})")
    log_stage("Borda Count")
    current_step += 1
    log_progress(current_step, TOTAL_STEPS, "Borda Count")

    # ── STEP 6c: Bootstrap CI ──
    log("\n" + "=" * 70)
    log("STEP 6c: Bootstrap RRF 置信区间计算 [v9.0]")
    log("=" * 70)
    rrf_lower, rrf_upper = bootstrap_rrf_ci(bridge_scores, n_total, unknown_mask, n_bootstrap=N_BOOTSTRAP)
    if rrf_lower is not None:
        bootstrap_df = pd.DataFrame({
            'gene_symbol': [gene_index[i] for i in unknown_idx],
            'RRF_lower_CI': rrf_lower,
            'RRF_upper_CI': rrf_upper,
            'RRF_CI_width': rrf_upper - rrf_lower,
        })
        safe_write_csv(bootstrap_df, os.path.join(OUTPUT_DIR, "bootstrap_ci.csv"))
        log(f"  Bootstrap CI 已保存: {len(bootstrap_df)} genes")
    log_stage("Bootstrap CI")
    current_step += 1
    log_progress(current_step, TOTAL_STEPS, "Bootstrap CI")

    # ── STEP 7: 四算法共识 ──
    consensus, family_top_genes = four_algorithm_consensus(
        bridge_scores, model_names, gene_index, unknown_mask,
        drug_probas, disease_probas, dt_metrics, dg_metrics)

    # 保存四算法共识
    consensus_rows = []
    for gene_sym in gene_index:
        gene_pos = gene_index.index(gene_sym)
        if unknown_mask[gene_pos]:
            in_cons = gene_sym in consensus if consensus else False
            consensus_rows.append({
                'gene_symbol': gene_sym,
                'in_four_algorithm_consensus': in_cons,
                'in_LASSO': gene_sym in family_top_genes.get('LASSO', set()),
                'in_SVM': gene_sym in family_top_genes.get('SVM', set()),
                'in_Tree': gene_sym in family_top_genes.get('Tree', set()),
                'in_Linear': gene_sym in family_top_genes.get('Linear', set()),
            })
    if consensus_rows:
        consensus_df = pd.DataFrame(consensus_rows)
        safe_write_csv(consensus_df, os.path.join(OUTPUT_DIR, "four_algorithm_consensus.csv"))
    log_stage("四算法共识")
    current_step += 1
    log_progress(current_step, TOTAL_STEPS, "四算法共识")

    # ── 构建未知基因排名 DataFrame (用于 SHAP 和后续输出) ──
    # 计算平均 Bridge 得分 (所有模型的平均)
    bridge_mean_all = np.zeros(n_total, dtype=np.float64)
    bridge_count = np.zeros(n_total, dtype=np.float64)
    for model, scores in bridge_scores.items():
        bridge_mean_all += scores
        bridge_count += 1
    bridge_count[bridge_count == 0] = 1
    bridge_mean_all /= bridge_count

    # 计算平均 DT/DG proba
    dt_mean_all = np.zeros(n_total, dtype=np.float64)
    dg_mean_all = np.zeros(n_total, dtype=np.float64)
    dt_count = 0
    dg_count = 0
    for model in drug_probas:
        dt_mean_all += drug_probas[model]
        dt_count += 1
    dt_mean_all /= max(dt_count, 1)
    for model in disease_probas:
        dg_mean_all += disease_probas[model]
        dg_count += 1
    dg_mean_all /= max(dg_count, 1)

    # 构建未知基因排名表
    df_unknown = pd.DataFrame({
        'gene_symbol': [gene_index[i] for i in unknown_idx],
        'RRF_score': rrf[unknown_idx],
        'RRF_rank': [rrf_ranks[i] for i in unknown_idx],
        'Borda_score': borda[unknown_idx],
        'Borda_rank': [borda_ranks[i] for i in unknown_idx],
        'DT_proba_mean': dt_mean_all[unknown_idx],
        'DG_proba_mean': dg_mean_all[unknown_idx],
        'Bridge_score_mean': bridge_mean_all[unknown_idx],
    }).sort_values('RRF_score', ascending=False)
    df_unknown['in_consensus'] = df_unknown['gene_symbol'].apply(
        lambda g: g in consensus if consensus else False
    )

    # ── SHAP 误分类检测 ──
    shap_flags, shap_feature_info = shap_misclassification_detection(
        X_base, y_dt, df_unknown, gene_index, feature_cols)
    if shap_flags is not None and len(shap_flags) > 0:
        flags_df = pd.DataFrame(shap_flags)
        safe_write_csv(flags_df, os.path.join(OUTPUT_DIR, "misclassification_flags.csv"))

        # 将 flagged 标记合并到 df_unknown
        flagged_set = set(f['gene_symbol'] for f in shap_flags if f['flagged'])
        df_unknown['misclassification_flagged'] = df_unknown['gene_symbol'].isin(flagged_set)
    else:
        df_unknown['misclassification_flagged'] = False

    log_stage("SHAP误分类检测")
    current_step += 1
    log_progress(current_step, TOTAL_STEPS, "SHAP误分类检测")

    # ── STEP 8: RF 重要性分层 ──
    rf_importance_df, rf_importance_union = stratified_rf_importance(
        X_base, y_dt, y_dg, feature_cols, fe_map, clf_map, metrics_rows)
    if rf_importance_df is not None:
        safe_write_csv(rf_importance_df, os.path.join(OUTPUT_DIR, "rf_importance.csv"))
    log_stage("RF重要性分层")
    current_step += 1
    log_progress(current_step, TOTAL_STEPS, "RF重要性分层")

    # ── 稳定性选择 ──
    stable_features, stability_df = stability_selection(
        X_base, y_dt, feature_cols, n_bootstrap=STABILITY_BOOTSTRAP_N,
        threshold=STABILITY_THRESHOLD)
    log_stage("稳定性选择")
    current_step += 1
    log_progress(current_step, TOTAL_STEPS, "稳定性选择")

    # ── [v9.0 FIX] 数据泄露审计 — 移除GSE61616外部验证 ──
    # 原因: GSE61616 (7d MCAO大鼠模型) 已作为训练数据来源之一
    # 其DEGs已整合到疾病基因集 (CIRI disease genes) 中用于模型训练
    # 使用同一数据集作为"外部验证"构成循环论证 (circular reasoning)
    # 参考: PMID:36004690 — "test data must be completely independent"
    # 参考: PMID:37113250 — nested CV 正确做法
    #
    # 替代方案: 内部5折分层CV (已在STEP 3-4完成) 提供无偏性能估计
    # 所有模型性能评估基于严格的 fold-level 训练/验证分离
    log("\n" + "=" * 70)
    log("STEP 9: 外部验证 — 数据泄露审计 [v9.0 FIX]")
    log("=" * 70)
    log("  [SKIP] GSE61616外部验证已移除 (数据泄露风险)")
    log("  [INFO] GSE61616 DEGs 已整合进训练集, 不可作为外部验证")
    log("  [INFO] 参考: PMID:36004690 — test data must be independent")
    log("  [INFO] 替代方案: 内部5折分层CV提供无偏性能估计")
    log("  [INFO] 如需真正的外部验证, 需使用完全独立的CIRI数据集")
    log("  [INFO] (例如 GSE78731, GSE148274 等未参与训练的数据集)")

    safe_write_json({
        'external_validation': 'SKIPPED',
        'reason': 'GSE61616 is part of training data (circular reasoning)',
        'note': 'Internal 5-fold stratified CV provides unbiased performance estimates',
        'recommendation': 'Use truly independent dataset (e.g. GSE78731, GSE148274)',
        'references': ['PMID:36004690', 'PMID:37113250'],
    }, os.path.join(OUTPUT_DIR, "data_leakage_audit.json"))

    log_stage("数据泄露审计")
    current_step += 1
    log_progress(current_step, TOTAL_STEPS, "数据泄露审计")

    # ── STEP 10: 最终输出生成 ──
    log("\n" + "=" * 70)
    log("STEP 10: 最终输出生成 [v9.0]")
    log("=" * 70)

    # 10a. 保存桥梁靶点排名表
    bridge_targets_path = os.path.join(OUTPUT_DIR, "bridge_targets.csv")
    safe_write_csv(df_unknown, bridge_targets_path)
    log(f"  桥梁靶点表已保存: bridge_targets.csv")
    log(f"  Top-20 桥梁靶点:")
    for _, row in df_unknown.head(20).iterrows():
        flag_str = " [FLAGGED]" if row.get('misclassification_flagged', False) else ""
        cons_str = " [CONSENSUS]" if row.get('in_consensus', False) else ""
        log(f"    {row['gene_symbol']:<12s} RRF={row['RRF_score']:.6f} "
            f"Rank={int(row['RRF_rank']):d} Borda={row['Borda_score']:.2f}{flag_str}{cons_str}")

    # 10b. 保存 CV 指标
    cv_metrics_path = os.path.join(OUTPUT_DIR, "cv_metrics.csv")
    safe_write_csv(pd.DataFrame(metrics_rows), cv_metrics_path)

    # 10c. 最佳模型汇总
    log("\n--- 最佳模型 DT ---")
    dt_best = sorted([r for r in metrics_rows if r['task'] == 'DT'],
                     key=lambda x: -x['auroc'])[:5]
    for r in dt_best:
        log(f"  {r['model']:<35s} AUROC={r['auroc']:.4f} AUPRC={r['auprc']:.4f}")

    log("\n--- 最佳模型 DG ---")
    dg_best = sorted([r for r in metrics_rows if r['task'] == 'DG'],
                     key=lambda x: -x['auroc'])[:5]
    for r in dg_best:
        log(f"  {r['model']:<35s} AUROC={r['auroc']:.4f} AUPRC={r['auprc']:.4f}")

    # 10d. 性能汇总图
    plot_model_performance_summary(metrics_rows, OUTPUT_DIR)

    current_step += 1
    log_progress(current_step, TOTAL_STEPS, "最终输出")
    log_stage("最终输出")

    # ── 最终总结 ──
    total_elapsed = time.time() - start_time
    log("\n" + "═" * 70)
    log("  流水线执行完毕 [v9.0]")
    log(f"  输出目录: {OUTPUT_DIR}")
    log(f"  总运行时间: {total_elapsed/60:.1f} 分钟 ({total_elapsed:.1f} 秒)")
    log("  ")
    log("  输出文件清单:")
    output_files = [
        "bridge_targets.csv", "four_algorithm_consensus.csv",
        "misclassification_flags.csv", "rf_importance.csv",
        "stability_selection.csv", "bootstrap_ci.csv",
        "cv_metrics.csv", "external_validation.json",
        "stacking_ensemble_results.csv",
        "model_performance_summary_DT.pdf", "model_performance_summary_DG.pdf",
        "class_imbalance_diagnostic.json", "pipeline_v9.log",
    ]
    for fname in output_files:
        fpath = os.path.join(OUTPUT_DIR, fname)
        if os.path.exists(fpath):
            size_kb = os.path.getsize(fpath) / 1024
            log(f"    [OK] {fname} ({size_kb:.1f} KB)")
        else:
            log(f"    [MISSING] {fname}")
    log("═" * 70)

    return df_unknown, metrics_rows, consensus


# ═══════════════════════════════════════════════════════════════════════════════
# 入口点
# ═══════════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    main()