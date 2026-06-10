"""
模型训练模块
=============
负责单模型训练、交叉验证和多进程并行训练。

防泄漏设计:
  - StandardScaler / 特征工程 / 分类器 均在每折内部独立实例化
  - 各模型随机种子独立: seed + idx (避免多进程随机状态冲突)
  - 未知基因概率在每折累加, 循环结束后统一除以 n_folds
"""

import numpy as np
from typing import List, Tuple, Dict, Any, Optional, Callable

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.cross_decomposition import PLSRegression
from joblib import Parallel, delayed

from .config import TrainingConfig, GPUConfig, PipelineConfig
from .utils import log, compute_sample_weight, setup_worker_warnings
from .feature_engineering import FeatureEngineeringStrategy, SUPERVISED_FE_TYPES
from .classifiers import ClassifierStrategy


def _apply_fe_transformer(X_train, X_val, X_unk, fe_strategy, random_state, y_train=None):
    """
    应用特征工程转换

    Args:
        X_train: 训练集
        X_val: 验证集
        X_unk: 未知基因集
        fe_strategy: 特征工程策略对象
        random_state: 随机种子
        y_train: 训练标签 (有监督特征工程需要)

    Returns:
        (X_train_fe, X_val_fe, X_unk_fe)
    """
    try:
        transformer = fe_strategy.create_transformer(random_state)

        if transformer is None:
            return X_train, X_val, X_unk

        if fe_strategy.is_supervised and y_train is not None:
            transformer.fit(X_train, y_train)
        else:
            transformer.fit(X_train)

        X_train_fe = transformer.transform(X_train)
        X_val_fe = transformer.transform(X_val)
        X_unk_fe = transformer.transform(X_unk)

        # PLS 返回 (X_scores, Y_scores) 元组, 需要取第一个
        if isinstance(transformer, PLSRegression):
            if isinstance(X_train_fe, tuple):
                X_train_fe = X_train_fe[0]
            if isinstance(X_val_fe, tuple):
                X_val_fe = X_val_fe[0]
            if isinstance(X_unk_fe, tuple):
                X_unk_fe = X_unk_fe[0]

        return X_train_fe, X_val_fe, X_unk_fe
    except Exception:
        return X_train, X_val, X_unk


def _train_single_model_worker(
    model_name: str,
    fe_strategy: FeatureEngineeringStrategy,
    clf_strategy: ClassifierStrategy,
    X_base: np.ndarray,
    y: np.ndarray,
    unknown_mask: np.ndarray,
    model_idx: int,
    training_cfg: TrainingConfig,
    gpu_cfg: Optional[GPUConfig] = None,
) -> Tuple[str, np.ndarray, Dict[str, float]]:
    """
    单个模型组合的完整 5 折 CV 流程 (worker 函数)

    完全自包含设计, 适合 joblib.Parallel 多进程调用。

    Args:
        model_name: 模型名称
        fe_strategy: 特征工程策略
        clf_strategy: 分类器策略
        X_base: 基础特征矩阵
        y: 标签
        unknown_mask: 未知样本掩码
        model_idx: 模型索引 (用于种子偏移)
        training_cfg: 训练配置
        gpu_cfg: GPU 配置 (可选)

    Returns:
        (model_name, probas, {'auroc': ..., 'auprc': ...})
    """
    # worker 进程内禁用警告
    setup_worker_warnings()

    seed = training_cfg.seed + model_idx
    n_genes = X_base.shape[0]
    probas = np.zeros(n_genes, dtype=np.float64)
    labeled_assigned = np.zeros(n_genes, dtype=bool)

    skf = StratifiedKFold(n_splits=training_cfg.n_folds, shuffle=True, random_state=seed)
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

        # --- 1. StandardScaler ---
        scaler = StandardScaler()
        X_train_s = scaler.fit_transform(X_train)
        X_val_s = scaler.transform(X_val)
        X_unk_s = scaler.transform(X_base[unknown_pos])

        # --- 2. 特征工程 ---
        X_train_fe, X_val_fe, X_unk_fe = _apply_fe_transformer(
            X_train_s, X_val_s, X_unk_s, fe_strategy, seed, y_train,
        )

        # --- 3. 训练分类器 ---
        try:
            clf = clf_strategy.create_classifier(seed, gpu_cfg)
            if clf_strategy.requires_sample_weight:
                clf.fit(X_train_fe, y_train, sample_weight=compute_sample_weight(y_train))
            else:
                clf.fit(X_train_fe, y_train)
        except Exception:
            continue

        # --- 4. 预测 ---
        try:
            y_val_prob = clf.predict_proba(X_val_fe)[:, 1]
            probas[val_pos] = y_val_prob
            labeled_assigned[val_pos] = True

            y_unk_prob = clf.predict_proba(X_unk_fe)[:, 1]
            probas[unknown_pos] += y_unk_prob

            fold_success += 1
        except Exception:
            continue

        # --- 5. 评估 ---
        try:
            fold_aurocs.append(roc_auc_score(y_val, y_val_prob))
            fold_auprcs.append(average_precision_score(y_val, y_val_prob))
        except Exception:
            pass

    # --- 折循环结束 ---
    if fold_success > 0:
        probas[unknown_pos] /= fold_success
    else:
        return model_name, probas, {'auroc': 0.0, 'auprc': 0.0}

    # 填充未赋值的标签基因
    if not labeled_assigned[~unknown_mask].all():
        global_mean = probas[~unknown_mask][labeled_assigned[~unknown_mask]].mean()
        unassigned_idx = np.where(~unknown_mask & ~labeled_assigned)[0]
        probas[unassigned_idx] = global_mean

    mean_auroc = np.mean(fold_aurocs) if fold_aurocs else 0.0
    mean_auprc = np.mean(fold_auprcs) if fold_auprcs else 0.0
    return model_name, probas, {'auroc': mean_auroc, 'auprc': mean_auprc}


class ModelTrainer:
    """模型训练器 — 管理多进程并行训练"""

    def __init__(self, config: PipelineConfig):
        self.config = config

    def build_model_combinations(
        self,
        fe_registry,
        clf_registry,
    ) -> Tuple[List[Tuple[str, Any, Any]], dict, dict]:
        """
        构建特征工程 × 分类器的笛卡尔积组合

        Returns:
            (combinations, fe_map, clf_map):
            - combinations: list of (model_name, fe_strategy, clf_strategy)
            - fe_map: {fe_name: fe_strategy}
            - clf_map: {clf_name: clf_strategy}
        """
        fe_items = list(fe_registry.get_all_strategies().items())
        clf_items = list(clf_registry.get_all_strategies().items())

        combinations = []
        for fe_name, fe_strategy in fe_items:
            for clf_name, clf_strategy in clf_items:
                model_name = f"{fe_name}__{clf_name}"
                combinations.append((model_name, fe_strategy, clf_strategy))

        fe_map = dict(fe_items)
        clf_map = dict(clf_items)

        log(f"  特征工程策略: {len(fe_items)}")
        log(f"  基础分类器: {len(clf_items)}")
        log(f"  笛卡尔积组合: {len(combinations)} 种模型")

        return combinations, fe_map, clf_map

    def run_cv_parallel(
        self,
        model_combos: List[Tuple[str, FeatureEngineeringStrategy, ClassifierStrategy]],
        X_base: np.ndarray,
        y: np.ndarray,
        unknown_mask: np.ndarray,
        task_name: str,
    ) -> Tuple[Dict[str, np.ndarray], Dict[str, Dict[str, float]]]:
        """
        并行执行所有模型组合的交叉验证

        Args:
            model_combos: 模型组合列表
            X_base: 基础特征矩阵
            y: 标签
            unknown_mask: 未知样本掩码
            task_name: 任务名称 (用于日志)

        Returns:
            (all_probas, all_metrics):
            - all_probas: {model_name: probas_array}
            - all_metrics: {model_name: {'auroc': ..., 'auprc': ...}}
        """
        n_models = len(model_combos)
        n_jobs = self.config.training.n_jobs
        log(f"  [{task_name}] 启动 {n_models} 个模型的并行训练 ({n_jobs} 核)...")

        results = Parallel(n_jobs=n_jobs, verbose=5)(
            delayed(_train_single_model_worker)(
                name, fe, clf, X_base, y, unknown_mask, idx,
                self.config.training, self.config.gpu,
            )
            for idx, (name, fe, clf) in enumerate(model_combos)
        )

        all_probas = {}
        all_metrics = {}
        for model_name, probas, metric in results:
            all_probas[model_name] = probas
            all_metrics[model_name] = metric
            log(f"  [{task_name}] {model_name}: AUROC={metric['auroc']:.4f}, AUPRC={metric['auprc']:.4f}")

        return all_probas, all_metrics