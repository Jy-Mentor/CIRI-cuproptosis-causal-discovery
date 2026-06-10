# -*- coding: utf-8 -*-
"""工具函数：日志系统、随机种子、特征标准化/降维。"""

import random
import logging
import sys
from pathlib import Path
from typing import Optional, Tuple, List

import numpy as np
import torch
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


def setup_logging(name: str = "hgt_pipeline", level: str = "INFO",
                  log_file: Optional[str] = None,
                  tensorboard_dir: Optional[str] = None,
                  enable_tensorboard: bool = False) -> logging.Logger:
    """配置日志系统，同时输出到控制台和文件。"""
    logger = logging.getLogger(name)
    logger.setLevel(getattr(logging, level.upper(), logging.INFO))
    logger.handlers.clear()

    fmt = logging.Formatter(
        "[%(asctime)s] %(levelname)s %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(fmt)
    logger.addHandler(ch)

    if log_file:
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(log_file, encoding="utf-8")
        fh.setFormatter(fmt)
        logger.addHandler(fh)

    return logger


def setup_tensorboard(log_dir: str) -> Optional[object]:
    """初始化 TensorBoard SummaryWriter。"""
    try:
        from torch.utils.tensorboard import SummaryWriter
        writer = SummaryWriter(log_dir=log_dir)
        return writer
    except ImportError:
        return None


def set_seed(seed: int) -> None:
    """固定所有随机种子以确保可复现性。"""
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


def zscore_normalize(arr: np.ndarray, eps: float = 1e-8) -> np.ndarray:
    """按特征维度执行 Z-score 标准化。"""
    mean = arr.mean(axis=0, keepdims=True)
    std = arr.std(axis=0, keepdims=True)
    return (arr - mean) / (std + eps)


def minmax_normalize(arr: np.ndarray) -> np.ndarray:
    """按特征维度执行 Min-Max 标准化到 [0, 1]。"""
    vmin = arr.min(axis=0, keepdims=True)
    vmax = arr.max(axis=0, keepdims=True)
    denom = vmax - vmin
    denom[denom == 0] = 1.0
    return (arr - vmin) / denom


def pca_reduce(arr: np.ndarray, target_dim: int, seed: int = 42) -> Tuple[np.ndarray, PCA]:
    """PCA 降维并返回降维后数组和拟合的 PCA 对象。"""
    if arr.shape[1] <= target_dim:
        return arr, None
    n_samples = arr.shape[0]
    if n_samples < target_dim:
        return arr[:, :min(target_dim, arr.shape[1])].copy(), None
    pca = PCA(n_components=target_dim, random_state=seed)
    result = pca.fit_transform(arr).astype(np.float32)
    return result, pca


def pca_then_umap(arr: np.ndarray, pca_dim: int, umap_dim: int,
                  seed: int = 42) -> np.ndarray:
    """先 PCA 降至中等维度，再 UMAP 保留局部结构。"""
    if arr.shape[1] <= pca_dim:
        pca_arr = arr
    else:
        pca = PCA(n_components=pca_dim, random_state=seed)
        pca_arr = pca.fit_transform(arr).astype(np.float64)

    try:
        import umap
        reducer = umap.UMAP(n_components=umap_dim, random_state=seed, n_jobs=1)
        result = reducer.fit_transform(pca_arr).astype(np.float32)
        return result
    except ImportError:
        import logging
        logging.getLogger("hgt_pipeline").warning("UMAP not installed, falling back to PCA")
        return pca_reduce(arr, umap_dim, seed)[0]


def normalize_features(arr: np.ndarray, method: str = "zscore") -> np.ndarray:
    """对特征矩阵执行标准化。"""
    if method == "zscore":
        return zscore_normalize(arr)
    elif method == "minmax":
        return minmax_normalize(arr)
    return arr


def try_compile(model: torch.nn.Module) -> torch.nn.Module:
    """若 PyTorch >= 2.0 且 CUDA 可用且 Triton 已安装，则使用 torch.compile 加速。

    注意: dynamic=True 是必须的，因为 HGT 包含大量动态形状。
    首次运行可能稍慢，后续 epoch 可提速 20-30%。
    Windows 上 Triton 通常不可用，此时回退到 eager 模式。
    可通过环境变量 TORCH_COMPILE=1 强制启用（需先安装 Triton）。
    """
    import os
    force_compile = os.environ.get("TORCH_COMPILE", "0") == "1"

    if not force_compile and not hasattr(torch, "compile"):
        return model

    if hasattr(torch, "compile") and torch.cuda.is_available():
        try:
            import triton  # noqa: F401
            return torch.compile(model, dynamic=True)
        except (ImportError, Exception) as e:
            import logging
            _logger = logging.getLogger("hgt_pipeline")
            _logger.warning(
                "torch.compile not available: %s. "
                "Set TORCH_COMPILE=1 and install triton to enable. "
                "Falling back to eager mode.",
                str(e) if str(e) else "Triton not installed",
            )
    return model


def setup_wandb(config: object, logger: logging.Logger) -> Optional[object]:
    """初始化 Weights & Biases 实验追踪。

    若 wandb 未安装或配置禁用，则返回 None。
    参考: Biewald, "Experiment Tracking with Weights and Biases", 2020
    """
    log_cfg = config.logging
    if not log_cfg.wandb:
        return None
    try:
        import wandb
    except ImportError:
        logger.warning("wandb not installed, skipping W&B integration")
        return None

    m_cfg = config.model
    t_cfg = config.training

    wandb.init(
        project=log_cfg.wandb_project,
        entity=log_cfg.wandb_entity,
        config={
            "model": {
                "hidden_dim": m_cfg.hidden_dim,
                "num_heads": m_cfg.num_heads,
                "num_layers": m_cfg.num_layers,
                "dropout": m_cfg.dropout,
                "initial_residual": m_cfg.initial_residual,
                "decoder_bias": m_cfg.decoder_bias,
                "decoder_factorization": m_cfg.decoder_factorization,
                "use_input_bn": m_cfg.use_input_bn,
            },
            "training": {
                "epochs": t_cfg.epochs,
                "lr": t_cfg.lr,
                "weight_decay": t_cfg.weight_decay,
                "focal_alpha": t_cfg.focal_alpha,
                "focal_gamma": t_cfg.focal_gamma,
                "adaptive_neg_ratio": t_cfg.adaptive_neg_ratio,
                "early_stop_metric": t_cfg.early_stop_metric,
                "composite_weight": t_cfg.composite_weight,
            },
            "seed": config.seed,
            "device": config.device,
        },
    )
    logger.info("W&B initialized: project=%s", log_cfg.wandb_project)
    return wandb