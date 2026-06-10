"""
GPU 检测模块
=============
负责在模块加载时检测 GPU 可用性，并缓存检测结果供后续使用。

设计原则:
  - 一次性检测, 模块级缓存 (全局变量)
  - 自动回退 CPU 机制
  - 显存安全模式检测 (多进程 GPU 并发 → OOM 风险)
"""

import os
import warnings
import numpy as np
from typing import Optional

from .config import GPUConfig
from .utils import log


def detect_gpu_environment(gpu_enabled: bool = True, n_jobs: int = -1) -> GPUConfig:
    """
    运行所有 GPU 检测, 自动检测显存冲突风险.

    Args:
        gpu_enabled: 是否启用 GPU
        n_jobs: 并行核心数 (-1 = 全部)

    Returns:
        GPUConfig 对象
    """
    gpu_cfg = GPUConfig(enabled=gpu_enabled)

    log("-" * 50)
    log("GPU 环境检测")
    log("-" * 50)

    # 计算有效并行数
    effective_jobs = os.cpu_count() if n_jobs == -1 else n_jobs
    log(f"  并行核心: {n_jobs} → 有效 {effective_jobs} 进程")

    cuda_avail, gpu_name, total_vram = _detect_pytorch_cuda()

    if not gpu_enabled:
        log("  [GPU] GPU_ENABLED=False, 所有模型使用 CPU")
        return gpu_cfg

    # ── 显存安全检测: 多进程 GPU 并发 → OOM 风险 ──
    gpu_parallel_threshold = 4  # 8GB VRAM 安全阈值
    if cuda_avail and effective_jobs > gpu_parallel_threshold and total_vram is not None:
        log(f"  ⚠️  显存安全模式: N_JOBS={effective_jobs} > {gpu_parallel_threshold}")
        log(f"     多进程 GPU 并发可能耗尽 {total_vram:.0f}GB VRAM")
        log(f"     → GPU 模型自动回退 CPU, 纯 CPU 模型仍可 24 进程并行 ✅")
        xgb_gpu_safe = False
        lgb_gpu_safe = False
    else:
        xgb_gpu_safe = gpu_enabled
        lgb_gpu_safe = gpu_enabled

    # XGBoost GPU 检测
    _detect_xgb_gpu(gpu_cfg, enabled=xgb_gpu_safe)

    # LightGBM GPU 检测
    _detect_lgb_gpu(gpu_cfg, enabled=lgb_gpu_safe)

    return gpu_cfg


def _detect_pytorch_cuda() -> tuple:
    """检测 PyTorch CUDA 可用性"""
    try:
        import torch
        cuda_avail = torch.cuda.is_available()
        cuda_ver = torch.version.cuda or "N/A"
        n_gpu = torch.cuda.device_count() if cuda_avail else 0
        total_vram = None
        gpu_name = "N/A"

        if cuda_avail and n_gpu > 0:
            gpu_name = torch.cuda.get_device_name(0)
            total_vram = torch.cuda.get_device_properties(0).total_memory / 1e9
            log(f"  PyTorch CUDA: 可用 ✅ (v{cuda_ver}, {n_gpu}× {gpu_name}, {total_vram:.1f}GB VRAM)")
        else:
            log(f"  PyTorch CUDA: 不可用")
        return cuda_avail, gpu_name, total_vram
    except ImportError:
        log("  PyTorch CUDA: 未安装")
        return False, "N/A", None


def _detect_xgb_gpu(gpu_cfg: GPUConfig, enabled: bool = True):
    """检测 XGBoost GPU 可用性"""
    if not enabled:
        gpu_cfg.xgb_available = False
        gpu_cfg.xgb_params = {'tree_method': 'hist', 'predictor': 'cpu_predictor'}
        gpu_cfg.xgb_info = "显存安全模式: 回退 CPU"
        log(f"  [GPU] XGBoost: {gpu_cfg.xgb_info}")
        return

    try:
        import xgboost as xgb
        params = {
            'tree_method': 'gpu_hist', 'predictor': 'gpu_predictor',
            'n_estimators': 1, 'verbosity': 0, 'n_jobs': 1,
        }
        X_dummy = np.random.randn(10, 5).astype(np.float32)
        y_dummy = np.random.randint(0, 2, 10)
        clf = xgb.XGBClassifier(**params)
        clf.fit(X_dummy, y_dummy)

        gpu_cfg.xgb_available = True
        gpu_cfg.xgb_params = {'tree_method': 'gpu_hist', 'predictor': 'gpu_predictor'}
        try:
            import subprocess
            result = subprocess.run(
                ['nvidia-smi', '--query-gpu=name,compute_cap,driver_version',
                 '--format=csv,noheader'],
                capture_output=True, text=True, timeout=10,
            )
            gpu_cfg.xgb_info = result.stdout.strip().replace(',', ' | ')
        except Exception:
            gpu_cfg.xgb_info = "GPU 可用 (nvidia-smi 信息获取失败)"
        log(f"  [GPU] XGBoost GPU 加速: 启用 ✅ | {gpu_cfg.xgb_info}")
    except Exception as e:
        gpu_cfg.xgb_available = False
        gpu_cfg.xgb_params = {'tree_method': 'hist', 'predictor': 'cpu_predictor'}
        err_msg = str(e).split('\n')[0]
        gpu_cfg.xgb_info = f"回退到 CPU ({err_msg})"
        log(f"  [GPU] XGBoost GPU 加速: 不可用 ⚠️ → {gpu_cfg.xgb_info}")


def _detect_lgb_gpu(gpu_cfg: GPUConfig, enabled: bool = True):
    """检测 LightGBM GPU 可用性"""
    if not enabled:
        gpu_cfg.lgb_available = False
        gpu_cfg.lgb_params = {}
        gpu_cfg.lgb_info = "显存安全模式: 回退 CPU"
        log(f"  [GPU] LightGBM: {gpu_cfg.lgb_info}")
        return

    try:
        import lightgbm as lgb
        params = {
            'device': 'gpu', 'gpu_platform_id': 0, 'gpu_device_id': 0,
            'n_estimators': 1, 'verbose': -1,
        }
        X_dummy = np.random.randn(10, 5).astype(np.float32)
        y_dummy = np.random.randint(0, 2, 10)
        clf = lgb.LGBMClassifier(**params)
        clf.fit(X_dummy, y_dummy)

        gpu_cfg.lgb_available = True
        gpu_cfg.lgb_params = {'device': 'gpu', 'gpu_platform_id': 0, 'gpu_device_id': 0}
        gpu_cfg.lgb_info = "GPU 加速启用"
        log(f"  [GPU] LightGBM GPU 加速: 启用 ✅")
    except Exception as e:
        gpu_cfg.lgb_available = False
        gpu_cfg.lgb_params = {}
        err_msg = str(e).split('\n')[0]
        gpu_cfg.lgb_info = f"回退到 CPU ({err_msg})"
        log(f"  [GPU] LightGBM GPU 加速: 不可用 ⚠️ → {gpu_cfg.lgb_info}")