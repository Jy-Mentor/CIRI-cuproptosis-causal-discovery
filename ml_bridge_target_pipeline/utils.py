"""
工具函数模块
=============
提供日志记录、样本权重计算等通用工具函数。
"""

import time
import warnings
import numpy as np
from sklearn.exceptions import ConvergenceWarning

# 全局警告过滤 (模块加载时执行一次)
warnings.filterwarnings("ignore")
warnings.filterwarnings("ignore", category=ConvergenceWarning)


def log(msg: str):
    """带时间戳的日志输出"""
    timestamp = time.strftime("%H:%M:%S")
    print(f"[{timestamp}] {msg}", flush=True)


def compute_sample_weight(y: np.ndarray) -> np.ndarray:
    """
    计算样本权重 (用于不支持 class_weight 的分类器)

    使用 balanced 策略: 每个类别的样本权重 = 总样本数 / (2 * 该类样本数)

    Args:
        y: 标签数组 (0/1)

    Returns:
        样本权重数组
    """
    n_pos = y.sum()
    n_neg = len(y) - n_pos
    return np.where(
        y == 1,
        len(y) / (2 * max(n_pos, 1)),
        len(y) / (2 * max(n_neg, 1))
    )


def setup_worker_warnings():
    """
    在 joblib worker 进程中配置警告过滤
    (joblib 不继承主进程的 warnings filter)
    """
    warnings.filterwarnings("ignore", category=ConvergenceWarning)
    warnings.filterwarnings("ignore", category=FutureWarning)
    warnings.filterwarnings("ignore", message="X does not have valid feature names")


def print_header(title: str, width: int = 70):
    """打印分隔标题"""
    log("=" * width)
    log(title)
    log("=" * width)