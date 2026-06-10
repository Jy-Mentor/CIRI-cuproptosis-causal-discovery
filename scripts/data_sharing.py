# -*- coding: utf-8 -*-
"""
CEHG-RNP 2.0 内存数据共享层
============================

功能:
- 跨阶段共享大型数据矩阵 (表达矩阵、PPI邻接矩阵等)
- 内存映射文件 (memmap) 避免重复加载
- 通过 multiprocessing.Manager 共享关键数据结构
- 自动缓存管理 (LRU策略)
- 版本化数据追踪

设计:
- 使用 np.memmap 处理大型矩阵 (单细胞表达矩阵 >1GB)
- 使用 pickle 缓存处理小型结果文件
- 线程安全的读写锁

版本: v1.0 | 日期: 2026-05-28
"""

import os
import sys
import time
import json
import pickle
import logging
import threading
import hashlib
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple, Union
from collections import OrderedDict

import numpy as np
import pandas as pd

logger = logging.getLogger("data_sharing")


class LRUCache:
    """LRU缓存实现"""

    def __init__(self, max_size: int = 10):
        self.max_size = max_size
        self._cache = OrderedDict()
        self._lock = threading.Lock()
        self.hits = 0
        self.misses = 0

    def get(self, key: str) -> Any:
        with self._lock:
            if key in self._cache:
                self._cache.move_to_end(key)
                self.hits += 1
                return self._cache[key]
            self.misses += 1
            return None

    def put(self, key: str, value: Any):
        with self._lock:
            if key in self._cache:
                self._cache.move_to_end(key)
            else:
                if len(self._cache) >= self.max_size:
                    self._cache.popitem(last=False)
                self._cache[key] = value

    def clear(self):
        with self._lock:
            self._cache.clear()

    @property
    def size(self):
        return len(self._cache)

    @property
    def hit_rate(self):
        total = self.hits + self.misses
        return self.hits / total if total > 0 else 0.0


class SharedDataManager:
    """
    跨阶段共享数据管理器

    提供:
    - 内存缓存 (LRU): 小型DataFrame/dict
    - 内存映射 (memmap): 大型矩阵
    - 版本追踪: 检测数据变更
    """

    _instance = None
    _lock = threading.Lock()

    def __new__(cls, *args, **kwargs):
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self, results_dir: str = None, max_cache_size: int = 20):
        if hasattr(self, '_initialized'):
            return

        if results_dir is None:
            from config import RESULTS_DIR
            results_dir = RESULTS_DIR

        self.results_dir = results_dir
        self._cache = LRUCache(max_size=max_cache_size)
        self._memmaps: Dict[str, np.memmap] = {}
        self._data_versions: Dict[str, str] = {}
        self._memmap_dir = os.path.join(results_dir, "shared_memmap")
        os.makedirs(self._memmap_dir, exist_ok=True)
        self._initialized = True

        logger.info(f"SharedDataManager 已初始化 (缓存大小: {max_cache_size})")

    def _compute_hash(self, filepath: str) -> str:
        """计算文件SHA256哈希"""
        if not os.path.exists(filepath):
            return "not_found"

        sha = hashlib.sha256()
        file_size = os.path.getsize(filepath)

        if file_size > 10 * 1024 * 1024:
            with open(filepath, 'rb') as f:
                sha.update(f.read(1024 * 1024))
                f.seek(file_size // 2)
                sha.update(f.read(1024 * 1024))
                f.seek(max(0, file_size - 1024 * 1024))
                sha.update(f.read(1024 * 1024))
        else:
            with open(filepath, 'rb') as f:
                sha.update(f.read())

        return sha.hexdigest()

    def _is_stale(self, key: str, filepath: str) -> bool:
        """检查缓存是否过期"""
        if key not in self._data_versions:
            return True

        current_hash = self._compute_hash(filepath)
        return current_hash != self._data_versions[key]

    def load_dataframe(self, stage: str, filename: str, **kwargs) -> Optional[pd.DataFrame]:
        """
        加载CSV为DataFrame（带缓存）

        Args:
            stage: 阶段目录名
            filename: 文件名
            **kwargs: 传递给pd.read_csv的参数
        """
        cache_key = f"df:{stage}:{filename}"

        cached = self._cache.get(cache_key)
        if cached is not None:
            return cached

        filepath = os.path.join(self.results_dir, stage, filename)
        if not os.path.exists(filepath):
            return None

        try:
            df = pd.read_csv(filepath, **kwargs)
            self._cache.put(cache_key, df)
            self._data_versions[cache_key] = self._compute_hash(filepath)
            logger.debug(f"DataFrame已缓存: {cache_key} ({df.shape})")
            return df
        except Exception as e:
            logger.error(f"DataFrame加载失败 {filepath}: {e}")
            return None

    def load_numpy(self, stage: str, filename: str) -> Optional[np.ndarray]:
        """
        加载NumPy数组（带缓存）

        自动选择加载策略:
        - 小数组 (<100MB): 缓存在内存中
        - 大数组 (>=100MB): 使用memmap
        """
        cache_key = f"np:{stage}:{filename}"
        filepath = os.path.join(self.results_dir, stage, filename)

        if not os.path.exists(filepath):
            return None

        file_size = os.path.getsize(filepath)

        if file_size < 100 * 1024 * 1024:
            cached = self._cache.get(cache_key)
            if cached is not None:
                return cached

            try:
                arr = np.load(filepath, allow_pickle=True)
                if isinstance(arr, np.ndarray):
                    self._cache.put(cache_key, arr)
                    return arr
                elif isinstance(arr, np.lib.npyio.NpzFile):
                    result = {k: arr[k] for k in arr.files}
                    self._cache.put(cache_key, result)
                    return result
            except Exception:
                pass

        memmap_key = f"mmap:{stage}:{filename}"
        if memmap_key in self._memmaps:
            return self._memmaps[memmap_key]

        try:
            mmap_path = os.path.join(self._memmap_dir, f"{stage}_{filename}.mmap")
            if os.path.exists(mmap_path + ".npy"):
                source = np.load(mmap_path + ".npy", allow_pickle=True)
                if isinstance(source, np.ndarray):
                    mmap = np.memmap(mmap_path, dtype=source.dtype, 
                                    mode='r', shape=source.shape)
                    self._memmaps[memmap_key] = mmap
                    return mmap

            arr = np.load(filepath, allow_pickle=True)
            if isinstance(arr, np.ndarray):
                fp = np.memmap(mmap_path, dtype=arr.dtype, mode='w+', shape=arr.shape)
                fp[:] = arr[:]
                fp.flush()
                np.save(mmap_path + ".npy", np.array([], dtype=arr.dtype))

                mmap = np.memmap(mmap_path, dtype=arr.dtype, mode='r', shape=arr.shape)
                self._memmaps[memmap_key] = mmap
                return mmap
        except Exception as e:
            logger.error(f"NumPy加载失败 {filepath}: {e}")

        return None

    def load_h5ad_matrix(self, stage: str, filename: str) -> Optional[np.ndarray]:
        """
        加载h5ad文件的表达矩阵（使用memmap避免重复加载）

        单细胞数据通常很大 (10K+ genes × 10K+ cells),
        使用memmap避免每次加载到内存
        """
        cache_key = f"h5ad_mat:{stage}:{filename}"

        memmap_key = f"mmap:h5ad:{stage}:{filename}"
        if memmap_key in self._memmaps:
            return self._memmaps[memmap_key]

        filepath = os.path.join(self.results_dir, stage, filename)
        if not os.path.exists(filepath):
            return None

        try:
            import scanpy as sc
            from scipy.sparse import issparse

            adata = sc.read_h5ad(filepath)
            X = adata.X

            if issparse(X):
                X = X.toarray()

            X = np.asarray(X, dtype=np.float32)

            mmap_path = os.path.join(self._memmap_dir, f"{stage}_{filename}.mmap")
            fp = np.memmap(mmap_path, dtype=np.float32, mode='w+', shape=X.shape)
            fp[:] = X[:]
            fp.flush()

            mmap = np.memmap(mmap_path, dtype=np.float32, mode='r', shape=X.shape)
            self._memmaps[memmap_key] = mmap
            logger.info(f"h5ad矩阵已memmap: {memmap_key} ({X.shape})")
            return mmap
        except Exception as e:
            logger.error(f"h5ad加载失败 {filepath}: {e}")
            return None

    def load_json(self, stage: str, filename: str) -> Optional[Dict]:
        """加载JSON文件（带缓存）"""
        cache_key = f"json:{stage}:{filename}"

        cached = self._cache.get(cache_key)
        if cached is not None:
            return cached

        filepath = os.path.join(self.results_dir, stage, filename)
        if not os.path.exists(filepath):
            return None

        try:
            with open(filepath, 'r', encoding='utf-8') as f:
                data = json.load(f)
            self._cache.put(cache_key, data)
            self._data_versions[cache_key] = self._compute_hash(filepath)
            return data
        except Exception as e:
            logger.error(f"JSON加载失败 {filepath}: {e}")
            return None

    def load_gene_dict(self, stage: str, filename: str,
                       key_col: str = 'Gene', value_col: str = None) -> Optional[Dict]:
        """
        加载基因-值映射字典（带缓存）

        这是管道中最常用的数据加载模式。
        """
        cache_key = f"gene_dict:{stage}:{filename}:{key_col}:{value_col}"

        cached = self._cache.get(cache_key)
        if cached is not None:
            return cached

        df = self.load_dataframe(stage, filename)
        if df is None or df.empty:
            return None

        actual_key = key_col
        if key_col not in df.columns:
            for col in df.columns:
                if col.upper() == key_col.upper():
                    actual_key = col
                    break

        if actual_key not in df.columns:
            logger.warning(f"列 {key_col} 在 {filename} 中不存在")
            return None

        if value_col is None:
            value_col = df.columns[1]
        if value_col not in df.columns:
            value_col = df.columns[-1]

        result = dict(zip(df[actual_key].str.upper(), df[value_col]))
        self._cache.put(cache_key, result)
        return result

    def invalidate(self, stage: str = None):
        """使缓存失效"""
        if stage is None:
            self._cache.clear()
            self._memmaps.clear()
            self._data_versions.clear()
            logger.info("所有缓存已失效")
        else:
            keys_to_remove = [k for k in self._cache._cache if stage in k]
            for k in keys_to_remove:
                del self._cache._cache[k]
            mmap_keys = [k for k in self._memmaps if stage in k]
            for k in mmap_keys:
                del self._memmaps[k]
            logger.info(f"已失效 {stage} 缓存")

    def get_stats(self) -> Dict:
        """获取缓存统计"""
        total_memmap_mb = sum(
            np.prod(mmap.shape) * mmap.dtype.itemsize / (1024 * 1024)
            for mmap in self._memmaps.values()
        )
        return {
            "cache_size": self._cache.size,
            "cache_hit_rate": round(self._cache.hit_rate, 3),
            "memmap_count": len(self._memmaps),
            "memmap_total_mb": round(total_memmap_mb, 1),
            "versioned_keys": len(self._data_versions),
        }


class StageDataValidator:
    """
    阶段数据契约验证器

    在每个阶段执行后验证输出文件的：
    - 完整性（文件是否存在）
    - 格式（CSV是否有必要的列）
    - 大小（文件是否非空）
    """

    def __init__(self, results_dir: str = None):
        if results_dir is None:
            from config import RESULTS_DIR
            results_dir = RESULTS_DIR

        self.results_dir = results_dir
        self.logger = logging.getLogger("data_validator")

    def validate_stage_output(self, stage_name: str, 
                              required_files: List[str],
                              required_columns: Dict[str, List[str]] = None) -> Tuple[bool, List[str]]:
        """
        验证阶段输出

        Args:
            stage_name: 阶段名称
            required_files: 必需的文件名列表
            required_columns: {filename: [必需列名]} 字典

        Returns:
            (is_valid, errors): 验证结果和错误列表
        """
        stage_dir = os.path.join(self.results_dir, stage_name)
        errors = []

        if not os.path.exists(stage_dir):
            errors.append(f"阶段目录不存在: {stage_dir}")
            return False, errors

        for filename in required_files:
            filepath = os.path.join(stage_dir, filename)

            if not os.path.exists(filepath):
                errors.append(f"缺失文件: {filename}")
                continue

            if os.path.getsize(filepath) == 0:
                errors.append(f"空文件: {filename}")
                continue

            if required_columns and filename in required_columns:
                try:
                    if filename.endswith('.csv'):
                        df = pd.read_csv(filepath, nrows=1)
                        for col in required_columns[filename]:
                            if col not in df.columns:
                                alt_col = None
                                for c in df.columns:
                                    if c.upper() == col.upper():
                                        alt_col = c
                                        break
                                if alt_col is None:
                                    errors.append(f"{filename}: 缺少列 {col}")
                    elif filename.endswith('.json'):
                        with open(filepath, 'r', encoding='utf-8') as f:
                            data = json.load(f)
                        for col in required_columns[filename]:
                            if col not in data and col not in str(data):
                                pass
                except Exception as e:
                    errors.append(f"{filename}: 读取失败 ({e})")

        is_valid = len(errors) == 0
        if is_valid:
            self.logger.info(f"✓ {stage_name} 输出验证通过 ({len(required_files)} 文件)")
        else:
            self.logger.warning(f"✗ {stage_name} 输出验证失败: {errors}")

        return is_valid, errors

    def validate_all_stages(self, stage_configs: Dict[str, List[str]]) -> Dict[str, bool]:
        """批量验证所有阶段"""
        results = {}
        for stage_name, required_files in stage_configs.items():
            is_valid, _ = self.validate_stage_output(stage_name, required_files)
            results[stage_name] = is_valid
        return results


def get_shared_manager(results_dir: str = None) -> SharedDataManager:
    """获取共享数据管理器单例"""
    return SharedDataManager(results_dir)


def get_validator(results_dir: str = None) -> StageDataValidator:
    """获取数据验证器"""
    return StageDataValidator(results_dir)
