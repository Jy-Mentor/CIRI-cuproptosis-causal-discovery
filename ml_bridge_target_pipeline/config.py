"""
配置管理模块
==============
封装所有配置项，消除全局状态变量。
通过 dataclass 实现类型安全，支持环境变量覆盖。
"""

import os
from dataclasses import dataclass, field
from typing import Dict, Any, Optional


@dataclass
class GPUConfig:
    """GPU 配置"""
    enabled: bool = True
    xgb_params: Dict[str, Any] = field(default_factory=lambda: {
        'tree_method': 'hist', 'predictor': 'cpu_predictor'
    })
    lgb_params: Dict[str, Any] = field(default_factory=dict)
    xgb_available: bool = False
    lgb_available: bool = False
    xgb_info: str = "未检测"
    lgb_info: str = "未检测"


@dataclass
class TrainingConfig:
    """训练配置"""
    seed: int = 42
    n_folds: int = 5
    rrf_k: int = 60
    n_jobs: int = -1  # -1 = 使用所有核心


@dataclass
class PathConfig:
    """路径配置"""
    data_dir: str = ""
    feature_path: str = ""
    drug_targets_path: str = ""
    disease_genes_path: str = ""
    subgraph_genes_path: str = ""
    gat_bridge_path: str = ""
    output_dir: str = ""

    @classmethod
    def create_default(cls) -> 'PathConfig':
        """创建默认路径配置 (基于原始 v3 硬编码路径)"""
        base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        return cls(
            data_dir=r"D:/反向网络药理学/GAT拓展维度/cache",
            feature_path=r"D:/反向网络药理学/GAT拓展维度/cache/enhanced_gene_features.csv",
            drug_targets_path=r"C:/Users/Jy-Mentor-7/Desktop/GAT/drug_targets.txt",
            disease_genes_path=r"C:/Users/Jy-Mentor-7/Desktop/GAT/disease_genes.txt",
            subgraph_genes_path=r"C:/Users/Jy-Mentor-7/Desktop/GAT/subgraph_genes.txt",
            gat_bridge_path=r"C:/Users/Jy-Mentor-7/Desktop/GAT/top20_bridge_genes.csv",
            output_dir=os.path.join(base_dir, "ml_output_v3"),
        )


@dataclass
class PipelineConfig:
    """管道全局配置"""
    gpu: GPUConfig = field(default_factory=GPUConfig)
    training: TrainingConfig = field(default_factory=TrainingConfig)
    paths: PathConfig = field(default_factory=PathConfig)

    @classmethod
    def create_default(cls) -> 'PipelineConfig':
        """创建默认配置"""
        return cls(
            gpu=GPUConfig(enabled=True),
            training=TrainingConfig(),
            paths=PathConfig.create_default(),
        )

    def ensure_output_dir(self):
        """确保输出目录存在"""
        os.makedirs(self.paths.output_dir, exist_ok=True)