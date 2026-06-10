# -*- coding: utf-8 -*-
"""配置加载器：从 YAML 读取所有超参数和路径，支持环境变量展开。"""

import os
import re
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional, Dict, Any, List

import yaml


def _resolve_env(val: str) -> str:
    """展开 $ENV 和 ${ENV} 环境变量"""
    if not isinstance(val, str):
        return val
    return os.path.expandvars(val)


def _resolve_refs(config: Dict[str, Any], root: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """递归解析 ${paths.xxx} 引用"""
    if root is None:
        root = config

    def _resolve_str(s: str) -> str:
        pattern = re.compile(r'\$\{([^}]+)\}')
        while True:
            m = pattern.search(s)
            if not m:
                break
            keys = m.group(1).split('.')
            val = root
            for k in keys:
                val = val[k] if isinstance(val, dict) else getattr(val, k, None)
                if val is None:
                    break
            if val is not None:
                s = s[:m.start()] + str(val) + s[m.end():]
            else:
                break
        return s

    for key, value in config.items():
        if isinstance(value, dict):
            _resolve_refs(value, root)
        elif isinstance(value, str):
            config[key] = _resolve_str(value)
        elif isinstance(value, list):
            config[key] = [_resolve_str(v) if isinstance(v, str) else v for v in value]
    return config


@dataclass
class PathsConfig:
    data_dir: str = ""
    gat_data_dir: str = ""
    cache_dir: str = ""
    project_dir: str = ""
    bridge_genes: str = ""
    enhanced_gene_features: str = ""
    drug_fingerprint: str = ""
    disease_features: str = ""
    pathway_features: str = ""
    ppi: str = ""
    coexp: str = ""
    tf_target: str = ""
    gene_pathway: str = ""
    subgraph_genes: str = ""
    methylation: str = ""
    mirna_target: Optional[str] = None
    reactome_gmt: str = ""
    trrust: str = ""
    pathway_list: str = ""
    pathway_pca_cache: str = ""
    model_save: str = ""
    folds_save_dir: str = ""
    log_dir: str = ""


@dataclass
class PreprocessingConfig:
    feature_dim: int = 256
    use_pca: bool = True
    use_umap: bool = False
    pca_intermediate_dim: int = 512
    normalize: str = "zscore"


@dataclass
class GraphConfig:
    subsample_ppi: bool = True
    ppi_max_edges: int = 80000
    ppi_score_threshold: int = 700
    ppi_split_by_score: bool = True
    ppi_strong_threshold: int = 900
    subsample_coexp: bool = True
    methylation_directed: bool = True
    add_pathway_hierarchy: bool = True
    add_mirna_edges: bool = False
    add_disease_pathway: bool = False


@dataclass
class ModelConfig:
    hidden_dim: int = 128
    num_heads: int = 4
    num_layers: int = 3
    dropout: float = 0.3
    initial_residual: bool = True
    drop_edge_p: float = 0.1
    decoder_bias: bool = True
    decoder_factorization: str = "distmult"
    use_input_bn: bool = True


@dataclass
class TrainingConfig:
    epochs: int = 200
    patience: int = 50
    lr: float = 0.001
    lr_patience: int = 20
    lr_factor: float = 0.5
    min_lr: float = 1e-6
    weight_decay: float = 0.0005
    grad_clip_norm: float = 1.0
    use_grad_scaling: bool = True
    eval_batch: int = 16384  # 评估批大小
    focal_alpha: float = 0.5
    focal_gamma: float = 2.0
    neg_sample_ratio: int = 3
    adaptive_neg_ratio: bool = True
    neg_sample_ratio_start: int = 1
    neg_sample_ratio_end: int = 10
    neg_sampling_mode: str = "degree"
    neg_degree_power: float = 0.75
    swa_epochs: int = 50
    early_stop_metric: str = "composite"
    composite_weight: float = 0.3  # 0.3*AUROC + 0.7*AUPRC
    val_ratio: float = 0.1
    val_neg_refresh_interval: int = 50
    use_minibatch: bool = False
    batch_size: int = 4096
    neighbor_sizes: List[int] = field(default_factory=lambda: [10, 10])
    num_workers: int = 0


@dataclass
class CVConfig:
    n_folds: int = 3
    save_fold_models: bool = True
    ensemble_inference: bool = True
    use_random_link_split: bool = False


@dataclass
class InferenceConfig:
    top_k: int = 10
    temperature: float = 1.5
    mc_dropout_samples: int = 10
    calibrate: bool = True
    eval_inductive: bool = True


@dataclass
class PretrainingConfig:
    enabled: bool = False
    epochs: int = 100
    lr: float = 0.001
    weight_decay: float = 0.0001
    mask_ratio: float = 0.3
    pretrain_save: str = ""


@dataclass
class LoggingConfig:
    level: str = "INFO"
    log_file: str = ""
    tensorboard: bool = True
    wandb: bool = False
    wandb_project: str = "hgt-pathway-discovery"
    wandb_entity: Optional[str] = None


@dataclass
class Config:
    paths: PathsConfig = field(default_factory=PathsConfig)
    device: str = "cuda"
    tf32: bool = True
    seed: int = 42
    preprocessing: PreprocessingConfig = field(default_factory=PreprocessingConfig)
    graph: GraphConfig = field(default_factory=GraphConfig)
    model: ModelConfig = field(default_factory=ModelConfig)
    training: TrainingConfig = field(default_factory=TrainingConfig)
    cv: CVConfig = field(default_factory=CVConfig)
    inference: InferenceConfig = field(default_factory=InferenceConfig)
    pretraining: PretrainingConfig = field(default_factory=PretrainingConfig)
    logging: LoggingConfig = field(default_factory=LoggingConfig)

    def __post_init__(self):
        if self.tf32 and self.device == "cuda":
            import torch
            torch.backends.cuda.matmul.allow_tf32 = True
            torch.backends.cudnn.allow_tf32 = True


def load_config(config_path: str = "config.yaml") -> Config:
    """从 YAML 文件加载配置并解析为 Config 数据类。"""
    base_dir = Path(config_path).parent

    if not Path(config_path).exists():
        config_path = str(base_dir / "hgt_pathway_discovery" / "config.yaml")

    with open(config_path, "r", encoding="utf-8") as f:
        raw = yaml.safe_load(f)

    raw = _resolve_refs(raw)
    raw = _resolve_env_recursive(raw)

    cfg = Config()

    if "paths" in raw:
        cfg.paths = PathsConfig(**raw["paths"])
    if "device" in raw:
        device_str = raw["device"]
        import torch
        if device_str == "cuda" and not torch.cuda.is_available():
            cfg.device = "cpu"
        else:
            cfg.device = device_str
    if "tf32" in raw:
        cfg.tf32 = raw["tf32"]
    if "seed" in raw:
        cfg.seed = raw["seed"]
    if "preprocessing" in raw:
        cfg.preprocessing = PreprocessingConfig(**raw["preprocessing"])
    if "graph" in raw:
        cfg.graph = GraphConfig(**raw["graph"])
    if "model" in raw:
        cfg.model = ModelConfig(**raw["model"])
    if "training" in raw:
        cfg.training = TrainingConfig(**raw["training"])
    if "cv" in raw:
        cfg.cv = CVConfig(**raw["cv"])
    if "inference" in raw:
        cfg.inference = InferenceConfig(**raw["inference"])
    if "pretraining" in raw:
        cfg.pretraining = PretrainingConfig(**raw["pretraining"])
    if "logging" in raw:
        cfg.logging = LoggingConfig(**raw["logging"])

    if cfg.tf32 and cfg.device == "cuda":
        import torch
        torch.backends.cuda.matmul.allow_tf32 = True
        torch.backends.cudnn.allow_tf32 = True

    return cfg


def _resolve_env_recursive(config: Dict[str, Any]) -> Dict[str, Any]:
    for key, value in config.items():
        if isinstance(value, dict):
            _resolve_env_recursive(value)
        elif isinstance(value, str):
            config[key] = _resolve_env(value)
    return config