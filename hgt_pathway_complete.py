# -*- coding: utf-8 -*-
"""
HGT (Heterogeneous Graph Transformer) 基因‑通路关联预测模型 v2.0
用于网络药理学第二阶段通路发现任务。

单文件完整版 — 包含所有模块: config / utils / data_loader / build_graph /
model / pretrain / train / inference / main

使用方法:
    python hgt_pathway_complete.py

参考文献:
  - HGT: Heterogeneous Graph Transformer, Hu et al., WWW 2020
  - GCNII: Simple and Deep GCN, Chen et al., ICML 2020
  - DropEdge: Rong et al., ICLR 2020
  - Feature Propagation for Missing Node Features, Rossi et al., NeurIPS 2021
  - Platt Scaling, Platt, 1999
  - Focal Loss, Lin et al., ICCV 2017
  - Deep Ensembles, Lakshminarayanan et al., NeurIPS 2017
  - Dropout as Bayesian Approximation, Gal & Ghahramani, ICML 2016
"""

import copy
import logging
import math
import os
import random
import re
import sys
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set, Any

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from torch import Tensor, nn
from torch.optim import AdamW
from torch.optim.lr_scheduler import ReduceLROnPlateau
from torch_geometric.data import HeteroData
from torch_geometric.nn import HGTConv, Linear
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import KFold

import yaml

warnings.filterwarnings("ignore")


# ============================================================================
# 内嵌 YAML 配置 (等价于 config.yaml)
# ============================================================================

_EMBEDDED_YAML = """
paths:
  data_dir: "D:/反向网络药理学/GAT拓展维度"
  gat_data_dir: "C:/Users/Jy-Mentor-7/Desktop/GAT"
  cache_dir: "${paths.data_dir}/cache"
  project_dir: "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙"
  bridge_genes: "${paths.project_dir}/data/bridge_genes.csv"
  enhanced_gene_features: "${paths.cache_dir}/enhanced_gene_features.csv"
  drug_fingerprint: "${paths.gat_data_dir}/drug_fingerprint.csv"
  disease_features: "${paths.data_dir}/disease_features.npy"
  pathway_features: "${paths.data_dir}/pathway_features.npy"
  ppi: "${paths.gat_data_dir}/ppi_subgraph.csv"
  coexp: "${paths.data_dir}/gene_coexp_edges.txt"
  tf_target: "${paths.data_dir}/tf_target_edges.txt"
  gene_pathway: "${paths.data_dir}/gene_pathway_edges.txt"
  subgraph_genes: "${paths.gat_data_dir}/subgraph_genes.txt"
  methylation: "${paths.data_dir}/gene_methylation_edges.txt"
  mirna_target: null
  reactome_gmt: "${paths.gat_data_dir}/ReactomePathways.gmt"
  reactome_relation: "${paths.data_dir}/ReactomePathwaysRelation.txt"
  trrust: "${paths.gat_data_dir}/TRRUST_Network_human.txt"
  complex2_pathway: "C:/Users/Jy-Mentor-7/Downloads/Complex_2_Pathway_human.txt"
  pathway_list: "${paths.cache_dir}/pathway_nodes.csv"
  pathway_pca_cache: "${paths.cache_dir}/pathway_features_pca256.npy"
  model_save: "${paths.project_dir}/hgt_pathway_model.pt"
  folds_save_dir: "${paths.project_dir}/folds"
  log_dir: "${paths.project_dir}/logs"

device: "cuda"
tf32: true
seed: 42

preprocessing:
  feature_dim: 256
  use_pca: true
  use_umap: false
  pca_intermediate_dim: 512
  normalize: "zscore"

graph:
  subsample_ppi: true
  ppi_max_edges: 80000
  ppi_score_threshold: 700
  ppi_split_by_score: true
  ppi_strong_threshold: 900
  subsample_coexp: true
  methylation_directed: true
  add_pathway_hierarchy: true
  add_mirna_edges: false
  add_disease_pathway: false

model:
  hidden_dim: 128
  num_heads: 4
  num_layers: 3
  dropout: 0.3
  initial_residual: true
  drop_edge_p: 0.1
  decoder_bias: true
  decoder_factorization: "distmult"
  use_input_bn: true

training:
  epochs: 200
  patience: 50
  lr: 0.001
  lr_patience: 20
  lr_factor: 0.5
  min_lr: 1.0e-6
  weight_decay: 0.0005
  grad_clip_norm: 1.0
  use_grad_scaling: true
  eval_batch: 16384
  focal_alpha: 0.5
  focal_gamma: 2.0
  neg_sample_ratio: 3
  adaptive_neg_ratio: true
  neg_sample_ratio_start: 1
  neg_sample_ratio_end: 10
  neg_sampling_mode: "degree"
  neg_degree_power: 0.75
  swa_epochs: 50
  early_stop_metric: "composite"
  composite_weight: 0.3
  val_ratio: 0.1
  val_neg_refresh_interval: 10
  use_minibatch: false
  batch_size: 4096
  neighbor_sizes: [10, 10]
  num_workers: 0

cv:
  n_folds: 3
  save_fold_models: true
  ensemble_inference: true
  use_random_link_split: false

inference:
  top_k: 10
  temperature: 1.5
  mc_dropout_samples: 0
  calibrate: true
  eval_inductive: true

pretraining:
  enabled: false
  epochs: 100
  lr: 0.001
  weight_decay: 0.0001
  mask_ratio: 0.3
  pretrain_save: "${paths.project_dir}/pretrain_model.pt"

logging:
  level: "INFO"
  log_file: "${paths.log_dir}/hgt_pipeline.log"
  tensorboard: true
  wandb: false
  wandb_project: "hgt-pathway-discovery"
  wandb_entity: null
"""


# ============================================================================
# config.py — 配置加载器
# ============================================================================

def _resolve_env(val: str) -> str:
    if not isinstance(val, str):
        return val
    return os.path.expandvars(val)


def _resolve_refs(config: Dict[str, Any], root: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
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


def _resolve_env_recursive(config: Dict[str, Any]) -> Dict[str, Any]:
    for key, value in config.items():
        if isinstance(value, dict):
            _resolve_env_recursive(value)
        elif isinstance(value, str):
            config[key] = _resolve_env(value)
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
    reactome_relation: str = ""
    trrust: str = ""
    complex2_pathway: str = ""
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
    eval_batch: int = 16384
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
    composite_weight: float = 0.3
    val_ratio: float = 0.1
    val_neg_refresh_interval: int = 10
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
            torch.backends.cuda.matmul.allow_tf32 = True
            torch.backends.cudnn.allow_tf32 = True


def load_config(config_path: Optional[str] = None) -> Config:
    if config_path and Path(config_path).exists():
        with open(config_path, "r", encoding="utf-8") as f:
            raw = yaml.safe_load(f)
    else:
        raw = yaml.safe_load(_EMBEDDED_YAML)

    raw = _resolve_refs(raw)
    raw = _resolve_env_recursive(raw)

    cfg = Config()

    if "paths" in raw:
        cfg.paths = PathsConfig(**raw["paths"])
    if "device" in raw:
        device_str = raw["device"]
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
        torch.backends.cuda.matmul.allow_tf32 = True
        torch.backends.cudnn.allow_tf32 = True

    return cfg


# ============================================================================
# utils.py — 工具函数
# ============================================================================

logger_utils = logging.getLogger("hgt_pipeline")


def setup_logging(name: str = "hgt_pipeline", level: str = "INFO",
                  log_file: Optional[str] = None,
                  tensorboard_dir: Optional[str] = None,
                  enable_tensorboard: bool = False) -> logging.Logger:
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
    try:
        from torch.utils.tensorboard import SummaryWriter
        writer = SummaryWriter(log_dir=log_dir)
        return writer
    except ImportError:
        return None


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


def zscore_normalize(arr: np.ndarray, eps: float = 1e-8) -> np.ndarray:
    mean = arr.mean(axis=0, keepdims=True)
    std = arr.std(axis=0, keepdims=True)
    return (arr - mean) / (std + eps)


def minmax_normalize(arr: np.ndarray) -> np.ndarray:
    vmin = arr.min(axis=0, keepdims=True)
    vmax = arr.max(axis=0, keepdims=True)
    denom = vmax - vmin
    denom[denom == 0] = 1.0
    return (arr - vmin) / denom


def pca_reduce(arr: np.ndarray, target_dim: int, seed: int = 42) -> Tuple[np.ndarray, Optional[PCA]]:
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
        logger_utils.warning("UMAP not installed, falling back to PCA")
        return pca_reduce(arr, umap_dim, seed)[0]


def normalize_features(arr: np.ndarray, method: str = "zscore") -> np.ndarray:
    if method == "zscore":
        return zscore_normalize(arr)
    elif method == "minmax":
        return minmax_normalize(arr)
    return arr


def try_compile(model: nn.Module) -> nn.Module:
    force_compile = os.environ.get("TORCH_COMPILE", "0") == "1"
    if not force_compile and not hasattr(torch, "compile"):
        return model
    if hasattr(torch, "compile") and torch.cuda.is_available():
        try:
            import triton
            return torch.compile(model, dynamic=True)
        except (ImportError, Exception) as e:
            logger_utils.warning(
                "torch.compile not available: %s. "
                "Set TORCH_COMPILE=1 and install triton to enable. "
                "Falling back to eager mode.",
                str(e) if str(e) else "Triton not installed",
            )
    return model


def setup_wandb(config: object, logger: logging.Logger) -> Optional[object]:
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


# ============================================================================
# data_loader.py — 数据加载
# ============================================================================


def load_txt(path: Path, skip_header_if_colname: bool = True) -> List[str]:
    if not path or not path.exists():
        return []
    for enc in ("utf-8", "gbk", "latin1"):
        try:
            with open(path, "r", encoding=enc) as f:
                lines = [ln.strip() for ln in f if ln.strip()]
            break
        except (UnicodeDecodeError, OSError):
            continue
    else:
        return []
    if skip_header_if_colname and lines and (
        lines[0].lower().startswith("gene") or lines[0].lower().startswith("symbol")
    ):
        lines = lines[1:]
    return list(dict.fromkeys(g.upper() for g in lines))


def load_edge_list(path: Optional[Path], sep: Optional[str] = None,
                   col_pair: Tuple[int, int] = (0, 1)) -> List[Tuple[str, str]]:
    if not path or not path.exists():
        return []
    if sep is None:
        sep = "\t" if path.suffix == ".txt" else ","
    raw: List[Tuple[str, str]] = []
    for enc in ("utf-8", "gbk", "latin1"):
        try:
            with open(path, "r", encoding=enc) as f:
                for ln in f:
                    ln = ln.strip()
                    if not ln or ln.startswith("#"):
                        continue
                    parts = ln.split(sep)
                    if len(parts) > max(col_pair):
                        raw.append((parts[col_pair[0]].upper().strip(),
                                    parts[col_pair[1]].strip()))
            break
        except (UnicodeDecodeError, OSError):
            continue
    return raw


def load_trrust_edges(path: Path) -> List[Tuple[str, str]]:
    """加载 TRRUST TF→Target 调控关系 (TF\tTarget\tRegulation\tPMID)"""
    if not path.exists():
        return []
    edges: List[Tuple[str, str]] = []
    for enc in ("utf-8", "latin1"):
        try:
            with open(path, "r", encoding=enc) as f:
                for ln in f:
                    ln = ln.strip()
                    if not ln or ln.startswith("#"):
                        continue
                    parts = ln.split("\t")
                    if len(parts) < 2:
                        continue
                    tf = parts[0].strip().upper()
                    target = parts[1].strip().upper()
                    if tf and target:
                        edges.append((tf, target))
        except (UnicodeDecodeError, OSError):
            continue
    edges = list(set(edges))
    print(f"[TRRUST] Loaded {len(edges)} TF→target edges from {path.name}")
    return edges


def load_gene_features(path: Path) -> Tuple[np.ndarray, List[str]]:
    df = pd.read_csv(path, index_col=0)
    df.index = df.index.astype(str).str.strip().str.upper()
    arr = df.apply(pd.to_numeric, errors="coerce").fillna(0.0).values.astype(np.float32)
    gene_names = df.index.tolist()
    return arr, gene_names


def load_drug_fingerprint(path: Path) -> np.ndarray:
    df = pd.read_csv(path)
    arr = df.apply(pd.to_numeric, errors="coerce").fillna(0).values.astype(np.float32)
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    return arr


def load_disease_features(path: Path) -> Optional[np.ndarray]:
    if not path.exists():
        return None
    arr = np.load(str(path))
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    elif arr.ndim == 2 and arr.shape[0] > 1:
        arr = np.mean(arr, axis=0, keepdims=True)
    return arr


def load_pathway_features(path: Path) -> np.ndarray:
    arr = np.load(str(path))
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    return arr


def load_ppi(path: Path, score_thresh: int = 700,
             bridge_set: Optional[Set[str]] = None,
             max_edges: int = 80000,
             subsample: bool = True) -> List[Tuple[str, str, float]]:
    df = pd.read_csv(path, sep="\t", encoding="utf-8")
    col_a, col_b = df.columns[0], df.columns[1]
    score_col = df.columns[2] if len(df.columns) >= 3 else None

    df = df.dropna(subset=[col_a, col_b])
    df[col_a] = df[col_a].astype(str).str.strip().str.upper()
    df[col_b] = df[col_b].astype(str).str.strip().str.upper()

    if score_col is not None:
        df[score_col] = pd.to_numeric(df[score_col], errors="coerce")
        df = df[df[score_col] >= score_thresh]
        df = df.sort_values(score_col, ascending=False)

    if subsample and bridge_set and len(bridge_set) > 0:
        mask = df[col_a].isin(bridge_set) | df[col_b].isin(bridge_set)
        df_bridge = df[mask]
        df_other = df[~mask]
        if len(df_bridge) > max_edges:
            df_bridge = df_bridge.head(max_edges)
        n_left = max(0, max_edges - len(df_bridge))
        df_other = df_other.head(n_left) if n_left > 0 else df_other.iloc[:0]
        df = pd.concat([df_bridge, df_other], ignore_index=True)
    elif subsample and len(df) > max_edges:
        df = df.head(max_edges)

    if score_col is not None:
        edges = list(zip(df[col_a], df[col_b], df[score_col]))
    else:
        edges = [(a, b, 0.0) for a, b in zip(df[col_a], df[col_b])]
    return edges


def load_pathway_list(path: Path) -> List[str]:
    for enc in ("utf-8", "gbk", "latin1"):
        try:
            df = pd.read_csv(path, encoding=enc)
            break
        except UnicodeDecodeError:
            continue
    name_col = "pathway_name" if "pathway_name" in df.columns else df.columns[0]
    return df[name_col].tolist()


def load_gene_pathway_edges(path: Path) -> List[Tuple[str, str]]:
    edges: List[Tuple[str, str]] = []
    for enc in ("utf-8", "gbk", "latin1"):
        try:
            with open(path, "r", encoding=enc) as f:
                for ln in f:
                    ln = ln.strip()
                    if not ln:
                        continue
                    parts = ln.split("\t")
                    if len(parts) >= 2:
                        edges.append((parts[0].upper(), parts[1].strip()))
            break
        except (UnicodeDecodeError, OSError):
            continue
    return edges


def parse_reactome_relation_file(relation_path: Path) -> List[Tuple[str, str]]:
    if not relation_path.exists():
        return []
    hierarchy: List[Tuple[str, str]] = []
    with open(relation_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            child_id, parent_id = parts[0].strip(), parts[1].strip()
            hierarchy.append((parent_id, child_id))
    hierarchy = list(set(hierarchy))
    return hierarchy


def parse_reactome_hierarchy(gmt_path: Path) -> List[Tuple[str, str]]:
    if not gmt_path.exists():
        return []

    hierarchy: List[Tuple[str, str]] = []
    pathway_data: Dict[str, Tuple[str, str, Set[str]]] = {}

    with open(gmt_path, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            name = parts[0].strip()
            reactome_id = parts[1].strip() if len(parts) > 1 else ""
            genes = set(g.strip().upper() for g in parts[2:] if g.strip())
            pathway_data[name] = (reactome_id, name, genes)

    pathway_names = sorted(pathway_data.keys())

    for i, name_a in enumerate(pathway_names):
        genes_a = pathway_data[name_a][2]
        if len(genes_a) < 3:
            continue
        for name_b in pathway_names[i + 1:]:
            genes_b = pathway_data[name_b][2]
            if len(genes_b) < 3:
                continue
            overlap = len(genes_a & genes_b)
            if overlap == 0:
                continue
            ratio_a = overlap / len(genes_a)
            ratio_b = overlap / len(genes_b)
            if ratio_a > 0.9 and ratio_b < 0.9:
                hierarchy.append((name_b, name_a))
            elif ratio_b > 0.9 and ratio_a < 0.9:
                hierarchy.append((name_a, name_b))

    hierarchy = list(set(hierarchy))
    return hierarchy


def load_complex_2_pathway_hierarchy(gmt_path: Path,
                                      complex_path: Path) -> Tuple[List[Tuple[str, str]], Dict[str, str]]:
    """利用 GMT (pathway_name↔R-HSA) + Complex_2_Pathway_human.txt (R-HSA hierarchy) 构建人源通路层级.

    Returns:
        hierarchy: List[(parent_name, child_name)] — 通路名称级层级边
        r_hsa_to_name: R-HSA-ID → pathway_name 映射字典
    """
    # 1. 加载 GMT → 建立 R-HSA ↔ name 映射
    r_hsa_to_name: Dict[str, str] = {}
    if gmt_path.exists():
        with open(gmt_path, "r", encoding="utf-8") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 2:
                    continue
                name = parts[0].strip()
                rhsa = parts[1].strip()
                if rhsa.startswith("R-HSA-"):
                    r_hsa_to_name[rhsa] = name
        print(f"[Complex2Pathway] GMT loaded: {len(r_hsa_to_name)} R-HSA→name mappings")

    # 2. 加载 Complex_2_Pathway → 提取 (top_level_pathway, pathway) 层级
    raw_hierarchy: List[Tuple[str, str]] = []
    if complex_path.exists():
        with open(complex_path, "r", encoding="utf-8") as f:
            header = next(f, "").strip().lower()
            has_header = "complex" in header and "pathway" in header
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                pathway_id = parts[1].strip()
                top_level_id = parts[2].strip()
                if pathway_id and top_level_id:
                    raw_hierarchy.append((top_level_id, pathway_id))
        raw_hierarchy = list(set(raw_hierarchy))
        print(f"[Complex2Pathway] Raw R-HSA hierarchy pairs: {len(raw_hierarchy)}")

    # 3. 转换为通路名称级
    name_hierarchy: List[Tuple[str, str]] = []
    for parent_id, child_id in raw_hierarchy:
        parent_name = r_hsa_to_name.get(parent_id)
        child_name = r_hsa_to_name.get(child_id)
        if parent_name and child_name and parent_name != child_name:
            name_hierarchy.append((parent_name, child_name))

    name_hierarchy = list(set(name_hierarchy))
    print(f"[Complex2Pathway] Name-level hierarchy edges: {len(name_hierarchy)}")
    return name_hierarchy, r_hsa_to_name


def resolve_pathway_hierarchy(config_paths: object,
                               valid_names: Optional[Set[str]] = None,
                               complex2_hierarchy: Optional[List[Tuple[str, str]]] = None) -> List[Tuple[str, str]]:
    """解析通路层级关系.

    优先级:
        1. Complex_2_Pathway_human.txt (R-HSA 官方层级, 需 GMT 映射)
        2. ReactomePathwaysRelation.txt (R-BTA ID, 仅当 ID 匹配通路名时)
        3. GMT 基因重叠启发式推断

    Returns:
        List[(parent_name, child_name)]
    """
    # 第 1 优先: Complex_2_Pathway 精确人源层级
    if complex2_hierarchy:
        if valid_names:
            matched = sum(1 for p, c in complex2_hierarchy
                          if p in valid_names and c in valid_names)
            if matched > 0:
                print(f"[Hierarchy] Complex_2_Pathway: {len(complex2_hierarchy)} edges, "
                      f"{matched} matched pathway names")
                # 合并 GMT 启发式层级 (补充 deeper hierarchy)
                gmt_path = Path(config_paths.reactome_gmt) if config_paths.reactome_gmt else None
                if gmt_path and gmt_path.exists():
                    gmt_edges = parse_reactome_hierarchy(gmt_path)
                    merged = list(set(complex2_hierarchy) | set(gmt_edges))
                    print(f"[Hierarchy] GMT heuristic: {len(gmt_edges)} edges, "
                          f"merged total: {len(merged)}")
                    return merged
                return complex2_hierarchy
            print(f"[Hierarchy] Complex_2_Pathway: {len(complex2_hierarchy)} edges, "
                  f"NONE matched pathway names")

    # 第 2 优先: ReactomePathwaysRelation (R-BTA ID)
    rel_path = Path(config_paths.reactome_relation) if config_paths.reactome_relation else None
    if rel_path and rel_path.exists():
        edges = parse_reactome_relation_file(rel_path)
        if edges:
            if valid_names:
                matched = sum(1 for p, c in edges if p in valid_names or c in valid_names)
                if matched > 0:
                    print(f"[Hierarchy] Relation file: {len(edges)} edges, {matched} matched")
                    return edges
                print(f"[Hierarchy] Relation file: {len(edges)} edges, NONE matched "
                      f"(ID mismatch, e.g. R-BTA vs R-HSA), falling back to GMT")
            else:
                return edges

    # 第 3 优先: GMT 启发式推断
    gmt_path = Path(config_paths.reactome_gmt) if config_paths.reactome_gmt else None
    if gmt_path and gmt_path.exists():
        gmt_edges = parse_reactome_hierarchy(gmt_path)
        print(f"[Hierarchy] GMT heuristic: {len(gmt_edges)} edges")
        return gmt_edges

    return []


def load_bridge_genes(path: Path, alt_path: Optional[Path] = None) -> List[str]:
    if path.exists():
        df = pd.read_csv(path)
        col = "gene_symbol" if "gene_symbol" in df.columns else df.columns[0]
        genes = df[col].astype(str).str.strip().str.upper().tolist()
        return genes
    if alt_path and alt_path.exists():
        df = pd.read_csv(alt_path)
        genes = df["gene_symbol"].astype(str).str.strip().str.upper().tolist()
        return genes
    return []


def load_mirna_edges(path: Optional[Path]) -> List[Tuple[str, str]]:
    if not path or not path.exists():
        return []
    return load_edge_list(path, sep="\t", col_pair=(0, 1))


def load_disease_pathway_edges(path: Optional[Path]) -> List[Tuple[str, str]]:
    if not path or not path.exists():
        return []
    return load_edge_list(path, sep="\t", col_pair=(0, 1))


def load_all_data(config: object) -> Dict:
    result = {}
    p = config.paths
    graph_cfg = config.graph

    result["gene_feat_arr"], result["gene_feat_names"] = load_gene_features(
        Path(p.enhanced_gene_features))

    result["pathway_feat_arr"] = load_pathway_features(Path(p.pathway_features))

    result["drug_fp_arr"] = load_drug_fingerprint(Path(p.drug_fingerprint))

    result["disease_feat_arr"] = load_disease_features(Path(p.disease_features))

    result["pathway_names"] = load_pathway_list(Path(p.pathway_list))

    result["bridge_genes"] = load_bridge_genes(
        Path(p.bridge_genes), alt_path=Path(p.gat_data_dir) / "all_bridge_genes.csv")

    bridge_set = set(result["bridge_genes"]) if result["bridge_genes"] else None
    result["ppi_edges"] = load_ppi(
        Path(p.ppi),
        score_thresh=graph_cfg.ppi_score_threshold,
        bridge_set=bridge_set,
        max_edges=graph_cfg.ppi_max_edges,
        subsample=graph_cfg.subsample_ppi,
    )

    result["coexp_edges"] = load_edge_list(Path(p.coexp), sep="\t")
    result["tf_edges"] = load_edge_list(Path(p.tf_target), sep="\t")

    # TRRUST 转录因子调控增强
    trrust_path = Path(p.trrust) if p.trrust else None
    if trrust_path and trrust_path.exists():
        trrust_edges = load_trrust_edges(trrust_path)
        existing_tf_set = set(result["tf_edges"])
        merged = list(existing_tf_set | set(trrust_edges))
        print(f"[TRRUST] Merged {len(trrust_edges)} TRRUST edges with {len(result['tf_edges'])} "
              f"existing TF edges → {len(merged)} total")
        result["tf_edges"] = merged
    else:
        print(f"[TRRUST] File not found at {trrust_path}, skipping")
    result["gene_pathway_edges"] = load_gene_pathway_edges(Path(p.gene_pathway))
    result["all_genes_list"] = load_txt(Path(p.subgraph_genes))

    meth_path = Path(p.methylation)
    result["methyl_edges"] = load_edge_list(meth_path, sep=",", col_pair=(0, 1)) if meth_path.exists() else None

    mirna_path = p.mirna_target
    if mirna_path:
        result["mirna_edges"] = load_mirna_edges(Path(mirna_path))
    else:
        result["mirna_edges"] = None

    result["pathway_reactome_ids"] = load_pathway_list(Path(p.pathway_list))

    # Complex_2_Pathway_human.txt: 基于 R-HSA ID 的精确人源通路层级
    gmt_path = Path(p.reactome_gmt) if p.reactome_gmt else None
    complex_path = Path(p.complex2_pathway) if p.complex2_pathway else None
    complex2_hierarchy: Optional[List[Tuple[str, str]]] = None
    if gmt_path and gmt_path.exists() and complex_path and complex_path.exists():
        complex2_hierarchy, _ = load_complex_2_pathway_hierarchy(gmt_path, complex_path)

    pathway_names_set = set(result["pathway_reactome_ids"])
    result["pathway_hierarchy"] = resolve_pathway_hierarchy(
        p, valid_names=pathway_names_set,
        complex2_hierarchy=complex2_hierarchy,
    )

    dp_path = Path(p.data_dir) / "disease_pathway_edges.txt"
    result["disease_pathway_edges"] = load_disease_pathway_edges(dp_path) if dp_path.exists() else []

    return result


# ============================================================================
# build_graph.py — 异构图构建
# ============================================================================


def normalize_node_features(data: HeteroData, method: str = "zscore",
                            eps: float = 1e-8) -> HeteroData:
    for nt in data.node_types:
        x = data[nt].x
        if x.size(0) == 0:
            continue
        if method == "zscore":
            mean = x.mean(dim=0)
            std = x.std(dim=0, unbiased=False)
            std[std < eps] = 1.0
            data[nt].x = (x - mean) / std
            data[nt].norm_stats = {"mean": mean, "std": std, "method": "zscore"}
        elif method == "minmax":
            vmin = x.min(dim=0)[0]
            vmax = x.max(dim=0)[0]
            denom = vmax - vmin
            denom[denom < eps] = 1.0
            data[nt].x = (x - vmin) / denom
            data[nt].norm_stats = {"min": vmin, "max": vmax, "method": "minmax"}
        else:
            data[nt].norm_stats = {"method": "none"}
    return data


def apply_cached_normalization(data: HeteroData) -> HeteroData:
    for nt in data.node_types:
        if hasattr(data[nt], "norm_stats") and data[nt].norm_stats.get("method") == "zscore":
            stats = data[nt].norm_stats
            data[nt].x = (data[nt].x - stats["mean"]) / stats["std"]
        elif hasattr(data[nt], "norm_stats") and data[nt].norm_stats.get("method") == "minmax":
            stats = data[nt].norm_stats
            data[nt].x = (data[nt].x - stats["min"]) / (stats["max"] - stats["min"])
    return data


def build_hetero_graph(
    gene_feat_arr: np.ndarray,
    gene_feat_names: List[str],
    drug_fp_arr: np.ndarray,
    disease_feat_arr: Optional[np.ndarray],
    pathway_feat_arr: np.ndarray,
    pathway_names: List[str],
    ppi_edges: List[Tuple[str, str, float]],
    coexp_edges: List[Tuple[str, str]],
    tf_edges: List[Tuple[str, str]],
    gene_pathway_edges: List[Tuple[str, str]],
    all_genes_list: Optional[List[str]] = None,
    bridge_genes: Optional[List[str]] = None,
    methyl_edges: Optional[List[Tuple[str, str]]] = None,
    mirna_edges: Optional[List[Tuple[str, str]]] = None,
    pathway_hierarchy: Optional[List[Tuple[str, str]]] = None,
    disease_pathway_edges: Optional[List[Tuple[str, str]]] = None,
    config: Optional[object] = None,
) -> Tuple[HeteroData, Dict[str, int], List[str], Dict[str, int]]:
    graph_cfg = config.graph if config else None
    preproc_cfg = config.preprocessing if config else None

    gene_set = set(gene_feat_names)
    if all_genes_list:
        gene_set &= set(all_genes_list)

    for a, b, _ in ppi_edges:
        gene_set.add(a); gene_set.add(b)
    for a, b in coexp_edges:
        gene_set.add(a); gene_set.add(b)
    for a, b in tf_edges:
        gene_set.add(a); gene_set.add(b)
    for a, b in gene_pathway_edges:
        gene_set.add(a)
    if methyl_edges:
        for a, b in methyl_edges:
            gene_set.add(a)
    if mirna_edges:
        for a, b in mirna_edges:
            gene_set.add(a); gene_set.add(b)

    gene_list = sorted(gene_set & set(gene_feat_names))
    if not gene_list:
        raise ValueError("No genes with features in the graph!")

    gene_to_idx = {g: i for i, g in enumerate(gene_list)}
    n_genes = len(gene_list)

    gene_feat_dict = dict(zip(gene_feat_names, gene_feat_arr))
    gene_feat = np.zeros((n_genes, gene_feat_arr.shape[1]), dtype=np.float32)
    for i, g in enumerate(gene_list):
        if g in gene_feat_dict:
            gene_feat[i] = gene_feat_dict[g]

    pathway_feat = pathway_feat_arr.astype(np.float32).copy()
    if pathway_feat.ndim == 1:
        pathway_feat = pathway_feat.reshape(1, -1)

    data = HeteroData()
    data["gene"].x = torch.from_numpy(gene_feat)
    data["drug"].x = torch.from_numpy(drug_fp_arr.reshape(1, -1).astype(np.float32))

    if disease_feat_arr is not None:
        dis_feat = disease_feat_arr.reshape(1, -1)
    else:
        dis_feat = np.mean(gene_feat, axis=0, keepdims=True)
    data["disease"].x = torch.from_numpy(dis_feat.astype(np.float32))
    data["pathway"].x = torch.from_numpy(pathway_feat)

    normalize_method = preproc_cfg.normalize if preproc_cfg else "zscore"
    if normalize_method != "none":
        data = normalize_node_features(data, normalize_method)
        print(f"[Build] Node features normalized: method={normalize_method}")

    pathway_name_to_idx = {name: i for i, name in enumerate(pathway_names)}

    ppi_strong_threshold = graph_cfg.ppi_strong_threshold if graph_cfg else 900
    ppi_split = graph_cfg.ppi_split_by_score if graph_cfg else True

    if ppi_split:
        strong_src, strong_dst = [], []
        weak_src, weak_dst = [], []
        for a, b, score in ppi_edges:
            if a not in gene_to_idx or b not in gene_to_idx:
                continue
            if score >= ppi_strong_threshold:
                strong_src.extend([gene_to_idx[a], gene_to_idx[b]])
                strong_dst.extend([gene_to_idx[b], gene_to_idx[a]])
            else:
                weak_src.extend([gene_to_idx[a], gene_to_idx[b]])
                weak_dst.extend([gene_to_idx[b], gene_to_idx[a]])

        data["gene", "strong_ppi", "gene"].edge_index = (
            torch.tensor([strong_src, strong_dst], dtype=torch.long)
            if strong_src else torch.zeros((2, 0), dtype=torch.long)
        )
        data["gene", "weak_ppi", "gene"].edge_index = (
            torch.tensor([weak_src, weak_dst], dtype=torch.long)
            if weak_src else torch.zeros((2, 0), dtype=torch.long)
        )
        print(f"[Build] strong_ppi: {len(strong_src)}, weak_ppi: {len(weak_src)}")
    else:
        ppi_src, ppi_dst = [], []
        for a, b, _ in ppi_edges:
            if a in gene_to_idx and b in gene_to_idx:
                ppi_src.extend([gene_to_idx[a], gene_to_idx[b]])
                ppi_dst.extend([gene_to_idx[b], gene_to_idx[a]])
        data["gene", "interacts", "gene"].edge_index = (
            torch.tensor([ppi_src, ppi_dst], dtype=torch.long)
            if ppi_src else torch.zeros((2, 0), dtype=torch.long)
        )
        print(f"[Build] PPI (unified): {len(ppi_src)}")

    if graph_cfg and graph_cfg.subsample_coexp and bridge_genes:
        bridge_set = set(bridge_genes)
        coe_filt = [(a, b) for a, b in coexp_edges if a in bridge_set or b in bridge_set]
        coe_used = coe_filt if coe_filt else coexp_edges
    else:
        coe_used = coexp_edges

    coe_src, coe_dst = [], []
    for a, b in coe_used:
        if a in gene_to_idx and b in gene_to_idx:
            coe_src.extend([gene_to_idx[a], gene_to_idx[b]])
            coe_dst.extend([gene_to_idx[b], gene_to_idx[a]])
    data["gene", "coexpressed", "gene"].edge_index = (
        torch.tensor([coe_src, coe_dst], dtype=torch.long)
        if coe_src else torch.zeros((2, 0), dtype=torch.long)
    )
    print(f"[Build] Coexp edges: {len(coe_src)} (used={len(coe_used)})")

    tf_src, tf_dst = [], []
    for a, b in tf_edges:
        if a in gene_to_idx and b in gene_to_idx:
            tf_src.append(gene_to_idx[a])
            tf_dst.append(gene_to_idx[b])
    data["gene", "regulates", "gene"].edge_index = (
        torch.tensor([tf_src, tf_dst], dtype=torch.long)
        if tf_src else torch.zeros((2, 0), dtype=torch.long)
    )
    print(f"[Build] TF edges: {len(tf_src)}")

    drug_targets = load_txt(Path(config.paths.gat_data_dir) / "drug_targets.txt")
    dt_src, dt_dst = [], []
    for g in drug_targets:
        if g in gene_to_idx:
            dt_src.append(0)
            dt_dst.append(gene_to_idx[g])
    data["drug", "targets", "gene"].edge_index = (
        torch.tensor([dt_src, dt_dst], dtype=torch.long)
        if dt_src else torch.zeros((2, 0), dtype=torch.long)
    )
    print(f"[Build] drug->gene edges: {len(dt_src)}")

    disease_genes = load_txt(Path(config.paths.gat_data_dir) / "disease_genes.txt")
    gd_src, gd_dst = [], []
    for g in disease_genes:
        if g in gene_to_idx:
            gd_src.append(gene_to_idx[g])
            gd_dst.append(0)
    data["gene", "assoc_with", "disease"].edge_index = (
        torch.tensor([gd_src, gd_dst], dtype=torch.long)
        if gd_src else torch.zeros((2, 0), dtype=torch.long)
    )
    print(f"[Build] gene->disease edges: {len(gd_src)}")

    gp_src, gp_dst = [], []
    for a, b in gene_pathway_edges:
        if a in gene_to_idx and b in pathway_name_to_idx:
            gp_src.append(gene_to_idx[a])
            gp_dst.append(pathway_name_to_idx[b])
    data["gene", "involved_in", "pathway"].edge_index = (
        torch.tensor([gp_src, gp_dst], dtype=torch.long)
        if gp_src else torch.zeros((2, 0), dtype=torch.long)
    )
    n_gp = len(gp_src)
    print(f"[Build] gene->pathway edges: {n_gp} (positive samples)")

    if pathway_hierarchy:
        ph_src, ph_dst = [], []
        for parent, child in pathway_hierarchy:
            if parent in pathway_name_to_idx and child in pathway_name_to_idx:
                ph_src.append(pathway_name_to_idx[parent])
                ph_dst.append(pathway_name_to_idx[child])
        data["pathway", "parent_of", "pathway"].edge_index = (
            torch.tensor([ph_src, ph_dst], dtype=torch.long)
            if ph_src else torch.zeros((2, 0), dtype=torch.long)
        )
        print(f"[Build] pathway hierarchy edges: {len(ph_src)}")

    if disease_pathway_edges:
        dp_src, dp_dst = [], []
        for disease_name, pathway_name in disease_pathway_edges:
            if pathway_name in pathway_name_to_idx:
                dp_src.append(0)
                dp_dst.append(pathway_name_to_idx[pathway_name])
        if dp_src:
            data["disease", "assoc_with", "pathway"].edge_index = torch.tensor(
                [dp_src, dp_dst], dtype=torch.long
            )
            print(f"[Build] disease->pathway edges: {len(dp_src)}")

    methylation_directed = graph_cfg.methylation_directed if graph_cfg else True

    if methyl_edges:
        cpg_set: Set[str] = set()
        for g, cpg in methyl_edges:
            if g in gene_to_idx:
                cpg_set.add(cpg)
        cpg_list = sorted(cpg_set)
        cpg_to_idx = {c: i for i, c in enumerate(cpg_list)}

        if cpg_list:
            cpg_feat = np.zeros((len(cpg_list), gene_feat.shape[1]), dtype=np.float32)
            cpg_count = np.zeros(len(cpg_list), dtype=np.float32)
            for g, cpg in methyl_edges:
                if g in gene_to_idx and cpg in cpg_to_idx:
                    cpg_feat[cpg_to_idx[cpg]] += gene_feat[gene_to_idx[g]]
                    cpg_count[cpg_to_idx[cpg]] += 1.0

            mask = cpg_count > 0
            cpg_feat[mask] /= cpg_count[mask, np.newaxis]
            cpg_feat[~mask] = gene_feat.mean(axis=0)

            quality_mask = np.zeros((len(cpg_list), 1), dtype=np.float32)
            quality_mask[mask] = 1.0
            cpg_feat_with_quality = np.concatenate([cpg_feat, quality_mask], axis=1)

            data["cpg"].x = torch.from_numpy(cpg_feat_with_quality).float()
            data["cpg"].propagation_mask = torch.from_numpy(mask)

            gm_src, gm_dst = [], []
            for g, cpg in methyl_edges:
                if g in gene_to_idx and cpg in cpg_to_idx:
                    gm_src.append(gene_to_idx[g])
                    gm_dst.append(cpg_to_idx[cpg])

            if methylation_directed:
                data["gene", "methylated_at", "cpg"].edge_index = torch.tensor(
                    [gm_src, gm_dst], dtype=torch.long
                )
            else:
                data["gene", "methylated_at", "cpg"].edge_index = torch.tensor(
                    [gm_src + gm_dst, gm_dst + gm_src], dtype=torch.long
                )

            n_isolated = int((cpg_count == 0).sum())
            edge_count = data["gene", "methylated_at", "cpg"].edge_index.size(1)
            print(f"[Build] Methylation edges: {edge_count} "
                  f"({'directed' if methylation_directed else 'undirected'}), "
                  f"CpG: {len(cpg_list)} (feature_propagation, isolated={n_isolated}, "
                  f"mean_degree={cpg_count.mean():.1f})")

    if mirna_edges:
        mirna_set: Set[str] = set()
        for m, g in mirna_edges:
            if g in gene_to_idx:
                mirna_set.add(m)
        mirna_list = sorted(mirna_set)
        mirna_to_idx = {m: i for i, m in enumerate(mirna_list)}

        if mirna_list:
            mirna_feat = np.zeros((len(mirna_list), gene_feat.shape[1]), dtype=np.float32)
            mirna_count = np.zeros(len(mirna_list), dtype=np.float32)
            for m, g in mirna_edges:
                if g in gene_to_idx and m in mirna_to_idx:
                    mirna_feat[mirna_to_idx[m]] += gene_feat[gene_to_idx[g]]
                    mirna_count[mirna_to_idx[m]] += 1.0
            mask_m = mirna_count > 0
            mirna_feat[mask_m] /= mirna_count[mask_m, np.newaxis]
            mirna_feat[~mask_m] = gene_feat.mean(axis=0)
            data["mirna"].x = torch.from_numpy(mirna_feat).float()

            mr_src, mr_dst = [], []
            for m, g in mirna_edges:
                if g in gene_to_idx and m in mirna_to_idx:
                    mr_src.append(mirna_to_idx[m])
                    mr_dst.append(gene_to_idx[g])
            data["mirna", "targets", "gene"].edge_index = (
                torch.tensor([mr_src, mr_dst], dtype=torch.long)
                if mr_src else torch.zeros((2, 0), dtype=torch.long)
            )
            print(f"[Build] miRNA->gene edges: {len(mr_src)}, miRNA nodes: {len(mirna_list)}")

    print(f"[Build] Node types: {list(data.node_types)}")
    print(f"[Build] Edge types: {list(data.edge_types)}")
    return data, gene_to_idx, gene_list, pathway_name_to_idx


# ============================================================================
# model.py — HGT 编码器-解码器模型
# ============================================================================


class HGTEncoder(nn.Module):
    def __init__(self, metadata: Tuple[List[str], List[Tuple[str, str, str]]],
                 hidden_dim: int, num_heads: int, num_layers: int,
                 dropout: float, initial_residual: bool = True,
                 drop_edge_p: float = 0.0):
        super().__init__()
        self.metadata = metadata
        self.num_layers = num_layers
        self.initial_residual = initial_residual
        self.drop_edge_p = drop_edge_p

        self.convs = nn.ModuleList()
        self.norms = nn.ModuleList()
        for _ in range(num_layers):
            self.convs.append(
                HGTConv(in_channels=hidden_dim, out_channels=hidden_dim,
                        metadata=metadata, heads=num_heads)
            )
            self.norms.append(nn.LayerNorm(hidden_dim))

        if initial_residual and num_layers > 1:
            self.skip_alphas = nn.ParameterList([
                nn.Parameter(torch.tensor(0.5)) for _ in range(num_layers)
            ])
        else:
            self.skip_alphas = None

        self.dropout = dropout

    @staticmethod
    def _drop_edges(edge_index_dict: Dict, drop_p: float) -> Dict:
        if drop_p <= 0:
            return edge_index_dict
        dropped = {}
        for et, ei in edge_index_dict.items():
            n_edges = ei.size(1)
            if n_edges <= 1:
                dropped[et] = ei
                continue
            keep_mask = torch.rand(n_edges, device=ei.device) > drop_p
            if keep_mask.sum() == 0:
                keep_mask[0] = True
            dropped[et] = ei[:, keep_mask]
        return dropped

    def forward(self, x_dict: Dict[str, Tensor],
                edge_index_dict: Dict,
                x0_dict: Optional[Dict[str, Tensor]] = None) -> Dict[str, Tensor]:
        for layer_idx, (conv, norm) in enumerate(zip(self.convs, self.norms)):
            if self.training and self.drop_edge_p > 0:
                ei_dict = self._drop_edges(edge_index_dict, self.drop_edge_p)
            else:
                ei_dict = edge_index_dict

            x_dict_new = conv(x_dict, ei_dict)
            for k in x_dict:
                if k not in x_dict_new:
                    x_dict_new[k] = x_dict[k]

            x_dict = {}
            for k in x_dict_new:
                residual = x_dict.get(k, x_dict_new[k])
                x_dict[k] = F.relu(norm(x_dict_new[k] + residual))

            if self.initial_residual and x0_dict is not None and self.skip_alphas is not None:
                alpha = torch.sigmoid(self.skip_alphas[layer_idx])
                for k in x_dict:
                    if k in x0_dict:
                        x_dict[k] = (1 - alpha) * x_dict[k] + alpha * x0_dict[k]

            x_dict = {
                k: F.dropout(v, p=self.dropout, training=self.training)
                for k, v in x_dict.items()
            }

        return x_dict


class DecoupledDecoder(nn.Module):
    def __init__(self, hidden_dim: int, dropout: float,
                 factorization: str = "distmult"):
        super().__init__()
        self.factorization = factorization

        self.gene_bias_proj = nn.Sequential(
            Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            Linear(hidden_dim, hidden_dim),
        )
        self.path_bias_proj = nn.Sequential(
            Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            Linear(hidden_dim, hidden_dim),
        )

        if factorization == "distmult":
            self.W = nn.Parameter(torch.randn(hidden_dim) * 0.1)
            self.output = nn.Linear(1, 1)
        else:
            self.decoder_mlp = nn.Sequential(
                Linear(hidden_dim * 2, hidden_dim),
                nn.BatchNorm1d(hidden_dim),
                nn.ReLU(),
                nn.Dropout(dropout),
                Linear(hidden_dim, hidden_dim // 2),
                nn.BatchNorm1d(hidden_dim // 2),
                nn.ReLU(),
                nn.Dropout(dropout * 0.5),
                Linear(hidden_dim // 2, 1),
            )

    def forward(self, z_gene: Tensor, z_path: Tensor) -> Tensor:
        gene_bias = self.gene_bias_proj(z_gene)
        path_bias = self.path_bias_proj(z_path)
        z_gene_biased = z_gene + 0.1 * gene_bias
        z_path_biased = z_path + 0.1 * path_bias

        if self.factorization == "distmult":
            scores = torch.sum(z_gene_biased * self.W * z_path_biased, dim=-1, keepdim=True)
            return self.output(scores).squeeze(-1)
        else:
            h = torch.cat([z_gene_biased, z_path_biased], dim=-1)
            return self.decoder_mlp(h).squeeze(-1)


class HGTModel(nn.Module):
    def __init__(self, metadata: Tuple[List[str], List[Tuple[str, str, str]]],
                 dim_dict: Dict[str, int], hidden_dim: int,
                 num_heads: int, num_layers: int,
                 dropout: float, initial_residual: bool = True,
                 decoder_bias: bool = True,
                 decoder_factorization: str = "distmult",
                 use_input_bn: bool = True,
                 drop_edge_p: float = 0.0):
        super().__init__()
        self.node_types = metadata[0]
        self.edge_types = metadata[1]
        self.hidden_dim = hidden_dim
        self.initial_residual = initial_residual
        self.use_input_bn = use_input_bn
        self.drop_edge_p = drop_edge_p

        self.proj = nn.ModuleDict()
        for nt, d_in in dim_dict.items():
            self.proj[nt] = Linear(d_in, hidden_dim)

        if use_input_bn:
            self.input_bn = nn.ModuleDict()
            for nt in dim_dict:
                self.input_bn[nt] = nn.BatchNorm1d(hidden_dim)
        else:
            self.input_bn = None

        self.encoder = HGTEncoder(
            metadata, hidden_dim, num_heads, num_layers,
            dropout, initial_residual=initial_residual,
            drop_edge_p=drop_edge_p,
        )

        self.decoder = DecoupledDecoder(
            hidden_dim, dropout, factorization=decoder_factorization,
        )

        self._init_weights()

    def _init_weights(self) -> None:
        for m in self.modules():
            if isinstance(m, (nn.Linear, Linear)):
                nn.init.xavier_uniform_(m.weight, gain=nn.init.calculate_gain("relu"))
                if m.bias is not None:
                    nn.init.zeros_(m.bias)
            elif isinstance(m, nn.BatchNorm1d):
                nn.init.ones_(m.weight)
                nn.init.zeros_(m.bias)
            elif isinstance(m, nn.LayerNorm):
                nn.init.ones_(m.weight)
                nn.init.zeros_(m.bias)

    def forward(self, x_dict: Dict[str, Tensor],
                edge_index_dict: Dict) -> Dict[str, Tensor]:
        if hasattr(self, "cpg_embed"):
            cpg_x = self.cpg_embed
            if hasattr(self, "cpg_quality_mask"):
                cpg_x = cpg_x * self.cpg_quality_mask
            if hasattr(self, "cpg_bias"):
                cpg_x = cpg_x + self.cpg_bias
            x_dict["cpg"] = cpg_x

        x_proj = {}
        for k, v in x_dict.items():
            if k in self.proj:
                x_proj[k] = self.proj[k](v)
            else:
                x_proj[k] = torch.zeros(v.size(0), self.hidden_dim, device=v.device)

        if self.input_bn is not None:
            for k in x_proj:
                if k in self.input_bn and x_proj[k].size(0) > 1:
                    x_proj[k] = self.input_bn[k](x_proj[k])

        x0_dict = copy.copy(x_proj) if self.initial_residual else None
        z_dict = self.encoder(x_proj, edge_index_dict, x0_dict=x0_dict)
        return z_dict

    def decode(self, z_dict: Dict[str, Tensor], edge_index: Tensor) -> Tensor:
        z_gene = z_dict["gene"][edge_index[0]]
        z_path = z_dict["pathway"][edge_index[1]]
        return self.decoder(z_gene, z_path)

    def decode_chunked(self, z_dict: Dict[str, Tensor],
                       edge_index: Tensor, chunk_size: int = 16384) -> Tensor:
        n_edges = edge_index.size(1)
        if n_edges <= chunk_size:
            return self.decode(z_dict, edge_index)
        scores_list = []
        for start in range(0, n_edges, chunk_size):
            end = min(start + chunk_size, n_edges)
            chunk_ei = edge_index[:, start:end]
            scores_list.append(self.decode(z_dict, chunk_ei))
        return torch.cat(scores_list, dim=0)

    def to_cpg_learnable(self, cpg_init_feat: Tensor,
                         quality_mask: Optional[Tensor] = None) -> None:
        self.cpg_embed = nn.Parameter(cpg_init_feat.clone())
        self.register_parameter("cpg_embed", self.cpg_embed)

        if quality_mask is not None:
            self.cpg_quality_mask = quality_mask.float().unsqueeze(1)
            self.register_buffer("cpg_quality_mask", self.cpg_quality_mask)

        self.cpg_bias = nn.Parameter(torch.zeros(cpg_init_feat.size(1)))
        self.register_parameter("cpg_bias", self.cpg_bias)


def focal_bce_loss(logits: Tensor, labels: Tensor,
                   alpha: float = 0.5, gamma: float = 2.0) -> Tensor:
    logits = torch.clamp(logits, -10, 10)
    bce = F.binary_cross_entropy_with_logits(logits, labels, reduction="none")
    pt = torch.where(labels == 1, torch.sigmoid(logits), 1 - torch.sigmoid(logits))
    focal_weight = (1 - pt) ** gamma
    alpha_weight = torch.where(labels == 1, alpha, 1 - alpha)
    return (alpha_weight * focal_weight * bce).mean()


class NegEdgeSampler:
    def __init__(self, pos_edges: List[Tuple[int, int]], n_src: int, n_dst: int,
                 seed: int = 42, mode: str = "degree",
                 src_degrees: Optional[np.ndarray] = None,
                 dst_degrees: Optional[np.ndarray] = None,
                 degree_power: float = 0.75):
        self.n_src = n_src
        self.n_dst = n_dst
        self.mode = mode
        self.degree_power = degree_power
        self.exclude = {(int(s), int(d)) for s, d in pos_edges}
        self._exclude_frozen = frozenset(self.exclude)
        self.rng = np.random.RandomState(seed)
        self._fixed: Optional[Tensor] = None
        self._fixed_n_pos: Optional[int] = None
        self._fixed_neg_ratio: Optional[int] = None

        if mode == "degree" and src_degrees is not None and dst_degrees is not None:
            self.src_prob = np.power(np.maximum(src_degrees, 1.0), degree_power)
            self.src_prob /= self.src_prob.sum()
            self.dst_prob = np.power(np.maximum(dst_degrees, 1.0), degree_power)
            self.dst_prob /= self.dst_prob.sum()
        else:
            self.src_prob = None
            self.dst_prob = None

    def sample(self, n_pos: int, neg_ratio: int = 3) -> Tensor:
        n_neg = n_pos * neg_ratio
        src: np.ndarray
        dst: np.ndarray
        if self.src_prob is not None and self.dst_prob is not None:
            src = self.rng.choice(self.n_src, size=n_neg * 2, p=self.src_prob)
            dst = self.rng.choice(self.n_dst, size=n_neg * 2, p=self.dst_prob)
        else:
            src = self.rng.randint(0, self.n_src, size=n_neg * 2)
            dst = self.rng.randint(0, self.n_dst, size=n_neg * 2)

        neg_edges: List[Tuple[int, int]] = []
        for s, d in zip(src, dst):
            if len(neg_edges) >= n_neg:
                break
            key = (int(s), int(d))
            if key not in self._exclude_frozen:
                neg_edges.append(key)

        while len(neg_edges) < n_neg:
            if self.src_prob is not None:
                s = int(self.rng.choice(self.n_src, p=self.src_prob))
                d = int(self.rng.choice(self.n_dst, p=self.dst_prob))
            else:
                s = int(self.rng.randint(0, self.n_src))
                d = int(self.rng.randint(0, self.n_dst))
            key = (s, d)
            if key not in self._exclude_frozen:
                neg_edges.append(key)

        return torch.tensor(
            [[e[0] for e in neg_edges], [e[1] for e in neg_edges]], dtype=torch.long
        )

    def sample_fixed(self, n_pos: int, neg_ratio: int = 3,
                     force_refresh: bool = False) -> Tensor:
        if force_refresh or self._fixed is None or self._fixed_n_pos != n_pos or self._fixed_neg_ratio != neg_ratio:
            self._fixed = self.sample(n_pos, neg_ratio)
            self._fixed_n_pos = n_pos
            self._fixed_neg_ratio = neg_ratio
        return self._fixed


# ============================================================================
# train.py — 训练模块
# ============================================================================


def remove_edges_from_data(data: HeteroData, edge_type: Tuple[str, str, str],
                            mask: Tensor) -> HeteroData:
    new_data = HeteroData()
    for nt in data.node_types:
        new_data[nt].x = data[nt].x
    for et in data.edge_types:
        if et == edge_type:
            new_data[et].edge_index = data[et].edge_index[:, mask]
        else:
            new_data[et].edge_index = data[et].edge_index
    return new_data


def cosine_neg_ratio(epoch: int, total_epochs: int,
                     start_ratio: int = 1, end_ratio: int = 10) -> int:
    progress = epoch / max(total_epochs - 1, 1)
    factor = (1 + math.cos(math.pi * progress)) / 2
    ratio = int(round(end_ratio + (start_ratio - end_ratio) * factor))
    return max(1, ratio)


class PlattScaler:
    def __init__(self):
        self._calibrator: Optional[LogisticRegression] = None

    def fit(self, logits: np.ndarray, labels: np.ndarray) -> None:
        logits = np.asarray(logits, dtype=np.float64).reshape(-1, 1)
        labels = np.asarray(labels, dtype=np.float64)
        self._calibrator = LogisticRegression(
            C=1.0, solver="lbfgs", max_iter=1000,
        )
        self._calibrator.fit(logits, labels)

    def predict_proba(self, logits: np.ndarray) -> np.ndarray:
        if self._calibrator is None:
            return 1.0 / (1.0 + np.exp(-np.asarray(logits)))
        logits = np.asarray(logits, dtype=np.float64).reshape(-1, 1)
        return self._calibrator.predict_proba(logits)[:, 1]

    @property
    def is_fitted(self) -> bool:
        return self._calibrator is not None


@torch.inference_mode()
def evaluate(model: nn.Module, data: HeteroData, gp_edge_index: Tensor,
             gp_pos_idx: Tensor, neg_sampler: NegEdgeSampler,
             neg_ratio: int = 3, eval_batch: int = 65536,
             force_refresh_val_neg: bool = False) -> Tuple[float, float, np.ndarray, np.ndarray]:
    model.eval()
    device = next(model.parameters()).device
    data_device = data if data["gene"].x.device == device else data.to(device)

    z_dict = model(data_device.x_dict, data_device.edge_index_dict)

    pos_ei = gp_edge_index[:, gp_pos_idx].to(device)
    neg_ei = neg_sampler.sample_fixed(
        gp_pos_idx.shape[0], neg_ratio,
        force_refresh=force_refresh_val_neg,
    ).to(device)

    eval_ei = torch.cat([pos_ei, neg_ei], dim=1)
    n_total = eval_ei.size(1)
    labels = torch.cat([
        torch.ones(pos_ei.size(1), device=device),
        torch.zeros(neg_ei.size(1), device=device),
    ]).cpu().numpy()

    logits_list: List[np.ndarray] = []
    for start in range(0, n_total, eval_batch):
        end = min(start + eval_batch, n_total)
        batch_logits = model.decode(z_dict, eval_ei[:, start:end])
        logits_list.append(batch_logits.detach().cpu().numpy())
    logits_arr = np.concatenate(logits_list)
    scores = 1.0 / (1.0 + np.exp(-logits_arr))

    auroc = roc_auc_score(labels, scores) if len(np.unique(labels)) > 1 else 0.5
    auprc = average_precision_score(labels, scores)
    return auroc, auprc, scores, labels, logits_arr


def compute_composite_metric(auroc: float, auprc: float,
                              weight: float = 0.7) -> float:
    return weight * auroc + (1 - weight) * auprc


def train_fold(model: nn.Module, data_train: HeteroData, gp_edge_index: Tensor,
               train_idx: Tensor, val_idx: Tensor,
               neg_sampler: NegEdgeSampler,
               config: object,
               tb_writer: Optional[object] = None,
               fold_idx: int = 0) -> Tuple[float, float, Optional[PlattScaler]]:
    device = next(model.parameters()).device
    cfg = config.training
    model_cfg = config.model

    optimizer = AdamW(model.parameters(), lr=cfg.lr, weight_decay=cfg.weight_decay)
    print(f"    [DEBUG] Optimizer created, starting training loop...", flush=True)
    scheduler = ReduceLROnPlateau(optimizer, mode="max", factor=cfg.lr_factor,
                                  patience=cfg.lr_patience, min_lr=cfg.min_lr)

    best_auroc = 0.0
    best_auprc = 0.0
    best_composite = 0.0
    best_state = None
    patience_cnt = 0
    neg_ratio = cfg.neg_sample_ratio

    train_ei = gp_edge_index[:, train_idx]
    data_train_device = data_train if data_train["gene"].x.device == device else data_train.to(device)

    use_scaler = cfg.use_grad_scaling and device.type == "cuda"
    scaler = torch.amp.GradScaler(device="cuda") if use_scaler else None

    swa_model: Optional[Dict] = None
    swa_n = 0
    swa_active = False

    for epoch in range(cfg.epochs):
        model.train()

        if epoch == 0:
            print(f"    [DEBUG] Starting training loop, device={device}, "
                  f"train_edges={train_ei.shape[1]}, neg_ratio={neg_ratio}",
                  flush=True)

        if cfg.adaptive_neg_ratio:
            neg_ratio = cosine_neg_ratio(
                epoch, cfg.epochs,
                cfg.neg_sample_ratio_start, cfg.neg_sample_ratio_end,
            )

        neg_ei = neg_sampler.sample(train_ei.shape[1], neg_ratio).to(device)
        batch_ei = torch.cat([train_ei.to(device), neg_ei], dim=1)
        batch_labels = torch.cat([
            torch.ones(train_ei.size(1), device=device),
            torch.zeros(neg_ei.size(1), device=device),
        ])

        perm = torch.randperm(batch_ei.size(1), device=device)
        batch_ei = batch_ei[:, perm]
        batch_labels = batch_labels[perm]

        if scaler is not None:
            optimizer.zero_grad()
            with torch.amp.autocast(device_type="cuda"):
                z_dict = model(data_train_device.x_dict, data_train_device.edge_index_dict)
                logits = model.decode_chunked(z_dict, batch_ei)
                loss = focal_bce_loss(logits, batch_labels, cfg.focal_alpha, cfg.focal_gamma)
            scaler.scale(loss).backward()
            scaler.unscale_(optimizer)
            nn.utils.clip_grad_norm_(model.parameters(), cfg.grad_clip_norm)
            scaler.step(optimizer)
            scaler.update()
        else:
            optimizer.zero_grad()
            z_dict = model(data_train_device.x_dict, data_train_device.edge_index_dict)
            logits = model.decode_chunked(z_dict, batch_ei)
            loss = focal_bce_loss(logits, batch_labels, cfg.focal_alpha, cfg.focal_gamma)
            loss.backward()
            nn.utils.clip_grad_norm_(model.parameters(), cfg.grad_clip_norm)
            optimizer.step()

        eval_interval = 5
        if (epoch + 1) % eval_interval == 0:
            force_refresh = (epoch + 1) % cfg.val_neg_refresh_interval == 0
            auroc, auprc, val_scores, val_labels, _ = evaluate(
                model, data_train_device, gp_edge_index, val_idx,
                neg_sampler, neg_ratio, cfg.eval_batch,
                force_refresh_val_neg=force_refresh,
            )

            composite = compute_composite_metric(auroc, auprc, cfg.composite_weight)

            if cfg.early_stop_metric == "composite":
                is_better = composite > best_composite
            elif cfg.early_stop_metric == "auprc":
                is_better = auprc > best_auprc
            else:
                is_better = auroc > best_auroc

            if is_better:
                best_auroc = auroc
                best_auprc = auprc
                best_composite = composite
                best_state = copy.deepcopy(model.state_dict())
                patience_cnt = 0
                swa_active = False
            else:
                patience_cnt += 1
                if not swa_active and patience_cnt >= cfg.patience // 2:
                    swa_active = True
                if patience_cnt >= cfg.patience:
                    break

            scheduler.step(auroc)

            if tb_writer is not None:
                step = epoch + 1
                tb_writer.add_scalar(f"Fold{fold_idx}/Loss", loss.item(), step)
                tb_writer.add_scalar(f"Fold{fold_idx}/AUROC", auroc, step)
                tb_writer.add_scalar(f"Fold{fold_idx}/AUPRC", auprc, step)
                tb_writer.add_scalar(f"Fold{fold_idx}/Composite", composite, step)
                tb_writer.add_scalar(f"Fold{fold_idx}/NegRatio", neg_ratio, step)
                tb_writer.add_scalar(f"Fold{fold_idx}/LR", optimizer.param_groups[0]["lr"], step)

            if (epoch + 1) % 10 == 0 or epoch == 0:
                current_lr = optimizer.param_groups[0]["lr"]
                print(f"    Epoch {epoch+1} | Loss: {loss:.4f} | "
                      f"AUROC: {auroc:.4f} AUPRC: {auprc:.4f} "
                      f"Comp: {composite:.4f} | NegRatio: {neg_ratio} | "
                      f"LR: {current_lr:.2e}", flush=True)

        if swa_active:
            if swa_model is None:
                swa_model = copy.deepcopy(model.state_dict())
                swa_n = 1
            else:
                swa_n += 1
                for key in swa_model:
                    if swa_model[key].dtype in (torch.float32, torch.float64):
                        swa_model[key].data += model.state_dict()[key].data

        if device.type == "cuda" and (epoch + 1) % 50 == 0:
            torch.cuda.empty_cache()

        sys.stdout.flush()

    if swa_model is not None and swa_n > 1:
        for key in swa_model:
            if swa_model[key].dtype in (torch.float32, torch.float64):
                swa_model[key].data.div_(swa_n)
        model.load_state_dict(swa_model)
        swa_auroc, swa_auprc, _, _, _ = evaluate(
            model, data_train_device, gp_edge_index, val_idx,
            neg_sampler, neg_ratio, cfg.eval_batch,
        )
        swa_composite = compute_composite_metric(swa_auroc, swa_auprc, cfg.composite_weight)

        swa_better = False
        if cfg.early_stop_metric == "composite":
            swa_better = swa_composite > best_composite
        elif cfg.early_stop_metric == "auprc":
            swa_better = swa_auprc > best_auprc
        else:
            swa_better = swa_auroc > best_auroc

        if swa_better:
            best_auroc, best_auprc = swa_auroc, swa_auprc
            print(f"  Using SWA (SWA: AUROC={swa_auroc:.4f}, Best: val)")
        else:
            model.load_state_dict(best_state)
            print(f"  Using best checkpoint (Best: AUROC={best_auroc:.4f})")
    elif best_state is not None:
        model.load_state_dict(best_state)

    platt_scaler: Optional[PlattScaler] = None
    if config.inference.calibrate:
        _, _, _, _, val_logits = evaluate(
            model, data_train_device, gp_edge_index, val_idx,
            neg_sampler, neg_ratio, cfg.eval_batch,
        )
        platt_scaler = PlattScaler()
        platt_scaler.fit(val_logits, val_labels)

    return best_auroc, best_auprc, platt_scaler


def train_fold_minibatch(model: nn.Module, data_train: HeteroData,
                          gp_edge_index: Tensor,
                          train_idx: Tensor, val_idx: Tensor,
                          neg_sampler: NegEdgeSampler,
                          config: object,
                          tb_writer: Optional[object] = None,
                          wandb_run: Optional[object] = None,
                          fold_idx: int = 0) -> Tuple[float, float, Optional[PlattScaler]]:
    from torch_geometric.loader import HGTLoader

    device = next(model.parameters()).device
    cfg = config.training
    model_cfg = config.model

    optimizer = AdamW(model.parameters(), lr=cfg.lr, weight_decay=cfg.weight_decay)
    scheduler = ReduceLROnPlateau(optimizer, mode="max", factor=cfg.lr_factor,
                                  patience=cfg.lr_patience, min_lr=cfg.min_lr)

    best_auroc = 0.0
    best_auprc = 0.0
    best_composite = 0.0
    best_state = None
    patience_cnt = 0
    neg_ratio = cfg.neg_sample_ratio

    train_ei = gp_edge_index[:, train_idx]
    data_train_device = data_train if data_train["gene"].x.device == device else data_train.to(device)

    neighbor_sizes_dict = {}
    for nt in data_train.node_types:
        neighbor_sizes_dict[nt] = cfg.neighbor_sizes[:model_cfg.num_layers]

    train_src_nodes = train_ei[0].tolist()
    train_dst_nodes = train_ei[1].tolist()

    use_scaler = cfg.use_grad_scaling and device.type == "cuda"
    scaler = torch.amp.GradScaler(device="cuda") if use_scaler else None

    swa_model: Optional[Dict] = None
    swa_n = 0
    swa_active = False

    for epoch in range(cfg.epochs):
        model.train()

        if cfg.adaptive_neg_ratio:
            neg_ratio = cosine_neg_ratio(
                epoch, cfg.epochs,
                cfg.neg_sample_ratio_start, cfg.neg_sample_ratio_end,
            )

        perm = torch.randperm(len(train_src_nodes))
        total_loss = 0.0
        n_batches_processed = 0

        for batch_start in range(0, len(train_src_nodes), cfg.batch_size):
            batch_end = min(batch_start + cfg.batch_size, len(train_src_nodes))
            batch_perm = perm[batch_start:batch_end]

            batch_input_nodes = {
                "gene": [train_src_nodes[i] for i in batch_perm.tolist()],
                "pathway": [train_dst_nodes[i] for i in batch_perm.tolist()],
            }

            try:
                loader = HGTLoader(
                    data_train_device,
                    num_samples=neighbor_sizes_dict,
                    input_nodes=batch_input_nodes,
                    batch_size=cfg.batch_size,
                    shuffle=False,
                    num_workers=cfg.num_workers,
                )
                batch_data = next(iter(loader))

                if scaler is not None:
                    optimizer.zero_grad()
                    with torch.amp.autocast(device_type="cuda"):
                        z_dict = model(batch_data.x_dict, batch_data.edge_index_dict)
                        pos_logits = model.decode(z_dict, batch_data["gene", "involved_in", "pathway"].edge_index)
                        neg_ei = neg_sampler.sample(pos_logits.size(0), neg_ratio).to(device)
                        neg_logits = model.decode(z_dict, neg_ei)
                        logits = torch.cat([pos_logits, neg_logits])
                        labels = torch.cat([
                            torch.ones(pos_logits.size(0), device=device),
                            torch.zeros(neg_logits.size(0), device=device),
                        ])
                        loss = focal_bce_loss(logits, labels, cfg.focal_alpha, cfg.focal_gamma)
                    scaler.scale(loss).backward()
                    scaler.unscale_(optimizer)
                    nn.utils.clip_grad_norm_(model.parameters(), cfg.grad_clip_norm)
                    scaler.step(optimizer)
                    scaler.update()
                else:
                    optimizer.zero_grad()
                    z_dict = model(batch_data.x_dict, batch_data.edge_index_dict)
                    pos_logits = model.decode(z_dict, batch_data["gene", "involved_in", "pathway"].edge_index)
                    neg_ei = neg_sampler.sample(pos_logits.size(0), neg_ratio).to(device)
                    neg_logits = model.decode(z_dict, neg_ei)
                    logits = torch.cat([pos_logits, neg_logits])
                    labels = torch.cat([
                        torch.ones(pos_logits.size(0), device=device),
                        torch.zeros(neg_logits.size(0), device=device),
                    ])
                    loss = focal_bce_loss(logits, labels, cfg.focal_alpha, cfg.focal_gamma)
                    loss.backward()
                    nn.utils.clip_grad_norm_(model.parameters(), cfg.grad_clip_norm)
                    optimizer.step()

                total_loss += loss.item()
                n_batches_processed += 1

            except (StopIteration, RuntimeError) as e:
                continue

        avg_loss = total_loss / max(n_batches_processed, 1)

        eval_interval = 5
        if (epoch + 1) % eval_interval == 0:
            force_refresh = (epoch + 1) % cfg.val_neg_refresh_interval == 0
            auroc, auprc, val_scores, val_labels, _ = evaluate(
                model, data_train_device, gp_edge_index, val_idx,
                neg_sampler, neg_ratio, cfg.eval_batch,
                force_refresh_val_neg=force_refresh,
            )

            composite = compute_composite_metric(auroc, auprc, cfg.composite_weight)

            if cfg.early_stop_metric == "composite":
                is_better = composite > best_composite
            elif cfg.early_stop_metric == "auprc":
                is_better = auprc > best_auprc
            else:
                is_better = auroc > best_auroc

            if is_better:
                best_auroc = auroc
                best_auprc = auprc
                best_composite = composite
                best_state = copy.deepcopy(model.state_dict())
                patience_cnt = 0
                swa_active = False
            else:
                patience_cnt += 1
                if not swa_active and patience_cnt >= cfg.patience // 2:
                    swa_active = True
                if patience_cnt >= cfg.patience:
                    break

            scheduler.step(auroc)

            if tb_writer is not None:
                step = epoch + 1
                tb_writer.add_scalar(f"Fold{fold_idx}/Loss", avg_loss, step)
                tb_writer.add_scalar(f"Fold{fold_idx}/AUROC", auroc, step)
                tb_writer.add_scalar(f"Fold{fold_idx}/AUPRC", auprc, step)
                tb_writer.add_scalar(f"Fold{fold_idx}/Composite", composite, step)
                tb_writer.add_scalar(f"Fold{fold_idx}/NegRatio", neg_ratio, step)
                tb_writer.add_scalar(f"Fold{fold_idx}/LR", optimizer.param_groups[0]["lr"], step)

            if wandb_run is not None:
                wandb_run.log({
                    f"fold{fold_idx}/loss": avg_loss,
                    f"fold{fold_idx}/auroc": auroc,
                    f"fold{fold_idx}/auprc": auprc,
                    f"fold{fold_idx}/composite": composite,
                    f"fold{fold_idx}/neg_ratio": neg_ratio,
                    f"fold{fold_idx}/lr": optimizer.param_groups[0]["lr"],
                }, step=epoch + 1)

            if (epoch + 1) % 10 == 0 or epoch == 0:
                current_lr = optimizer.param_groups[0]["lr"]
                print(f"    Epoch {epoch+1} | Loss: {avg_loss:.4f} | "
                      f"AUROC: {auroc:.4f} AUPRC: {auprc:.4f} "
                      f"Comp: {composite:.4f} | NegRatio: {neg_ratio} | "
                      f"LR: {current_lr:.2e}", flush=True)

        if swa_active:
            if swa_model is None:
                swa_model = copy.deepcopy(model.state_dict())
                swa_n = 1
            else:
                swa_n += 1
                for key in swa_model:
                    if swa_model[key].dtype in (torch.float32, torch.float64):
                        swa_model[key].data += model.state_dict()[key].data

        if device.type == "cuda" and (epoch + 1) % 50 == 0:
            torch.cuda.empty_cache()

        sys.stdout.flush()

    if swa_model is not None and swa_n > 1:
        for key in swa_model:
            if swa_model[key].dtype in (torch.float32, torch.float64):
                swa_model[key].data.div_(swa_n)
        model.load_state_dict(swa_model)
        swa_auroc, swa_auprc, _, _, _ = evaluate(
            model, data_train_device, gp_edge_index, val_idx,
            neg_sampler, neg_ratio, cfg.eval_batch,
        )
        swa_composite = compute_composite_metric(swa_auroc, swa_auprc, cfg.composite_weight)

        swa_better = False
        if cfg.early_stop_metric == "composite":
            swa_better = swa_composite > best_composite
        elif cfg.early_stop_metric == "auprc":
            swa_better = swa_auprc > best_auprc
        else:
            swa_better = swa_auroc > best_auroc

        if swa_better:
            best_auroc, best_auprc = swa_auroc, swa_auprc
            print(f"  Using SWA (SWA: AUROC={swa_auroc:.4f}, Best: val)")
        else:
            model.load_state_dict(best_state)
            print(f"  Using best checkpoint (Best: AUROC={best_auroc:.4f})")
    elif best_state is not None:
        model.load_state_dict(best_state)

    platt_scaler: Optional[PlattScaler] = None
    if config.inference.calibrate:
        _, _, _, _, val_logits = evaluate(
            model, data_train_device, gp_edge_index, val_idx,
            neg_sampler, neg_ratio, cfg.eval_batch,
        )
        platt_scaler = PlattScaler()
        platt_scaler.fit(val_logits, val_labels)

    return best_auroc, best_auprc, platt_scaler


def cross_validate(data: HeteroData, gp_edge_index: Tensor,
                   n_genes: int, n_pathways: int,
                   gene_degrees: Optional[np.ndarray] = None,
                   pathway_degrees: Optional[np.ndarray] = None,
                   config: Optional[object] = None,
                   tb_writer: Optional[object] = None,
                   wandb_run: Optional[object] = None) -> Tuple[List[Dict], List[HGTModel], List[PlattScaler], Tensor]:
    cv_cfg = config.cv
    n_edges = gp_edge_index.size(1)
    if n_edges < cv_cfg.n_folds:
        return [], [], [], torch.zeros((2, 0), dtype=torch.long)

    kf = KFold(n_splits=cv_cfg.n_folds, shuffle=True, random_state=config.seed)
    cv_scores: List[Dict] = []
    cv_models: List[HGTModel] = []
    cv_scalers: List[PlattScaler] = []
    cv_heldout_edges_list: List[Tensor] = []

    folds_dir = Path(config.paths.folds_save_dir)
    folds_dir.mkdir(parents=True, exist_ok=True)

    for fold, (train_idx, val_idx) in enumerate(kf.split(range(n_edges))):
        print(f"\n{'='*50}")
        print(f"[CV] Fold {fold+1}/{cv_cfg.n_folds}")
        print(f"{'='*50}")

        train_mask_np = np.ones(n_edges, dtype=bool)
        train_mask_np[val_idx] = False

        data_train = remove_edges_from_data(
            data, ("gene", "involved_in", "pathway"),
            torch.from_numpy(train_mask_np),
        )

        all_pos_idx = np.concatenate([train_idx, val_idx])
        neg_sampler = NegEdgeSampler(
            pos_edges=gp_edge_index[:, all_pos_idx].t().tolist(),
            n_src=n_genes, n_dst=n_pathways,
            seed=config.seed + fold,
            mode=config.training.neg_sampling_mode,
            src_degrees=gene_degrees,
            dst_degrees=pathway_degrees,
            degree_power=config.training.neg_degree_power,
        )
        print(f"    [DEBUG] NegSampler created, creating model...", flush=True)

        model_cfg = config.model
        model = HGTModel(
            metadata=data.metadata(),
            dim_dict={nt: data[nt].x.size(-1) for nt in data.node_types},
            hidden_dim=model_cfg.hidden_dim,
            num_heads=model_cfg.num_heads,
            num_layers=model_cfg.num_layers,
            dropout=model_cfg.dropout,
            initial_residual=model_cfg.initial_residual,
            drop_edge_p=getattr(model_cfg, "drop_edge_p", 0.0),
            decoder_bias=model_cfg.decoder_bias,
            decoder_factorization=model_cfg.decoder_factorization,
            use_input_bn=getattr(model_cfg, "use_input_bn", True),
        ).to(config.device)

        if "cpg" in data.node_types:
            model.to_cpg_learnable(
                data["cpg"].x,
                quality_mask=data["cpg"].propagation_mask if hasattr(data["cpg"], "propagation_mask") else None,
            )

        print(f"    [DEBUG] Model created and on {config.device}, creating sampler...", flush=True)

        if hasattr(torch, "compile") and config.device == "cuda":
            try:
                import triton
                model = torch.compile(model, dynamic=True)
                print(f"    [DEBUG] torch.compile enabled (dynamic=True)", flush=True)
            except (ImportError, Exception) as e:
                print(f"    [DEBUG] torch.compile not available: {e}, using eager mode", flush=True)

        gp_ei_device = gp_edge_index.to(config.device)

        print(f"    [DEBUG] Calling train_fold...", flush=True)
        t_cfg = config.training
        if t_cfg.use_minibatch:
            auroc, auprc, platt_scaler = train_fold_minibatch(
                model, data_train.to(config.device), gp_ei_device,
                torch.from_numpy(train_idx).long(),
                torch.from_numpy(val_idx).long(),
                neg_sampler, config, tb_writer, wandb_run, fold,
            )
        else:
            auroc, auprc, platt_scaler = train_fold(
                model, data_train.to(config.device), gp_ei_device,
                torch.from_numpy(train_idx).long(),
                torch.from_numpy(val_idx).long(),
                neg_sampler, config, tb_writer, fold,
            )

        cv_scores.append({"fold": fold + 1, "auroc": auroc, "auprc": auprc})
        cv_models.append(model)
        cv_scalers.append(platt_scaler)
        cv_heldout_edges_list.append(gp_edge_index[:, val_idx])

        print(f"[CV] Fold {fold+1}: AUROC={auroc:.4f}, AUPRC={auprc:.4f}")

        if cv_cfg.save_fold_models:
            torch.save(model.state_dict(), str(folds_dir / f"fold_{fold+1}.pt"))

        if config.device == "cuda":
            torch.cuda.empty_cache()

    cv_heldout_edges = torch.cat(cv_heldout_edges_list, dim=1) if cv_heldout_edges_list else torch.zeros((2, 0), dtype=torch.long)
    return cv_scores, cv_models, cv_scalers, cv_heldout_edges


def train_final(model: nn.Module, data: HeteroData,
                gp_edge_index: Tensor,
                neg_sampler: NegEdgeSampler,
                config: object,
                tb_writer: Optional[object] = None) -> Tuple[float, float, Optional[PlattScaler]]:
    device = next(model.parameters()).device
    cfg = config.training

    n_edges = gp_edge_index.size(1)

    data_train = data.to(device)

    auroc, auprc, platt_scaler = train_fold(
        model, data_train, gp_edge_index.to(device),
        torch.arange(n_edges, device=device), torch.arange(n_edges, device=device),
        neg_sampler, config, tb_writer, fold_idx=999,
    )

    print(f"  Final training complete (full graph, no validation split)")
    return auroc, auprc, platt_scaler


# ============================================================================
# inference.py — 推理模块
# ============================================================================


@torch.inference_mode()
def _encode_once(model: nn.Module, data_device: HeteroData) -> Dict[str, Tensor]:
    return model(data_device.x_dict, data_device.edge_index_dict)


def mc_dropout_predict_per_gene(model: nn.Module, data_device: HeteroData,
                                gene_idx: int, pathway_idx: Tensor,
                                n_samples: int = 10) -> Tuple[np.ndarray, np.ndarray]:
    z_dict = _encode_once(model, data_device)

    edge_idx = torch.stack([
        torch.full((len(pathway_idx),), gene_idx, dtype=torch.long, device=data_device["gene"].x.device),
        pathway_idx,
    ])

    model.train()
    all_logits = []
    for _ in range(n_samples):
        logits = model.decode(z_dict, edge_idx)
        all_logits.append(logits.detach().cpu().numpy())

    model.eval()
    all_logits = np.stack(all_logits, axis=0)
    mean_logits = all_logits.mean(axis=0)
    std_logits = all_logits.std(axis=0)
    return mean_logits, std_logits


@torch.inference_mode()
def predict_bridge_pathways(model: nn.Module, data: HeteroData,
                             bridge_genes: List[str],
                             gene_to_idx: Dict[str, int],
                             pathway_names: List[str],
                             gene_idx_to_name: Optional[Dict[int, str]] = None,
                             gp_edge_index: Optional[Tensor] = None,
                             config: Optional[object] = None,
                             platt_scaler: Optional[PlattScaler] = None) -> pd.DataFrame:
    inf_cfg = config.inference if config else None
    top_k = inf_cfg.top_k if inf_cfg else 10
    temperature = inf_cfg.temperature if inf_cfg else 1.5
    mc_samples = inf_cfg.mc_dropout_samples if inf_cfg else 0
    calibrate = inf_cfg.calibrate if inf_cfg else True

    model.eval()
    device = next(model.parameters()).device
    data_device = data if data["gene"].x.device == device else data.to(device)

    known_gp_genes: Set[str] = set()
    inductive_ground_truth: Dict[str, Set[int]] = {}
    if gp_edge_index is not None:
        if gene_idx_to_name is not None:
            for i in range(gp_edge_index.size(1)):
                known_gp_genes.add(gene_idx_to_name[int(gp_edge_index[0, i])])
        else:
            idx_to_gene = {v: k for k, v in gene_to_idx.items()}
            for i in range(gp_edge_index.size(1)):
                known_gp_genes.add(idx_to_gene[int(gp_edge_index[0, i])])

    z_dict = _encode_once(model, data_device)

    results: List[Dict] = []
    n_pathways = len(pathway_names)
    pathway_idx_tensor = torch.arange(n_pathways, dtype=torch.long, device=device)

    for gene in bridge_genes:
        if gene not in gene_to_idx:
            continue
        gi = gene_to_idx[gene]

        if mc_samples > 0:
            mean_logits, std_logits = mc_dropout_predict_per_gene(
                model, data_device, gi, pathway_idx_tensor,
                n_samples=mc_samples,
            )
            raw_logits = mean_logits
            has_uncertainty = True
        else:
            edge_idx = torch.stack([
                torch.full((n_pathways,), gi, dtype=torch.long, device=device),
                pathway_idx_tensor,
            ])
            raw_logits = model.decode(z_dict, edge_idx).detach().cpu().numpy()
            has_uncertainty = False

        if calibrate and platt_scaler is not None and platt_scaler.is_fitted:
            scores = platt_scaler.predict_proba(raw_logits)
        else:
            scores = 1.0 / (1.0 + np.exp(-raw_logits / temperature))

        is_inductive = gene not in known_gp_genes if known_gp_genes else None

        for pi, pname in enumerate(pathway_names):
            row = {
                "gene_symbol": gene,
                "pathway_name": pname,
                "score": float(scores[pi]),
                "eval_mode": "inductive" if is_inductive else "transductive",
            }
            if has_uncertainty:
                row["uncertainty"] = float(std_logits[pi])
            results.append(row)

    df = pd.DataFrame(results)
    if df.empty:
        print("[Predict] No bridge genes found in graph!")
        return df

    df["rank"] = df.groupby("gene_symbol")["score"].rank(ascending=False, method="dense")
    df = df.sort_values(["gene_symbol", "rank"]).reset_index(drop=True)

    save_path = Path(config.paths.project_dir) / "bridge_pathway_scores.csv"
    df.to_csv(save_path, index=False)
    calib_tag = "+calibrated" if (calibrate and platt_scaler is not None and platt_scaler.is_fitted) else "+temp_scale"
    print(f"[Output] Saved bridge_pathway_scores.csv ({len(df)} rows{calib_tag})")

    topk = df[df["rank"] <= top_k].copy()
    print(f"\n{'='*60}")
    print(f"  Top-{top_k} Pathways per Bridge Gene")
    print(f"{'='*60}")
    for gene in bridge_genes:
        gene_rows = topk[topk["gene_symbol"] == gene]
        if gene_rows.empty:
            continue
        mode_tag = gene_rows["eval_mode"].iloc[0]
        print(f"\n  {gene} [{mode_tag}]:")
        for _, row in gene_rows.iterrows():
            pname = row["pathway_name"]
            pname_trim = pname[:60] + "..." if len(pname) > 60 else pname
            uncert_str = f" ±{row['uncertainty']:.4f}" if "uncertainty" in row else ""
            print(f"    #{int(row['rank']):2d} {pname_trim:60s} {row['score']:.4f}{uncert_str}")

    if known_gp_genes:
        n_trans = sum(1 for g in bridge_genes if g in known_gp_genes)
        n_ind = len(bridge_genes) - n_trans
        print(f"\n[Eval Mode] Transductive: {n_trans} genes, Inductive: {n_ind} genes")
        if n_ind > 0:
            induct_mask = df["gene_symbol"].apply(lambda g: g not in known_gp_genes)
            induct_df = df[induct_mask]
            if not induct_df.empty:
                gene_stats = induct_df.groupby("gene_symbol")["score"].agg(["mean", "max"])
                gene_stats.columns = ["mean_score", "max_score"]
                gene_stats = gene_stats.sort_values("max_score", ascending=False)
                print(f"  Inductive genes (top 20 by max_score):")
                for gene, row in gene_stats.head(20).iterrows():
                    print(f"    {gene}: mean={row['mean_score']:.4f}, max={row['max_score']:.4f}")
                if len(gene_stats) > 20:
                    print(f"    ... and {len(gene_stats) - 20} more")
                print(f"  Inductive overall: mean(mean_score)={gene_stats['mean_score'].mean():.4f}, "
                      f"mean(max_score)={gene_stats['max_score'].mean():.4f}")

    return df


def evaluate_inductive_subset(model: nn.Module, data: HeteroData,
                               gp_edge_index: Tensor,
                               bridge_genes: List[str],
                               gene_to_idx: Dict[str, int],
                               pathway_names: List[str],
                               gene_idx_to_name: Dict[int, str],
                               neg_sampler, config,
                               heldout_gp_edges: Optional[Tensor] = None) -> Tuple[float, float]:
    device = next(model.parameters()).device
    data_device = data if data["gene"].x.device == device else data.to(device)

    known_gp_genes: Set[str] = set()
    for i in range(gp_edge_index.size(1)):
        known_gp_genes.add(gene_idx_to_name[int(gp_edge_index[0, i])])

    inductive_genes = [g for g in bridge_genes if g in gene_to_idx and g not in known_gp_genes]
    if not inductive_genes:
        print("[Inductive Eval] No inductive genes to evaluate")
        return -1.0, -1.0

    n_pathways = len(pathway_names)
    model.eval()
    z_dict = _encode_once(model, data_device)

    all_scores = []
    for gene in inductive_genes:
        gi = gene_to_idx[gene]
        edge_idx = torch.stack([
            torch.full((n_pathways,), gi, dtype=torch.long, device=device),
            torch.arange(n_pathways, dtype=torch.long, device=device),
        ])
        logits = model.decode(z_dict, edge_idx)
        scores = torch.sigmoid(logits).detach().cpu().numpy()
        all_scores.append(scores)

    all_scores = np.concatenate(all_scores)

    print(f"[Inductive Eval] {len(inductive_genes)} inductive genes, "
          f"mean_score={all_scores.mean():.4f}, "
          f"max_score={all_scores.max():.4f}, "
          f"score_std={all_scores.std():.4f}")

    if heldout_gp_edges is not None and heldout_gp_edges.size(1) > 0:
        heldout_pos: Set[Tuple[int, int]] = set()
        for i in range(heldout_gp_edges.size(1)):
            g_idx = int(heldout_gp_edges[0, i])
            p_idx = int(heldout_gp_edges[1, i])
            if g_idx < len(gene_idx_to_name):
                g_name = gene_idx_to_name[g_idx]
                if g_name in inductive_genes:
                    heldout_pos.add((g_idx, p_idx))

        if heldout_pos:
            n_pos = len(heldout_pos)
            neg_ei = neg_sampler.sample(n_pos, neg_ratio=3).to(device)

            pos_scores_list = []
            for g_idx, p_idx in heldout_pos:
                edge_idx = torch.tensor([[g_idx], [p_idx]], dtype=torch.long, device=device)
                logits = model.decode(z_dict, edge_idx)
                pos_scores_list.append(torch.sigmoid(logits).cpu().numpy())
            pos_scores = np.concatenate(pos_scores_list)

            neg_scores_list = []
            for start in range(0, neg_ei.size(1), 16384):
                end = min(start + 16384, neg_ei.size(1))
                batch_ei = neg_ei[:, start:end]
                neg_logits = model.decode(z_dict, batch_ei)
                neg_scores_list.append(torch.sigmoid(neg_logits).cpu().numpy())
            neg_scores = np.concatenate(neg_scores_list)

            labels = np.concatenate([np.ones(len(pos_scores)), np.zeros(len(neg_scores))])
            scores = np.concatenate([pos_scores, neg_scores])

            auroc = roc_auc_score(labels, scores) if len(np.unique(labels)) > 1 else 0.5
            auprc = average_precision_score(labels, scores)
            print(f"[Inductive Eval] With {n_pos} heldout positive edges: "
                  f"AUROC={auroc:.4f}, AUPRC={auprc:.4f}")
            return auroc, auprc
        else:
            print("[Inductive Eval] No heldout edges match inductive genes, "
                  "cannot compute AUROC/AUPRC")

    return -1.0, -1.0


def ensemble_predict_bridge_pathways(
    models: List[nn.Module], data: HeteroData,
    bridge_genes: List[str], gene_to_idx: Dict[str, int],
    pathway_names: List[str],
    gp_edge_index: Optional[Tensor] = None,
    config: Optional[object] = None,
    platt_scalers: Optional[List[PlattScaler]] = None,
) -> pd.DataFrame:
    if not models:
        print("[Ensemble] No models provided!")
        return pd.DataFrame()

    inf_cfg = config.inference if config else None
    top_k = inf_cfg.top_k if inf_cfg else 10
    temperature = inf_cfg.temperature if inf_cfg else 1.5

    device = next(models[0].parameters()).device
    data_device = data if data["gene"].x.device == device else data.to(device)

    gene_idx_to_name = {v: k for k, v in gene_to_idx.items()}
    known_gp_genes: Set[str] = set()
    if gp_edge_index is not None:
        for i in range(gp_edge_index.size(1)):
            known_gp_genes.add(gene_idx_to_name[int(gp_edge_index[0, i])])

    n_pathways = len(pathway_names)
    n_genes_valid = sum(1 for g in bridge_genes if g in gene_to_idx)
    all_logits = np.zeros((n_genes_valid, n_pathways), dtype=np.float32)
    gene_order: List[str] = []

    for gene in bridge_genes:
        if gene in gene_to_idx:
            gene_order.append(gene)

    for model in models:
        model.eval()

    with torch.inference_mode():
        for model in models:
            z_dict = _encode_once(model, data_device)
            for idx, gene in enumerate(gene_order):
                gi = gene_to_idx[gene]
                edge_idx = torch.stack([
                    torch.full((n_pathways,), gi, dtype=torch.long, device=device),
                    torch.arange(n_pathways, dtype=torch.long, device=device),
                ])
                logits = model.decode(z_dict, edge_idx)
                all_logits[idx] += logits.cpu().numpy()

    all_logits /= len(models)

    results: List[Dict] = []
    for idx, gene in enumerate(gene_order):
        scores = 1.0 / (1.0 + np.exp(-all_logits[idx] / temperature))
        is_inductive = gene not in known_gp_genes if known_gp_genes else None
        for pi, pname in enumerate(pathway_names):
            results.append({
                "gene_symbol": gene,
                "pathway_name": pname,
                "score": float(scores[pi]),
                "eval_mode": "inductive" if is_inductive else "transductive",
            })

    df = pd.DataFrame(results)
    if df.empty:
        return df

    df["rank"] = df.groupby("gene_symbol")["score"].rank(ascending=False, method="dense")
    df = df.sort_values(["gene_symbol", "rank"]).reset_index(drop=True)

    save_path = Path(config.paths.project_dir) / "bridge_pathway_scores_ensemble.csv"
    df.to_csv(save_path, index=False)
    print(f"[Output] Saved bridge_pathway_scores_ensemble.csv ({len(df)} rows, {len(models)} models)")

    topk = df[df["rank"] <= top_k].copy()
    print(f"\n{'='*60}")
    print(f"  Top-{top_k} Pathways per Bridge Gene (Ensemble of {len(models)} models)")
    print(f"{'='*60}")
    for gene in gene_order[:30]:
        gene_rows = topk[topk["gene_symbol"] == gene]
        if gene_rows.empty:
            continue
        mode_tag = gene_rows["eval_mode"].iloc[0]
        print(f"\n  {gene} [{mode_tag}]:")
        for _, row in gene_rows.iterrows():
            pname = str(row["pathway_name"])
            pname_trim = pname[:60] + "..." if len(pname) > 60 else pname
            print(f"    #{int(row['rank']):2d} {pname_trim:60s} {row['score']:.4f}")

    return df


def filter_and_report_brain_ischemia(df: pd.DataFrame,
                                     bridge_genes: List[str],
                                     top_k: int = 5,
                                     score_thresh: float = 0.8) -> pd.DataFrame:
    keywords = [
        "apoptosis", "necroptosis", "autophagy", "ferroptosis",
        "NF-kappa B", "TNF", "Toll-like receptor", "MAPK",
        "PI3K-Akt", "HIF-1", "VEGF", "p53", "JAK-STAT",
        "oxidative stress", "reactive oxygen", "calcium",
        "cGMP-PKG", "neurotrophin", "glutamate", "GABA",
        "focal adhesion", "tight junction", "leukocyte",
        "complement", "NOD-like receptor", "inflammasome",
        "interleukin", "chemokine", "integrin", "endocytosis",
        "autophagy", "mitophagy", "ubiquitin", "proteasome",
    ]
    pattern = '|'.join(keywords)
    mask = df["pathway_name"].str.contains(pattern, case=False, na=False)
    candidate = df[mask & (df["rank"] <= top_k) & (df["score"] >= score_thresh)]

    if not candidate.empty:
        save_path = Path.cwd() / "bridge_pathway_cirI_candidates.csv"
        try:
            candidate.to_csv(save_path, index=False)
            saved_ok = True
        except PermissionError:
            alt_path = Path.cwd() / f"bridge_pathway_cirI_candidates_{pd.Timestamp.now():%Y%m%d_%H%M%S}.csv"
            candidate.to_csv(alt_path, index=False)
            save_path = alt_path
            saved_ok = True
            print(f"  [Warning] Primary file locked, saved to {alt_path.name}")
        print(f"\n{'='*60}")
        print(f"  CIRI-relevant pathway candidates (Top-{top_k}, score>={score_thresh})")
        print(f"{'='*60}")
        for gene in candidate["gene_symbol"].unique():
            gene_rows = candidate[candidate["gene_symbol"] == gene]
            print(f"\n  {gene}:")
            for _, row in gene_rows.iterrows():
                pname = str(row["pathway_name"])
                pname_trim = pname[:60] + "..." if len(pname) > 60 else pname
                print(f"    #{int(row['rank']):2d} {pname_trim:60s} {row['score']:.4f}")
        print(f"\n[Output] Saved bridge_pathway_cirI_candidates.csv ({len(candidate)} rows)")
    else:
        print("\n[Filter] No CIRI-relevant pathways found at current threshold")

    return candidate


# ============================================================================
# pretrain.py — GAE 预训练
# ============================================================================


class HeteroGAE(nn.Module):
    def __init__(self, metadata: Tuple[List[str], List[Tuple[str, str, str]]],
                 dim_dict: Dict[str, int], hidden_dim: int,
                 num_heads: int, num_layers: int, dropout: float):
        super().__init__()
        self.metadata = metadata
        self.hidden_dim = hidden_dim

        self.proj = nn.ModuleDict()
        for nt, d_in in dim_dict.items():
            self.proj[nt] = Linear(d_in, hidden_dim)

        self.encoder = HGTEncoder(
            metadata, hidden_dim, num_heads, num_layers,
            dropout, initial_residual=True,
        )

        self._init_weights()

    def _init_weights(self) -> None:
        for m in self.modules():
            if isinstance(m, (nn.Linear, Linear)):
                nn.init.xavier_uniform_(m.weight, gain=nn.init.calculate_gain("relu"))
                if m.bias is not None:
                    nn.init.zeros_(m.bias)
            elif isinstance(m, nn.LayerNorm):
                nn.init.ones_(m.weight)
                nn.init.zeros_(m.bias)

    def encode(self, data: HeteroData) -> Dict[str, Tensor]:
        x_proj = {}
        for k, v in data.x_dict.items():
            if k in self.proj:
                x_proj[k] = self.proj[k](v)
            else:
                x_proj[k] = torch.zeros(v.size(0), self.hidden_dim, device=v.device)

        x0_dict = copy.copy(x_proj)
        z_dict = self.encoder(x_proj, data.edge_index_dict, x0_dict=x0_dict)
        return z_dict

    def decode(self, z_dict: Dict[str, Tensor],
               edge_index: Tensor, src_type: str = "gene",
               dst_type: str = "gene") -> Tensor:
        z_src = z_dict[src_type][edge_index[0]]
        z_dst = z_dict[dst_type][edge_index[1]]
        scores = (z_src * z_dst).sum(dim=-1)
        return scores

    def forward(self, data: HeteroData,
                pos_edge_index: Tensor,
                neg_edge_index: Tensor) -> Tuple[Tensor, Tensor]:
        z_dict = self.encode(data)
        pos_scores = self.decode(z_dict, pos_edge_index)
        neg_scores = self.decode(z_dict, neg_edge_index)
        return pos_scores, neg_scores

    def transfer_weights_to(self, target_model: HGTModel) -> None:
        for nt in self.proj:
            if nt in target_model.proj:
                target_model.proj[nt].load_state_dict(self.proj[nt].state_dict())

        target_enc = target_model.encoder
        pretrain_enc = self.encoder

        for i, (t_conv, p_conv) in enumerate(zip(target_enc.convs, pretrain_enc.convs)):
            t_conv.load_state_dict(p_conv.state_dict())

        for i, (t_norm, p_norm) in enumerate(zip(target_enc.norms, pretrain_enc.norms)):
            t_norm.load_state_dict(p_norm.state_dict())

        if (target_enc.skip_alphas is not None and
                pretrain_enc.skip_alphas is not None):
            for i, (t_alpha, p_alpha) in enumerate(
                zip(target_enc.skip_alphas, pretrain_enc.skip_alphas)
            ):
                t_alpha.data.copy_(p_alpha.data)


def mask_edges(edge_index: Tensor, mask_ratio: float = 0.3,
               seed: int = 42) -> Tuple[Tensor, Tensor, Tensor]:
    rng = np.random.RandomState(seed)
    n_edges = edge_index.size(1)
    n_mask = max(1, int(n_edges * mask_ratio))
    perm = rng.permutation(n_edges)
    mask_idx = perm[:n_mask]
    train_idx = perm[n_mask:]

    mask_ei = edge_index[:, mask_idx]
    train_ei = edge_index[:, train_idx]

    n_src = int(edge_index[0].max().item()) + 1
    n_dst = int(edge_index[1].max().item()) + 1

    pos_set = set()
    for i in range(n_edges):
        pos_set.add((int(edge_index[0, i]), int(edge_index[1, i])))

    neg_edges_list = []
    while len(neg_edges_list) < n_mask:
        s = rng.randint(0, n_src)
        d = rng.randint(0, n_dst)
        if (s, d) not in pos_set:
            neg_edges_list.append((s, d))

    neg_ei = torch.tensor(
        [[e[0] for e in neg_edges_list], [e[1] for e in neg_edges_list]],
        dtype=torch.long,
    )
    return train_ei, mask_ei, neg_ei


def pretrain_gae(data: HeteroData, config: object,
                 device: torch.device) -> Optional[HeteroGAE]:
    pt_cfg = config.pretraining
    if not pt_cfg.enabled:
        logger_utils.info("Pretraining disabled, skipping...")
        return None

    m_cfg = config.model

    logger_utils.info("=" * 60)
    logger_utils.info("[Pretrain] Heterogeneous GAE Pre-training")
    logger_utils.info(f"  Epochs: {pt_cfg.epochs}, LR: {pt_cfg.lr}, "
                      f"Mask Ratio: {pt_cfg.mask_ratio}")
    logger_utils.info("=" * 60)

    ppi_edge_types = [
        ("gene", "strong_ppi", "gene"),
        ("gene", "weak_ppi", "gene"),
        ("gene", "interacts", "gene"),
    ]
    ppi_et = None
    for et in ppi_edge_types:
        if et in data.edge_types:
            ppi_et = et
            break

    if ppi_et is None:
        logger_utils.warning("No PPI edges found for pretraining, skipping...")
        return None

    ppi_ei = data[ppi_et].edge_index
    logger_utils.info(f"  Using edge type: {ppi_et} ({ppi_ei.size(1)} edges)")

    train_ei, mask_ei, neg_ei = mask_edges(
        ppi_ei, pt_cfg.mask_ratio, config.seed,
    )
    logger_utils.info(f"  Train edges: {train_ei.size(1)}, "
                      f"Masked edges: {mask_ei.size(1)}, "
                      f"Neg edges: {neg_ei.size(1)}")

    mask_np = np.zeros(ppi_ei.size(1), dtype=bool)
    rng = np.random.RandomState(config.seed)
    n_mask = mask_ei.size(1)
    n_total = ppi_ei.size(1)
    perm = rng.permutation(n_total)
    mask_np[perm[:n_mask]] = True
    data_train = remove_edges_from_data(
        data, ppi_et, torch.from_numpy(~mask_np),
    ).to(device)

    gae = HeteroGAE(
        metadata=data.metadata(),
        dim_dict={nt: data[nt].x.size(-1) for nt in data.node_types},
        hidden_dim=m_cfg.hidden_dim,
        num_heads=m_cfg.num_heads,
        num_layers=m_cfg.num_layers,
        dropout=m_cfg.dropout,
    ).to(device)

    optimizer = torch.optim.AdamW(
        gae.parameters(), lr=pt_cfg.lr,
        weight_decay=pt_cfg.weight_decay,
    )

    mask_ei_dev = mask_ei.to(device)
    neg_ei_dev = neg_ei.to(device)

    best_loss = float("inf")
    best_state = None

    for epoch in range(pt_cfg.epochs):
        gae.train()
        optimizer.zero_grad()

        pos_scores, neg_scores = gae(data_train, mask_ei_dev, neg_ei_dev)

        pos_loss = F.binary_cross_entropy_with_logits(
            pos_scores, torch.ones_like(pos_scores),
        )
        neg_loss = F.binary_cross_entropy_with_logits(
            neg_scores, torch.zeros_like(neg_scores),
        )
        loss = pos_loss + neg_loss

        loss.backward()
        optimizer.step()

        with torch.no_grad():
            pos_acc = (torch.sigmoid(pos_scores) > 0.5).float().mean().item()
            neg_acc = (torch.sigmoid(neg_scores) < 0.5).float().mean().item()

        if loss.item() < best_loss:
            best_loss = loss.item()
            best_state = copy.deepcopy(gae.state_dict())

        if (epoch + 1) % 10 == 0 or epoch == 0:
            logger_utils.info(
                f"  Pretrain Epoch {epoch+1:3d}/{pt_cfg.epochs} | "
                f"Loss: {loss:.4f} | PosAcc: {pos_acc:.3f} | "
                f"NegAcc: {neg_acc:.3f}"
            )

    if best_state is not None:
        gae.load_state_dict(best_state)
        logger_utils.info(f"  Pretrain complete: Best Loss = {best_loss:.4f}")

    save_path = config.paths.pretrain_save if hasattr(config.paths, "pretrain_save") else None
    if save_path:
        torch.save(gae.state_dict(), save_path)
        logger_utils.info(f"  Pretrained model saved to {save_path}")

    return gae


# ============================================================================
# main.py — 主入口
# ============================================================================


def main() -> None:
    config_path = Path(__file__).resolve().parent / "config.yaml"
    config = load_config(str(config_path) if config_path.exists() else None)

    log_cfg = config.logging
    log_file = log_cfg.log_file if log_cfg.log_file else None
    logger = setup_logging(
        "hgt_pipeline", log_cfg.level, log_file,
        enable_tensorboard=log_cfg.tensorboard,
    )

    tb_writer = None
    if log_cfg.tensorboard:
        log_dir = Path(config.paths.log_dir)
        tb_writer = setup_tensorboard(str(log_dir))

    wandb_run = setup_wandb(config, logger)

    set_seed(config.seed)
    logger.info(f"Device: {config.device}, Seed: {config.seed}")
    logger.info(f"Config loaded from embedded YAML")

    m_cfg = config.model
    t_cfg = config.training
    logger.info(f"Model: layers={m_cfg.num_layers}, hidden={m_cfg.hidden_dim}, "
                f"heads={m_cfg.num_heads}, dropout={m_cfg.dropout}")
    logger.info(f"Model: initial_residual={m_cfg.initial_residual}, "
                f"decoder_bias={m_cfg.decoder_bias}, factorization={m_cfg.decoder_factorization}")
    logger.info(f"Training: epochs={t_cfg.epochs}, patience={t_cfg.patience}, "
                f"lr={t_cfg.lr}, wd={t_cfg.weight_decay}")
    logger.info(f"Training: neg_mode={t_cfg.neg_sampling_mode}, "
                f"adaptive_neg={t_cfg.adaptive_neg_ratio}, "
                f"early_stop={t_cfg.early_stop_metric}")
    logger.info(f"Graph: ppi_split={config.graph.ppi_split_by_score}, "
                f"methyl_directed={config.graph.methylation_directed}")

    # [1/6] 数据加载
    logger.info("=" * 60)
    logger.info("[1/6] Loading data...")

    data_dict = load_all_data(config)

    gene_feat_arr = data_dict["gene_feat_arr"]
    gene_feat_names = data_dict["gene_feat_names"]
    pathway_feat_arr = data_dict["pathway_feat_arr"]
    drug_fp_arr = data_dict["drug_fp_arr"]
    disease_feat_arr = data_dict["disease_feat_arr"]
    pathway_names = data_dict["pathway_names"]
    bridge_genes = data_dict["bridge_genes"]
    ppi_edges = data_dict["ppi_edges"]
    coexp_edges = data_dict["coexp_edges"]
    tf_edges = data_dict["tf_edges"]
    gene_pathway_edges = data_dict["gene_pathway_edges"]
    all_genes_list = data_dict["all_genes_list"]
    methyl_edges = data_dict.get("methyl_edges")
    mirna_edges = data_dict.get("mirna_edges")
    pathway_hierarchy = data_dict.get("pathway_hierarchy", [])
    disease_pathway_edges = data_dict.get("disease_pathway_edges", [])

    logger.info(f"  Gene features: {gene_feat_arr.shape}")
    logger.info(f"  Pathway features: {pathway_feat_arr.shape}")
    logger.info(f"  Bridge genes: {len(bridge_genes)}")
    logger.info(f"  PPI edges: {len(ppi_edges)}")
    logger.info(f"  Gene-pathway edges: {len(gene_pathway_edges)}")
    logger.info(f"  Pathway hierarchy: {len(pathway_hierarchy)}")

    # [2/6] 特征预处理
    logger.info("=" * 60)
    logger.info("[2/6] Preprocessing features...")

    preproc_cfg = config.preprocessing
    feature_dim = preproc_cfg.feature_dim

    if preproc_cfg.use_pca and gene_feat_arr.shape[1] > feature_dim:
        if preproc_cfg.use_umap and gene_feat_arr.shape[1] > preproc_cfg.pca_intermediate_dim:
            logger.info(f"  Gene features: PCA({preproc_cfg.pca_intermediate_dim}) -> UMAP({feature_dim})")
            gene_feat_arr = pca_then_umap(
                gene_feat_arr,
                preproc_cfg.pca_intermediate_dim,
                feature_dim,
                seed=config.seed,
            )
        else:
            gene_feat_arr, _ = pca_reduce(gene_feat_arr, feature_dim, config.seed)
            logger.info(f"  Gene features PCA: {gene_feat_arr.shape}")

    pathway_pca_cache = Path(config.paths.pathway_pca_cache)
    if preproc_cfg.use_pca and pathway_feat_arr.shape[1] > feature_dim:
        if pathway_pca_cache.exists():
            pathway_feat_arr = np.load(str(pathway_pca_cache)).astype(np.float32)
            logger.info(f"  Pathway PCA cache loaded: {pathway_feat_arr.shape}")
        else:
            pathway_feat_arr, _ = pca_reduce(pathway_feat_arr, feature_dim, config.seed)
            np.save(str(pathway_pca_cache), pathway_feat_arr)
            logger.info(f"  Pathway PCA: {pathway_feat_arr.shape} (cached)")

    if len(pathway_names) != pathway_feat_arr.shape[0]:
        logger.warning(f"  pathway_names ({len(pathway_names)}) != features ({pathway_feat_arr.shape[0]})")
        pathway_names = [f"pathway_{i}" for i in range(pathway_feat_arr.shape[0])]

    # [3/6] 构建异构图
    logger.info("=" * 60)
    logger.info("[3/6] Building heterogeneous graph...")

    data, gene_to_idx, gene_list, pathway_name_to_idx = build_hetero_graph(
        gene_feat_arr, gene_feat_names,
        drug_fp_arr, disease_feat_arr,
        pathway_feat_arr, pathway_names,
        ppi_edges, coexp_edges, tf_edges,
        gene_pathway_edges,
        all_genes_list=all_genes_list,
        bridge_genes=bridge_genes,
        methyl_edges=methyl_edges,
        mirna_edges=mirna_edges,
        pathway_hierarchy=pathway_hierarchy,
        disease_pathway_edges=disease_pathway_edges,
        config=config,
    )

    gp_edge_index = data["gene", "involved_in", "pathway"].edge_index
    n_genes = data["gene"].x.size(0)
    n_pathways = data["pathway"].x.size(0)

    if gp_edge_index.size(1) == 0:
        logger.error("No gene-pathway edges found! Aborting.")
        return

    logger.info(f"  Gene nodes: {n_genes}, Pathway nodes: {n_pathways}")
    logger.info(f"  Gene-pathway edges (positive): {gp_edge_index.size(1)}")
    logger.info(f"  Edge types: {list(data.edge_types)}")

    logger.info("  Normalization stats stored per node type for future reload.")
    logger.info("  To reload a saved HeteroData object, call:")
    logger.info("    apply_cached_normalization(data)")

    # [4/6] 节点度计算
    logger.info("[4/6] Computing node degrees...")

    gp_ei_np = gp_edge_index.cpu().numpy()
    gene_degrees = np.ones(n_genes, dtype=np.float64)
    pathway_degrees = np.ones(n_pathways, dtype=np.float64)
    for i in range(gp_ei_np.shape[1]):
        gene_degrees[gp_ei_np[0, i]] += 1
        pathway_degrees[gp_ei_np[1, i]] += 1
    logger.info(f"  Gene degree: max={gene_degrees.max():.0f}, mean={gene_degrees.mean():.1f}")
    logger.info(f"  Pathway degree: max={pathway_degrees.max():.0f}, mean={pathway_degrees.mean():.1f}")

    # 4.5 GAE 预训练
    pretrained_gae = None
    if config.pretraining.enabled:
        pretrained_gae = pretrain_gae(data, config, torch.device(config.device))

    # [5/6] 交叉验证
    logger.info("=" * 60)
    logger.info(f"[5/6] {config.cv.n_folds}-fold Cross Validation...")

    cv_results, cv_models, cv_scalers, cv_heldout_edges = cross_validate(
        data, gp_edge_index, n_genes, n_pathways,
        gene_degrees, pathway_degrees, config, tb_writer, wandb_run,
    )

    if cv_results:
        aurocs = [m["auroc"] for m in cv_results]
        auprcs = [m["auprc"] for m in cv_results]
        logger.info(f"  CV AUROC: {np.mean(aurocs):.4f} +/- {np.std(aurocs):.4f}")
        logger.info(f"  CV AUPRC: {np.mean(auprcs):.4f} +/- {np.std(auprcs):.4f}")

    # [6/6] 最终训练 + 推理
    logger.info("=" * 60)
    logger.info("[6/6] Final training & inference...")

    final_neg_sampler = NegEdgeSampler(
        pos_edges=gp_edge_index.t().tolist(),
        n_src=n_genes, n_dst=n_pathways,
        seed=config.seed + 999,
        mode=config.training.neg_sampling_mode,
        src_degrees=gene_degrees,
        dst_degrees=pathway_degrees,
        degree_power=config.training.neg_degree_power,
    )

    final_model = HGTModel(
        metadata=data.metadata(),
        dim_dict={nt: data[nt].x.size(-1) for nt in data.node_types},
        hidden_dim=m_cfg.hidden_dim, num_heads=m_cfg.num_heads,
        num_layers=m_cfg.num_layers, dropout=m_cfg.dropout,
        initial_residual=m_cfg.initial_residual,
        drop_edge_p=getattr(m_cfg, "drop_edge_p", 0.0),
        decoder_bias=m_cfg.decoder_bias,
        decoder_factorization=m_cfg.decoder_factorization,
        use_input_bn=getattr(m_cfg, "use_input_bn", True),
    ).to(config.device)

    if "cpg" in data.node_types:
        final_model.to_cpg_learnable(
            data["cpg"].x,
            quality_mask=data["cpg"].propagation_mask if hasattr(data["cpg"], "propagation_mask") else None,
        )
        logger.info("  CpG learnable parameters registered")

    final_model = try_compile(final_model)

    if pretrained_gae is not None:
        logger.info("  Transferring pretrained GAE weights to HGTModel...")
        pretrained_gae.transfer_weights_to(final_model)
        logger.info("  Pretrained weights transferred successfully")

    final_auroc, final_auprc, final_scaler = train_final(
        final_model, data, gp_edge_index.to(config.device),
        final_neg_sampler, config, tb_writer,
    )
    logger.info(f"  Final training: full graph, no validation split (CV provides generalization)")

    model_path = Path(config.paths.model_save)
    torch.save(final_model.state_dict(), str(model_path))
    logger.info(f"  Model saved to {model_path}")

    gene_idx_to_name = {v: k for k, v in gene_to_idx.items()}

    if bridge_genes:
        logger.info("  Running final inference...")
        mc_samples = config.inference.mc_dropout_samples
        if mc_samples > 0:
            logger.info(f"  MC Dropout enabled: {mc_samples} samples per gene")

        predict_bridge_pathways(
            final_model, data.to(config.device),
            bridge_genes, gene_to_idx, pathway_names,
            gene_idx_to_name=gene_idx_to_name,
            gp_edge_index=gp_edge_index,
            config=config,
            platt_scaler=final_scaler,
        )

        if config.inference.eval_inductive:
            evaluate_inductive_subset(
                final_model, data.to(config.device),
                gp_edge_index, bridge_genes,
                gene_to_idx, pathway_names,
                gene_idx_to_name,
                final_neg_sampler, config,
                heldout_gp_edges=cv_heldout_edges,
            )

    if config.cv.ensemble_inference and cv_models:
        logger.info("  Running ensemble inference...")
        ensemble_predict_bridge_pathways(
            cv_models, data.to(config.device),
            bridge_genes, gene_to_idx, pathway_names,
            gp_edge_index=gp_edge_index,
            config=config,
            platt_scalers=cv_scalers,
        )

    logger.info("  Filtering CIRI-relevant pathways...")
    scores_path = Path(config.paths.project_dir) / "bridge_pathway_scores.csv"
    if scores_path.exists():
        try:
            df_scores = pd.read_csv(str(scores_path))
            filter_and_report_brain_ischemia(
                df_scores, bridge_genes,
                top_k=config.inference.top_k,
            )
        except PermissionError as e:
            logger.warning(f"  Could not read/write CIRI filter files: {e}")
    else:
        logger.warning("  bridge_pathway_scores.csv not found, skipping CIRI filter")

    if tb_writer is not None:
        tb_writer.close()
    if wandb_run is not None:
        wandb_run.finish()
    if config.device == "cuda":
        torch.cuda.empty_cache()

    logger.info("=" * 60)
    logger.info("  Pipeline Complete!")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()