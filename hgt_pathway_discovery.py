# -*- coding: utf-8 -*-
"""
HGT (Heterogeneous Graph Transformer) 基因‑通路关联预测模型
用于网络药理学第二阶段通路发现任务。

参考文献:
  - HGT: Heterogeneous Graph Transformer, Hu et al., WWW 2020
  - MHGTMDA: Molecular Heterogeneous Graph Transformer for miRNA-disease, Zou et al., 2024
  - HGTDR: Advancing Drug Repurposing with Heterogeneous Graph Transformers, Gharizadeh et al., 2024
  - PyG official hgt_dblp example (https://github.com/pyg-team/pytorch_geometric)
  - Heterogeneous GNNs for Link Prediction in Biomedical Networks, Hu et al., 2025

核心设计:
  1. 多组学异构图构建 (gene/drug/disease/pathway + PPI/coexp/TF/methylation)
  2. HGT编码器 + MLP解码器 链路预测
  3. Focal Loss处理类别不平衡 + SWA提升泛化
  4. 严格数据隔离: CV中剔除验证集gene-pathway边
  5. 内存优化: PPI采样 / PCA降维 / GradScaler
"""

import os
import sys
import copy
import random
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set

import numpy as np
import pandas as pd

import torch
import torch.nn.functional as F
from torch import Tensor, nn
from torch.optim import AdamW
from torch.optim.lr_scheduler import ReduceLROnPlateau

import torch_geometric
from torch_geometric.data import HeteroData
from torch_geometric.nn import HGTConv, Linear

from sklearn.model_selection import KFold
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.decomposition import PCA

import warnings
warnings.filterwarnings("ignore")


# ============================================================================
# 0. 全局配置
# ============================================================================

DATA_DIR = Path(r"D:\反向网络药理学\GAT拓展维度")
GAT_DATA_DIR = Path(r"C:\Users\Jy-Mentor-7\Desktop\GAT")
CACHE_DIR = DATA_DIR / "cache"
PROJECT_DIR = Path(r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙")

BRIDGE_GENES_PATH = PROJECT_DIR / "data" / "bridge_genes.csv"
ENHANCED_GENE_FEATURES_PATH = CACHE_DIR / "enhanced_gene_features.csv"
DRUG_FINGERPRINT_PATH = GAT_DATA_DIR / "drug_fingerprint.csv"
DISEASE_FEATURES_PATH = DATA_DIR / "disease_features.npy"
PATHWAY_FEATURES_PATH = DATA_DIR / "pathway_features.npy"
PPI_PATH = GAT_DATA_DIR / "ppi_subgraph.csv"
COEXP_PATH = DATA_DIR / "gene_coexp_edges.txt"
TF_TARGET_PATH = DATA_DIR / "tf_target_edges.txt"
GENE_PATHWAY_PATH = DATA_DIR / "gene_pathway_edges.txt"
SUBGRAPH_GENES_PATH = GAT_DATA_DIR / "subgraph_genes.txt"
METHYLATION_PATH = DATA_DIR / "gene_methylation_edges.txt"
MIRNA_PATH = None
PATHWAY_LIST_PATH = CACHE_DIR / "pathway_nodes.csv"

HIDDEN_DIM = 128
NUM_HEADS = 4
NUM_HGT_LAYERS = 2
DROPOUT = 0.3

LR = 1e-3
LR_PATIENCE = 20
LR_FACTOR = 0.5
WEIGHT_DECAY = 5e-4
EPOCHS = 200
PATIENCE = 50
NEG_SAMPLE_RATIO = 3
FOCAL_ALPHA = 0.5
FOCAL_GAMMA = 2.0
GRAD_CLIP_NORM = 1.0

N_FOLDS = 3
CV_RANDOM_STATE = 42
SWA_EPOCHS = 50
TOP_K = 10

SUBSAMPLE_PPI = True
PPI_MAX_EDGES = 80000
SUBSAMPLE_COEXP = True
USE_PCA_REDUCTION = True
GENE_FEATURE_DIM = 256
USE_GRAD_SCALING = True

DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

print(f"[Config] Device: {DEVICE}")
print(f"[Config] Hidden: {HIDDEN_DIM}, Heads: {NUM_HEADS}, Layers: {NUM_HGT_LAYERS}")
print(f"[Config] Epochs: {EPOCHS}, Patience: {PATIENCE}, NegRatio: {NEG_SAMPLE_RATIO}")
print(f"[Config] Focal: alpha={FOCAL_ALPHA}, gamma={FOCAL_GAMMA}")
print(f"[Config] LR={LR}, WD={WEIGHT_DECAY}, LR patience={LR_PATIENCE}, factor={LR_FACTOR}")
print(f"[Config] Memory: SUBSAMPLE_PPI={SUBSAMPLE_PPI} PPI_MAX={PPI_MAX_EDGES}")
print(f"[Config] PCA={USE_PCA_REDUCTION}(->{GENE_FEATURE_DIM})")


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


set_seed(CV_RANDOM_STATE)


# ============================================================================
# 特征降维
# ============================================================================

def maybe_reduce_features(arr: np.ndarray, target_dim: int) -> np.ndarray:
    if arr.shape[1] <= target_dim:
        return arr
    n_samples = arr.shape[0]
    if n_samples < target_dim:
        return arr[:, :n_samples].copy()
    pca = PCA(n_components=target_dim, random_state=CV_RANDOM_STATE)
    arr_r = pca.fit_transform(arr).astype(np.float32)
    print(f"  [PCA] {arr.shape[1]} -> {target_dim} (var={pca.explained_variance_ratio_.sum():.3f})")
    return arr_r


# ============================================================================
# 1. 数据加载
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


def load_gene_features(path: Path) -> Tuple[np.ndarray, List[str]]:
    df = pd.read_csv(path, index_col=0)
    df.index = df.index.astype(str).str.strip().str.upper()
    arr = df.apply(pd.to_numeric, errors="coerce").fillna(0.0).values.astype(np.float32)
    genes = df.index.tolist()
    print(f"[Load] gene features: {arr.shape}")
    return arr, genes


def load_drug_fingerprint(path: Path) -> np.ndarray:
    df = pd.read_csv(path)
    arr = df.apply(pd.to_numeric, errors="coerce").fillna(0).values.astype(np.float32)
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    print(f"[Load] drug fingerprint: {arr.shape}")
    return arr


def load_disease_features(path: Path) -> Optional[np.ndarray]:
    if not path.exists():
        print("[Load] disease_features.npy not found, will use gene mean")
        return None
    arr = np.load(str(path))
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    elif arr.ndim == 2 and arr.shape[0] > 1:
        arr = np.mean(arr, axis=0, keepdims=True)
    print(f"[Load] disease features: {arr.shape}")
    return arr


def load_pathway_features(path: Path) -> np.ndarray:
    arr = np.load(str(path))
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    print(f"[Load] pathway features: {arr.shape}")
    return arr


def load_ppi(
    path: Path,
    score_thresh: int = 700,
    bridge_set: Optional[Set[str]] = None,
    max_edges: int = 80000,
    subsample: bool = True,
) -> List[Tuple[str, str]]:
    df = pd.read_csv(path, sep="\t", encoding="utf-8")
    col_a, col_b = df.columns[0], df.columns[1]
    score_col = df.columns[2] if len(df.columns) >= 3 else None
    df = df.dropna(subset=[col_a, col_b])
    df[col_a] = df[col_a].astype(str).str.strip().str.upper()
    df[col_b] = df[col_b].astype(str).str.strip().str.upper()
    if score_col is not None:
        df[score_col] = pd.to_numeric(df[score_col], errors="coerce")
        df = df[df[score_col] >= score_thresh]
    if score_col is not None:
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
        print(f"  [PPI Subsample] bridge_edges={len(df_bridge)}, other_edges={len(df_other)}")
    elif subsample and len(df) > max_edges:
        df = df.head(max_edges)
    edges = list(zip(df[col_a], df[col_b]))
    print(f"[Load] PPI edges: {len(edges)} (score>={score_thresh})")
    return edges


def load_edge_list(path: Path, sep: Optional[str] = None, col_pair: Tuple[int, int] = (0, 1)) -> List[Tuple[str, str]]:
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
                        raw.append((parts[col_pair[0]].upper(), parts[col_pair[1]].upper()))
            break
        except (UnicodeDecodeError, OSError):
            continue
    print(f"[Load] edge list ({path.name}): {len(raw)}")
    return raw


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
    print(f"[Load] gene-pathway edges: {len(edges)}")
    return edges


def load_bridge_genes(path: Path) -> List[str]:
    if path.exists():
        df = pd.read_csv(path)
        col = "gene_symbol" if "gene_symbol" in df.columns else df.columns[0]
        genes = df[col].astype(str).str.strip().str.upper().tolist()
        print(f"[Load] bridge genes: {len(genes)}")
        return genes
    print(f"[Warn] bridge_genes.csv not found, trying all_bridge_genes.csv")
    alt = GAT_DATA_DIR / "all_bridge_genes.csv"
    if alt.exists():
        df = pd.read_csv(alt)
        genes = df["gene_symbol"].astype(str).str.strip().str.upper().tolist()
        print(f"[Load] bridge genes: {len(genes)} (all_bridge_genes)")
        return genes
    return []


# ============================================================================
# 2. 异构图构建
# ============================================================================

def build_hetero_graph(
    gene_feat_arr: np.ndarray,
    gene_feat_names: List[str],
    drug_fp_arr: np.ndarray,
    disease_feat_arr: Optional[np.ndarray],
    pathway_feat_arr: np.ndarray,
    pathway_names: List[str],
    ppi_edges: List[Tuple[str, str]],
    coexp_edges: List[Tuple[str, str]],
    tf_edges: List[Tuple[str, str]],
    gene_pathway_edges: List[Tuple[str, str]],
    methyl_edges: Optional[List[Tuple[str, str]]] = None,
    mirna_edges: Optional[List[Tuple[str, str]]] = None,
    all_genes_list: Optional[List[str]] = None,
    bridge_genes: Optional[List[str]] = None,
) -> Tuple[HeteroData, Dict[str, int], List[str], Dict[str, int]]:
    """构建多组学异构图，包含 gene/drug/disease/pathway 节点及多种边类型。"""
    gene_set = set(gene_feat_names)
    if all_genes_list:
        gene_set &= set(all_genes_list)

    for a, b in ppi_edges:
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
    print(f"[Build] Gene nodes: {n_genes}")

    gene_feat_dict = dict(zip(gene_feat_names, gene_feat_arr))
    gene_feat = np.zeros((n_genes, gene_feat_arr.shape[1]), dtype=np.float32)
    for i, g in enumerate(gene_list):
        if g in gene_feat_dict:
            gene_feat[i] = gene_feat_dict[g]

    data = HeteroData()
    data["gene"].x = torch.from_numpy(gene_feat)
    data["drug"].x = torch.from_numpy(drug_fp_arr.reshape(1, -1))

    if disease_feat_arr is not None:
        dis_feat = disease_feat_arr.reshape(1, -1)
    else:
        dis_feat = np.mean(gene_feat, axis=0, keepdims=True)
    data["disease"].x = torch.from_numpy(dis_feat.astype(np.float32))
    data["pathway"].x = torch.from_numpy(pathway_feat_arr.astype(np.float32))

    # PPI (无向)
    ppi_src, ppi_dst = [], []
    for a, b in ppi_edges:
        if a in gene_to_idx and b in gene_to_idx:
            ppi_src.append(gene_to_idx[a])
            ppi_dst.append(gene_to_idx[b])
    data["gene", "interacts", "gene"].edge_index = torch.tensor(
        [ppi_src + ppi_dst, ppi_dst + ppi_src], dtype=torch.long
    )
    print(f"[Build] PPI edges: {data['gene','interacts','gene'].edge_index.size(1)}")

    # 共表达 (无向)
    if SUBSAMPLE_COEXP and bridge_genes:
        bridge_set = set(bridge_genes)
        coe_filt = [(a, b) for a, b in coexp_edges if a in bridge_set or b in bridge_set]
        coe_edges_used = coe_filt if coe_filt else coexp_edges
    else:
        coe_edges_used = coexp_edges
    coe_src, coe_dst = [], []
    for a, b in coe_edges_used:
        if a in gene_to_idx and b in gene_to_idx:
            coe_src.append(gene_to_idx[a])
            coe_dst.append(gene_to_idx[b])
    if coe_src:
        data["gene", "coexpressed", "gene"].edge_index = torch.tensor(
            [coe_src + coe_dst, coe_dst + coe_src], dtype=torch.long
        )
    else:
        data["gene", "coexpressed", "gene"].edge_index = torch.zeros((2, 0), dtype=torch.long)
    print(f"[Build] Coexp edges: {data['gene','coexpressed','gene'].edge_index.size(1)} (filtered={len(coe_edges_used)})")

    # TF→靶基因 (有向)
    tf_src, tf_dst = [], []
    for a, b in tf_edges:
        if a in gene_to_idx and b in gene_to_idx:
            tf_src.append(gene_to_idx[a])
            tf_dst.append(gene_to_idx[b])
    if tf_src:
        data["gene", "regulates", "gene"].edge_index = torch.tensor([tf_src, tf_dst], dtype=torch.long)
    else:
        data["gene", "regulates", "gene"].edge_index = torch.zeros((2, 0), dtype=torch.long)
    print(f"[Build] TF edges: {len(tf_src)}")

    # 药物→基因 (有向)
    drug_targets = load_txt(GAT_DATA_DIR / "drug_targets.txt")
    dt_src, dt_dst = [], []
    for g in drug_targets:
        if g in gene_to_idx:
            dt_src.append(0)
            dt_dst.append(gene_to_idx[g])
    data["drug", "targets", "gene"].edge_index = (
        torch.tensor([dt_src, dt_dst], dtype=torch.long) if dt_src else torch.zeros((2, 0), dtype=torch.long)
    )
    print(f"[Build] drug->gene edges: {len(dt_src)}")

    # 基因→疾病 (有向)
    disease_genes = load_txt(GAT_DATA_DIR / "disease_genes.txt")
    gd_src, gd_dst = [], []
    for g in disease_genes:
        if g in gene_to_idx:
            gd_src.append(gene_to_idx[g])
            gd_dst.append(0)
    data["gene", "assoc_with", "disease"].edge_index = (
        torch.tensor([gd_src, gd_dst], dtype=torch.long) if gd_src else torch.zeros((2, 0), dtype=torch.long)
    )
    print(f"[Build] gene->disease edges: {len(gd_src)}")

    # 基因→通路 (有向, 监督任务)
    gp_src, gp_dst = [], []
    pathway_name_to_idx = {name: i for i, name in enumerate(pathway_names)}
    for a, b in gene_pathway_edges:
        if a in gene_to_idx and b in pathway_name_to_idx:
            gp_src.append(gene_to_idx[a])
            gp_dst.append(pathway_name_to_idx[b])
    data["gene", "involved_in", "pathway"].edge_index = (
        torch.tensor([gp_src, gp_dst], dtype=torch.long) if gp_src else torch.zeros((2, 0), dtype=torch.long)
    )
    print(f"[Build] gene->pathway edges: {len(gp_src)} (positive samples)")

    # 甲基化 (可选)
    if methyl_edges:
        cpg_set: Set[str] = set()
        for g, cpg in methyl_edges:
            if g in gene_to_idx:
                cpg_set.add(cpg)
        cpg_list = sorted(cpg_set)
        cpg_to_idx = {c: i for i, c in enumerate(cpg_list)}
        if cpg_list:
            data["cpg"].x = torch.randn(len(cpg_list), gene_feat.shape[1])
            gm_src, gm_dst = [], []
            for g, cpg in methyl_edges:
                if g in gene_to_idx and cpg in cpg_to_idx:
                    gm_src.append(gene_to_idx[g])
                    gm_dst.append(cpg_to_idx[cpg])
            data["gene", "methylated_at", "cpg"].edge_index = torch.tensor(
                [gm_src + gm_dst, gm_dst + gm_src], dtype=torch.long
            )
            print(f"[Build] Methylation edges: {data['gene','methylated_at','cpg'].edge_index.size(1)}, CpG: {len(cpg_list)}")

    print(f"[Build] Edge types: {list(data.edge_types)}")
    return data, gene_to_idx, gene_list, pathway_name_to_idx


# ============================================================================
# 3. 负采样 (参考 PyG negative_sampling 策略)
# ============================================================================

class NegEdgeSampler:
    def __init__(self, pos_edges: List[Tuple[int, int]], n_src: int, n_dst: int, seed: int = 42):
        self.n_src = n_src
        self.n_dst = n_dst
        self.exclude = {(int(s), int(d)) for s, d in pos_edges}
        self.rng = np.random.RandomState(seed)

    def sample(self, n_pos: int) -> Tensor:
        n_neg = n_pos * NEG_SAMPLE_RATIO
        src = self.rng.randint(0, self.n_src, size=n_neg * 2)
        dst = self.rng.randint(0, self.n_dst, size=n_neg * 2)
        neg_edges: List[Tuple[int, int]] = []
        for s, d in zip(src, dst):
            if len(neg_edges) >= n_neg:
                break
            key = (int(s), int(d))
            if key not in self.exclude:
                neg_edges.append(key)
                self.exclude.add(key)
        while len(neg_edges) < n_neg:
            s = int(self.rng.randint(0, self.n_src))
            d = int(self.rng.randint(0, self.n_dst))
            key = (s, d)
            if key not in self.exclude:
                neg_edges.append(key)
                self.exclude.add(key)
        return torch.tensor(
            [[e[0] for e in neg_edges], [e[1] for e in neg_edges]], dtype=torch.long
        )


# ============================================================================
# 4. HGT 模型 (参考 HGT WWW 2020 + PyG official hgt_dblp)
# ============================================================================

class HGTEncoder(nn.Module):
    """多层 HGTConv 编码器: 残差连接 + LayerNorm + ReLU + Dropout"""

    def __init__(self, metadata: Tuple[List[str], List[Tuple[str, str, str]]],
                 hidden_dim: int, num_heads: int, num_layers: int, dropout: float):
        super().__init__()
        self.metadata = metadata
        self.convs = nn.ModuleList()
        self.norms = nn.ModuleList()
        for _ in range(num_layers):
            self.convs.append(
                HGTConv(in_channels=hidden_dim, out_channels=hidden_dim,
                        metadata=metadata, heads=num_heads)
            )
            self.norms.append(nn.LayerNorm(hidden_dim))
        self.dropout = dropout

    def forward(self, x_dict: Dict[str, Tensor], edge_index_dict: Dict) -> Dict[str, Tensor]:
        for conv, norm in zip(self.convs, self.norms):
            x_dict_new = conv(x_dict, edge_index_dict)
            for k in x_dict:
                if k not in x_dict_new:
                    x_dict_new[k] = x_dict[k]
            x_dict = {
                k: F.relu(norm(x_dict_new[k]) + x_dict[k])
                for k in x_dict_new
            }
            x_dict = {
                k: F.dropout(v, p=self.dropout, training=self.training)
                for k, v in x_dict.items()
            }
        return x_dict


class HGTModel(nn.Module):
    """HGT 链路预测模型: 类型感知投影 + HGT编码器 + MLP解码器"""

    def __init__(self, metadata: Tuple[List[str], List[Tuple[str, str, str]]],
                 dim_dict: Dict[str, int], hidden_dim: int,
                 num_heads: int, num_layers: int, dropout: float):
        super().__init__()
        self.node_types = metadata[0]
        self.edge_types = metadata[1]

        self.proj = nn.ModuleDict()
        for nt, d_in in dim_dict.items():
            self.proj[nt] = Linear(d_in, hidden_dim)

        self.encoder = HGTEncoder(metadata, hidden_dim, num_heads, num_layers, dropout)

        self.decoder = nn.Sequential(
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

    def forward(self, x_dict: Dict[str, Tensor],
                edge_index_dict: Dict) -> Dict[str, Tensor]:
        x_proj = {}
        for k, v in x_dict.items():
            x_proj[k] = self.proj[k](v) if k in self.proj else v
        return self.encoder(x_proj, edge_index_dict)

    def decode(self, z_dict: Dict[str, Tensor], edge_index: Tensor) -> Tensor:
        z_gene = z_dict["gene"][edge_index[0]]
        z_path = z_dict["pathway"][edge_index[1]]
        return self.decoder(torch.cat([z_gene, z_path], dim=-1)).squeeze(-1)


# ============================================================================
# 5. Focal Loss (Lin et al., ICCV 2017)
# ============================================================================

def focal_bce_loss(logits: Tensor, labels: Tensor,
                   alpha: float = FOCAL_ALPHA, gamma: float = FOCAL_GAMMA) -> Tensor:
    logits = torch.clamp(logits, -10, 10)
    bce = F.binary_cross_entropy_with_logits(logits, labels, reduction="none")
    pt = torch.where(labels == 1, torch.sigmoid(logits), 1 - torch.sigmoid(logits))
    focal_weight = (1 - pt) ** gamma
    alpha_weight = torch.where(labels == 1, alpha, 1 - alpha)
    return (alpha_weight * focal_weight * bce).mean()


# ============================================================================
# 6. 子图工具: 从消息传递图中剔除指定边
# ============================================================================

def remove_edges_from_data(data: HeteroData, edge_type: Tuple[str, str, str],
                            mask: Tensor) -> HeteroData:
    data_out = data.clone()
    data_out[edge_type].edge_index = data_out[edge_type].edge_index[:, mask]
    return data_out


# ============================================================================
# 7. 训练与评估
# ============================================================================

@torch.no_grad()
def evaluate(model: nn.Module, data: HeteroData, gp_edge_index: Tensor,
             gp_pos_idx: Tensor, neg_sampler: NegEdgeSampler) -> Tuple[float, float, np.ndarray, np.ndarray]:
    model.eval()
    device = next(model.parameters()).device
    data_device = data if data.x_dict["gene"].device == device else data.to(device)

    z_dict = model(data_device.x_dict, data_device.edge_index_dict)

    pos_ei = gp_edge_index[:, gp_pos_idx].to(device)
    neg_ei = neg_sampler.sample(gp_pos_idx.shape[0]).to(device)

    eval_ei = torch.cat([pos_ei, neg_ei], dim=1)
    labels = np.concatenate([np.ones(pos_ei.size(1)), np.zeros(neg_ei.size(1))])
    logits = model.decode(z_dict, eval_ei)
    scores = torch.sigmoid(logits).cpu().numpy()

    auroc = roc_auc_score(labels, scores) if len(np.unique(labels)) > 1 else 0.5
    auprc = average_precision_score(labels, scores)
    return auroc, auprc, scores, labels


def train_fold(model: nn.Module, data: HeteroData, gp_edge_index: Tensor,
               train_idx: Tensor, val_idx: Tensor,
               neg_sampler: NegEdgeSampler) -> Tuple[float, float]:
    optimizer = AdamW(model.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)
    scheduler = ReduceLROnPlateau(optimizer, mode="max", factor=LR_FACTOR,
                                  patience=LR_PATIENCE, min_lr=1e-6)

    best_auroc = 0.0
    best_auprc = 0.0
    best_state = None
    patience_cnt = 0

    train_ei = gp_edge_index[:, train_idx]
    device = next(model.parameters()).device
    use_scaler = USE_GRAD_SCALING and DEVICE.type == "cuda"
    scaler = torch.amp.GradScaler(device="cuda") if use_scaler else None

    for epoch in range(EPOCHS):
        model.train()

        neg_ei = neg_sampler.sample(train_ei.shape[1]).to(device)
        batch_ei = torch.cat([train_ei.to(device), neg_ei], dim=1)
        batch_labels = torch.cat([
            torch.ones(train_ei.size(1), device=device),
            torch.zeros(neg_ei.size(1), device=device),
        ])

        perm = torch.randperm(batch_ei.size(1), device=device)
        batch_ei = batch_ei[:, perm]
        batch_labels = batch_labels[perm]

        data_device = data if data.x_dict["gene"].device == device else data.to(device)

        if scaler is not None:
            optimizer.zero_grad()
            with torch.amp.autocast(device_type="cuda"):
                z_dict = model(data_device.x_dict, data_device.edge_index_dict)
                logits = model.decode(z_dict, batch_ei)
                loss = focal_bce_loss(logits, batch_labels)
            scaler.scale(loss).backward()
            scaler.unscale_(optimizer)
            nn.utils.clip_grad_norm_(model.parameters(), GRAD_CLIP_NORM)
            scaler.step(optimizer)
            scaler.update()
        else:
            optimizer.zero_grad()
            z_dict = model(data_device.x_dict, data_device.edge_index_dict)
            logits = model.decode(z_dict, batch_ei)
            loss = focal_bce_loss(logits, batch_labels)
            loss.backward()
            nn.utils.clip_grad_norm_(model.parameters(), GRAD_CLIP_NORM)
            optimizer.step()

        if (epoch + 1) % 20 == 0:
            auroc, auprc, _, _ = evaluate(model, data, gp_edge_index, val_idx, neg_sampler)
            scheduler.step(auroc)

            if auroc > best_auroc:
                best_auroc = auroc
                best_auprc = auprc
                best_state = copy.deepcopy(model.state_dict())
                patience_cnt = 0
            else:
                patience_cnt += 1
                if patience_cnt >= PATIENCE:
                    if (epoch + 1) % 100 != 0:
                        print(f"    Early stop at epoch {epoch+1}, best AUROC={best_auroc:.4f}")
                    break

            if (epoch + 1) % 100 == 0:
                current_lr = optimizer.param_groups[0]["lr"]
                print(f"    Epoch {epoch+1} | Loss: {loss:.4f} | "
                      f"Val AUROC: {auroc:.4f} AUPRC: {auprc:.4f} | LR: {current_lr:.2e}")

        if DEVICE.type == "cuda" and (epoch + 1) % 50 == 0:
            torch.cuda.empty_cache()

    if best_state is not None:
        model.load_state_dict(best_state)
    return best_auroc, best_auprc


# ============================================================================
# 8. 交叉验证 (严格数据隔离)
# ============================================================================

def cross_validate(data: HeteroData, gp_edge_index: Tensor,
                   n_genes: int, n_pathways: int) -> List[Dict]:
    n_edges = gp_edge_index.size(1)
    if n_edges < N_FOLDS:
        print(f"[CV] Too few edges ({n_edges}) for {N_FOLDS}-fold CV")
        return []

    kf = KFold(n_splits=N_FOLDS, shuffle=True, random_state=CV_RANDOM_STATE)
    cv_scores: List[Dict] = []

    for fold, (train_idx, val_idx) in enumerate(kf.split(range(n_edges))):
        print(f"\n{'='*50}")
        print(f"[CV] Fold {fold+1}/{N_FOLDS}")
        print(f"{'='*50}")

        train_mask_np = np.ones(n_edges, dtype=bool)
        train_mask_np[val_idx] = False

        data_train = remove_edges_from_data(
            data, ("gene", "involved_in", "pathway"),
            torch.from_numpy(train_mask_np),
        )

        neg_sampler = NegEdgeSampler(
            pos_edges=gp_edge_index[:, train_idx].t().tolist(),
            n_src=n_genes, n_dst=n_pathways,
            seed=CV_RANDOM_STATE + fold,
        )

        model = HGTModel(
            metadata=data.metadata(),
            dim_dict={nt: data[nt].x.size(-1) for nt in data.node_types},
            hidden_dim=HIDDEN_DIM, num_heads=NUM_HEADS,
            num_layers=NUM_HGT_LAYERS, dropout=DROPOUT,
        ).to(DEVICE)

        gp_ei_device = gp_edge_index.to(DEVICE)

        auroc, auprc = train_fold(
            model, data_train.to(DEVICE), gp_ei_device,
            torch.from_numpy(train_idx).long(),
            torch.from_numpy(val_idx).long(),
            neg_sampler,
        )
        cv_scores.append({"fold": fold + 1, "auroc": auroc, "auprc": auprc})
        print(f"[CV] Fold {fold+1}: AUROC={auroc:.4f}, AUPRC={auprc:.4f}")

        if DEVICE.type == "cuda":
            torch.cuda.empty_cache()

    return cv_scores


# ============================================================================
# 9. 最终训练 (全量数据 + SWA)
# ============================================================================

def train_final(model: nn.Module, data: HeteroData, gp_edge_index: Tensor,
                neg_sampler: NegEdgeSampler) -> nn.Module:
    optimizer = AdamW(model.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)
    scheduler = ReduceLROnPlateau(optimizer, mode="max", factor=LR_FACTOR,
                                  patience=LR_PATIENCE, min_lr=1e-6)

    best_auroc = 0.0
    best_auprc = 0.0
    best_state = None
    patience_cnt = 0

    n_edges = gp_edge_index.size(1)
    val_size = max(1, int(n_edges * 0.1))
    perm = torch.randperm(n_edges, generator=torch.Generator().manual_seed(CV_RANDOM_STATE))
    train_idx = perm[val_size:]
    val_idx = perm[:val_size]

    val_mask = torch.zeros(n_edges, dtype=torch.bool)
    val_mask[val_idx] = True
    data_train = remove_edges_from_data(
        data, ("gene", "involved_in", "pathway"), ~val_mask,
    )

    device = next(model.parameters()).device
    data_train_device = data_train.to(device)
    gp_ei_device = gp_edge_index.to(device)

    swa_model: Optional[Dict] = None
    swa_n = 0

    use_scaler = USE_GRAD_SCALING and DEVICE.type == "cuda"
    scaler = torch.amp.GradScaler(device="cuda") if use_scaler else None

    for epoch in range(EPOCHS):
        model.train()

        train_ei = gp_ei_device[:, train_idx]
        neg_ei = neg_sampler.sample(train_ei.shape[1]).to(device)
        batch_ei = torch.cat([train_ei, neg_ei], dim=1)
        batch_labels = torch.cat([
            torch.ones(train_ei.size(1), device=device),
            torch.zeros(neg_ei.size(1), device=device),
        ])

        perm_i = torch.randperm(batch_ei.size(1), device=device)
        batch_ei = batch_ei[:, perm_i]
        batch_labels = batch_labels[perm_i]

        if scaler is not None:
            optimizer.zero_grad()
            with torch.amp.autocast(device_type="cuda"):
                z_dict = model(data_train_device.x_dict, data_train_device.edge_index_dict)
                logits = model.decode(z_dict, batch_ei)
                loss = focal_bce_loss(logits, batch_labels)
            scaler.scale(loss).backward()
            scaler.unscale_(optimizer)
            nn.utils.clip_grad_norm_(model.parameters(), GRAD_CLIP_NORM)
            scaler.step(optimizer)
            scaler.update()
        else:
            optimizer.zero_grad()
            z_dict = model(data_train_device.x_dict, data_train_device.edge_index_dict)
            logits = model.decode(z_dict, batch_ei)
            loss = focal_bce_loss(logits, batch_labels)
            loss.backward()
            nn.utils.clip_grad_norm_(model.parameters(), GRAD_CLIP_NORM)
            optimizer.step()

        if epoch >= EPOCHS - SWA_EPOCHS:
            if swa_model is None:
                swa_model = copy.deepcopy(model.state_dict())
                swa_n = 1
            else:
                swa_n += 1
                for key in swa_model:
                    if swa_model[key].dtype in (torch.float32, torch.float64):
                        swa_model[key].data += model.state_dict()[key].data

        if (epoch + 1) % 20 == 0:
            auroc, auprc, _, _ = evaluate(model, data, gp_ei_device, val_idx, neg_sampler)
            scheduler.step(auroc)
            if auroc > best_auroc:
                best_auroc = auroc
                best_auprc = auprc
                best_state = copy.deepcopy(model.state_dict())
                patience_cnt = 0
            else:
                patience_cnt += 1
                if patience_cnt >= PATIENCE:
                    if (epoch + 1) % 100 != 0:
                        print(f"  Early stop at epoch {epoch+1}")
                    break
            if (epoch + 1) % 100 == 0:
                current_lr = optimizer.param_groups[0]["lr"]
                print(f"  Epoch {epoch+1} | Loss: {loss:.4f} | "
                      f"Val AUROC: {auroc:.4f} AUPRC: {auprc:.4f} | LR: {current_lr:.2e}")

        if DEVICE.type == "cuda" and (epoch + 1) % 50 == 0:
            torch.cuda.empty_cache()

    if swa_model is not None and swa_n > 1:
        for key in swa_model:
            if swa_model[key].dtype in (torch.float32, torch.float64):
                swa_model[key].data.div_(swa_n)
        model.load_state_dict(swa_model)
        swa_auroc, swa_auprc, _, _ = evaluate(model, data, gp_ei_device, val_idx, neg_sampler)
        print(f"  SWA ({swa_n} epochs): AUROC={swa_auroc:.4f}, AUPRC={swa_auprc:.4f}")
        if swa_auroc > best_auroc:
            print(f"  Using SWA (best={best_auroc:.4f} vs SWA={swa_auroc:.4f})")
            best_auroc, best_auprc = swa_auroc, swa_auprc
        else:
            model.load_state_dict(best_state)
            print(f"  Using best checkpoint (best={best_auroc:.4f} vs SWA={swa_auroc:.4f})")
    elif best_state is not None:
        model.load_state_dict(best_state)

    print(f"  Final val AUROC={best_auroc:.4f}, AUPRC={best_auprc:.4f}")
    return model


# ============================================================================
# 10. 推理
# ============================================================================

@torch.no_grad()
def predict_bridge_pathways(model: nn.Module, data: HeteroData,
                            bridge_genes: List[str], gene_to_idx: Dict[str, int],
                            pathway_names: List[str], top_k: int = TOP_K) -> pd.DataFrame:
    model.eval()
    device = next(model.parameters()).device
    data_device = data if data.x_dict["gene"].device == device else data.to(device)

    z_dict = model(data_device.x_dict, data_device.edge_index_dict)

    results: List[Dict] = []
    n_pathways = len(pathway_names)
    for gene in bridge_genes:
        if gene not in gene_to_idx:
            continue
        gi = gene_to_idx[gene]
        edge_idx = torch.stack([
            torch.full((n_pathways,), gi, dtype=torch.long, device=device),
            torch.arange(n_pathways, dtype=torch.long, device=device),
        ])
        logits = model.decode(z_dict, edge_idx)
        scores = torch.sigmoid(logits).cpu().numpy()
        for pi, pname in enumerate(pathway_names):
            results.append({"gene_symbol": gene, "pathway_name": pname, "score": float(scores[pi])})

    df = pd.DataFrame(results)
    if df.empty:
        print("[Predict] No bridge genes found in graph!")
        return df

    df["rank"] = df.groupby("gene_symbol")["score"].rank(ascending=False, method="dense")
    df = df.sort_values(["gene_symbol", "rank"]).reset_index(drop=True)

    save_path = PROJECT_DIR / "bridge_pathway_scores.csv"
    df.to_csv(save_path, index=False)
    print(f"[Output] Saved bridge_pathway_scores.csv ({len(df)} rows)")

    topk = df[df["rank"] <= top_k].copy()
    print(f"\n{'='*60}")
    print(f"  Top-{top_k} Pathways per Bridge Gene")
    print(f"{'='*60}")
    for gene in bridge_genes:
        gene_rows = topk[topk["gene_symbol"] == gene]
        if gene_rows.empty:
            continue
        print(f"\n  {gene}:")
        for _, row in gene_rows.iterrows():
            pname = row["pathway_name"]
            pname_trim = pname[:60] + "..." if len(pname) > 60 else pname
            print(f"    #{int(row['rank']):2d} {pname_trim:60s} {row['score']:.4f}")

    return df


# ============================================================================
# 11. 主流程
# ============================================================================

def main() -> None:
    print("=" * 60)
    print("  HGT Gene-Pathway Association Prediction")
    print("  Network Pharmacology Stage 2: Pathway Discovery")
    print("=" * 60)

    print("\n[1] Loading data...")
    gene_feat_arr, gene_feat_names = load_gene_features(ENHANCED_GENE_FEATURES_PATH)
    if USE_PCA_REDUCTION and gene_feat_arr.shape[1] > GENE_FEATURE_DIM:
        gene_feat_arr = maybe_reduce_features(gene_feat_arr, GENE_FEATURE_DIM)

    pathway_feat_arr = load_pathway_features(PATHWAY_FEATURES_PATH)
    if USE_PCA_REDUCTION and pathway_feat_arr.shape[1] > GENE_FEATURE_DIM:
        pathway_feat_arr = maybe_reduce_features(pathway_feat_arr, GENE_FEATURE_DIM)

    drug_fp_arr = load_drug_fingerprint(DRUG_FINGERPRINT_PATH)
    disease_feat_arr = load_disease_features(DISEASE_FEATURES_PATH)
    pathway_names = load_pathway_list(PATHWAY_LIST_PATH)

    if len(pathway_names) != pathway_feat_arr.shape[0]:
        print(f"[Warn] pathway_names ({len(pathway_names)}) != features ({pathway_feat_arr.shape[0]}), using indices")
        pathway_names = [f"pathway_{i}" for i in range(pathway_feat_arr.shape[0])]

    bridge_genes = load_bridge_genes(BRIDGE_GENES_PATH)

    ppi_edges = load_ppi(PPI_PATH, bridge_set=set(bridge_genes) if bridge_genes else None)
    coexp_edges = load_edge_list(COEXP_PATH, sep="\t")
    tf_edges = load_edge_list(TF_TARGET_PATH, sep="\t")
    gene_pathway_edges = load_gene_pathway_edges(GENE_PATHWAY_PATH)
    all_genes_list = load_txt(SUBGRAPH_GENES_PATH)

    methyl_edges = load_edge_list(METHYLATION_PATH, sep=",", col_pair=(0, 1)) if METHYLATION_PATH and METHYLATION_PATH.exists() else None
    mirna_edges = load_edge_list(MIRNA_PATH) if MIRNA_PATH else None

    print(f"[Load] Bridge genes: {len(bridge_genes)}")

    print("\n[2] Building heterogeneous graph...")
    data, gene_to_idx, gene_list, pathway_name_to_idx = build_hetero_graph(
        gene_feat_arr, gene_feat_names,
        drug_fp_arr, disease_feat_arr,
        pathway_feat_arr, pathway_names,
        ppi_edges, coexp_edges, tf_edges,
        gene_pathway_edges,
        methyl_edges=methyl_edges,
        mirna_edges=mirna_edges,
        all_genes_list=all_genes_list,
        bridge_genes=bridge_genes,
    )

    gp_edge_index = data["gene", "involved_in", "pathway"].edge_index
    n_genes = data["gene"].x.size(0)
    n_pathways = data["pathway"].x.size(0)
    print(f"\n[Train] Gene-pathway edges: {gp_edge_index.size(1)}")
    print(f"[Train] Gene nodes: {n_genes}, Pathway nodes: {n_pathways}")

    if gp_edge_index.size(1) == 0:
        print("[FATAL] No gene-pathway edges found!")
        return

    print(f"\n[3] {N_FOLDS}-fold Cross Validation...")
    cv_results = cross_validate(data, gp_edge_index, n_genes, n_pathways)

    if cv_results:
        aurocs = [m["auroc"] for m in cv_results]
        auprcs = [m["auprc"] for m in cv_results]
        print(f"\n{'='*60}")
        print(f"  CV Results ({N_FOLDS}-fold, no leak)")
        print(f"  AUROC: {np.mean(aurocs):.4f} +/- {np.std(aurocs):.4f}")
        print(f"  AUPRC: {np.mean(auprcs):.4f} +/- {np.std(auprcs):.4f}")
        print(f"{'='*60}")
    else:
        print("[CV] Skipped (insufficient edges)")

    print(f"\n[4] Final training on all edges...")
    final_neg_sampler = NegEdgeSampler(
        pos_edges=gp_edge_index.t().tolist(),
        n_src=n_genes, n_dst=n_pathways,
        seed=CV_RANDOM_STATE + 999,
    )

    final_model = HGTModel(
        metadata=data.metadata(),
        dim_dict={nt: data[nt].x.size(-1) for nt in data.node_types},
        hidden_dim=HIDDEN_DIM, num_heads=NUM_HEADS,
        num_layers=NUM_HGT_LAYERS, dropout=DROPOUT,
    ).to(DEVICE)

    final_model = train_final(final_model, data, gp_edge_index.to(DEVICE), final_neg_sampler)

    model_path = PROJECT_DIR / "hgt_pathway_model.pt"
    torch.save(final_model.state_dict(), str(model_path))
    print(f"[Model] Saved to {model_path}")

    print(f"\n[5] Predicting pathway scores for bridge genes...")
    if bridge_genes:
        predict_bridge_pathways(
            final_model, data.to(DEVICE), bridge_genes,
            gene_to_idx, pathway_names, top_k=TOP_K,
        )
    else:
        print("[Skip] No bridge genes to score")

    print(f"\n{'='*60}")
    print("  Pipeline Complete!")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()