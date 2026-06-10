# -*- coding: utf-8 -*-
"""
HeCo-Style 异构图对比预训练模块

在链路预测微调之前，使用对比学习（HeCo, WWW 2022）预训练 MRHormer 编码器。

架构:
  - View 1 (Schema View): 标准前向传播 (原始图结构 + 原始特征)
  - View 2 (Network Schema View): 增强图 (边丢弃 + 特征掩码 + 节点丢弃)
  - 投影头: 每种节点类型的 MLP → 对比嵌入空间
  - 损失: 每节点类型 InfoNCE + 跨类型联合优化

用法:
  python hgt_pretrain_contrastive.py                       # 默认配置
  python hgt_pretrain_contrastive.py --epochs 300 --lr 1e-3 --tau 0.5

参考:
  - HeCo: Heterogeneous Graph Contrastive Learning, Wang et al., WWW 2022
  - SimCLR: A Simple Framework for Contrastive Learning, Chen et al., ICML 2020
  - MRHormer: TRLA + MRGA encoder from hgt_pathway_discovery.py
  - GCC: Generative & Contrastive Graph Learning, Cen et al., KDD 2022
"""

import sys
import copy
import math
import random
import os
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set, Union
import argparse
import json
import gc

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from torch import Tensor, nn
from torch.optim import AdamW
from torch.cuda.amp import GradScaler, autocast

from torch_geometric.data import HeteroData

import warnings
warnings.filterwarnings("ignore")

# ============================================================================
# 导入 MRHormer 模型及数据加载函数 (从单片模块导入)
# ============================================================================
# 注: 使用 importlib 加载 hgt_pathway_discovery.py 而非 hgt_pathway_discovery/
#     包，因为单片文件中包含 MRHormerModel (TRLA+MRGA) 架构。
#     hgt_pathway_discovery/ 目录包包含不同的 HGTModel 架构。
# ============================================================================
import importlib.util as _import_util

_HGT_MODULE_PATH = os.path.join(os.path.dirname(__file__), "hgt_pathway_discovery.py")
_HGT_SPEC = _import_util.spec_from_file_location("_hgt_pathway_discovery_mod", _HGT_MODULE_PATH)
_HGT_MODULE = _import_util.module_from_spec(_HGT_SPEC)
_HGT_SPEC.loader.exec_module(_HGT_MODULE)

# MRHormer 模型类
MRHormerModel = _HGT_MODULE.MRHormerModel
TRLAEncoder = _HGT_MODULE.TRLAEncoder
TRLAConv = _HGT_MODULE.TRLAConv

# 数据加载函数
load_gene_features = _HGT_MODULE.load_gene_features
load_drug_fingerprint = _HGT_MODULE.load_drug_fingerprint
load_ppi = _HGT_MODULE.load_ppi
load_txt = _HGT_MODULE.load_txt
load_edge_list = _HGT_MODULE.load_edge_list
load_pathway_list = _HGT_MODULE.load_pathway_list
maybe_reduce_features = _HGT_MODULE.maybe_reduce_features

# 路径和配置
GAT_DATA_DIR = _HGT_MODULE.GAT_DATA_DIR
DATA_DIR = _HGT_MODULE.DATA_DIR
BASE_DEVICE = _HGT_MODULE.DEVICE

# ============================================================================
# 0. 对比学习配置
# ============================================================================

# ---- 训练超参数 ----
CONTRASTIVE_EPOCHS = 200
CONTRASTIVE_LR = 1e-3
CONTRASTIVE_TAU = 0.5              # InfoNCE 温度
CONTRASTIVE_EDGE_DROPOUT = 0.3     # View 2 边丢弃概率
CONTRASTIVE_FEAT_MASK = 0.2        # View 2 特征掩码概率
CONTRASTIVE_HIDDEN_DIM = 64        # 编码器隐藏维度
CONTRASTIVE_NUM_HEADS = 2          # 注意力头数
CONTRASTIVE_NUM_LAYERS = 2         # TRLA 层数
CONTRASTIVE_BATCH_SIZE = 4096      # 对比损失节点批大小
CONTRASTIVE_GRAD_ACCUM = 4         # 梯度累积步数
CONTRASTIVE_WEIGHT_DECAY = 5e-4
CONTRASTIVE_DROPOUT = 0.2
CONTRASTIVE_WARMUP_EPOCHS = 10     # 学习率预热轮数
CONTRASTIVE_PROJ_DIM = 128         # 投影头输出维度
CONTRASTIVE_TAU_LEARNABLE = True   # 每节点类型可学习温度
CONTRASTIVE_NODE_DROPOUT = 0.0     # 节点丢弃概率 (默认关闭)
CONTRASTIVE_LOSS_WEIGHTS = None    # 每类型损失权重, None=等权

# ---- 路径 ----
PRETRAINED_WEIGHTS_PATH = Path(__file__).parent / "pretrained_contrastive_weights.pt"
PRETRAINED_CONFIG_PATH = Path(__file__).parent / "pretrained_contrastive_config.json"

# ---- 设备 ----
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
if DEVICE.type == "cuda":
    torch.backends.cuda.matmul.allow_tf32 = True
    torch.backends.cudnn.allow_tf32 = True

print(f"[Contrastive] Device: {DEVICE}")
print(f"[Contrastive] Epochs={CONTRASTIVE_EPOCHS}, LR={CONTRASTIVE_LR}, "
      f"Tau={CONTRASTIVE_TAU}, EdgeDrop={CONTRASTIVE_EDGE_DROPOUT}, "
      f"FeatMask={CONTRASTIVE_FEAT_MASK}, Hidden={CONTRASTIVE_HIDDEN_DIM}, "
      f"Heads={CONTRASTIVE_NUM_HEADS}, Layers={CONTRASTIVE_NUM_LAYERS}")


def set_seed(seed: int = 42) -> None:
    """全局确定性种子。"""
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


set_seed(42)


# ============================================================================
# 1. 异构图数据增强
# ============================================================================

class HeteroGraphAugmentation:
    """异构图数据增强模块（对比学习 View 2 使用）。

    提供三种增强策略:
      1. 边丢弃 (Edge Dropout): 每种边类型独立随机丢弃
      2. 特征掩码 (Feature Masking): 随机将部分特征置零
      3. 节点丢弃 (Node Dropout): 随机将部分节点特征全部置零

    参考:
      - HeCo (WWW 2022): 网络模式视图的增强策略
      - GraphCL (ICLR 2021): 图级对比学习的增强框架
    """

    def __init__(
        self,
        edge_dropout_p: float = 0.3,
        feat_mask_p: float = 0.2,
        node_dropout_p: float = 0.0,
        mask_token: float = 0.0,
        # 节点类型差异化增强: 每种类型独立概率
        # None 表示使用全局默认值
        per_type_edge_dropout: Optional[Dict[str, float]] = None,
        per_type_feat_mask: Optional[Dict[str, float]] = None,
    ):
        self.edge_dropout_p = edge_dropout_p
        self.feat_mask_p = feat_mask_p
        self.node_dropout_p = node_dropout_p
        self.mask_token = mask_token
        # 按节点类型差异化增强 (HeCo, WWW 2022):
        # 不同类型的拓扑和特征重要性不同, 应使用不同增强强度
        # e.g. gene 节点度大 → 高边丢弃; cpg 节点少 → 低特征掩码
        self.per_type_edge_dropout = per_type_edge_dropout or {}
        self.per_type_feat_mask = per_type_feat_mask or {}

    def edge_dropout(self, edge_index_dict: Dict) -> Dict:
        """对每种边类型随机丢弃边，保留至少一条边。
        
        支持按边类型差异化增强: e.g. PPI 边密度高可丢弃更多,
        gene_pathway 边稀疏应少丢弃。
        """
        dropped = {}
        for etype, ei in edge_index_dict.items():
            # 获取该边类型的丢弃概率 (per_type 优先, 否则全局默认)
            p = self.per_type_edge_dropout.get(etype, self.edge_dropout_p)
            if p <= 0:
                dropped[etype] = ei
                continue
            n_edges = ei.size(1)
            if n_edges <= 1:
                dropped[etype] = ei
                continue
            keep_mask = torch.rand(n_edges, device=ei.device) > p
            if keep_mask.sum() == 0:
                keep_mask[0] = True
            dropped[etype] = ei[:, keep_mask]
        return dropped

    def feature_masking(self, x_dict: Dict[str, Tensor]) -> Dict[str, Tensor]:
        """对每种节点类型的特征随机掩码。
        
        支持按节点类型差异化: e.g. drug/FP 特征稀疏应少掩码,
        gene/RNA-seq 特征冗余可多掩码。
        """
        masked = {}
        for nt, x in x_dict.items():
            p = self.per_type_feat_mask.get(nt, self.feat_mask_p)
            if p <= 0:
                masked[nt] = x
                continue
            mask = torch.rand_like(x) < p
            x_aug = x.clone()
            x_aug[mask] = self.mask_token
            masked[nt] = x_aug
        return masked

    def node_dropout(self, x_dict: Dict[str, Tensor]) -> Dict[str, Tensor]:
        """对每种节点类型随机丢弃节点（特征全部置零）。"""
        if self.node_dropout_p <= 0:
            return x_dict
        dropped = {}
        for nt, x in x_dict.items():
            mask = torch.rand(x.size(0), 1, device=x.device) < self.node_dropout_p
            x_aug = x.clone()
            x_aug[mask.expand_as(x)] = 0.0
            dropped[nt] = x_aug
        return dropped

    def __call__(
        self,
        x_dict: Dict[str, Tensor],
        edge_index_dict: Dict,
    ) -> Tuple[Dict[str, Tensor], Dict]:
        """对 View 2 应用所有增强。

        Args:
            x_dict: 原始节点特征
            edge_index_dict: 原始边索引

        Returns:
            x_aug: 增强后的节点特征
            edge_aug: 增强后的边索引
        """
        x_aug = self.feature_masking(x_dict)
        x_aug = self.node_dropout(x_aug)
        edge_aug = self.edge_dropout(edge_index_dict)
        return x_aug, edge_aug


# ============================================================================
# 2. 对比学习模型
# ============================================================================

class ContrastiveMRHormer(nn.Module):
    """HeCo 风格异构图对比学习模块。

    架构:
      - encoder: MRHormerModel (TRLA + MRGA, 权重在两个视图间共享)
      - proj_heads: 每种节点类型的 MLP 投影头 (hidden_dim → proj_dim)
      - logit_scale: 每节点类型可学习温度参数 (或固定温度)

    输入:
      - x_dict, edge_index_dict: View 1 (Schema View, 原始图)
      - x_aug, edge_aug:         View 2 (Network Schema View, 增强图)

    输出:
      - z1_dict, z2_dict: 两个视图的对比嵌入 (投影后, L2 归一化)

    参考:
      - HeCo: 双视图异构图对比学习 (WWW 2022)
      - CLIP: 可学习 logit scale 用于温度自适应 (ICML 2021)
    """

    def __init__(
        self,
        encoder: MRHormerModel,
        node_types: List[str],
        proj_dim: int = CONTRASTIVE_PROJ_DIM,
        tau: float = CONTRASTIVE_TAU,
        tau_learnable: bool = CONTRASTIVE_TAU_LEARNABLE,
    ):
        super().__init__()
        self.encoder = encoder
        self.node_types = [nt for nt in node_types if nt in dict(encoder.named_children()) or True]
        self.hidden_dim = encoder.hidden_dim
        self.proj_dim = proj_dim
        self.tau_fixed = tau

        # 投影头: 每种节点类型独立 MLP (hidden_dim → hidden_dim → proj_dim)
        self.proj_heads = nn.ModuleDict()
        for nt in node_types:
            self.proj_heads[nt] = nn.Sequential(
                nn.Linear(self.hidden_dim, self.hidden_dim),
                nn.BatchNorm1d(self.hidden_dim),
                nn.ReLU(inplace=True),
                nn.Linear(self.hidden_dim, proj_dim),
            )

        # 可学习温度参数 (每节点类型独立)
        if tau_learnable:
            self.logit_scale = nn.ParameterDict()
            for nt in node_types:
                # 初始化为 log(1/tau), 使初始温度接近 tau
                init_log_scale = math.log(1.0 / max(tau, 0.1))
                self.logit_scale[nt] = nn.Parameter(torch.tensor(init_log_scale))
        else:
            self.logit_scale = None
        self.tau_learnable = tau_learnable

        self._init_weights()

    def _init_weights(self) -> None:
        """投影头权重初始化。"""
        for name, module in self.proj_heads.named_modules():
            if isinstance(module, nn.Linear):
                nn.init.xavier_uniform_(module.weight, gain=nn.init.calculate_gain("relu"))
                if module.bias is not None:
                    nn.init.zeros_(module.bias)
            elif isinstance(module, nn.BatchNorm1d):
                nn.init.ones_(module.weight)
                nn.init.zeros_(module.bias)

    def forward(
        self,
        x_dict: Dict[str, Tensor],
        edge_index_dict: Dict,
        x_aug: Dict[str, Tensor],
        edge_aug: Dict,
        node_types_to_use: Optional[List[str]] = None,
    ) -> Tuple[Dict[str, Tensor], Dict[str, Tensor]]:
        """前向传播：计算两个视图的对比嵌入。

        Args:
            x_dict:          View 1 节点特征
            edge_index_dict: View 1 边索引
            x_aug:           View 2 节点特征 (增强后)
            edge_aug:        View 2 边索引 (增强后)
            node_types_to_use: 仅计算的节点类型子集 (None=全部)

        Returns:
            z1_dict: View 1 投影嵌入 {node_type: Tensor[N, proj_dim]}
            z2_dict: View 2 投影嵌入 {node_type: Tensor[N, proj_dim]}
        """
        # View 1: Schema view (原始图 → MRHormer 编码)
        z1_raw = self.encoder(x_dict, edge_index_dict)

        # View 2: Network schema view (增强图 → MRHormer 编码)
        z2_raw = self.encoder(x_aug, edge_aug)

        # 投影头
        types_to_use = node_types_to_use if node_types_to_use else self.node_types
        z1_dict = {}
        z2_dict = {}
        for nt in types_to_use:
            if nt in z1_raw and nt in z2_raw:
                z1_proj = self.proj_heads[nt](z1_raw[nt])
                z2_proj = self.proj_heads[nt](z2_raw[nt])
                # L2 归一化 (对比学习中至关重要)
                z1_dict[nt] = F.normalize(z1_proj, dim=-1)
                z2_dict[nt] = F.normalize(z2_proj, dim=-1)

        return z1_dict, z2_dict

    def get_tau(self, node_type: str) -> float:
        """获取指定节点类型的有效温度参数。"""
        if self.tau_learnable and self.logit_scale is not None and node_type in self.logit_scale:
            # logit_scale = log(1/tau) → tau = 1/exp(logit_scale)
            # 钳制 logit_scale ≥ 0 (即 tau ≤ 1), 避免温度过大导致梯度消失
            return 1.0 / torch.exp(self.logit_scale[node_type].clamp(min=0.0))
        return self.tau_fixed

    @torch.no_grad()
    def get_encoder_weights(self) -> Dict[str, Tensor]:
        """提取编码器权重，用于加载到下游模型。"""
        return self.encoder.state_dict()

    def update_encoder_drop_edge(self, drop_edge_p: float) -> None:
        """更新编码器的 DropEdge 概率（微调时使用）。"""
        if hasattr(self.encoder, 'encoder') and hasattr(self.encoder.encoder, 'drop_edge_p'):
            self.encoder.encoder.drop_edge_p = drop_edge_p

    def get_pretrained_config(self) -> Dict:
        """返回预训练配置字典。"""
        return {
            "hidden_dim": self.hidden_dim,
            "proj_dim": self.proj_dim,
            "node_types": self.node_types,
            "tau_fixed": self.tau_fixed,
            "tau_learnable": self.tau_learnable,
            "model_architecture": "ContrastiveMRHormer",
        }


# ============================================================================
# 3. InfoNCE 对比损失
# ============================================================================

def contrastive_info_nce(
    z1: Tensor,
    z2: Tensor,
    tau: float = 0.5,
) -> Tensor:
    """InfoNCE 对比损失（每节点类型独立计算）。

    公式:
      L = -1/N * Σ_i log( exp(sim(z1_i, z2_i)/τ) / Σ_j exp(sim(z1_i, z2_j)/τ) )

    其中 sim(a, b) = a · b (余弦相似度, 输入已 L2 归一化).

    Args:
        z1: View 1 嵌入, shape (N, D), 已 L2 归一化
        z2: View 2 嵌入, shape (N, D), 已 L2 归一化
        tau: 温度参数

    Returns:
        标量损失

    参考:
      - SimCLR InfoNCE (Chen et al., ICML 2020)
      - HeCo 使用相同公式 (Wang et al., WWW 2022)
    """
    N = z1.size(0)
    if N <= 1:
        return torch.tensor(0.0, device=z1.device, requires_grad=True)

    # 相似度矩阵: S[i][j] = z1[i] · z2[j] (已 L2 归一化)
    sim = torch.mm(z1, z2.t()) / tau  # (N, N)

    # 数值稳定性: 减去每行最大值
    sim_max = sim.max(dim=1, keepdim=True)[0].detach()
    sim = sim - sim_max

    # InfoNCE: L_i = -log(exp(S[i,i]) / Σ_j exp(S[i,j]))
    #         = -(S[i,i] - log(Σ_j exp(S[i,j])))
    pos_sim = sim.diag()  # (N,)
    log_sum_exp = torch.logsumexp(sim, dim=1)  # (N,)
    loss_per_node = -(pos_sim - log_sum_exp)

    return loss_per_node.mean()


def contrastive_info_nce_symmetric(
    z1: Tensor,
    z2: Tensor,
    tau: float = 0.5,
) -> Tensor:
    """对称 InfoNCE: L = (L(z1→z2) + L(z2→z1)) / 2.

    参考: HeCo 使用对称对比损失增强双向对比信号.
    """
    loss_fwd = contrastive_info_nce(z1, z2, tau)
    loss_bwd = contrastive_info_nce(z2, z1, tau)
    return (loss_fwd + loss_bwd) / 2.0


def compute_contrastive_loss(
    z1_dict: Dict[str, Tensor],
    z2_dict: Dict[str, Tensor],
    model: ContrastiveMRHormer,
    loss_weights: Optional[Dict[str, float]] = None,
    symmetric: bool = True,
) -> Tuple[Tensor, Dict[str, Tensor]]:
    """计算整个异构图的对比损失（所有节点类型加权求和）。

    Args:
        z1_dict: View 1 嵌入
        z2_dict: View 2 嵌入
        model: 对比学习模型 (用于获取每类型的温度)
        loss_weights: 每类型损失权重, None=等权
        symmetric: 是否使用对称 InfoNCE

    Returns:
        total_loss: 加权总损失
        per_type_loss: 每类型损失字典
    """
    per_type_loss = {}
    valid_types = []

    for nt in z1_dict:
        if nt not in z2_dict:
            continue
        z1, z2 = z1_dict[nt], z2_dict[nt]
        if z1.size(0) <= 1 or z2.size(0) <= 1:
            continue

        tau = model.get_tau(nt)
        if symmetric:
            loss_nt = contrastive_info_nce_symmetric(z1, z2, tau)
        else:
            loss_nt = contrastive_info_nce(z1, z2, tau)

        per_type_loss[nt] = loss_nt.detach()
        valid_types.append(nt)

    if not valid_types:
        return torch.tensor(0.0, device=DEVICE, requires_grad=True), {}

    # 加权求和
    if loss_weights is None:
        total_loss = torch.stack([per_type_loss[nt] for nt in valid_types]).mean()
    else:
        weights = torch.tensor(
            [loss_weights.get(nt, 1.0) for nt in valid_types],
            device=DEVICE,
        )
        weights = weights / weights.sum()
        losses = torch.stack([per_type_loss[nt] for nt in valid_types])
        total_loss = (weights * losses).sum()

    return total_loss, per_type_loss


# ============================================================================
# 4. 数据加载与图构建
# ============================================================================

def load_pretrain_data(
    use_full_graph: bool = True,
) -> Tuple[HeteroData, List[str], List[Tuple[str, str, str]], Dict[str, int]]:
    """加载预训练数据：构建异构图用于对比学习。

    使用与 hgt_pathway_discovery.py 相同的数据加载函数，
    构建包含 gene/drug/disease/pathway/cpg 的多组学异构图。

    Args:
        use_full_graph: 是否加载完整图 (含 pathway/cpg)
                        False 则仅使用 drug/gene/disease (快速调试)

    Returns:
        data: HeteroData 异构图对象
        node_types: 节点类型列表
        edge_types: 边类型列表
        dim_dict: 每节点类型的输入维度
    """
    print(f"\n{'='*60}")
    print(f"  [Data] 加载预训练数据...")
    print(f"{'='*60}")

    # ---- 加载基因特征 ----
    gene_path = GAT_DATA_DIR / "subgraph_embeddings.csv"
    gene_feat_arr, gene_feat_names = load_gene_features(gene_path)

    # ---- PCA 降维 ----
    gene_feat_arr = maybe_reduce_features(gene_feat_arr, CONTRASTIVE_HIDDEN_DIM)

    # ---- 药物指纹 ----
    drug_fp_arr = load_drug_fingerprint(GAT_DATA_DIR / "drug_fingerprint.csv")

    # ---- PPI 边 ----
    ppi_edges = load_ppi(
        GAT_DATA_DIR / "ppi_subgraph.csv",
        score_thresh=700,
        max_edges=40000,
        subsample=True,
    )

    # ---- 药物靶点和疾病基因 ----
    drug_targets = load_txt(GAT_DATA_DIR / "drug_targets.txt")
    disease_genes = load_txt(GAT_DATA_DIR / "disease_genes.txt")
    subgraph_genes = load_txt(GAT_DATA_DIR / "subgraph_genes.txt")
    all_gene_set = set(gene_feat_names)

    # ---- 确定基因全集 ----
    if subgraph_genes:
        all_gene_set &= set(subgraph_genes)
    for a, b in ppi_edges:
        all_gene_set.add(a)
        all_gene_set.add(b)
    for g in drug_targets:
        all_gene_set.add(g)
    for g in disease_genes:
        all_gene_set.add(g)

    gene_list = sorted(all_gene_set & set(gene_feat_names))
    if not gene_list:
        raise ValueError("No genes with features in the graph!")
    gene_to_idx = {g: i for i, g in enumerate(gene_list)}
    n_genes = len(gene_list)
    print(f"[Data] Gene nodes: {n_genes}")

    # ---- 基因特征矩阵 ----
    gene_feat_dict_map = dict(zip(gene_feat_names, gene_feat_arr))
    gene_feat = np.zeros((n_genes, gene_feat_arr.shape[1]), dtype=np.float32)
    for i, g in enumerate(gene_list):
        if g in gene_feat_dict_map:
            gene_feat[i] = gene_feat_dict_map[g]

    # ---- 构建 HeteroData ----
    data = HeteroData()
    data["gene"].x = torch.from_numpy(gene_feat)
    data["drug"].x = torch.from_numpy(drug_fp_arr.reshape(1, -1))

    # 疾病特征: 使用疾病基因的均值
    disease_in_graph = [g for g in disease_genes if g in gene_to_idx]
    if disease_in_graph:
        dis_indices = [gene_to_idx[g] for g in disease_in_graph]
        disease_feat = np.mean(gene_feat[dis_indices], axis=0, keepdims=True)
    else:
        disease_feat = np.mean(gene_feat, axis=0, keepdims=True)
    data["disease"].x = torch.from_numpy(disease_feat.astype(np.float32))

    node_types = ["gene", "drug", "disease"]

    # ---- 共表达边 ----
    coexp_edges = load_edge_list(DATA_DIR / "gene_coexp_edges.txt")
    # TF-靶基因边
    tf_edges = load_edge_list(DATA_DIR / "tf_target_edges.txt")
    # 基因-通路边
    gene_pathway_edges = load_edge_list(DATA_DIR / "gene_pathway_edges.txt")

    if use_full_graph:
        # ---- 通路特征 ----
        pathway_path = DATA_DIR / "pathway_features.npy"
        if pathway_path.exists():
            pathway_feat_arr = np.load(str(pathway_path))
            if pathway_feat_arr.ndim == 1:
                pathway_feat_arr = pathway_feat_arr.reshape(1, -1)
            pathway_feat_arr = maybe_reduce_features(
                pathway_feat_arr, CONTRASTIVE_HIDDEN_DIM
            )
        else:
            # 回退: 使用基因特征均值
            n_pathways = len(load_pathway_list(DATA_DIR / "pathway_nodes.csv"))
            pathway_feat_arr = np.zeros(
                (max(n_pathways, 1), gene_feat.shape[1]), dtype=np.float32
            )
            pathway_feat_arr[:] = gene_feat.mean(axis=0)

        # 通路名称
        pathway_names = load_pathway_list(DATA_DIR / "pathway_nodes.csv")
        if not pathway_names:
            pathway_names = [f"pathway_{i}" for i in range(pathway_feat_arr.shape[0])]

        data["pathway"].x = torch.from_numpy(pathway_feat_arr.astype(np.float32))
        node_types.append("pathway")

        # ---- 甲基化边 (可选) ----
        methyl_edges = load_edge_list(DATA_DIR / "gene_methylation_edges.txt")

        # CpG 节点 (特征传播初始化)
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
                data["cpg"].x = torch.from_numpy(cpg_feat.astype(np.float32))
                node_types.append("cpg")

    # ---- 边构建 ----
    # PPI (无向)
    ppi_src, ppi_dst = [], []
    for a, b in ppi_edges:
        if a in gene_to_idx and b in gene_to_idx:
            ppi_src.append(gene_to_idx[a])
            ppi_dst.append(gene_to_idx[b])
    # 双向 PPI
    data["gene", "interacts", "gene"].edge_index = torch.tensor(
        [ppi_src + ppi_dst, ppi_dst + ppi_src], dtype=torch.long
    )
    print(f"[Data] PPI edges: {data['gene','interacts','gene'].edge_index.size(1)}")

    # 共表达 (无向)
    if coexp_edges:
        coe_src, coe_dst = [], []
        for a, b in coexp_edges:
            if a in gene_to_idx and b in gene_to_idx:
                coe_src.append(gene_to_idx[a])
                coe_dst.append(gene_to_idx[b])
        if coe_src:
            data["gene", "coexpressed", "gene"].edge_index = torch.tensor(
                [coe_src + coe_dst, coe_dst + coe_src], dtype=torch.long
            )
    if ("gene", "coexpressed", "gene") not in data.edge_types:
        data["gene", "coexpressed", "gene"].edge_index = torch.zeros((2, 0), dtype=torch.long)

    # TF 调控 (有向)
    if tf_edges:
        tf_src, tf_dst = [], []
        for a, b in tf_edges:
            if a in gene_to_idx and b in gene_to_idx:
                tf_src.append(gene_to_idx[a])
                tf_dst.append(gene_to_idx[b])
        if tf_src:
            data["gene", "regulates", "gene"].edge_index = torch.tensor(
                [tf_src, tf_dst], dtype=torch.long
            )
    if ("gene", "regulates", "gene") not in data.edge_types:
        data["gene", "regulates", "gene"].edge_index = torch.zeros((2, 0), dtype=torch.long)

    # 药物→基因 (有向)
    dt_src, dt_dst = [], []
    for g in drug_targets:
        if g in gene_to_idx:
            dt_src.append(0)
            dt_dst.append(gene_to_idx[g])
    data["drug", "targets", "gene"].edge_index = (
        torch.tensor([dt_src, dt_dst], dtype=torch.long) if dt_src
        else torch.zeros((2, 0), dtype=torch.long)
    )
    print(f"[Data] drug→gene edges: {len(dt_src)}")

    # 基因→疾病 (有向)
    gd_src, gd_dst = [], []
    for g in disease_genes:
        if g in gene_to_idx:
            gd_src.append(gene_to_idx[g])
            gd_dst.append(0)
    data["gene", "assoc_with", "disease"].edge_index = (
        torch.tensor([gd_src, gd_dst], dtype=torch.long) if gd_src
        else torch.zeros((2, 0), dtype=torch.long)
    )
    print(f"[Data] gene→disease edges: {len(gd_src)}")

    # 反向边: gene → targeted_by → drug
    if dt_src:
        data["gene", "targeted_by", "drug"].edge_index = torch.tensor(
            [dt_dst, dt_src], dtype=torch.long
        )

    if use_full_graph:
        # 基因→通路 (有向)
        if gene_pathway_edges:
            pathway_name_to_idx = {name: i for i, name in enumerate(pathway_names)}
            gp_src, gp_dst = [], []
            for a, b in gene_pathway_edges:
                if a in gene_to_idx and b in pathway_name_to_idx:
                    gp_src.append(gene_to_idx[a])
                    gp_dst.append(pathway_name_to_idx[b])
            data["gene", "involved_in", "pathway"].edge_index = (
                torch.tensor([gp_src, gp_dst], dtype=torch.long) if gp_src
                else torch.zeros((2, 0), dtype=torch.long)
            )
            print(f"[Data] gene→pathway edges: {len(gp_src)}")
        if ("gene", "involved_in", "pathway") not in data.edge_types:
            data["gene", "involved_in", "pathway"].edge_index = torch.zeros((2, 0), dtype=torch.long)

        # 甲基化边 (双向)
        if methyl_edges and cpg_list:
            gm_src, gm_dst = [], []
            for g, cpg in methyl_edges:
                if g in gene_to_idx and cpg in cpg_to_idx:
                    gm_src.append(gene_to_idx[g])
                    gm_dst.append(cpg_to_idx[cpg])
            if gm_src:
                data["gene", "methylated_at", "cpg"].edge_index = torch.tensor(
                    [gm_src + gm_dst, gm_dst + gm_src], dtype=torch.long
                )
                print(f"[Data] methylation edges: {data['gene','methylated_at','cpg'].edge_index.size(1)}")
        if ("gene", "methylated_at", "cpg") not in data.edge_types:
            data["gene", "methylated_at", "cpg"].edge_index = torch.zeros((2, 0), dtype=torch.long)

    # ---- 边类型列表 ----
    edge_types = list(data.edge_types)
    print(f"[Data] Edge types: {edge_types}")

    # ---- 每类型维度 ----
    dim_dict = {}
    for nt in data.node_types:
        dim_dict[nt] = data[nt].x.size(1)
    print(f"[Data] Node types: {list(data.node_types)}, Dims: {dim_dict}")

    print(f"{'='*60}\n")
    return data, node_types, edge_types, dim_dict


# ============================================================================
# 5. 预训练入口
# ============================================================================

def pretrain(
    data: HeteroData,
    node_types: List[str],
    edge_types: List[Tuple[str, str, str]],
    dim_dict: Dict[str, int],
    epochs: int = CONTRASTIVE_EPOCHS,
    lr: float = CONTRASTIVE_LR,
    tau: float = CONTRASTIVE_TAU,
    edge_dropout_p: float = CONTRASTIVE_EDGE_DROPOUT,
    feat_mask_p: float = CONTRASTIVE_FEAT_MASK,
    hidden_dim: int = CONTRASTIVE_HIDDEN_DIM,
    num_heads: int = CONTRASTIVE_NUM_HEADS,
    num_layers: int = CONTRASTIVE_NUM_LAYERS,
    batch_size: int = CONTRASTIVE_BATCH_SIZE,
    grad_accum: int = CONTRASTIVE_GRAD_ACCUM,
    weight_decay: float = CONTRASTIVE_WEIGHT_DECAY,
    dropout: float = CONTRASTIVE_DROPOUT,
    warmup_epochs: int = CONTRASTIVE_WARMUP_EPOCHS,
    proj_dim: int = CONTRASTIVE_PROJ_DIM,
    tau_learnable: bool = CONTRASTIVE_TAU_LEARNABLE,
    loss_weights: Optional[Dict[str, float]] = None,
    save_path: Path = PRETRAINED_WEIGHTS_PATH,
    config_save_path: Path = PRETRAINED_CONFIG_PATH,
    verbose: bool = True,
) -> ContrastiveMRHormer:
    """执行 HeCo 风格对比预训练。

    流程:
      1. 创建 MRHormer 编码器 + 投影头
      2. 每轮: 对 View 2 随机增强 → 双视图编码 → InfoNCE 损失
      3. 保存最佳检查点 (最低对比损失)

    Args:
        data: HeteroData 异构图
        node_types: 节点类型列表
        edge_types: 边类型列表
        dim_dict: 每节点类型的特征维度
        epochs: 训练轮数
        lr: 学习率
        tau: InfoNCE 温度
        edge_dropout_p: View 2 边丢弃概率
        feat_mask_p: View 2 特征掩码概率
        hidden_dim: 编码器隐藏维度
        num_heads: 注意力头数
        num_layers: TRLA 层数
        batch_size: 对比损失节点批大小 (当前为全图, 保留参数)
        grad_accum: 梯度累积步数
        weight_decay: 权重衰减
        dropout: Dropout 率
        warmup_epochs: 学习率预热轮数
        proj_dim: 投影头输出维度
        tau_learnable: 是否使用可学习温度
        loss_weights: 每节点类型损失权重
        save_path: 权重保存路径
        config_save_path: 配置保存路径
        verbose: 是否打印训练日志

    Returns:
        训练好的 ContrastiveMRHormer 模型
    """
    print(f"\n{'='*60}")
    print(f"  HeCo 风格对比预训练")
    print(f"{'='*60}")
    print(f"  编码器: MRHormer (TRLA+MRGA)")
    print(f"  节点类型: {node_types}")
    print(f"  边类型: {edge_types}")
    print(f"  轮数: {epochs}, LR: {lr}, Tau: {tau}")
    print(f"  增强: EdgeDrop={edge_dropout_p}, FeatMask={feat_mask_p}")
    print(f"  维度: Hidden={hidden_dim}, Proj={proj_dim}, Heads={num_heads}")
    print(f"{'='*60}\n")

    # ---- 初始化编码器 ----
    encoder = MRHormerModel(
        node_types=node_types,
        edge_types=edge_types,
        dim_dict=dim_dict,
        hidden_dim=hidden_dim,
        num_heads=num_heads,
        num_layers=num_layers,
        dropout=dropout,
        initial_residual=True,
        drop_edge_p=0.0,  # 对比学习自行处理边丢弃
        use_mrga=True,
        use_input_bn=True,
    ).to(DEVICE)

    # ---- 初始化对比模型 ----
    model = ContrastiveMRHormer(
        encoder=encoder,
        node_types=node_types,
        proj_dim=proj_dim,
        tau=tau,
        tau_learnable=tau_learnable,
    ).to(DEVICE)

    # ---- 优化器 & 调度器 ----
    optimizer = AdamW(model.parameters(), lr=lr, weight_decay=weight_decay)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
        optimizer, T_max=epochs, eta_min=lr * 0.01,
    )

    # ---- 数据增强器 ----
    augmentor = HeteroGraphAugmentation(
        edge_dropout_p=edge_dropout_p,
        feat_mask_p=feat_mask_p,
        node_dropout_p=0.0,
    )

    # ---- 准备数据 ----
    data = data.to(DEVICE)
    x_dict = {nt: data[nt].x for nt in data.node_types}
    edge_index_dict = {et: data[et].edge_index for et in data.edge_types}

    # 自动混合精度
    scaler = GradScaler() if DEVICE.type == "cuda" else None
    use_amp = scaler is not None

    # ---- 训练循环 ----
    best_loss = float("inf")
    best_state_dict = None
    patience_counter = 0
    patience = max(30, epochs // 5)

    for epoch in range(1, epochs + 1):
        model.train()

        # ---- 学习率预热 ----
        if epoch <= warmup_epochs:
            warmup_lr = lr * epoch / warmup_epochs
            for param_group in optimizer.param_groups:
                param_group["lr"] = warmup_lr

        # ---- View 2: 数据增强 ----
        with torch.no_grad() if not model.training else torch.enable_grad():
            x_aug, edge_aug = augmentor(x_dict, edge_index_dict)

        # ---- 前向 + 损失 ----
        if use_amp:
            with autocast():
                z1_dict, z2_dict = model(x_dict, edge_index_dict, x_aug, edge_aug)
                total_loss, per_type_loss = compute_contrastive_loss(
                    z1_dict, z2_dict, model, loss_weights=loss_weights,
                )
            # 梯度累积
            loss = total_loss / grad_accum
            scaler.scale(loss).backward()
        else:
            z1_dict, z2_dict = model(x_dict, edge_index_dict, x_aug, edge_aug)
            total_loss, per_type_loss = compute_contrastive_loss(
                z1_dict, z2_dict, model, loss_weights=loss_weights,
            )
            loss = total_loss / grad_accum
            loss.backward()

        # ---- 梯度累积 + 更新 ----
        if (epoch % grad_accum == 0) or (epoch == epochs):
            if use_amp:
                scaler.unscale_(optimizer)
                torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
                scaler.step(optimizer)
                scaler.update()
            else:
                torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
                optimizer.step()
            optimizer.zero_grad()

        # ---- 学习率调度 ----
        if epoch > warmup_epochs:
            scheduler.step()

        # ---- 日志 ----
        if verbose and (epoch == 1 or epoch % 20 == 0 or epoch == epochs):
            loss_str = " | ".join(
                f"{nt}={v.item():.4f}" for nt, v in per_type_loss.items()
            )
            current_lr = optimizer.param_groups[0]["lr"]
            print(f"  Epoch {epoch:4d}/{epochs} | Loss: {total_loss.item():.6f} | "
                  f"PerType: [{loss_str}] | LR: {current_lr:.2e}")

        # ---- 早停检查 ----
        if total_loss.item() < best_loss:
            best_loss = total_loss.item()
            best_state_dict = copy.deepcopy(model.state_dict())
            patience_counter = 0
        else:
            patience_counter += 1
            if patience_counter >= patience:
                if verbose:
                    print(f"  [Early stop] Epoch {epoch} — "
                          f"loss not improved for {patience} epochs")
                break

    # ---- 加载最佳模型 ----
    if best_state_dict is not None:
        model.load_state_dict(best_state_dict)
        print(f"\n[Pretrain] 加载最佳模型 (loss={best_loss:.6f})")

    # ---- 保存 ----
    torch.save(model.state_dict(), save_path)
    print(f"[Pretrain] 权重保存至: {save_path}")

    config = model.get_pretrained_config()
    config["best_loss"] = best_loss
    config["epochs_trained"] = min(epoch, epochs)
    config["final_loss"] = total_loss.item()
    with open(config_save_path, "w") as f:
        json.dump(config, f, indent=2)
    print(f"[Pretrain] 配置保存至: {config_save_path}")

    return model


# ============================================================================
# 6. 权重加载 (微调集成)
# ============================================================================

def load_pretrained_weights(
    target_model: nn.Module,
    pretrained_path: Path = PRETRAINED_WEIGHTS_PATH,
    strict: bool = False,
    verbose: bool = True,
) -> int:
    """将预训练对比学习权重加载到下游 MRHormerModel。

    匹配编码器部分的 state_dict 键名:
    - ContrastiveMRHormer.encoder.* → MRHormerModel.*
    - 跳过投影头 (proj_heads) 和温度参数 (logit_scale) 等非编码器参数

    Args:
        target_model: 目标 MRHormerModel (待微调) 或任意 nn.Module
        pretrained_path: 预训练权重路径
        strict: 是否严格匹配所有键 (默认 False, 跳过不匹配的键)
        verbose: 是否打印加载信息

    Returns:
        loaded_count: 成功加载的参数数量
    """
    if not pretrained_path.exists():
        if verbose:
            print(f"[Load] 预训练权重不存在: {pretrained_path}, 跳过加载")
        return 0

    # 加载预训练 state_dict
    pretrained_sd = torch.load(pretrained_path, map_location="cpu")
    target_sd = target_model.state_dict()

    # 提取 encoder 部分的键: strip "encoder." 前缀
    encoder_prefix = "encoder."
    loaded_count = 0
    skipped_count = 0
    mismatched_count = 0

    for pt_key, pt_val in pretrained_sd.items():
        # 跳过非编码器参数 (投影头、温度参数)
        if not pt_key.startswith(encoder_prefix):
            continue

        # 映射到目标模型键名
        mapped_key = pt_key[len(encoder_prefix):]  # 去掉 "encoder." 前缀

        if mapped_key in target_sd:
            if target_sd[mapped_key].shape == pt_val.shape:
                target_sd[mapped_key] = pt_val.to(target_sd[mapped_key].dtype)
                loaded_count += 1
            else:
                mismatched_count += 1
                if verbose:
                    print(f"  [Warn] 形状不匹配: {mapped_key} "
                          f"预训练={tuple(pt_val.shape)}, "
                          f"目标={tuple(target_sd[mapped_key].shape)}")
        else:
            skipped_count += 1

    # 加载到目标模型
    target_model.load_state_dict(target_sd, strict=strict)

    if verbose:
        print(f"\n[Load] 预训练权重: {pretrained_path}")
        print(f"  [Load] 成功加载: {loaded_count} 参数")
        if skipped_count > 0:
            print(f"  [Load] 跳过 (无匹配键): {skipped_count}")
        if mismatched_count > 0:
            print(f"  [Load] 形状不匹配: {mismatched_count}")
        print(f"  [Load] 编码器权重已注入目标模型")

    return loaded_count


def save_encoder_only(
    contrastive_model: ContrastiveMRHormer,
    save_path: Path = None,
) -> Path:
    """仅保存编码器权重（不含投影头），用于更灵活的微调加载。

    Args:
        contrastive_model: 训练好的对比模型
        save_path: 保存路径 (默认: 同目录 encoder_weights.pt)

    Returns:
        保存路径
    """
    if save_path is None:
        save_path = Path(__file__).parent / "pretrained_encoder_weights.pt"

    encoder_sd = contrastive_model.encoder.state_dict()
    torch.save(encoder_sd, save_path)
    print(f"[Save] 编码器权重保存至: {save_path} (不含投影头)")
    return save_path


# ============================================================================
# 7. 快速评估: 对比损失计算 (用于验证集监控)
# ============================================================================

@torch.no_grad()
def evaluate_contrastive_loss(
    model: ContrastiveMRHormer,
    data: HeteroData,
    augmentor: HeteroGraphAugmentation,
    loss_weights: Optional[Dict[str, float]] = None,
) -> float:
    """在验证集上计算平均对比损失（用于训练监控）。"""
    model.eval()
    data = data.to(DEVICE)
    x_dict = {nt: data[nt].x for nt in data.node_types}
    edge_index_dict = {et: data[et].edge_index for et in data.edge_types}

    x_aug, edge_aug = augmentor(x_dict, edge_index_dict)
    z1_dict, z2_dict = model(x_dict, edge_index_dict, x_aug, edge_aug)
    total_loss, _ = compute_contrastive_loss(
        z1_dict, z2_dict, model, loss_weights=loss_weights,
    )
    return total_loss.item()


# ============================================================================
# 8. 参数解析 & 主入口
# ============================================================================

def parse_args() -> argparse.Namespace:
    """解析命令行参数。"""
    parser = argparse.ArgumentParser(
        description="HeCo-Style Heterogeneous Graph Contrastive Pretraining",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # 训练参数
    parser.add_argument("--epochs", type=int, default=CONTRASTIVE_EPOCHS,
                        help="训练轮数")
    parser.add_argument("--lr", type=float, default=CONTRASTIVE_LR,
                        help="学习率")
    parser.add_argument("--tau", type=float, default=CONTRASTIVE_TAU,
                        help="InfoNCE 温度")
    parser.add_argument("--edge-dropout", type=float, default=CONTRASTIVE_EDGE_DROPOUT,
                        help="View 2 边丢弃概率")
    parser.add_argument("--feat-mask", type=float, default=CONTRASTIVE_FEAT_MASK,
                        help="View 2 特征掩码概率")
    parser.add_argument("--hidden-dim", type=int, default=CONTRASTIVE_HIDDEN_DIM,
                        help="编码器隐藏维度")
    parser.add_argument("--num-heads", type=int, default=CONTRASTIVE_NUM_HEADS,
                        help="注意力头数")
    parser.add_argument("--num-layers", type=int, default=CONTRASTIVE_NUM_LAYERS,
                        help="TRLA 层数")
    parser.add_argument("--proj-dim", type=int, default=CONTRASTIVE_PROJ_DIM,
                        help="投影头输出维度")
    parser.add_argument("--grad-accum", type=int, default=CONTRASTIVE_GRAD_ACCUM,
                        help="梯度累积步数")
    parser.add_argument("--weight-decay", type=float, default=CONTRASTIVE_WEIGHT_DECAY,
                        help="权重衰减")
    parser.add_argument("--seed", type=int, default=42,
                        help="随机种子")

    # 数据 & 保存
    parser.add_argument("--use-full-graph", action="store_true", default=True,
                        help="使用完整图 (含 pathway/cpg)")
    parser.add_argument("--no-full-graph", action="store_false", dest="use_full_graph",
                        help="仅使用 drug/gene/disease 子图")
    parser.add_argument("--save-path", type=str,
                        default=str(PRETRAINED_WEIGHTS_PATH),
                        help="权重保存路径")
    parser.add_argument("--resume", type=str, default=None,
                        help="从检查点恢复训练")

    return parser.parse_args()


def main() -> None:
    """主入口: 数据加载 → 对比预训练 → 权重保存。"""
    args = parse_args()

    # 设置种子
    set_seed(args.seed)

    print(f"\n{'#'*60}")
    print(f"  HeCo-Style 异构图对比预训练")
    print(f"  BCP × Cuproptosis × CIRI Target Screening System")
    print(f"{'#'*60}\n")

    # ---- 加载数据 ----
    data, node_types, edge_types, dim_dict = load_pretrain_data(
        use_full_graph=args.use_full_graph,
    )

    # ---- 执行预训练 ----
    model = pretrain(
        data=data,
        node_types=node_types,
        edge_types=edge_types,
        dim_dict=dim_dict,
        epochs=args.epochs,
        lr=args.lr,
        tau=args.tau,
        edge_dropout_p=args.edge_dropout,
        feat_mask_p=args.feat_mask,
        hidden_dim=args.hidden_dim,
        num_heads=args.num_heads,
        num_layers=args.num_layers,
        proj_dim=args.proj_dim,
        grad_accum=args.grad_accum,
        weight_decay=args.weight_decay,
        save_path=Path(args.save_path),
    )

    # ---- 额外保存编码器权重 (无需投影头) ----
    encoder_path = Path(args.save_path).parent / "pretrained_encoder_weights.pt"
    save_encoder_only(model, encoder_path)

    # ---- 打印总结 ----
    print(f"\n{'='*60}")
    print(f"  预训练完成")
    print(f"  - 编码器: MRHormer (TRLA+MRGA)")
    print(f"  - 节点类型: {node_types}")
    print(f"  - 权重: {args.save_path}")
    print(f"  - 编码器权重: {encoder_path}")
    print(f"  - 使用: load_pretrained_weights(model, '{args.save_path}')")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()