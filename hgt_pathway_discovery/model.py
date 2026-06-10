# -*- coding: utf-8 -*-
"""HGT 编码器-解码器模型：增强版架构。

改进:
  1. 3 层 HGT + GCNII 式初始残差 — 缓解深层过平滑
  2. 解耦解码器 (gene_bias + pathway_bias) — 关系特异性预测
  3. DistMult 因子分解 — 增强可解释性
  4. CpG 可学习参数 — 初始特征注入 nn.Parameter，训练中自适应优化
"""

import copy
from typing import Dict, List, Tuple, Optional

import numpy as np
import torch
import torch.nn.functional as F
from torch import Tensor, nn
from torch_geometric.nn import HGTConv, Linear


class HGTEncoder(nn.Module):
    """多层 HGTConv 编码器：初始残差 + LayerNorm + ReLU + Dropout + DropEdge。

    参考:
      - HGT: Hu et al., WWW 2020
      - GCNII: Chen et al., ICML 2020 (初始残差思想)
      - DropEdge: Rong et al., ICLR 2020 (随机丢弃边防止过平滑)
    """

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
            # GCNII: 每层以初始嵌入为残差参考
            self.skip_alphas = nn.ParameterList([
                nn.Parameter(torch.tensor(0.5)) for _ in range(num_layers)
            ])
        else:
            self.skip_alphas = None

        self.dropout = dropout

    @staticmethod
    def _drop_edges(edge_index_dict: Dict, drop_p: float) -> Dict:
        """对每种边类型独立执行 DropEdge。

        随机丢弃比例为 drop_p 的边，防止深层 GNN 过平滑。
        参考: Rong et al., "DropEdge: Towards Deep Graph Convolutional
        Networks on Node Classification", ICLR 2020.

        Args:
            edge_index_dict: 原始边索引字典 {edge_type: (2, N)}
            drop_p: 丢弃概率
        Returns:
            随机丢弃后的边索引字典 (原地返回 dict，值被替换)
        """
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
                keep_mask[0] = True  # 至少保留一条边
            dropped[et] = ei[:, keep_mask]
        return dropped

    def forward(self, x_dict: Dict[str, Tensor],
                edge_index_dict: Dict,
                x0_dict: Optional[Dict[str, Tensor]] = None) -> Dict[str, Tensor]:
        """前向传播，支持初始残差和 DropEdge。

        Args:
            x_dict: 当前层节点特征
            edge_index_dict: 边索引字典
            x0_dict: 初始投影特征 (用于 GCNII 初始残差)
        """
        for layer_idx, (conv, norm) in enumerate(zip(self.convs, self.norms)):
            # DropEdge: 训练时随机丢弃边，每层独立采样
            if self.training and self.drop_edge_p > 0:
                ei_dict = self._drop_edges(edge_index_dict, self.drop_edge_p)
            else:
                ei_dict = edge_index_dict

            x_dict_new = conv(x_dict, ei_dict)

            # 补全未参与消息传递的节点类型
            for k in x_dict:
                if k not in x_dict_new:
                    x_dict_new[k] = x_dict[k]

            # 残差连接 + LayerNorm + ReLU
            x_dict = {}
            for k in x_dict_new:
                residual = x_dict.get(k, x_dict_new[k])
                x_dict[k] = F.relu(norm(x_dict_new[k] + residual))

            # GCNII 初始残差
            if self.initial_residual and x0_dict is not None and self.skip_alphas is not None:
                alpha = torch.sigmoid(self.skip_alphas[layer_idx])
                for k in x_dict:
                    if k in x0_dict:
                        x_dict[k] = (1 - alpha) * x_dict[k] + alpha * x0_dict[k]

            # Dropout
            x_dict = {
                k: F.dropout(v, p=self.dropout, training=self.training)
                for k, v in x_dict.items()
            }

        return x_dict


class DecoupledDecoder(nn.Module):
    """解耦解码器：为基因和通路分别添加可学习偏置，增强关系特异性。

    结构: gene_bias(z_gene) || path_bias(z_path) → MLP → Output

    可选: DistMult 因子分解 (z_gene * W * z_path).sum()
    """

    def __init__(self, hidden_dim: int, dropout: float,
                 factorization: str = "distmult"):
        super().__init__()
        self.factorization = factorization

        # 类型特异性偏置投影
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
            # DistMult: 对角矩阵 W 将基因和通路嵌入映射到标量分数
            self.W = nn.Parameter(torch.randn(hidden_dim) * 0.1)
            self.output = nn.Linear(1, 1)
        else:
            # 标准 MLP 解码器
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
        """解码基因-通路对分数。

        Args:
            z_gene: (N, hidden_dim) 基因节点嵌入
            z_path: (N, hidden_dim) 通路节点嵌入
        Returns:
            scores: (N,) 预测分数 (logits)
        """
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
    """HGT 链路预测模型：类型感知投影 + BatchNorm + HGT编码器 + 解耦解码器。

    支持 CpG 可学习参数：将特征传播初始值注入 nn.Parameter，
    训练中通过梯度下降自适应优化，同时保留质量掩码引导。

    改进:
      - BatchNorm1d 作为投影后第一层：统一多模态特征分布，增强梯度稳定性
      - 参考: Ioffe & Szegedy, "Batch Normalization", ICML 2015
    """

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

        # 类型感知投影层
        self.proj = nn.ModuleDict()
        for nt, d_in in dim_dict.items():
            self.proj[nt] = Linear(d_in, hidden_dim)

        # BatchNorm 作为投影后第一层 (统一多模态特征分布)
        if use_input_bn:
            self.input_bn = nn.ModuleDict()
            for nt in dim_dict:
                self.input_bn[nt] = nn.BatchNorm1d(hidden_dim)
        else:
            self.input_bn = None

        # HGT 编码器
        self.encoder = HGTEncoder(
            metadata, hidden_dim, num_heads, num_layers,
            dropout, initial_residual=initial_residual,
            drop_edge_p=drop_edge_p,
        )

        # 解耦解码器
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
        """编码：类型投影 → BatchNorm → HGT 编码。

        若 CpG 可学习参数已注册，在前向传播中优先使用。
        质量掩码作为条件信息：cpg_x = cpg_embed * quality_mask + cpg_bias

        Returns:
            z_dict: 各节点类型的最终嵌入
        """
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

        # BatchNorm 统一多模态特征分布 (跳过单节点类型)
        if self.input_bn is not None:
            for k in x_proj:
                if k in self.input_bn and x_proj[k].size(0) > 1:
                    x_proj[k] = self.input_bn[k](x_proj[k])

        # 初始残差: 保存投影+BN后的初始嵌入
        x0_dict = copy.copy(x_proj) if self.initial_residual else None

        z_dict = self.encoder(x_proj, edge_index_dict, x0_dict=x0_dict)
        return z_dict

    def decode(self, z_dict: Dict[str, Tensor], edge_index: Tensor) -> Tensor:
        """解码：提取基因-通路对嵌入并计算分数。

        Args:
            z_dict: 编码器输出的节点嵌入
            edge_index: (2, N) 边索引 [gene_idx, pathway_idx]
        """
        z_gene = z_dict["gene"][edge_index[0]]
        z_path = z_dict["pathway"][edge_index[1]]
        return self.decoder(z_gene, z_path)

    def decode_chunked(self, z_dict: Dict[str, Tensor],
                       edge_index: Tensor, chunk_size: int = 16384) -> Tensor:
        """分块解码：防止大批次 OOM。

        Args:
            z_dict: 编码器输出的节点嵌入
            edge_index: (2, N) 边索引
            chunk_size: 每块大小
        """
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
        """将 CpG 初始特征注册为可学习参数，并集成质量掩码作为条件信息。

        调用时机: 图构建后、训练前。
        原理: 特征传播提供合理的初始值，训练中通过梯度下降进一步优化。
        质量掩码标记哪些 CpG 获得了真实传播特征（HGNN-IMA, Li et al., 2025），
        未被传播到的 CpG 其掩码值为 0，模型通过 bias 获得共享先验。

        Args:
            cpg_init_feat: 特征传播初始化的 CpG 特征 (N_cpg, D_features)
            quality_mask: 原始传播二进制掩码 (N_cpg,) — True 表示该 CpG 有邻居基因
        """
        self.cpg_embed = nn.Parameter(cpg_init_feat.clone())
        self.register_parameter("cpg_embed", self.cpg_embed)

        if quality_mask is not None:
            self.cpg_quality_mask = quality_mask.float().unsqueeze(1)
            self.register_buffer("cpg_quality_mask", self.cpg_quality_mask)

        self.cpg_bias = nn.Parameter(torch.zeros(cpg_init_feat.size(1)))
        self.register_parameter("cpg_bias", self.cpg_bias)


# ============================================================================
# Focal BCE Loss (Lin et al., ICCV 2017)
# ============================================================================

def focal_bce_loss(logits: Tensor, labels: Tensor,
                   alpha: float = 0.5, gamma: float = 2.0) -> Tensor:
    """Focal BCE Loss：通过调制因子聚焦难样本。"""
    logits = torch.clamp(logits, -10, 10)
    bce = F.binary_cross_entropy_with_logits(logits, labels, reduction="none")
    pt = torch.where(labels == 1, torch.sigmoid(logits), 1 - torch.sigmoid(logits))
    focal_weight = (1 - pt) ** gamma
    alpha_weight = torch.where(labels == 1, alpha, 1 - alpha)
    return (alpha_weight * focal_weight * bce).mean()


# ============================================================================
# 负边采样器 (增强版：自适应负采样比 + 定时刷新验证集)
# ============================================================================

class NegEdgeSampler:
    """负边采样器：支持均匀随机和度加权采样，动态排除已知正边。

    增强:
      - 自适应负采样比: 训练后期逐步提高负样本比例
      - 固定验证集 + 定时刷新: 每 refresh_interval epoch 重采样
    """

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
        self._exclude_frozen = frozenset(self.exclude)  # O(1) 成员检查 (不可变)
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
        """采样负边。

        Args:
            n_pos: 正样本数量
            neg_ratio: 负正比 (1:n)
        """
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
        """采样固定负边集（用于验证集，确保不同 epoch 可比）。

        Args:
            force_refresh: 强制刷新缓存
        """
        if force_refresh or self._fixed is None or self._fixed_n_pos != n_pos or self._fixed_neg_ratio != neg_ratio:
            self._fixed = self.sample(n_pos, neg_ratio)
            self._fixed_n_pos = n_pos
            self._fixed_neg_ratio = neg_ratio
        return self._fixed