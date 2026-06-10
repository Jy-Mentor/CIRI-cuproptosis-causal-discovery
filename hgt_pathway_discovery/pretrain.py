# -*- coding: utf-8 -*-
"""GAE 预训练模块：异构图自编码器预训练，获得信息更丰富的初始嵌入。

改进:
  1. 在异构图上游使用 GAE 做无监督预训练
  2. 重构边掩码任务：随机遮蔽部分 PPI 边，训练 HGT 编码器恢复它们
  3. 预训练权重注入下游 HGTModel，提升收敛速度和最终性能

参考:
  - Kipf & Welling, "Variational Graph Auto-Encoders", NeurIPS 2016
  - Hu et al., "Strategies for Pre-training Graph Neural Networks", ICLR 2020
  - You et al., "Graph Contrastive Learning", NeurIPS 2020
"""

import copy
import logging
from typing import Dict, Tuple, List, Optional

import numpy as np
import torch
import torch.nn.functional as F
from torch import Tensor, nn
from torch_geometric.data import HeteroData
from torch_geometric.nn import HGTConv, Linear

from .model import HGTModel, NegEdgeSampler
from .train import remove_edges_from_data


logger = logging.getLogger("hgt_pipeline")


# ============================================================================
# 异构图自编码器 (GAE on Heterogeneous Graph)
# ============================================================================

class HeteroGAE(nn.Module):
    """异构图自编码器：基于 HGT 编码器 + 内积解码器的无监督预训练。

    在 PPI 图上训练边重构任务，学习基因节点的结构嵌入。
    预训练完成后，编码器权重可迁移到下游 HGTModel。

    训练方式:
      随机遮蔽 PPI 边 → HGT 编码全图 → 内积解码被遮蔽边 → BCE Loss
    """

    def __init__(self, metadata: Tuple[List[str], List[Tuple[str, str, str]]],
                 dim_dict: Dict[str, int], hidden_dim: int,
                 num_heads: int, num_layers: int, dropout: float):
        super().__init__()
        self.metadata = metadata
        self.hidden_dim = hidden_dim

        # 类型感知投影层 — 与 HGTModel 结构一致，便于权重迁移
        self.proj = nn.ModuleDict()
        for nt, d_in in dim_dict.items():
            self.proj[nt] = Linear(d_in, hidden_dim)

        # HGT 编码器 (复用 HGTModel 的编码器结构)
        from .model import HGTEncoder
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
        """编码异构图，返回各节点类型的嵌入。"""
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
        """内积解码器：score = sigmoid(z_src @ z_dst.T)。"""
        z_src = z_dict[src_type][edge_index[0]]
        z_dst = z_dict[dst_type][edge_index[1]]
        scores = (z_src * z_dst).sum(dim=-1)
        return scores

    def forward(self, data: HeteroData,
                pos_edge_index: Tensor,
                neg_edge_index: Tensor) -> Tuple[Tensor, Tensor]:
        """前向传播：编码 → 解码正边和负边。

        Returns:
            pos_scores, neg_scores
        """
        z_dict = self.encode(data)
        pos_scores = self.decode(z_dict, pos_edge_index)
        neg_scores = self.decode(z_dict, neg_edge_index)
        return pos_scores, neg_scores

    def transfer_weights_to(self, target_model: HGTModel) -> None:
        """将预训练的投影层和编码器权重迁移到下游 HGTModel。

        仅迁移共有的层：proj (类型感知投影) 和 encoder。
        解码器权重保持不变（下游任务不同）。
        """
        # 迁移投影层
        for nt in self.proj:
            if nt in target_model.proj:
                target_model.proj[nt].load_state_dict(self.proj[nt].state_dict())
                logger.debug(f"  Pretrain → HGTModel: proj[{nt}] transferred")

        # 迁移编码器 (HGTEncoder)
        target_enc = target_model.encoder
        pretrain_enc = self.encoder

        for i, (t_conv, p_conv) in enumerate(zip(target_enc.convs, pretrain_enc.convs)):
            t_conv.load_state_dict(p_conv.state_dict())
            logger.debug(f"  Pretrain → HGTModel: conv[{i}] transferred")

        for i, (t_norm, p_norm) in enumerate(zip(target_enc.norms, pretrain_enc.norms)):
            t_norm.load_state_dict(p_norm.state_dict())
            logger.debug(f"  Pretrain → HGTModel: norm[{i}] transferred")

        if (target_enc.skip_alphas is not None and
                pretrain_enc.skip_alphas is not None):
            for i, (t_alpha, p_alpha) in enumerate(
                zip(target_enc.skip_alphas, pretrain_enc.skip_alphas)
            ):
                t_alpha.data.copy_(p_alpha.data)
                logger.debug(f"  Pretrain → HGTModel: skip_alpha[{i}] transferred")

        logger.info("Pretrained weights successfully transferred to HGTModel")


# ============================================================================
# 边掩码预训练
# ============================================================================

def mask_edges(edge_index: Tensor, mask_ratio: float = 0.3,
               seed: int = 42) -> Tuple[Tensor, Tensor, Tensor]:
    """随机遮蔽部分边用于自监督预训练。

    Args:
        edge_index: (2, n_edges) 边索引
        mask_ratio: 遮蔽比例
        seed: 随机种子

    Returns:
        train_ei: 用于训练的可见边
        mask_ei: 被遮蔽的边 (正样本，需要重构)
        neg_ei: 随机采样的负边
    """
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
    """执行异构图自编码器预训练。

    在 PPI 边上训练边重构任务，获得基因节点的结构嵌入。

    Args:
        data: 完整的异构图
        config: 全局配置
        device: 计算设备

    Returns:
        预训练好的 HeteroGAE 模型，或 None（若预训练禁用）
    """
    pt_cfg = config.pretraining
    if not pt_cfg.enabled:
        logger.info("Pretraining disabled, skipping...")
        return None

    m_cfg = config.model

    logger.info("=" * 60)
    logger.info("[Pretrain] Heterogeneous GAE Pre-training")
    logger.info(f"  Epochs: {pt_cfg.epochs}, LR: {pt_cfg.lr}, "
                f"Mask Ratio: {pt_cfg.mask_ratio}")
    logger.info("=" * 60)

    # 选择 PPI 边类型作为预训练任务
    # 优先级: strong_ppi > weak_ppi > interacts
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
        logger.warning("No PPI edges found for pretraining, skipping...")
        return None

    ppi_ei = data[ppi_et].edge_index
    logger.info(f"  Using edge type: {ppi_et} ({ppi_ei.size(1)} edges)")

    # 遮蔽边
    train_ei, mask_ei, neg_ei = mask_edges(
        ppi_ei, pt_cfg.mask_ratio, config.seed,
    )
    logger.info(f"  Train edges: {train_ei.size(1)}, "
                f"Masked edges: {mask_ei.size(1)}, "
                f"Neg edges: {neg_ei.size(1)}")

    # 构建训练图（移除被遮蔽的边）
    mask_np = np.zeros(ppi_ei.size(1), dtype=bool)
    rng = np.random.RandomState(config.seed)
    n_mask = mask_ei.size(1)
    n_total = ppi_ei.size(1)
    perm = rng.permutation(n_total)
    mask_np[perm[:n_mask]] = True
    data_train = remove_edges_from_data(
        data, ppi_et, torch.from_numpy(~mask_np),
    ).to(device)

    # 创建 GAE 模型
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
            logger.info(
                f"  Pretrain Epoch {epoch+1:3d}/{pt_cfg.epochs} | "
                f"Loss: {loss:.4f} | PosAcc: {pos_acc:.3f} | "
                f"NegAcc: {neg_acc:.3f}"
            )

    if best_state is not None:
        gae.load_state_dict(best_state)
        logger.info(f"  Pretrain complete: Best Loss = {best_loss:.4f}")

    # 保存预训练模型
    save_path = config.paths.pretrain_save if hasattr(config.paths, "pretrain_save") else None
    if save_path:
        torch.save(gae.state_dict(), save_path)
        logger.info(f"  Pretrained model saved to {save_path}")

    return gae