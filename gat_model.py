"""
石竹烯-CIRI GAT 模型定义（标准范式修正版）
使用 PyTorch Geometric 的 GATv2Conv 实现图注意力网络

核心修正：
1. 使用 GATv2Conv + edge_dim 将 STRING combined_score 作为边特征输入，
   让模型自主学习边权重与注意力的非线性映射（Hetero-KGraphDTI范式）。
2. 移除自定义 message() 方法，避免破坏 softmax 归一化。
3. add_self_loops=True 由 PyG 自动处理，无需手动干预 edge_weight 长度。
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GATv2Conv


class GAT(nn.Module):
    """
    图注意力网络 (Graph Attention Network) - GATv2 标准范式

    架构:
      - 输入层: in_channels=35
      - GATv2Conv 第1层: hidden_channels=64, heads=4, concat=True → 256维
      - GATv2Conv 第2层: out_channels=32, heads=4, concat=True → 128维
      - 输出层: Linear(128 → num_classes)

    边特征处理:
      - edge_attr: (E, 1) 一维边特征，即 STRING combined_score 的 z-score/CDF 归一化值
      - edge_dim=1 让 GATv2Conv 的 lin_edge 学习边特征到注意力空间的映射
    """

    def __init__(
        self,
        in_channels: int = 35,
        hidden_channels: int = 64,
        out_channels: int = 32,
        num_heads: int = 4,
        num_classes: int = 2,
        dropout: float = 0.5,
        attention_dropout: float = 0.5,
        use_edge_attr: bool = True,
        use_batch_norm: bool = False,
    ):
        super(GAT, self).__init__()

        self.use_edge_attr = use_edge_attr
        self.dropout = dropout
        self.use_batch_norm = use_batch_norm

        # 第1层 GATv2Conv（支持 edge_dim=1）
        self.gat1 = GATv2Conv(
            in_channels=in_channels,
            out_channels=hidden_channels,
            heads=num_heads,
            concat=True,
            dropout=attention_dropout,
            add_self_loops=True,
            edge_dim=1 if use_edge_attr else None,
        )

        # 第2层 GATv2Conv
        self.gat2 = GATv2Conv(
            in_channels=hidden_channels * num_heads,
            out_channels=out_channels,
            heads=num_heads,
            concat=True,
            dropout=attention_dropout,
            add_self_loops=True,
            edge_dim=1 if use_edge_attr else None,
        )

        # 可选 BatchNorm
        if use_batch_norm:
            self.bn1 = nn.BatchNorm1d(hidden_channels * num_heads)
            self.bn2 = nn.BatchNorm1d(out_channels * num_heads)

        # 输出分类层
        self.classifier = nn.Linear(out_channels * num_heads, num_classes)

        # 只初始化自己添加的层，不覆盖 PyG GATv2Conv 的默认初始化
        nn.init.xavier_uniform_(self.classifier.weight)
        nn.init.zeros_(self.classifier.bias)

    def forward(
        self,
        x: torch.Tensor,
        edge_index: torch.Tensor,
        edge_attr: torch.Tensor = None,
    ) -> torch.Tensor:
        """
        前向传播

        Args:
            x: 节点特征 (N, in_channels)
            edge_index: 边索引 (2, E)
            edge_attr: 边特征 (E, 1)，可选。当 use_edge_attr=True 时必须提供

        Returns:
            logits: (N, num_classes)
        """
        # 第1层 - 使用 edge_attr 参数
        x = self.gat1(
            x, edge_index,
            edge_attr=edge_attr if self.use_edge_attr else None
        )
        if self.use_batch_norm:
            x = self.bn1(x)
        x = F.elu(x)
        x = F.dropout(x, p=self.dropout, training=self.training)

        # 第2层
        x = self.gat2(
            x, edge_index,
            edge_attr=edge_attr if self.use_edge_attr else None
        )
        if self.use_batch_norm:
            x = self.bn2(x)
        x = F.elu(x)
        x = F.dropout(x, p=self.dropout, training=self.training)

        # 输出层
        logits = self.classifier(x)

        return logits

    def get_attention_weights(self, x, edge_index, edge_attr=None):
        """获取注意力权重（用于可视化分析）

        GATv2Conv 在 return_attention_weights=True 时返回:
        - (x, edge_index, attn_weights) 三元组
        attn_weights 形状: (E, heads)
        """
        result = self.gat1(
            x, edge_index,
            edge_attr=edge_attr if self.use_edge_attr else None,
            return_attention_weights=True
        )
        if isinstance(result, tuple):
            if len(result) == 3:
                _, _, attn_weights = result
            elif len(result) == 2:
                _, attn_weights = result
            else:
                attn_weights = result[-1]
            return attn_weights
        else:
            return result
