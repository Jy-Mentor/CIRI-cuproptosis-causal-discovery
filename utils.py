"""
石竹烯-CIRI GAT 训练工具模块
包含数据加载、mask 生成、评估指标计算等通用函数
"""

import os
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any

import numpy as np
import pandas as pd
import torch
import yaml
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score, f1_score,
    roc_auc_score, average_precision_score, confusion_matrix
)

logger = logging.getLogger(__name__)


def load_config(config_path: str = "config.yaml") -> Dict[str, Any]:
    """加载 YAML 配置文件"""
    with open(config_path, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)
    return config


def set_seed(seed: int = 42):
    """设置随机种子确保可复现性"""
    torch.manual_seed(seed)
    np.random.seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


def load_graph_data(processed_dir: str) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor, List[str]]:
    """
    加载图数据（GATv2 标准范式：返回 edge_attr 而非 edge_weight）
    
    Returns:
        x: 节点特征矩阵 (N, F)
        edge_index: 边索引 (2, E)
        edge_attr: 边特征 (E, 1)，即 STRING combined_score 的 z-score/CDF 归一化值
        gene_symbols: 基因符号列表 (N,)
    """
    processed_path = Path(processed_dir)
    
    # 加载节点特征
    node_df = pd.read_csv(processed_path / "node_features.csv")
    gene_symbols = node_df["GeneSymbol"].tolist()
    feature_cols = [c for c in node_df.columns if c != "GeneSymbol"]

    # 拓扑特征合法性检查：聚类系数必须在 [0, 1] 范围内
    if "ClusteringCoefficient" in node_df.columns:
        cc_col = node_df["ClusteringCoefficient"]
        if not cc_col.between(0, 1).all():
            logger.warning("聚类系数越界 [0, 1]，自动进行 Min-Max 归一化")
            cc_min = cc_col.min()
            cc_max = cc_col.max()
            if cc_max > cc_min:
                node_df["ClusteringCoefficient"] = (cc_col - cc_min) / (cc_max - cc_min)
            else:
                node_df["ClusteringCoefficient"] = 0.0
            # 保存修正后的文件（备份原文件）
            backup_path = processed_path / "node_features_backup.csv"
            if not backup_path.exists():
                pd.read_csv(processed_path / "node_features.csv").to_csv(backup_path, index=False)
            node_df.to_csv(processed_path / "node_features.csv", index=False)

    x = torch.tensor(node_df[feature_cols].values, dtype=torch.float32)

    # 加载边
    edge_df = pd.read_csv(processed_path / "edge_index.csv")
    
    # 创建基因到索引的映射
    gene_to_idx = {g: i for i, g in enumerate(gene_symbols)}
    
    # 过滤基因池外的边
    valid_edges = []
    valid_weights = []
    for _, row in edge_df.iterrows():
        src, tgt = row["Source"], row["Target"]
        if src in gene_to_idx and tgt in gene_to_idx:
            valid_edges.append((gene_to_idx[src], gene_to_idx[tgt]))
            valid_weights.append(float(row["Weight"]))
    
    if len(valid_edges) < len(edge_df):
        logger.warning(f"过滤了 {len(edge_df) - len(valid_edges)} 条基因池外边")
    
    # 构建无向图的 edge_index（双向）
    edge_list = []
    weight_list = []
    for (src, tgt), w in zip(valid_edges, valid_weights):
        edge_list.append([src, tgt])
        edge_list.append([tgt, src])
        weight_list.append(w)
        weight_list.append(w)
    
    edge_index = torch.tensor(edge_list, dtype=torch.long).t().contiguous()
    edge_weight = torch.tensor(weight_list, dtype=torch.float32)
    
    # --- GATv2 标准范式：将 edge_weight 转换为 edge_attr (E, 1) ---
    # 对 combined_score 做 z-score 归一化，使其均值为0，方差为1
    # 这样 lin_edge 可以学习边特征到注意力空间的非线性映射
    ew_mean = edge_weight.mean()
    ew_std = edge_weight.std()
    if ew_std > 0:
        edge_attr = ((edge_weight - ew_mean) / ew_std).unsqueeze(-1)
    else:
        edge_attr = edge_weight.unsqueeze(-1)
    
    logger.info(f"节点数: {x.shape[0]}, 特征维度: {x.shape[1]}")
    logger.info(f"边数: {edge_index.shape[1]}")
    logger.info(f"边特征 edge_attr 形状: {edge_attr.shape}")
    
    return x, edge_index, edge_attr, gene_symbols


def load_labels(processed_dir: str, gene_symbols: List[str]) -> torch.Tensor:
    """
    加载标签，确保与节点特征顺序一致
    
    Returns:
        y: 标签张量 (N,)，-1 表示未知
    """
    processed_path = Path(processed_dir)
    label_df = pd.read_csv(processed_path / "labels.csv")
    
    # 创建基因到标签的映射
    gene_to_label = dict(zip(label_df["GeneSymbol"], label_df["Label"]))
    
    # 按 gene_symbols 顺序提取标签
    labels = [gene_to_label.get(g, -1) for g in gene_symbols]
    y = torch.tensor(labels, dtype=torch.long)
    
    n_pos = (y == 1).sum().item()
    n_neg = (y == 0).sum().item()
    n_unk = (y == -1).sum().item()
    logger.info(f"标签分布 — 阳性: {n_pos}, 阴性: {n_neg}, 未知: {n_unk}")
    
    return y


def create_masks(y: torch.Tensor, train_ratio: float = 0.6,
                 val_ratio: float = 0.2, test_ratio: float = 0.2,
                 random_state: int = 42) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """
    创建训练/验证/测试 mask
    仅对有标签节点（0 和 1）进行划分，未知节点（-1）不参与
    
    Returns:
        train_mask, val_mask, test_mask: (N,) 的 bool 张量
    """
    N = len(y)
    labeled_indices = (y != -1).nonzero(as_tuple=True)[0].numpy()
    labeled_labels = y[labeled_indices].numpy()
    
    # Stratified split
    sss = StratifiedShuffleSplit(n_splits=1, test_size=test_ratio, random_state=random_state)
    train_val_idx, test_idx = next(sss.split(labeled_indices, labeled_labels))
    
    val_size = val_ratio / (train_ratio + val_ratio)
    sss_val = StratifiedShuffleSplit(n_splits=1, test_size=val_size, random_state=random_state)
    train_idx, val_idx = next(sss_val.split(train_val_idx, labeled_labels[train_val_idx]))
    
    train_indices = labeled_indices[train_val_idx[train_idx]]
    val_indices = labeled_indices[train_val_idx[val_idx]]
    test_indices = labeled_indices[test_idx]
    
    train_mask = torch.zeros(N, dtype=torch.bool)
    val_mask = torch.zeros(N, dtype=torch.bool)
    test_mask = torch.zeros(N, dtype=torch.bool)
    
    train_mask[train_indices] = True
    val_mask[val_indices] = True
    test_mask[test_indices] = True
    
    logger.info(f"训练集: {train_mask.sum().item()}, 验证集: {val_mask.sum().item()}, 测试集: {test_mask.sum().item()}")
    
    return train_mask, val_mask, test_mask


def compute_metrics(y_true: np.ndarray, y_pred: np.ndarray, y_prob: np.ndarray) -> Dict[str, float]:
    """计算分类指标"""
    metrics = {}
    
    metrics["accuracy"] = accuracy_score(y_true, y_pred)
    metrics["precision"] = precision_score(y_true, y_pred, zero_division=0)
    metrics["recall"] = recall_score(y_true, y_pred, zero_division=0)
    metrics["f1"] = f1_score(y_true, y_pred, zero_division=0)
    
    # ROC-AUC 和 PR-AUC
    if len(np.unique(y_true)) > 1:
        metrics["roc_auc"] = roc_auc_score(y_true, y_prob)
        metrics["pr_auc"] = average_precision_score(y_true, y_prob)
    else:
        metrics["roc_auc"] = 0.0
        metrics["pr_auc"] = 0.0
    
    return metrics


def plot_confusion_matrix(y_true: np.ndarray, y_pred: np.ndarray, save_path: str):
    """绘制混淆矩阵"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    cm = confusion_matrix(y_true, y_pred)
    plt.figure(figsize=(6, 5))
    sns.heatmap(cm, annot=True, fmt="d", cmap="Blues",
                xticklabels=["Non-target", "Target"],
                yticklabels=["Non-target", "Target"])
    plt.xlabel("Predicted")
    plt.ylabel("True")
    plt.title("Confusion Matrix")
    plt.tight_layout()
    plt.savefig(save_path, dpi=150)
    plt.close()
    logger.info(f"混淆矩阵已保存: {save_path}")


def plot_training_curves(history: Dict[str, List[float]], save_path: str):
    """绘制训练曲线"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    
    # Loss
    axes[0].plot(history["train_loss"], label="Train Loss")
    axes[0].plot(history["val_loss"], label="Val Loss")
    axes[0].set_xlabel("Epoch")
    axes[0].set_ylabel("Loss")
    axes[0].set_title("Loss Curve")
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    # F1
    axes[1].plot(history["val_f1"], label="Val F1", color="green")
    axes[1].set_xlabel("Epoch")
    axes[1].set_ylabel("F1 Score")
    axes[1].set_title("Validation F1")
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    # LR
    if "lr" in history and history["lr"]:
        axes[2].plot(history["lr"], label="Learning Rate", color="red")
        axes[2].set_xlabel("Epoch")
        axes[2].set_ylabel("LR")
        axes[2].set_title("Learning Rate")
        axes[2].legend()
        axes[2].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=150)
    plt.close()
    logger.info(f"训练曲线已保存: {save_path}")
