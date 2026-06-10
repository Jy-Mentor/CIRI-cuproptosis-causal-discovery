#!/usr/bin/env python3
"""
石竹烯-CIRI GATv2 靶点推断与排序脚本（标准范式修正版）
加载训练好的模型，对全部节点进行预测，输出 Top-K 候选靶点

修正点：
1. 使用 edge_attr 替代 edge_weight，适配 GATv2Conv 的 edge_dim 参数
"""

import os
import sys
import logging
import argparse
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F

from gat_model import GAT
from utils import load_config, load_graph_data, load_labels

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
logger = logging.getLogger(__name__)


def load_trained_model(model_path: str, config: Dict, device: torch.device) -> GAT:
    """加载训练好的模型"""
    model_cfg = config["model"]
    model = GAT(
        in_channels=model_cfg["in_channels"],
        hidden_channels=model_cfg["hidden_channels"],
        out_channels=model_cfg["out_channels"],
        num_heads=model_cfg["num_heads"],
        num_classes=model_cfg["num_classes"],
        dropout=0.0,  # 推断时关闭 dropout
        attention_dropout=0.0,
        use_edge_attr=model_cfg["use_edge_attr"],
        use_batch_norm=model_cfg["use_batch_norm"],
    ).to(device)

    checkpoint = torch.load(model_path, map_location=device)
    model.load_state_dict(checkpoint["model_state_dict"])
    model.eval()
    logger.info(f"加载模型: {model_path} (epoch {checkpoint.get('epoch', 'unknown')})")
    return model


@torch.no_grad()
def predict(model: GAT, x, edge_index, edge_attr) -> Tuple[np.ndarray, np.ndarray]:
    """
    对全部节点进行预测

    Returns:
        probs: (N, 2) 每个类别的概率
        logits: (N, 2) 原始 logits
    """
    logits = model(x, edge_index, edge_attr)
    probs = F.softmax(logits, dim=1)
    return probs.cpu().numpy(), logits.cpu().numpy()


def rank_targets(
    gene_symbols: List[str],
    probs: np.ndarray,
    labels: np.ndarray,
    node_features_df: pd.DataFrame,
    top_k_list: List[int] = [50, 100, 200],
    results_dir: Path = Path("./results"),
) -> pd.DataFrame:
    """
    对未知节点按 P(靶点) 降序排序，输出 Top-K 列表

    Args:
        gene_symbols: 基因符号列表
        probs: 预测概率 (N, 2)
        labels: 原始标签 (N,)，-1 为未知
        node_features_df: 节点特征 DataFrame（用于提取拓扑和通路信息）
        top_k_list: 输出 Top-K 列表
        results_dir: 结果保存目录

    Returns:
        full_results: 全部未知节点的排序结果 DataFrame
    """
    results_dir.mkdir(parents=True, exist_ok=True)

    # 只取未知节点
    unknown_mask = labels == -1
    unknown_indices = np.where(unknown_mask)[0]

    unknown_genes = [gene_symbols[i] for i in unknown_indices]
    unknown_probs = probs[unknown_indices, 1]  # P(target)

    # 构建结果 DataFrame
    results = pd.DataFrame({
        "GeneSymbol": unknown_genes,
        "P_target": unknown_probs,
        "Rank": np.arange(1, len(unknown_genes) + 1),
    })

    # 添加网络拓扑特征
    topo_cols = ["Degree", "PageRank", "ClusteringCoefficient"]
    for col in topo_cols:
        if col in node_features_df.columns:
            results[col] = node_features_df.loc[unknown_indices, col].values

    # 添加铜死亡标记
    cupro_genes = {"FDX1", "LIAS", "LIPT1", "DLAT", "PDHB", "PDHX",
                   "SLC31A1", "ATP7A", "ATP7B", "ATOX1", "NFE2L2",
                   "HIF1A", "MTOR", "NFKB1", "GPX4"}
    results["is_cuproptosis"] = results["GeneSymbol"].isin(cupro_genes).astype(int)

    # 添加炎症通路标记
    inflam_cols = ["toll_like_flag", "tnf_flag", "nfkb_flag"]
    results["is_inflammatory"] = 0
    for col in inflam_cols:
        if col in node_features_df.columns:
            results["is_inflammatory"] |= node_features_df.loc[unknown_indices, col].values
    results["is_inflammatory"] = results["is_inflammatory"].clip(0, 1).astype(int)

    # 按 P(target) 降序排序
    results = results.sort_values("P_target", ascending=False).reset_index(drop=True)
    results["Rank"] = np.arange(1, len(results) + 1)

    # 保存 Top-K
    for k in top_k_list:
        top_k = results.head(k).copy()
        save_path = results_dir / f"top_targets_{k}.csv"
        top_k.to_csv(save_path, index=False)
        logger.info(f"Top-{k} 候选靶点已保存: {save_path}")

    # 保存全部未知节点结果
    full_path = results_dir / "all_unknown_predictions.csv"
    results.to_csv(full_path, index=False)
    logger.info(f"全部未知节点预测已保存: {full_path}")

    return results


def print_top_summary(results: pd.DataFrame, top_n: int = 20):
    """打印 Top-N 候选靶点摘要"""
    logger.info("=" * 60)
    logger.info(f"Top-{top_n} 候选靶点摘要")
    logger.info("=" * 60)

    top = results.head(top_n)
    for _, row in top.iterrows():
        markers = []
        if row["is_cuproptosis"] == 1:
            markers.append("[铜死亡]")
        if row["is_inflammatory"] == 1:
            markers.append("[炎症]")
        marker_str = " ".join(markers) if markers else ""
        logger.info(
            f"Rank {int(row['Rank']):03d}: {row['GeneSymbol']:12s} "
            f"P={row['P_target']:.4f} "
            f"Degree={row.get('Degree', 0):.1f} "
            f"PageRank={row.get('PageRank', 0):.6f} "
            f"{marker_str}"
        )


def main():
    parser = argparse.ArgumentParser(description="Predict and rank potential targets")
    parser.add_argument("--config", type=str, default="config.yaml", help="Path to config file")
    parser.add_argument("--model", type=str, default=None, help="Path to model checkpoint")
    args = parser.parse_args()

    config = load_config(args.config)
    data_cfg = config["data"]
    out_cfg = config["output"]

    # 设备
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    logger.info(f"使用设备: {device}")

    # 加载数据（GATv2 标准范式：返回 edge_attr）
    processed_dir = data_cfg["processed_dir"]
    x, edge_index, edge_attr, gene_symbols = load_graph_data(processed_dir)
    y = load_labels(processed_dir, gene_symbols)

    x = x.to(device)
    edge_index = edge_index.to(device)
    edge_attr = edge_attr.to(device)

    # 加载节点特征（用于附加信息）
    node_features_df = pd.read_csv(Path(processed_dir) / "node_features.csv")

    # 加载模型
    model_path = args.model or (Path(out_cfg["model_dir"]) / out_cfg["best_model"])
    model = load_trained_model(str(model_path), config, device)

    # 推断
    logger.info("开始推断...")
    probs, logits = predict(model, x, edge_index, edge_attr)
    logger.info(f"推断完成，共 {len(probs)} 个节点")

    # 排序
    results_dir = Path(out_cfg["results_dir"])
    results = rank_targets(
        gene_symbols=gene_symbols,
        probs=probs,
        labels=y.numpy(),
        node_features_df=node_features_df,
        top_k_list=out_cfg["top_k"],
        results_dir=results_dir,
    )

    # 打印摘要
    print_top_summary(results, top_n=20)

    logger.info("=" * 60)
    logger.info("靶点排序完成")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
