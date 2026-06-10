#!/usr/bin/env python3
"""
石竹烯-CIRI GATv2 主训练脚本（标准范式修正版）
半监督学习在 3,918 个基因节点中排序潜在治疗靶点

修正点：
1. 使用 edge_attr 替代 edge_weight，适配 GATv2Conv 的 edge_dim 参数
2. 铜死亡先验损失仅作用于有标签节点（训练集内），避免污染未知集
"""

import os
import sys
import logging
import argparse
from pathlib import Path
from typing import Dict, List

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.optim import AdamW
from torch.optim.lr_scheduler import CosineAnnealingLR, ReduceLROnPlateau

from gat_model import GAT
from utils import (
    load_config, set_seed, load_graph_data, load_labels,
    create_masks, compute_metrics, plot_confusion_matrix, plot_training_curves
)

# ---------------------------------------------------------------------------
# 日志配置
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
    ],
)
logger = logging.getLogger(__name__)


def setup_logging(log_dir: Path):
    """设置日志文件"""
    log_dir.mkdir(parents=True, exist_ok=True)
    file_handler = logging.FileHandler(log_dir / "training_log.txt", encoding="utf-8")
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
    logger.addHandler(file_handler)


def train_epoch(model, optimizer, x, edge_index, edge_attr, y, train_mask, pos_weight,
                cupro_mask=None, cupro_weight=0.3) -> float:
    """训练一个 epoch

    Args:
        cupro_mask: 铜死亡基因的 mask (N,), bool 张量
        cupro_weight: 铜死亡先验损失权重
    """
    model.train()
    optimizer.zero_grad()

    logits = model(x, edge_index, edge_attr)

    # 仅对有标签节点计算损失
    train_logits = logits[train_mask]
    train_labels = y[train_mask]

    # 类别权重
    weight = torch.tensor([1.0, pos_weight], device=logits.device)
    loss = F.cross_entropy(train_logits, train_labels, weight=weight)

    # 铜死亡基因先验损失：仅作用于有标签节点（训练集内），避免污染未知集
    if cupro_mask is not None and cupro_mask.any():
        # 只取训练集中的铜死亡基因
        train_cupro_mask = cupro_mask & train_mask
        if train_cupro_mask.any():
            cupro_logits = logits[train_cupro_mask]
            # 软标签：目标概率 0.6 而非 1.0，避免过强先验
            cupro_targets = torch.full(
                (cupro_logits.shape[0],), 0.6,
                device=logits.device, dtype=torch.float32
            )
            # 使用 BCEWithLogitsLoss 配合软标签
            cupro_probs = F.softmax(cupro_logits, dim=1)[:, 1]
            cupro_loss = F.binary_cross_entropy(cupro_probs, cupro_targets)
            loss = loss + cupro_weight * cupro_loss

    # 检查 NaN
    if torch.isnan(loss):
        logger.warning("Loss 为 NaN，跳过此 epoch")
        return float('nan')

    loss.backward()
    torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
    optimizer.step()

    return loss.item()


@torch.no_grad()
def evaluate(model, x, edge_index, edge_attr, y, mask, pos_weight=1.0) -> Dict[str, float]:
    """评估模型，返回指标和损失"""
    model.eval()
    logits = model(x, edge_index, edge_attr)
    probs = F.softmax(logits, dim=1)

    mask_logits = logits[mask]
    mask_probs = probs[mask]
    mask_labels = y[mask]

    preds = mask_logits.argmax(dim=1).cpu().numpy()
    labels = mask_labels.cpu().numpy()
    prob_pos = mask_probs[:, 1].cpu().numpy()

    metrics = compute_metrics(labels, preds, prob_pos)

    # 计算验证集/测试集损失
    weight = torch.tensor([1.0, pos_weight], device=logits.device)
    loss = F.cross_entropy(mask_logits, mask_labels, weight=weight)
    metrics["loss"] = loss.item()

    return metrics


def main():
    parser = argparse.ArgumentParser(description="Train GATv2 for target prediction")
    parser.add_argument("--config", type=str, default="config.yaml", help="Path to config file")
    args = parser.parse_args()

    # 加载配置
    config = load_config(args.config)
    data_cfg = config["data"]
    model_cfg = config["model"]
    train_cfg = config["training"]
    split_cfg = config["split"]
    repro_cfg = config["reproducibility"]
    out_cfg = config["output"]
    eval_cfg = config["evaluation"]

    # 设置随机种子
    set_seed(repro_cfg["torch_seed"])

    # 创建输出目录
    model_dir = Path(out_cfg["model_dir"])
    results_dir = Path(out_cfg["results_dir"])
    logs_dir = Path(out_cfg["logs_dir"])
    for d in [model_dir, results_dir, logs_dir]:
        d.mkdir(parents=True, exist_ok=True)

    setup_logging(logs_dir)

    # 设备
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    logger.info(f"使用设备: {device}")
    if device.type == "cuda":
        logger.info(f"GPU: {torch.cuda.get_device_name(0)}")
        logger.info(f"显存: {torch.cuda.get_device_properties(0).total_memory / 1024**3:.1f} GB")

    # 加载数据（GATv2 标准范式：返回 edge_attr）
    processed_dir = data_cfg["processed_dir"]
    x, edge_index, edge_attr, gene_symbols = load_graph_data(processed_dir)
    y = load_labels(processed_dir, gene_symbols)

    # 创建 mask
    train_mask, val_mask, test_mask = create_masks(
        y,
        train_ratio=split_cfg["train_ratio"],
        val_ratio=split_cfg["val_ratio"],
        test_ratio=split_cfg["test_ratio"],
        random_state=split_cfg["random_state"],
    )

    # 移动到设备
    x = x.to(device)
    edge_index = edge_index.to(device)
    edge_attr = edge_attr.to(device)
    y = y.to(device)
    train_mask = train_mask.to(device)
    val_mask = val_mask.to(device)
    test_mask = test_mask.to(device)

    # 计算 pos_weight
    n_pos = (y[train_mask] == 1).sum().item()
    n_neg = (y[train_mask] == 0).sum().item()
    pos_weight = n_neg / max(n_pos, 1)
    logger.info(f"训练集 阳性: {n_pos}, 阴性: {n_neg}, pos_weight: {pos_weight:.2f}")

    # 铜死亡基因 mask
    cupro_genes = set(train_cfg.get("cupro_genes", [
        "FDX1", "LIAS", "LIPT1", "DLAT", "PDHB", "PDHX",
        "SLC31A1", "ATP7A", "ATP7B", "ATOX1", "NFE2L2",
        "HIF1A", "MTOR", "NFKB1", "GPX4",
    ]))
    cupro_mask = torch.tensor([g in cupro_genes for g in gene_symbols], dtype=torch.bool, device=device)
    cupro_weight = train_cfg.get("cupro_weight", 0.3)
    logger.info(f"铜死亡基因数: {cupro_mask.sum().item()}, 先验损失权重: {cupro_weight}")

    # 初始化模型（GATv2 标准范式）
    model = GAT(
        in_channels=model_cfg["in_channels"],
        hidden_channels=model_cfg["hidden_channels"],
        out_channels=model_cfg["out_channels"],
        num_heads=model_cfg["num_heads"],
        num_classes=model_cfg["num_classes"],
        dropout=model_cfg["dropout"],
        attention_dropout=model_cfg["attention_dropout"],
        use_edge_attr=model_cfg["use_edge_attr"],
        use_batch_norm=model_cfg["use_batch_norm"],
    ).to(device)

    logger.info(f"模型参数: {sum(p.numel() for p in model.parameters()):,}")

    # 优化器
    optimizer = AdamW(
        model.parameters(),
        lr=train_cfg["learning_rate"],
        weight_decay=train_cfg["weight_decay"],
    )

    # 学习率调度器
    scheduler_type = train_cfg["scheduler"]
    if scheduler_type == "cosine":
        scheduler = CosineAnnealingLR(optimizer, T_max=train_cfg["scheduler_t_max"])
    elif scheduler_type == "plateau":
        scheduler = ReduceLROnPlateau(optimizer, mode="max", factor=train_cfg["scheduler_factor"],
                                       patience=train_cfg["scheduler_patience"], verbose=True)
    else:
        scheduler = None

    # 训练历史
    history = {
        "train_loss": [],
        "val_loss": [],
        "val_f1": [],
        "val_roc_auc": [],
        "lr": [],
    }

    best_val_roc_auc = 0.0
    best_epoch = 0
    patience_counter = 0
    epochs = train_cfg["epochs"]
    patience = train_cfg["patience"]

    logger.info("=" * 60)
    logger.info("开始训练（早停指标: ROC-AUC）")
    logger.info("=" * 60)

    for epoch in range(1, epochs + 1):
        # 训练
        loss = train_epoch(model, optimizer, x, edge_index, edge_attr, y, train_mask, pos_weight,
                           cupro_mask=cupro_mask, cupro_weight=cupro_weight)

        if np.isnan(loss):
            # 自动降低学习率
            for param_group in optimizer.param_groups:
                param_group["lr"] *= 0.5
            logger.warning(f"Epoch {epoch}: Loss NaN，学习率降至 {optimizer.param_groups[0]['lr']:.6f}")
            continue

        history["train_loss"].append(loss)

        # 验证
        val_metrics = evaluate(model, x, edge_index, edge_attr, y, val_mask, pos_weight=pos_weight)
        val_f1 = val_metrics["f1"]
        val_roc_auc = val_metrics["roc_auc"]
        val_loss = val_metrics["loss"]
        history["val_f1"].append(val_f1)
        history["val_roc_auc"].append(val_roc_auc)
        history["val_loss"].append(val_loss)
        history["lr"].append(optimizer.param_groups[0]["lr"])

        # 学习率调度
        if scheduler_type == "cosine":
            scheduler.step()
        elif scheduler_type == "plateau":
            scheduler.step(val_roc_auc)

        # 保存最佳模型（基于 ROC-AUC）
        if val_roc_auc > best_val_roc_auc:
            best_val_roc_auc = val_roc_auc
            best_epoch = epoch
            patience_counter = 0
            torch.save({
                "epoch": epoch,
                "model_state_dict": model.state_dict(),
                "optimizer_state_dict": optimizer.state_dict(),
                "val_roc_auc": val_roc_auc,
                "val_f1": val_f1,
                "config": config,
            }, model_dir / out_cfg["best_model"])
            logger.info(f"Epoch {epoch:03d}: 新最佳模型，Val ROC-AUC={val_roc_auc:.4f}, Val F1={val_f1:.4f}")
        else:
            patience_counter += 1

        # 每 10 epoch 打印日志
        if epoch % 10 == 0 or epoch == 1:
            logger.info(
                f"Epoch {epoch:03d} | Loss: {loss:.4f} | Val Loss: {val_loss:.4f} | "
                f"Val ROC-AUC: {val_roc_auc:.4f} | Val F1: {val_f1:.4f} | "
                f"LR: {optimizer.param_groups[0]['lr']:.6f}"
            )

        # 早停
        if patience_counter >= patience:
            logger.info(f"早停触发，最佳 epoch: {best_epoch}, 最佳 Val ROC-AUC: {best_val_roc_auc:.4f}")
            break

    logger.info("=" * 60)
    logger.info("训练完成")
    logger.info("=" * 60)

    # 加载最佳模型进行测试
    checkpoint = torch.load(model_dir / out_cfg["best_model"], map_location=device)
    model.load_state_dict(checkpoint["model_state_dict"])
    logger.info(f"加载最佳模型 (epoch {checkpoint['epoch']}, val_f1={checkpoint['val_f1']:.4f})")

    # 测试集评估
    test_metrics = evaluate(model, x, edge_index, edge_attr, y, test_mask)
    logger.info("测试集指标:")
    for k, v in test_metrics.items():
        logger.info(f"  {k}: {v:.4f}")

    # 混淆矩阵
    if eval_cfg["plot_confusion_matrix"]:
        model.eval()
        with torch.no_grad():
            test_logits = model(x, edge_index, edge_attr)
            test_preds = test_logits[test_mask].argmax(dim=1).cpu().numpy()
            test_labels = y[test_mask].cpu().numpy()
        plot_confusion_matrix(test_labels, test_preds, str(logs_dir / "confusion_matrix.png"))

    # 训练曲线
    if eval_cfg["plot_training_curves"]:
        plot_training_curves(history, str(logs_dir / "training_curves.png"))

    logger.info("=" * 60)
    logger.info("GATv2 训练流程完成")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
