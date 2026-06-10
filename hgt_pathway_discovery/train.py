# -*- coding: utf-8 -*-
"""训练模块：增强版训练策略。

改进:
  1. 复合早停指标: 0.7*AUROC + 0.3*AUPRC
  2. 自适应负采样比: 余弦退火 1:1 → 10:1
  3. Platt Scaling 校准: 验证集学习校准参数
  4. 验证负边定时刷新: 每 N epoch 重采样
  5. TensorBoard 日志记录
  6. SWA (Stochastic Weight Averaging) 增强泛化
  7. NeighborLoader mini-batch: 适用于大规模图的子图采样训练
  8. WandB 实验追踪
"""

import copy
import sys
import math
from typing import Dict, Tuple, Optional, List
from pathlib import Path

import numpy as np
import torch
import torch.nn.functional as F
from torch import Tensor, nn
from torch.optim import AdamW
from torch.optim.lr_scheduler import ReduceLROnPlateau
from torch_geometric.data import HeteroData
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.linear_model import LogisticRegression

from .model import HGTModel, focal_bce_loss, NegEdgeSampler


# ============================================================================
# 子图工具: 从消息传递图中剔除指定边
# ============================================================================

def remove_edges_from_data(data: HeteroData, edge_type: Tuple[str, str, str],
                            mask: Tensor) -> HeteroData:
    """轻量级图复制：只复制边索引，节点特征张量共享内存。"""
    new_data = HeteroData()
    for nt in data.node_types:
        new_data[nt].x = data[nt].x
    for et in data.edge_types:
        if et == edge_type:
            new_data[et].edge_index = data[et].edge_index[:, mask]
        else:
            new_data[et].edge_index = data[et].edge_index
    return new_data


# ============================================================================
# 自适应负采样比
# ============================================================================

def cosine_neg_ratio(epoch: int, total_epochs: int,
                     start_ratio: int = 1, end_ratio: int = 10) -> int:
    """余弦退火负采样比: 训练初期低比例，后期逐步提高。

    公式: ratio = end + (start - end) * (1 + cos(pi * epoch / total)) / 2
    """
    progress = epoch / max(total_epochs - 1, 1)
    factor = (1 + math.cos(math.pi * progress)) / 2
    ratio = int(round(end_ratio + (start_ratio - end_ratio) * factor))
    return max(1, ratio)


# ============================================================================
# Platt Scaling 校准
# ============================================================================

class PlattScaler:
    """Platt Scaling: 使用逻辑回归将模型 logits 校准为概率。

    参考: Platt, "Probabilistic Outputs for SVMs", 1999
    """

    def __init__(self):
        self._calibrator: Optional[LogisticRegression] = None

    def fit(self, logits: np.ndarray, labels: np.ndarray) -> None:
        """在验证集 logits 上拟合校准器。"""
        logits = np.asarray(logits, dtype=np.float64).reshape(-1, 1)
        labels = np.asarray(labels, dtype=np.float64)
        self._calibrator = LogisticRegression(
            C=1.0, solver="lbfgs", max_iter=1000,
        )
        self._calibrator.fit(logits, labels)

    def predict_proba(self, logits: np.ndarray) -> np.ndarray:
        """校准后的概率。"""
        if self._calibrator is None:
            return 1.0 / (1.0 + np.exp(-np.asarray(logits)))
        logits = np.asarray(logits, dtype=np.float64).reshape(-1, 1)
        return self._calibrator.predict_proba(logits)[:, 1]

    @property
    def is_fitted(self) -> bool:
        return self._calibrator is not None


# ============================================================================
# 评估函数
# ============================================================================

@torch.inference_mode()
def evaluate(model: nn.Module, data: HeteroData, gp_edge_index: Tensor,
             gp_pos_idx: Tensor, neg_sampler: NegEdgeSampler,
             neg_ratio: int = 3, eval_batch: int = 65536,
             force_refresh_val_neg: bool = False) -> Tuple[float, float, np.ndarray, np.ndarray]:
    """评估模型在验证集上的 AUROC 和 AUPRC。

    Args:
        force_refresh_val_neg: 强制刷新验证负边集
    """
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

    scores_list: List[np.ndarray] = []
    for start in range(0, n_total, eval_batch):
        end = min(start + eval_batch, n_total)
        batch_logits = model.decode(z_dict, eval_ei[:, start:end])
        scores_list.append(torch.sigmoid(batch_logits).cpu().numpy())
    scores = np.concatenate(scores_list)

    auroc = roc_auc_score(labels, scores) if len(np.unique(labels)) > 1 else 0.5
    auprc = average_precision_score(labels, scores)
    return auroc, auprc, scores, labels


def compute_composite_metric(auroc: float, auprc: float,
                              weight: float = 0.7) -> float:
    """复合指标: weight*AUROC + (1-weight)*AUPRC"""
    return weight * auroc + (1 - weight) * auprc


# ============================================================================
# 单折训练
# ============================================================================

def train_fold(model: nn.Module, data_train: HeteroData, gp_edge_index: Tensor,
               train_idx: Tensor, val_idx: Tensor,
               neg_sampler: NegEdgeSampler,
               config: object,
               tb_writer: Optional[object] = None,
               fold_idx: int = 0) -> Tuple[float, float, Optional[PlattScaler]]:
    """单折训练循环。

    Returns:
        best_auroc, best_auprc, platt_scaler
    """
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

    # SWA: 从 patience_cnt 稳定期开始累积，而非固定最后N轮
    # 当 patience_cnt 达到 patience//2 且未被重置时，模型已进入性能平台期
    swa_model: Optional[Dict] = None
    swa_n = 0
    swa_active = False

    for epoch in range(cfg.epochs):
        model.train()

        if epoch == 0:
            print(f"    [DEBUG] Starting training loop, device={device}, "
                  f"train_edges={train_ei.shape[1]}, neg_ratio={neg_ratio}",
                  flush=True)

        # 自适应负采样比
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

        # 评估
        eval_interval = 5
        if (epoch + 1) % eval_interval == 0:
            force_refresh = (epoch + 1) % cfg.val_neg_refresh_interval == 0
            auroc, auprc, val_scores, val_labels = evaluate(
                model, data_train_device, gp_edge_index, val_idx,
                neg_sampler, neg_ratio, cfg.eval_batch,
                force_refresh_val_neg=force_refresh,
            )

            composite = compute_composite_metric(auroc, auprc, cfg.composite_weight)

            # 早停判断
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
                swa_active = False  # 性能提升时重置 SWA
            else:
                patience_cnt += 1
                # SWA 激活: patience 过半且未触发早停 → 模型已进入稳定平台期
                if not swa_active and patience_cnt >= cfg.patience // 2:
                    swa_active = True
                if patience_cnt >= cfg.patience:
                    break

            scheduler.step(auroc)

            # TensorBoard
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

        # SWA 记录：仅在稳定平台期累积权重
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

    # SWA vs Best 比较
    if swa_model is not None and swa_n > 1:
        for key in swa_model:
            if swa_model[key].dtype in (torch.float32, torch.float64):
                swa_model[key].data.div_(swa_n)
        model.load_state_dict(swa_model)
        swa_auroc, swa_auprc, _, _ = evaluate(
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

    # Platt Scaling 校准
    platt_scaler: Optional[PlattScaler] = None
    if config.inference.calibrate:
        _, _, val_scores, val_labels = evaluate(
            model, data_train_device, gp_edge_index, val_idx,
            neg_sampler, neg_ratio, cfg.eval_batch,
        )
        platt_scaler = PlattScaler()
        platt_scaler.fit(val_scores, val_labels)

    return best_auroc, best_auprc, platt_scaler


# ============================================================================
# Mini-batch 训练 (NeighborLoader)
# ============================================================================

def train_fold_minibatch(model: nn.Module, data_train: HeteroData,
                          gp_edge_index: Tensor,
                          train_idx: Tensor, val_idx: Tensor,
                          neg_sampler: NegEdgeSampler,
                          config: object,
                          tb_writer: Optional[object] = None,
                          wandb_run: Optional[object] = None,
                          fold_idx: int = 0) -> Tuple[float, float, Optional[PlattScaler]]:
    """单折训练循环 (NeighborLoader 子图采样版本)。

    适用于大规模图 (>100万边) 的 mini-batch 训练。
    使用 PyG 的 HGTLoader 按 batch 采样子图进行消息传递。

    参考:
      - Hamilton et al., "Inductive Representation Learning on Large Graphs", NeurIPS 2017
      - PyG HGTLoader: https://pytorch-geometric.readthedocs.io/en/latest/modules/loader.html
    """
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

    # HGTLoader: 对每种节点类型采样指定数量的邻居
    neighbor_sizes_dict = {}
    for nt in data_train.node_types:
        neighbor_sizes_dict[nt] = cfg.neighbor_sizes[:model_cfg.num_layers]

    # 构建训练边的输入节点列表（基因和通路节点对）
    train_src_nodes = train_ei[0].tolist()
    train_dst_nodes = train_ei[1].tolist()

    use_scaler = cfg.use_grad_scaling and device.type == "cuda"
    scaler = torch.amp.GradScaler(device="cuda") if use_scaler else None

    swa_model: Optional[Dict] = None
    swa_n = 0
    swa_active = False

    n_batches = max(1, len(train_src_nodes) // cfg.batch_size + 1)

    for epoch in range(cfg.epochs):
        model.train()

        # 自适应负采样比
        if cfg.adaptive_neg_ratio:
            neg_ratio = cosine_neg_ratio(
                epoch, cfg.epochs,
                cfg.neg_sample_ratio_start, cfg.neg_sample_ratio_end,
            )

        # 打乱训练样本顺序
        perm = torch.randperm(len(train_src_nodes))
        total_loss = 0.0
        n_batches_processed = 0

        for batch_start in range(0, len(train_src_nodes), cfg.batch_size):
            batch_end = min(batch_start + cfg.batch_size, len(train_src_nodes))
            batch_perm = perm[batch_start:batch_end]

            # 采样子图
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

        # 评估 (使用全图)
        eval_interval = 5
        if (epoch + 1) % eval_interval == 0:
            force_refresh = (epoch + 1) % cfg.val_neg_refresh_interval == 0
            auroc, auprc, val_scores, val_labels = evaluate(
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

        # SWA 记录：仅在稳定平台期累积权重
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

    # SWA vs Best 比较
    if swa_model is not None and swa_n > 1:
        for key in swa_model:
            if swa_model[key].dtype in (torch.float32, torch.float64):
                swa_model[key].data.div_(swa_n)
        model.load_state_dict(swa_model)
        swa_auroc, swa_auprc, _, _ = evaluate(
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
        _, _, val_scores, val_labels = evaluate(
            model, data_train_device, gp_edge_index, val_idx,
            neg_sampler, neg_ratio, cfg.eval_batch,
        )
        platt_scaler = PlattScaler()
        platt_scaler.fit(val_scores, val_labels)

    return best_auroc, best_auprc, platt_scaler

def cross_validate(data: HeteroData, gp_edge_index: Tensor,
                   n_genes: int, n_pathways: int,
                   gene_degrees: Optional[np.ndarray] = None,
                   pathway_degrees: Optional[np.ndarray] = None,
                   config: Optional[object] = None,
                   tb_writer: Optional[object] = None,
                   wandb_run: Optional[object] = None) -> Tuple[List[Dict], List[HGTModel], List[PlattScaler], Tensor]:
    """K-Fold 交叉验证，严格数据隔离。

    Returns:
        cv_scores: 每折评估结果
        cv_models: 训练好的模型列表
        cv_scalers: 每折的 Platt 校准器
        cv_heldout_edges: 所有折的验证边索引拼接 (用于归纳评估)
    """
    from sklearn.model_selection import KFold

    cv_cfg = config.cv
    n_edges = gp_edge_index.size(1)
    if n_edges < cv_cfg.n_folds:
        return [], [], []

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

        from .model import HGTModel as HGTModelCls
        model_cfg = config.model
        model = HGTModelCls(
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

        # CpG 可学习参数注册
        if "cpg" in data.node_types:
            model.to_cpg_learnable(
                data["cpg"].x,
                quality_mask=data["cpg"].propagation_mask if hasattr(data["cpg"], "propagation_mask") else None,
            )

        print(f"    [DEBUG] Model created and on {config.device}, creating sampler...", flush=True)

        # torch.compile 加速 (PyTorch >= 2.0, 需要 Triton)
        if hasattr(torch, "compile") and config.device == "cuda":
            try:
                import triton  # noqa: F401
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
        cv_heldout_edges_list.append(gp_edge_index[:, val_idx])  # 保留验证边用于归纳评估

        print(f"[CV] Fold {fold+1}: AUROC={auroc:.4f}, AUPRC={auprc:.4f}")

        # 保存每折模型
        if cv_cfg.save_fold_models:
            torch.save(model.state_dict(), str(folds_dir / f"fold_{fold+1}.pt"))

        if config.device == "cuda":
            torch.cuda.empty_cache()

    cv_heldout_edges = torch.cat(cv_heldout_edges_list, dim=1) if cv_heldout_edges_list else torch.zeros((2, 0), dtype=torch.long)
    return cv_scores, cv_models, cv_scalers, cv_heldout_edges


# ============================================================================
# 最终训练（全量数据）
# ============================================================================

def train_final(model: nn.Module, data: HeteroData, gp_edge_index: Tensor,
                neg_sampler: NegEdgeSampler,
                config: object,
                tb_writer: Optional[object] = None) -> Tuple[float, float, Optional[PlattScaler]]:
    """使用全量数据训练最终模型（保留 10% 作为监控验证集）。

    Returns:
        best_auroc, best_auprc, platt_scaler
    """
    device = next(model.parameters()).device
    cfg = config.training

    n_edges = gp_edge_index.size(1)
    val_size = max(1, int(n_edges * cfg.val_ratio))
    perm = torch.randperm(n_edges, generator=torch.Generator().manual_seed(config.seed + 999))
    train_idx = perm[val_size:]
    val_idx = perm[:val_size]

    val_mask = torch.zeros(n_edges, dtype=torch.bool)
    val_mask[val_idx] = True
    data_train = remove_edges_from_data(
        data, ("gene", "involved_in", "pathway"), ~val_mask,
    )

    auroc, auprc, platt_scaler = train_fold(
        model, data_train.to(device), gp_edge_index.to(device),
        train_idx.long(), val_idx.long(),
        neg_sampler, config, tb_writer, fold_idx=999,
    )

    print(f"  Final val: AUROC={auroc:.4f}, AUPRC={auprc:.4f}")
    return auroc, auprc, platt_scaler