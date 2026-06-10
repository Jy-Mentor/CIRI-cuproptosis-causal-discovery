# -*- coding: utf-8 -*-
"""推理模块：转导推理 + 归纳评估 + MC Dropout 不确定性 + 集成预测。

改进:
  1. 归纳子集独立评估: 计算 AUROC/AUPRC (若有真实标签)
  2. MC Dropout: 多次前向传播估计预测方差
  3. Platt Scaling 校准概率
  4. 集成预测: 平均多折 CV 模型 logits
"""

from typing import Dict, List, Tuple, Optional, Set
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from torch import Tensor, nn
from torch_geometric.data import HeteroData
from sklearn.metrics import roc_auc_score, average_precision_score

from .train import PlattScaler


@torch.inference_mode()
def _encode_once(model: nn.Module, data_device: HeteroData) -> Dict[str, Tensor]:
    """编码全图一次（复用嵌入）。"""
    return model(data_device.x_dict, data_device.edge_index_dict)


def mc_dropout_predict_per_gene(model: nn.Module, data_device: HeteroData,
                                gene_idx: int, pathway_idx: Tensor,
                                n_samples: int = 10, temperature: float = 1.0) -> Tuple[np.ndarray, np.ndarray]:
    """MC Dropout 预测（解码器级）：一次编码 + 多次解码器前向传播。

    注意：此函数仅对解码器执行 MC Dropout，编码器仅编码一次。
    若需编码器级 MC Dropout，请使用 mc_dropout_predict_full。

    参考: Gal & Ghahramani, "Dropout as a Bayesian Approximation",
          ICML 2016

    Args:
        model: 训练好的模型
        data_device: 图数据（已在目标设备）
        gene_idx: 基因索引
        pathway_idx: 通路索引列表
        n_samples: MC 采样次数
        temperature: 温度校准

    Returns:
        mean_scores: 平均预测分数
        std_scores: 预测标准差（不确定性估计）
    """
    z_dict = _encode_once(model, data_device)

    edge_idx = torch.stack([
        torch.full((len(pathway_idx),), gene_idx, dtype=torch.long, device=data_device["gene"].x.device),
        pathway_idx,
    ])

    model.train()
    all_scores = []
    for _ in range(n_samples):
        logits = model.decode(z_dict, edge_idx)
        scores = torch.sigmoid(logits / temperature).cpu().numpy()
        all_scores.append(scores)

    model.eval()
    all_scores = np.stack(all_scores, axis=0)
    mean_scores = all_scores.mean(axis=0)
    std_scores = all_scores.std(axis=0)
    return mean_scores, std_scores


@torch.inference_mode()
def predict_bridge_pathways(model: nn.Module, data: HeteroData,
                             bridge_genes: List[str],
                             gene_to_idx: Dict[str, int],
                             pathway_names: List[str],
                             gene_idx_to_name: Optional[Dict[int, str]] = None,
                             gp_edge_index: Optional[Tensor] = None,
                             config: Optional[object] = None,
                             platt_scaler: Optional[PlattScaler] = None) -> pd.DataFrame:
    """转导推理: 使用全图预测所有桥梁基因的通路关联分数。

    Args:
        model: 训练好的 HGTModel
        data: 全图 HeteroData
        bridge_genes: 桥梁基因列表
        gene_to_idx: 基因→索引映射
        pathway_names: 通路名列表
        gene_idx_to_name: 基因索引→名称映射 (用于归纳评估)
        gp_edge_index: 基因-通路边索引 (用于划分转导/归纳)
        config: 全局配置
        platt_scaler: Platt 校准器 (可选)

    Returns:
        df: 包含 gene_symbol, pathway_name, score, eval_mode, rank 的 DataFrame
    """
    inf_cfg = config.inference if config else None
    top_k = inf_cfg.top_k if inf_cfg else 10
    temperature = inf_cfg.temperature if inf_cfg else 1.5
    mc_samples = inf_cfg.mc_dropout_samples if inf_cfg else 0
    calibrate = inf_cfg.calibrate if inf_cfg else True

    model.eval()
    device = next(model.parameters()).device
    data_device = data if data["gene"].x.device == device else data.to(device)

    # 确定已知通路关联基因 (转导 vs 归纳)
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

    # 编码全图
    z_dict = _encode_once(model, data_device)

    results: List[Dict] = []
    n_pathways = len(pathway_names)
    pathway_idx_tensor = torch.arange(n_pathways, dtype=torch.long, device=device)

    for gene in bridge_genes:
        if gene not in gene_to_idx:
            continue
        gi = gene_to_idx[gene]

        # MC Dropout 预测 (优先)
        if mc_samples > 0:
            scores, stds = mc_dropout_predict_per_gene(
                model, data_device, gi, pathway_idx_tensor,
                n_samples=mc_samples, temperature=temperature,
            )
            has_uncertainty = True

            # MC Dropout 后可选 Platt 校准
            if calibrate and platt_scaler is not None and platt_scaler.is_fitted:
                scores = platt_scaler.predict_proba(scores)
        else:
            edge_idx = torch.stack([
                torch.full((n_pathways,), gi, dtype=torch.long, device=device),
                pathway_idx_tensor,
            ])
            logits = model.decode(z_dict, edge_idx)

            # Platt 校准
            if calibrate and platt_scaler is not None and platt_scaler.is_fitted:
                scores = platt_scaler.predict_proba(logits.cpu().numpy())
            else:
                scores = torch.sigmoid(logits / temperature).cpu().numpy()
            has_uncertainty = False

        is_inductive = gene not in known_gp_genes if known_gp_genes else None

        for pi, pname in enumerate(pathway_names):
            row = {
                "gene_symbol": gene,
                "pathway_name": pname,
                "score": float(scores[pi]),
                "eval_mode": "inductive" if is_inductive else "transductive",
            }
            if has_uncertainty:
                row["uncertainty"] = float(stds[pi])
            results.append(row)

    df = pd.DataFrame(results)
    if df.empty:
        print("[Predict] No bridge genes found in graph!")
        return df

    df["rank"] = df.groupby("gene_symbol")["score"].rank(ascending=False, method="dense")
    df = df.sort_values(["gene_symbol", "rank"]).reset_index(drop=True)

    # 保存
    save_path = Path(config.paths.project_dir) / "bridge_pathway_scores.csv"
    df.to_csv(save_path, index=False)
    calib_tag = "+calibrated" if (calibrate and platt_scaler is not None and platt_scaler.is_fitted) else ""
    print(f"[Output] Saved bridge_pathway_scores.csv ({len(df)} rows, T={temperature}{calib_tag})")

    # 打印 Top-K
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

    # 归纳统计
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
    """对归纳子集（无已知通路关联的 bridge 基因）计算 AUROC/AUPRC。

    若提供 heldout_gp_edges（从 CV 保留的验证边），可构建真实正负标签
    并计算定量指标。否则仅报告预测统计。

    Args:
        heldout_gp_edges: (2, M) 保留的基因-通路验证边索引，用于构建
                          归纳评估的真实标签。若为 None，则无法计算 AUROC/AUPRC。

    Returns:
        induct_auroc, induct_auprc (若无法计算则返回 -1, -1)
    """
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

    # 若提供了 heldout 边，构建真实标签进行定量评估
    if heldout_gp_edges is not None and heldout_gp_edges.size(1) > 0:
        # 提取属于归纳基因的 heldout 边
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
    """集成多模型预测：平均所有 CV 模型 logits。

    参考: Lakshminarayanan et al., "Simple and Scalable Predictive
          Uncertainty Estimation using Deep Ensembles", NeurIPS 2017
    """
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

    # 第一遍：确定基因顺序
    for gene in bridge_genes:
        if gene in gene_to_idx:
            gene_order.append(gene)

    # 多模型平均
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

    # Top-K 打印
    topk = df[df["rank"] <= top_k].copy()
    print(f"\n{'='*60}")
    print(f"  Top-{top_k} Pathways per Bridge Gene (Ensemble of {len(models)} models)")
    print(f"{'='*60}")
    for gene in gene_order[:30]:  # 仅打印前30个基因
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