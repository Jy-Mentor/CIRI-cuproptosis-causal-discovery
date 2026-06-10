# -*- coding: utf-8 -*-
"""异构图构建模块：增强版，支持有向甲基化、PPI权重拆分、pathway层级边。

改进:
  1. normalize_node_features: 在图构建后对每种节点特征做 Z-score 归一化，
     并存储 mean/std 到 HeteroData 以便推理时复用
  2. 防止 HGT 注意力系数被量级偏差主导 (gene/drug/disease/pathway)
"""

from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set

import numpy as np
import torch
from torch_geometric.data import HeteroData

from .data_loader import load_txt
from .utils import normalize_features, zscore_normalize


def normalize_node_features(data: HeteroData, method: str = "zscore",
                             eps: float = 1e-8) -> HeteroData:
    """对 HeteroData 中每种节点特征做归一化，并存储 mean/std 以便推理复用。

    基因、药物、疾病、通路等节点特征数值范围差异大（指纹 0-1，表达量波动，
    疾病特征是基因均值），直接输入 HGT 会导致注意力系数被量级偏差主导。

    改进: 在图构建后立即归一化，并将统计量保存在 data[node_type].norm_stats
    中，推理时可精确复用。跳过单节点类型 (drug/disease) 的全零标准差。

    Returns:
        data: 归一化后的 HeteroData (原地修改)
    """
    for nt in data.node_types:
        x = data[nt].x
        if x.size(0) == 0:
            continue

        if method == "zscore":
            mean = x.mean(dim=0)
            std = x.std(dim=0, unbiased=False)
            std[std < eps] = 1.0
            data[nt].x = (x - mean) / std
            data[nt].norm_stats = {"mean": mean, "std": std, "method": "zscore"}
        elif method == "minmax":
            vmin = x.min(dim=0)[0]
            vmax = x.max(dim=0)[0]
            denom = vmax - vmin
            denom[denom < eps] = 1.0
            data[nt].x = (x - vmin) / denom
            data[nt].norm_stats = {"min": vmin, "max": vmax, "method": "minmax"}
        else:
            data[nt].norm_stats = {"method": "none"}

    return data


def apply_cached_normalization(data: HeteroData) -> HeteroData:
    """使用缓存的归一化统计量对推理数据重新归一化。

    用于推理阶段，当输入数据可能与训练时分布不同时，
    复用训练时的 mean/std 进行归一化。
    """
    for nt in data.node_types:
        if hasattr(data[nt], "norm_stats") and data[nt].norm_stats.get("method") == "zscore":
            stats = data[nt].norm_stats
            data[nt].x = (data[nt].x - stats["mean"]) / stats["std"]
        elif hasattr(data[nt], "norm_stats") and data[nt].norm_stats.get("method") == "minmax":
            stats = data[nt].norm_stats
            data[nt].x = (data[nt].x - stats["min"]) / (stats["max"] - stats["min"])
    return data


def build_hetero_graph(
    gene_feat_arr: np.ndarray,
    gene_feat_names: List[str],
    drug_fp_arr: np.ndarray,
    disease_feat_arr: Optional[np.ndarray],
    pathway_feat_arr: np.ndarray,
    pathway_names: List[str],
    ppi_edges: List[Tuple[str, str, float]],
    coexp_edges: List[Tuple[str, str]],
    tf_edges: List[Tuple[str, str]],
    gene_pathway_edges: List[Tuple[str, str]],
    all_genes_list: Optional[List[str]] = None,
    bridge_genes: Optional[List[str]] = None,
    methyl_edges: Optional[List[Tuple[str, str]]] = None,
    mirna_edges: Optional[List[Tuple[str, str]]] = None,
    pathway_hierarchy: Optional[List[Tuple[str, str]]] = None,
    disease_pathway_edges: Optional[List[Tuple[str, str]]] = None,
    config: Optional[object] = None,
) -> Tuple[HeteroData, Dict[str, int], List[str], Dict[str, int]]:
    """构建增强型多组学异构图。

    Returns:
        data: 构建好的 HeteroData
        gene_to_idx: 基因名→索引映射
        gene_list: 基因名列表
        pathway_name_to_idx: 通路名→索引映射
    """
    graph_cfg = config.graph if config else None
    preproc_cfg = config.preprocessing if config else None

    # --- 确定图节点集合 ---
    gene_set = set(gene_feat_names)
    if all_genes_list:
        gene_set &= set(all_genes_list)

    for a, b, _ in ppi_edges:
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

    # --- 基因特征提取与标准化 ---
    gene_feat_dict = dict(zip(gene_feat_names, gene_feat_arr))
    gene_feat = np.zeros((n_genes, gene_feat_arr.shape[1]), dtype=np.float32)
    for i, g in enumerate(gene_list):
        if g in gene_feat_dict:
            gene_feat[i] = gene_feat_dict[g]

    # --- 通路特征标准化 ---
    pathway_feat = pathway_feat_arr.astype(np.float32).copy()
    if pathway_feat.ndim == 1:
        pathway_feat = pathway_feat.reshape(1, -1)

    # --- 构建 HeteroData（先填充原始特征，再统一归一化） ---
    data = HeteroData()

    # 基因节点
    data["gene"].x = torch.from_numpy(gene_feat)

    # 药物节点
    data["drug"].x = torch.from_numpy(drug_fp_arr.reshape(1, -1).astype(np.float32))

    # 疾病节点
    if disease_feat_arr is not None:
        dis_feat = disease_feat_arr.reshape(1, -1)
    else:
        dis_feat = np.mean(gene_feat, axis=0, keepdims=True)
    data["disease"].x = torch.from_numpy(dis_feat.astype(np.float32))

    # 通路节点
    data["pathway"].x = torch.from_numpy(pathway_feat)

    # --- Z-score 标准化（在 tensor 上执行，存储 mean/std） ---
    normalize_method = preproc_cfg.normalize if preproc_cfg else "zscore"
    if normalize_method != "none":
        data = normalize_node_features(data, normalize_method)
        print(f"[Build] Node features normalized: method={normalize_method}")

    pathway_name_to_idx = {name: i for i, name in enumerate(pathway_names)}

    # ========================================================================
    # PPI 边 - 按置信度拆分为 strong_ppi / weak_ppi
    # ========================================================================
    ppi_strong_threshold = graph_cfg.ppi_strong_threshold if graph_cfg else 900
    ppi_split = graph_cfg.ppi_split_by_score if graph_cfg else True

    if ppi_split:
        strong_src, strong_dst = [], []
        weak_src, weak_dst = [], []
        for a, b, score in ppi_edges:
            if a not in gene_to_idx or b not in gene_to_idx:
                continue
            if score >= ppi_strong_threshold:
                strong_src.extend([gene_to_idx[a], gene_to_idx[b]])
                strong_dst.extend([gene_to_idx[b], gene_to_idx[a]])
            else:
                weak_src.extend([gene_to_idx[a], gene_to_idx[b]])
                weak_dst.extend([gene_to_idx[b], gene_to_idx[a]])

        data["gene", "strong_ppi", "gene"].edge_index = (
            torch.tensor([strong_src, strong_dst], dtype=torch.long)
            if strong_src else torch.zeros((2, 0), dtype=torch.long)
        )
        data["gene", "weak_ppi", "gene"].edge_index = (
            torch.tensor([weak_src, weak_dst], dtype=torch.long)
            if weak_src else torch.zeros((2, 0), dtype=torch.long)
        )
        print(f"[Build] strong_ppi: {len(strong_src)}, weak_ppi: {len(weak_src)}")
    else:
        ppi_src, ppi_dst = [], []
        for a, b, _ in ppi_edges:
            if a in gene_to_idx and b in gene_to_idx:
                ppi_src.extend([gene_to_idx[a], gene_to_idx[b]])
                ppi_dst.extend([gene_to_idx[b], gene_to_idx[a]])
        data["gene", "interacts", "gene"].edge_index = (
            torch.tensor([ppi_src, ppi_dst], dtype=torch.long)
            if ppi_src else torch.zeros((2, 0), dtype=torch.long)
        )
        print(f"[Build] PPI (unified): {len(ppi_src)}")

    # ========================================================================
    # 共表达边 (无向，bridge基因优先)
    # ========================================================================
    if graph_cfg and graph_cfg.subsample_coexp and bridge_genes:
        bridge_set = set(bridge_genes)
        coe_filt = [(a, b) for a, b in coexp_edges if a in bridge_set or b in bridge_set]
        coe_used = coe_filt if coe_filt else coexp_edges
    else:
        coe_used = coexp_edges

    coe_src, coe_dst = [], []
    for a, b in coe_used:
        if a in gene_to_idx and b in gene_to_idx:
            coe_src.extend([gene_to_idx[a], gene_to_idx[b]])
            coe_dst.extend([gene_to_idx[b], gene_to_idx[a]])
    data["gene", "coexpressed", "gene"].edge_index = (
        torch.tensor([coe_src, coe_dst], dtype=torch.long)
        if coe_src else torch.zeros((2, 0), dtype=torch.long)
    )
    print(f"[Build] Coexp edges: {len(coe_src)} (used={len(coe_used)})")

    # ========================================================================
    # TF→靶基因 (有向)
    # ========================================================================
    tf_src, tf_dst = [], []
    for a, b in tf_edges:
        if a in gene_to_idx and b in gene_to_idx:
            tf_src.append(gene_to_idx[a])
            tf_dst.append(gene_to_idx[b])
    data["gene", "regulates", "gene"].edge_index = (
        torch.tensor([tf_src, tf_dst], dtype=torch.long)
        if tf_src else torch.zeros((2, 0), dtype=torch.long)
    )
    print(f"[Build] TF edges: {len(tf_src)}")

    # ========================================================================
    # 药物→基因 (有向)
    # ========================================================================
    drug_targets = load_txt(Path(config.paths.gat_data_dir) / "drug_targets.txt")
    dt_src, dt_dst = [], []
    for g in drug_targets:
        if g in gene_to_idx:
            dt_src.append(0)
            dt_dst.append(gene_to_idx[g])
    data["drug", "targets", "gene"].edge_index = (
        torch.tensor([dt_src, dt_dst], dtype=torch.long)
        if dt_src else torch.zeros((2, 0), dtype=torch.long)
    )
    print(f"[Build] drug->gene edges: {len(dt_src)}")

    # ========================================================================
    # 基因→疾病 (有向)
    # ========================================================================
    disease_genes = load_txt(Path(config.paths.gat_data_dir) / "disease_genes.txt")
    gd_src, gd_dst = [], []
    for g in disease_genes:
        if g in gene_to_idx:
            gd_src.append(gene_to_idx[g])
            gd_dst.append(0)
    data["gene", "assoc_with", "disease"].edge_index = (
        torch.tensor([gd_src, gd_dst], dtype=torch.long)
        if gd_src else torch.zeros((2, 0), dtype=torch.long)
    )
    print(f"[Build] gene->disease edges: {len(gd_src)}")

    # ========================================================================
    # 基因→通路 (有向, 监督任务)
    # ========================================================================
    gp_src, gp_dst = [], []
    for a, b in gene_pathway_edges:
        if a in gene_to_idx and b in pathway_name_to_idx:
            gp_src.append(gene_to_idx[a])
            gp_dst.append(pathway_name_to_idx[b])
    data["gene", "involved_in", "pathway"].edge_index = (
        torch.tensor([gp_src, gp_dst], dtype=torch.long)
        if gp_src else torch.zeros((2, 0), dtype=torch.long)
    )
    n_gp = len(gp_src)
    print(f"[Build] gene->pathway edges: {n_gp} (positive samples)")

    # ========================================================================
    # 通路层级边 (pathway → pathway, 有向 parent→child)
    # ========================================================================
    if pathway_hierarchy:
        ph_src, ph_dst = [], []
        for parent, child in pathway_hierarchy:
            if parent in pathway_name_to_idx and child in pathway_name_to_idx:
                ph_src.append(pathway_name_to_idx[parent])
                ph_dst.append(pathway_name_to_idx[child])
        data["pathway", "parent_of", "pathway"].edge_index = (
            torch.tensor([ph_src, ph_dst], dtype=torch.long)
            if ph_src else torch.zeros((2, 0), dtype=torch.long)
        )
        print(f"[Build] pathway hierarchy edges: {len(ph_src)}")

    # ========================================================================
    # 疾病→通路 (有向, 可选 - DisGeNET/KEGG)
    # ========================================================================
    if disease_pathway_edges:
        dp_src, dp_dst = [], []
        for disease_name, pathway_name in disease_pathway_edges:
            if pathway_name in pathway_name_to_idx:
                dp_src.append(0)  # disease node index = 0
                dp_dst.append(pathway_name_to_idx[pathway_name])
        if dp_src:
            data["disease", "assoc_with", "pathway"].edge_index = torch.tensor(
                [dp_src, dp_dst], dtype=torch.long
            )
            print(f"[Build] disease->pathway edges: {len(dp_src)}")

    # ========================================================================
    # 甲基化边 — 有向 (gene → CpG)
    #   参考: 甲基化是基因调控 CpG 位点的过程，反向传播无生物学依据
    #   CpG 节点初始化仍使用 Feature Propagation (Rossi et al., NeurIPS 2021)，
    #   但在模型中将初始特征注入 nn.Parameter 使其可学习。
    # ========================================================================
    methylation_directed = graph_cfg.methylation_directed if graph_cfg else True

    if methyl_edges:
        cpg_set: Set[str] = set()
        for g, cpg in methyl_edges:
            if g in gene_to_idx:
                cpg_set.add(cpg)
        cpg_list = sorted(cpg_set)
        cpg_to_idx = {c: i for i, c in enumerate(cpg_list)}

        if cpg_list:
            # Feature Propagation 初始化 CpG 特征
            cpg_feat = np.zeros((len(cpg_list), gene_feat.shape[1]), dtype=np.float32)
            cpg_count = np.zeros(len(cpg_list), dtype=np.float32)
            for g, cpg in methyl_edges:
                if g in gene_to_idx and cpg in cpg_to_idx:
                    cpg_feat[cpg_to_idx[cpg]] += gene_feat[gene_to_idx[g]]
                    cpg_count[cpg_to_idx[cpg]] += 1.0

            mask = cpg_count > 0
            cpg_feat[mask] /= cpg_count[mask, np.newaxis]
            cpg_feat[~mask] = gene_feat.mean(axis=0)

            # 质量掩码：标记哪些 CpG 获得了真实传播特征 (HGNN-IMA, Li et al., 2025)
            quality_mask = np.zeros((len(cpg_list), 1), dtype=np.float32)
            quality_mask[mask] = 1.0
            cpg_feat_with_quality = np.concatenate([cpg_feat, quality_mask], axis=1)

            # 存入 HeteroData — 后续模型中将此作为 CpG 可学习参数的初始值
            data["cpg"].x = torch.from_numpy(cpg_feat_with_quality).float()
            data["cpg"].propagation_mask = torch.from_numpy(mask)  # 保留原始传播标记

            gm_src, gm_dst = [], []
            for g, cpg in methyl_edges:
                if g in gene_to_idx and cpg in cpg_to_idx:
                    gm_src.append(gene_to_idx[g])
                    gm_dst.append(cpg_to_idx[cpg])

            if methylation_directed:
                # 仅保留 gene → cpg 方向
                data["gene", "methylated_at", "cpg"].edge_index = torch.tensor(
                    [gm_src, gm_dst], dtype=torch.long
                )
            else:
                # 双向 (旧行为)
                data["gene", "methylated_at", "cpg"].edge_index = torch.tensor(
                    [gm_src + gm_dst, gm_dst + gm_src], dtype=torch.long
                )

            n_isolated = int((cpg_count == 0).sum())
            edge_count = data["gene", "methylated_at", "cpg"].edge_index.size(1)
            print(f"[Build] Methylation edges: {edge_count} "
                  f"({'directed' if methylation_directed else 'undirected'}), "
                  f"CpG: {len(cpg_list)} (feature_propagation, isolated={n_isolated}, "
                  f"mean_degree={cpg_count.mean():.1f})")

    # ========================================================================
    # miRNA→基因 (有向, 可选)
    # ========================================================================
    if mirna_edges:
        mirna_set: Set[str] = set()
        for m, g in mirna_edges:
            if g in gene_to_idx:
                mirna_set.add(m)
        mirna_list = sorted(mirna_set)
        mirna_to_idx = {m: i for i, m in enumerate(mirna_list)}

        if mirna_list:
            mirna_feat = np.zeros((len(mirna_list), gene_feat.shape[1]), dtype=np.float32)
            mirna_count = np.zeros(len(mirna_list), dtype=np.float32)
            for m, g in mirna_edges:
                if g in gene_to_idx and m in mirna_to_idx:
                    mirna_feat[mirna_to_idx[m]] += gene_feat[gene_to_idx[g]]
                    mirna_count[mirna_to_idx[m]] += 1.0
            mask_m = mirna_count > 0
            mirna_feat[mask_m] /= mirna_count[mask_m, np.newaxis]
            mirna_feat[~mask_m] = gene_feat.mean(axis=0)
            data["mirna"].x = torch.from_numpy(mirna_feat).float()

            mr_src, mr_dst = [], []
            for m, g in mirna_edges:
                if g in gene_to_idx and m in mirna_to_idx:
                    mr_src.append(mirna_to_idx[m])
                    mr_dst.append(gene_to_idx[g])
            data["mirna", "targets", "gene"].edge_index = (
                torch.tensor([mr_src, mr_dst], dtype=torch.long)
                if mr_src else torch.zeros((2, 0), dtype=torch.long)
            )
            print(f"[Build] miRNA->gene edges: {len(mr_src)}, miRNA nodes: {len(mirna_list)}")

    print(f"[Build] Node types: {list(data.node_types)}")
    print(f"[Build] Edge types: {list(data.edge_types)}")
    return data, gene_to_idx, gene_list, pathway_name_to_idx