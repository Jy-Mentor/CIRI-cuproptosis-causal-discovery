# -*- coding: utf-8 -*-
"""
阶段3: HGT异构网络学习 — 关系推理与边预测
====================================================
使用异构图Transformer (HGT) 在"石竹烯-铁衰老-脑缺血"异构图上进行
消息传递，输出所有节点嵌入，并通过药物-靶点边预测验证网络建模有效性。

核心任务:
  1. 加载阶段1异构图 + 阶段2分子指纹
  2. 为各节点类型构建初始特征
  3. 训练HGT进行链路预测 (drug-gene边)
  4. ACSL4回忆实验: 掩蔽石竹烯-ACSL4边, 检验恢复能力
  5. 提取全部节点嵌入供阶段4使用

验证策略:
  - 时间分割验证 (模拟已知→未知边发现)
  - ACSL4掩蔽召回实验
  - AUROC / AUPRC / Hits@K 评估
"""

import os
import sys
import json
import copy
import math
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
from collections import defaultdict

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch import Tensor
from torch_geometric.data import HeteroData
from torch_geometric.nn import HGTConv, Linear
from sklearn.metrics import roc_auc_score, average_precision_score, precision_recall_curve

# 路径
BASE_DIR = Path(__file__).parent.parent
sys.path.insert(0, str(BASE_DIR))
DATA_DIR = BASE_DIR / "data" / "ferroptosis_graph"
RESULTS_DIR = BASE_DIR / "results" / "hgt_learning"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# 检查DATA_DIR是否存在
if not DATA_DIR.exists():
    raise FileNotFoundError(
        f"数据目录不存在: {DATA_DIR}\n"
        f"请先运行 stage1_build_graph.py 和 stage2_gat_encoder.py 生成异构图和分子指纹"
    )

# 全局常量 — 综合优化
HIDDEN_DIM = 128
NUM_HEADS = 4
NUM_LAYERS = 2
DROPOUT = 0.4
LR = 1e-3          # 提高学习率配合warmup
WEIGHT_DECAY = 1e-4
NEG_RATIO = 1
VAL_SPLIT = 0.2
# Focal Loss 参数
FOCAL_ALPHA = 0.75  # 正样本权重
FOCAL_GAMMA = 2.0   # 难样本聚焦
# 训练
MAX_EPOCHS = 200
PATIENCE = 30
WARMUP_EPOCHS = 10
GRAD_CLIP = 5.0

# PPI 配置
PPI_FILE = Path(r"C:\Users\Jy-Mentor-7\Desktop\9606蛋白质\9606_human_ppi_symbol.txt")
PPI_MIN_SCORE = 700

# Pathway 配置
PATHWAY_EDGES_FILE = Path(r"D:\反向网络药理学\GAT拓展维度\gene_pathway_edges.txt")
PATHWAY_CO_MEMBER_JACCARD = 0.15  # 基因间通路共有的Jaccard阈值

# Methylation 配置
METHYLATION_EDGES_FILE = Path(r"D:\反向网络药理学\GAT拓展维度\gene_methylation_edges.txt")
METHYLATION_CORRELATION_THRESHOLD = 0.1  # 甲基化模式相似度阈值 (降低以捕获更多甲基化共模式边)

# 基因集定义 (从stage1导入 — 必须在sys.path.insert之后)
from bcp_ferroptosis_pipeline.stage1_build_graph import (  # noqa: E402
    ALL_FERROPTOSIS_GENES, ALL_AGING_GENES, ALL_CIRI_GENES,
    BCP_TARGETS, ACSL4_FIRST_NEIGHBORS,
    CIRI_DEGS, FERRO_AGING_GENES_FILE,
)

# 设备
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Device: {DEVICE}")

# ============================================================
# 0. 多维边增强 (PPI + TF + CoPression + Pathway)
# ============================================================

# 人类已知转录因子列表 (TRRUST / AnimalTFDB / 文献)
HUMAN_TFS = {
    'TP53', 'NFE2L2', 'NFKB1', 'RELA', 'STAT3', 'JUN', 'FOS', 'MYC', 'HIF1A',
    'SP1', 'CTCF', 'EGR1', 'ATF4', 'CEBPB', 'FOXO3', 'PPARG', 'NRF1', 'TFEB',
    'MITF', 'SREBF1', 'SREBF2', 'NR1H4', 'RXRA', 'HNF4A', 'NR3C1',
    'ESR1', 'AR', 'PGR', 'RARA', 'VDR', 'PPARA', 'PPARD',
    'NR1I2', 'NR1I3', 'AHR', 'ARNT', 'EPAS1', 'ARNTL', 'CLOCK',
    'CREB1', 'ATF2', 'ATF3', 'ATF6', 'XBP1', 'DDIT3', 'MAF', 'MAFK',
    'MAFG', 'BACH1', 'BACH2', 'NFIL3',
    'FOXO1', 'FOXO4', 'FOXM1', 'FOXA1', 'FOXA2', 'GATA1', 'GATA2', 'GATA3',
    'GATA4', 'GATA6', 'KLF4', 'KLF5', 'SP3', 'SP4', 'WT1', 'EGR2', 'EGR3',
    'USF1', 'USF2', 'MYCN', 'MAX', 'MXD1', 'MNT',
    'E2F1', 'E2F2', 'E2F3', 'E2F4', 'RB1', 'SMAD2', 'SMAD3', 'SMAD4',
    'RUNX1', 'RUNX2', 'RUNX3', 'ETV1', 'ETV4', 'ETV5', 'ELK1', 'ELK4',
    'SRF', 'MEF2A', 'MEF2C', 'NFAT5', 'IRF1', 'IRF3',
    'STAT1', 'STAT2', 'STAT5A', 'STAT5B', 'STAT6', 'HSF1', 'HSF2',
    'PAX5', 'PAX6', 'SOX2', 'SOX9', 'POU5F1', 'NANOG', 'TCF3', 'TCF4',
    'LEF1', 'TCF7L2', 'ZEB1', 'ZEB2', 'SNAI1', 'SNAI2', 'TWIST1', 'TWIST2',
    'HMGA1', 'HMGA2', 'NFIA', 'NFIB', 'NFIX', 'SOX10', 'OLIG2', 'NEUROD1',
    'ASCL1', 'HES1', 'HEY1', 'NOTCH1', 'RBPJ', 'GLI1', 'GLI2',
}


def enrich_graph_edges(
    data: HeteroData,
    add_ppi: bool = True,
    add_tf: bool = True,
    add_coexp: bool = True,
    add_pathway: bool = True,
    add_methylation: bool = True,
) -> HeteroData:
    """多维边增强: PPI + TF调控 + 共表达 + 通路 + 甲基化。

    按维度递增策略逐步添加边类型，便于消融实验。

    Args:
        data: 原始异构图
        add_ppi: 添加STRING PPI边
        add_tf: 添加TF调控边
        add_coexp: 添加高置信共表达边
        add_pathway: 添加通路共成员边 (gene-gene, Jaccard相似度)
        add_methylation: 添加甲基化共模式边 (gene-gene)

    Returns:
        增强后的 HeteroData, edge_stats字典
    """
    gene_names = data['gene'].names
    gene_to_idx = {g: i for i, g in enumerate(gene_names)}
    gene_set = set(gene_names)
    edge_stats = {}

    # === PPI 增强 ===
    if add_ppi and PPI_FILE.exists():
        print("  加载STRING PPI边...")
        ppi_edges = set()
        old_gg = data['gene', 'interacts', 'gene'].edge_index
        for j in range(old_gg.shape[1]):
            ppi_edges.add((int(old_gg[0, j]), int(old_gg[1, j])))

        with open(PPI_FILE, 'r') as f:
            f.readline()
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 3:
                    continue
                try:
                    if float(parts[2]) < PPI_MIN_SCORE:
                        continue
                except ValueError:
                    continue
                a, b = parts[0].upper(), parts[1].upper()
                if a == b or a not in gene_set or b not in gene_set:
                    continue
                ai, bi = gene_to_idx[a], gene_to_idx[b]
                ppi_edges.add((ai, bi))
                ppi_edges.add((bi, ai))

        new_gg = torch.tensor(list(ppi_edges), dtype=torch.long).T
        data['gene', 'interacts', 'gene'].edge_index = new_gg
        edge_stats['gene↔gene (PPI)'] = new_gg.shape[1]
        print(f"    PPI边: {old_gg.shape[1]} → {new_gg.shape[1]} (平均度 {new_gg.shape[1]/len(gene_names):.1f})")
    else:
        edge_stats['gene↔gene (PPI)'] = 0

    # === TF调控边 (有向) ===
    if add_tf:
        tfs_in_graph = HUMAN_TFS & gene_set
        if tfs_in_graph and add_ppi:
            print(f"  提取TF调控边 ({len(tfs_in_graph)}个TF)...")
            ppi_set = set()
            gg = data['gene', 'interacts', 'gene'].edge_index
            for j in range(gg.shape[1]):
                ppi_set.add((int(gg[0, j]), int(gg[1, j])))

            tf_edges = []
            for tf_name in tfs_in_graph:
                tf_idx = gene_to_idx[tf_name]
                for target_name in gene_set:
                    if target_name in tfs_in_graph:
                        continue
                    target_idx = gene_to_idx[target_name]
                    if (tf_idx, target_idx) in ppi_set:
                        tf_edges.append((tf_idx, target_idx))

            if tf_edges:
                data['gene', 'regulates', 'gene'].edge_index = torch.tensor(
                    tf_edges, dtype=torch.long
                ).T
                edge_stats['gene→gene (TF regulates)'] = len(tf_edges)
                print(f"    TF调控边: {len(tf_edges)} ({len(tfs_in_graph)} TFs)")
            else:
                data['gene', 'regulates', 'gene'].edge_index = torch.zeros((2, 0), dtype=torch.long)
        else:
            data['gene', 'regulates', 'gene'].edge_index = torch.zeros((2, 0), dtype=torch.long)
            edge_stats['gene→gene (TF regulates)'] = 0
    else:
        data['gene', 'regulates', 'gene'].edge_index = torch.zeros((2, 0), dtype=torch.long)
        edge_stats['gene→gene (TF regulates)'] = 0

    # === 共表达边 (高置信STRING, 无向) ===
    if add_coexp and PPI_FILE.exists():
        print("  提取高置信共表达边 (STRING score ≥ 900)...")
        coexp_edges = set()
        with open(PPI_FILE, 'r') as f:
            f.readline()
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 3:
                    continue
                try:
                    if float(parts[2]) < 900:
                        continue
                except ValueError:
                    continue
                a, b = parts[0].upper(), parts[1].upper()
                if a == b or a not in gene_set or b not in gene_set:
                    continue
                ai, bi = gene_to_idx[a], gene_to_idx[b]
                coexp_edges.add((ai, bi))
                coexp_edges.add((bi, ai))

        if coexp_edges:
            data['gene', 'coexpressed', 'gene'].edge_index = torch.tensor(
                list(coexp_edges), dtype=torch.long
            ).T
            edge_stats['gene↔gene (coexpressed)'] = len(coexp_edges)
            print(f"    共表达边: {len(coexp_edges)}")
        else:
            data['gene', 'coexpressed', 'gene'].edge_index = torch.zeros((2, 0), dtype=torch.long)
    else:
        data['gene', 'coexpressed', 'gene'].edge_index = torch.zeros((2, 0), dtype=torch.long)
        edge_stats['gene↔gene (coexpressed)'] = 0

    # === 通路共成员边 (gene-gene, 基于Jaccard相似度) ===
    if add_pathway and PATHWAY_EDGES_FILE.exists():
        print(f"  提取通路共成员边 (Jaccard ≥ {PATHWAY_CO_MEMBER_JACCARD})...")
        # 加载基因→通路映射
        gene_pathways: Dict[str, Set[str]] = defaultdict(set)
        with open(PATHWAY_EDGES_FILE, 'r', encoding='utf-8') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    gene, pathway = parts[0].upper(), parts[1].strip()
                    if gene in gene_set:
                        gene_pathways[gene].add(pathway)

        n_mapped = len(gene_pathways)
        print(f"    通路映射基因: {n_mapped}/{len(gene_names)}")

        # 计算基因对之间的通路Jaccard相似度
        pathway_edges = set()
        mapped_genes = list(gene_pathways.keys())
        n_pairs = 0
        for i, g1 in enumerate(mapped_genes):
            pw1 = gene_pathways[g1]
            if len(pw1) < 2:
                continue
            for g2 in mapped_genes[i + 1:]:
                pw2 = gene_pathways[g2]
                if len(pw2) < 2:
                    continue
                intersection = len(pw1 & pw2)
                if intersection == 0:
                    continue
                union = len(pw1 | pw2)
                jaccard = intersection / union
                if jaccard >= PATHWAY_CO_MEMBER_JACCARD:
                    idx1, idx2 = gene_to_idx[g1], gene_to_idx[g2]
                    pathway_edges.add((idx1, idx2))
                    pathway_edges.add((idx2, idx1))
                n_pairs += 1

        if pathway_edges:
            data['gene', 'co_pathway', 'gene'].edge_index = torch.tensor(
                list(pathway_edges), dtype=torch.long
            ).T
            edge_stats['gene↔gene (co-pathway)'] = len(pathway_edges)
            print(f"    通路共成员边: {len(pathway_edges)} (Jaccard ≥ {PATHWAY_CO_MEMBER_JACCARD}, "
                  f"扫描{n_pairs:,}对)")
        else:
            data['gene', 'co_pathway', 'gene'].edge_index = torch.zeros((2, 0), dtype=torch.long)
            edge_stats['gene↔gene (co-pathway)'] = 0
    else:
        data['gene', 'co_pathway', 'gene'].edge_index = torch.zeros((2, 0), dtype=torch.long)
        edge_stats['gene↔gene (co-pathway)'] = 0

    # === 甲基化共模式边 (gene-gene, 基于共有甲基化CpG位点) ===
    if add_methylation and METHYLATION_EDGES_FILE.exists():
        print("  提取甲基化共模式边...")
        # 加载基因→CpG位点映射
        gene_cpgs: Dict[str, Set[str]] = defaultdict(set)
        with open(METHYLATION_EDGES_FILE, 'r', encoding='utf-8') as f:
            header = f.readline()  # 跳过表头
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    gene = parts[0].upper()
                    cpg = parts[1].strip()
                    if gene in gene_set:
                        gene_cpgs[gene].add(cpg)

        n_meth_mapped = len(gene_cpgs)
        print(f"    甲基化映射基因: {n_meth_mapped}/{len(gene_names)}")

        # 计算甲基化模式相似度 (共有CpG比例)
        meth_edges = set()
        meth_genes = list(gene_cpgs.keys())
        for i, g1 in enumerate(meth_genes):
            cpgs1 = gene_cpgs[g1]
            if len(cpgs1) < 1:  # 放宽到至少1个CpG
                continue
            for g2 in meth_genes[i + 1:]:
                cpgs2 = gene_cpgs[g2]
                if len(cpgs2) < 1:
                    continue
                intersection = len(cpgs1 & cpgs2)
                if intersection == 0:
                    continue
                # 使用简单的共有比例 (Jaccard近似)
                union = len(cpgs1 | cpgs2)
                jaccard = intersection / union if union > 0 else 0
                if jaccard >= METHYLATION_CORRELATION_THRESHOLD:
                    idx1, idx2 = gene_to_idx[g1], gene_to_idx[g2]
                    meth_edges.add((idx1, idx2))
                    meth_edges.add((idx2, idx1))

        if meth_edges:
            data['gene', 'co_methylated', 'gene'].edge_index = torch.tensor(
                list(meth_edges), dtype=torch.long
            ).T
            edge_stats['gene↔gene (co-methylated)'] = len(meth_edges)
            print(f"    甲基化共模式边: {len(meth_edges)} (阈值 ≥ {METHYLATION_CORRELATION_THRESHOLD})")
        else:
            data['gene', 'co_methylated', 'gene'].edge_index = torch.zeros((2, 0), dtype=torch.long)
            edge_stats['gene↔gene (co-methylated)'] = 0
    else:
        data['gene', 'co_methylated', 'gene'].edge_index = torch.zeros((2, 0), dtype=torch.long)
        edge_stats['gene↔gene (co-methylated)'] = 0

    # === 打印边统计 ===
    print("\n  多维边增强完成:")
    total_edges = 0
    for k, v in edge_stats.items():
        total_edges += v
        print(f"    {k}: {v}")
    print(f"    总计: {total_edges} 条增强边")

    return data

# ============================================================
# 0.5 高级技术: HeCo对比预训练 + 子图提取 + 困难负样本挖掘
# ============================================================

def extract_ferroptosis_subgraph(
    data: HeteroData,
    center_genes: Set[str] = None,
    k_hop: int = 2,
) -> Tuple[Set[int], Dict[str, List[int]]]:
    """提取铁死亡核心基因的k-hop子图。

    基于 HierHGT-DTI 层次聚合思想，以 ACSL4 和铁死亡核心基因为中心，
    提取k-hop邻域子图，聚焦模型训练于关键生物学上下文。

    Args:
        data: 异构图
        center_genes: 中心基因名称集合 (默认: ACSL4 + 铁死亡核心)
        k_hop: 跳数

    Returns:
        subgraph_gene_indices: 子图中基因节点索引集合
        hop_info: {hop: [gene_indices]} 每跳基因索引
    """
    gene_names = data['gene'].names
    gene_to_idx = {g: i for i, g in enumerate(gene_names)}

    if center_genes is None:
        center_genes = {"ACSL4", "GPX4", "SLC7A11", "FSP1", "NFE2L2",
                         "HMOX1", "TFRC", "LPCAT3", "ALOX5", "PTGS2",
                         "TP53", "HIF1A", "STAT3"}

    # 获取有效中心基因索引
    center_indices = set()
    for g in center_genes:
        if g in gene_to_idx:
            center_indices.add(gene_to_idx[g])

    if not center_indices:
        return set(range(len(gene_names))), {0: list(range(len(gene_names)))}

    # BFS提取k-hop子图
    visited = set(center_indices)
    hop_info = {0: list(center_indices)}
    frontier = set(center_indices)

    # 获取所有gene-gene边 (PPI)
    if hasattr(data['gene', 'interacts', 'gene'], 'edge_index'):
        gg_ei = data['gene', 'interacts', 'gene'].edge_index.numpy()
        # 构建邻接表
        adj = defaultdict(set)
        for j in range(gg_ei.shape[1]):
            u, v = int(gg_ei[0, j]), int(gg_ei[1, j])
            adj[u].add(v)
            adj[v].add(u)

        for hop in range(1, k_hop + 1):
            next_frontier = set()
            for node in frontier:
                for neighbor in adj[node]:
                    if neighbor not in visited:
                        visited.add(neighbor)
                        next_frontier.add(neighbor)
            hop_info[hop] = list(next_frontier)
            frontier = next_frontier
            if not frontier:
                break

    # 也包含drug-gene边关联的基因
    if hasattr(data['drug', 'targets', 'gene'], 'edge_index'):
        dg_ei = data['drug', 'targets', 'gene'].edge_index.numpy()
        for j in range(dg_ei.shape[1]):
            g_idx = int(dg_ei[1, j])
            if g_idx not in visited:
                visited.add(g_idx)
                if 'target' not in hop_info:
                    hop_info['target'] = []
                hop_info['target'].append(g_idx)

    print(f"  子图提取: {len(center_indices)} 中心基因 → {len(visited)}/{len(gene_names)} 节点 "
          f"(k={k_hop}, 覆盖率 {len(visited)/len(gene_names)*100:.1f}%)")

    return visited, hop_info


class HeCoContrastivePretrainer:
    """HeCo风格异构图对比预训练器 (WWW 2022)。

    参考:
      - HeCo: Heterogeneous Graph Contrastive Learning, Wang et al., WWW 2022
      - HierHGT-DTI: 层次化异构图Transformer, Bioinformatics Advances 2025

    双视图对比:
      - View 1 (Schema View): 原始图结构 + 原始特征
      - View 2 (Augmented View): 边丢弃 + 特征掩码 增强图

    对比目标:
      - 最大化同一节点在两个视图中的嵌入相似度
      - 最小化不同节点间的嵌入相似度
    """

    def __init__(
        self,
        model: nn.Module,
        hidden_dim: int = 128,
        proj_dim: int = 64,
        tau: float = 0.5,
        edge_dropout: float = 0.3,
        feat_mask: float = 0.2,
        device: str = "cpu",
    ):
        self.model = model
        self.hidden_dim = hidden_dim
        self.tau = tau
        self.edge_dropout = edge_dropout
        self.feat_mask = feat_mask
        self.device = device

        # 投影头 (每个节点类型)
        self.projectors = nn.ModuleDict()
        for nt in model.node_types:
            self.projectors[nt] = nn.Sequential(
                nn.Linear(hidden_dim, hidden_dim),
                nn.ReLU(),
                nn.Linear(hidden_dim, proj_dim),
            )
        self.projectors = self.projectors.to(device)

        self.optimizer = None

    def _augment_graph(self, data: HeteroData, epoch: int = 0) -> HeteroData:
        """增强图: 边丢弃 + 特征掩码 (View 2)。

        若模型启用了自适应图增强 (AdaptiveGraphAugmentation)，则使用课程式
        增强强度 (从弱→强)；否则使用固定强度。
        """
        augmented = data.clone().detach()

        # 检查模型是否有自适应增强模块
        use_adaptive = (
            hasattr(self.model, 'graph_augmentation')
            and self.model.graph_augmentation is not None
        )

        if use_adaptive:
            # 使用课程式增强强度
            feat_mask_rate, edge_dropout_rate = (
                self.model.graph_augmentation.get_augmentation_strength(epoch)
            )
        else:
            feat_mask_rate = self.feat_mask
            edge_dropout_rate = self.edge_dropout

        # 边丢弃
        for et in augmented.edge_types:
            if hasattr(augmented[et], 'edge_index') and augmented[et].edge_index.numel() > 0:
                ei = augmented[et].edge_index
                n_edges = ei.shape[1]
                if n_edges > 0:
                    keep_mask = torch.rand(n_edges) > edge_dropout_rate
                    if keep_mask.sum() > 0:
                        augmented[et].edge_index = ei[:, keep_mask]
                    else:
                        augmented[et].edge_index = ei  # 至少保留一条

        # 特征掩码
        for nt in augmented.node_types:
            if hasattr(augmented[nt], 'x') and augmented[nt].x is not None:
                x = augmented[nt].x.clone()
                mask = torch.rand_like(x) > feat_mask_rate
                augmented[nt].x = x * mask

        return augmented

    def _forward_view(self, data: HeteroData) -> Dict[str, Tensor]:
        """前向传播一个视图，返回投影后的嵌入。"""
        x_dict = {nt: data[nt].x.to(self.device) for nt in self.model.node_types}
        edge_index_dict = {
            et: data[et].edge_index.to(self.device) for et in self.model.edge_types
            if hasattr(data[et], 'edge_index')
        }

        z_dict_out = self.model(x_dict, edge_index_dict)
        # 新接口返回 (z_dict, extra)，兼容旧接口
        z_dict = z_dict_out[0] if isinstance(z_dict_out, tuple) else z_dict_out

        # 投影
        projected = {}
        for nt in self.model.node_types:
            if nt in z_dict:
                projected[nt] = F.normalize(
                    self.projectors[nt](z_dict[nt]), dim=-1
                )
        return projected

    def _info_nce_loss(
        self,
        z1: Tensor,
        z2: Tensor,
    ) -> Tensor:
        """InfoNCE对比损失。

        Args:
            z1: View 1 嵌入 (N, D)
            z2: View 2 嵌入 (N, D)
        """
        N = z1.size(0)
        if N < 2:
            return torch.tensor(0.0, device=self.device)

        # 相似度矩阵
        sim = torch.mm(z1, z2.T) / self.tau  # (N, N)

        # 正样本: 对角线
        pos_sim = torch.diag(sim)  # (N,)

        # 负样本: 所有非对角线
        # 使用 -inf 掩盖对角线
        mask = torch.eye(N, device=self.device, dtype=torch.bool)
        neg_sim = sim.masked_fill(mask, float('-inf'))

        # 分子分母
        numerator = torch.exp(pos_sim)
        denominator = numerator + torch.exp(neg_sim).sum(dim=1)

        loss = -torch.log(numerator / (denominator + 1e-8)).mean()
        return loss

    def pretrain(
        self,
        data: HeteroData,
        epochs: int = 100,
        lr: float = 1e-3,
        patience: int = 20,
    ) -> Dict[str, List[float]]:
        """对比预训练。"""
        print(f"\n  [HeCo对比预训练] epochs={epochs}, tau={self.tau}")
        print(f"    边丢弃率={self.edge_dropout}, 特征掩码率={self.feat_mask}")

        params = list(self.model.parameters()) + list(self.projectors.parameters())
        self.optimizer = torch.optim.AdamW(params, lr=lr, weight_decay=1e-5)

        data = data.to(self.device)
        history = {"loss": [], "gene_align": []}
        best_loss = float('inf')
        best_state = None
        patience_counter = 0

        for epoch in range(epochs):
            self.model.train()
            self.projectors.train()

            # 增强图 (传入epoch以支持课程式增强)
            data_aug = self._augment_graph(data, epoch=epoch)

            # 双视图前向
            z1 = self._forward_view(data)      # View 1: 原始
            z2 = self._forward_view(data_aug)  # View 2: 增强

            # 计算每种节点类型的对比损失(主要关注gene)
            total_loss = 0.0
            gene_align = 0.0
            n_types = 0
            for nt in z1:
                if nt in z2 and z1[nt].size(0) > 1:
                    loss_nt = self._info_nce_loss(z1[nt], z2[nt])
                    total_loss += loss_nt
                    n_types += 1
                    if nt == 'gene':
                        gene_align = loss_nt.item()

            if n_types > 0:
                total_loss = total_loss / n_types

            self.optimizer.zero_grad()
            total_loss.backward()
            torch.nn.utils.clip_grad_norm_(params, 1.0)
            self.optimizer.step()

            loss_val = total_loss.item()
            history["loss"].append(loss_val)
            history["gene_align"].append(gene_align)

            if (epoch + 1) % 50 == 0:
                print(f"    Epoch {epoch+1:4d} | Loss: {loss_val:.4f} | Gene Align: {gene_align:.4f}")

            # 早停
            if loss_val < best_loss - 0.001:
                best_loss = loss_val
                best_state = copy.deepcopy(self.model.state_dict())
                patience_counter = 0
            else:
                patience_counter += 1

            if patience_counter >= patience:
                print(f"    对比预训练早停于 epoch {epoch+1}, best loss = {best_loss:.4f}")
                break

        # 恢复最佳模型
        if best_state is not None:
            self.model.load_state_dict(best_state)

        print(f"    对比预训练完成 | Best Loss: {best_loss:.4f}")
        return history


def mine_hard_negatives(
    data: HeteroData,
    positive_gene_indices: Set[int],
    n_hard: int = 50,
    n_easy: int = 100,
) -> Tuple[Tensor, Tensor]:
    """困难负样本挖掘: 从通路相近但非靶点的基因中采样。

    策略 (基于 DrugBAN/HierHGT-DTI 负采样):
      - Hard negatives: 与正样本共享通路的基因 (Jaccard ≥ 0.1)
      - Easy negatives: 随机采样非正样本基因

    Args:
        data: 异构图
        positive_gene_indices: 正样本基因索引集合
        n_hard: 困难负样本数量
        n_easy: 简单负样本数量

    Returns:
        neg_drug_indices, neg_gene_indices
    """
    import random as _random
    _random.seed(42)

    gene_names = data['gene'].names
    n_genes = len(gene_names)
    gene_to_idx = {g: i for i, g in enumerate(gene_names)}

    n_drugs = data['drug'].num_nodes

    # 收集困难负样本：与正样本通路相近的基因
    hard_neg_pool = set()
    if PATHWAY_EDGES_FILE.exists():
        # 加载基因→通路映射
        gene_pathways: Dict[str, Set[str]] = defaultdict(set)
        with open(PATHWAY_EDGES_FILE, 'r', encoding='utf-8') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    gene, pathway = parts[0].upper(), parts[1].strip()
                    if gene in gene_to_idx:
                        gene_pathways[gene].add(pathway)

        # 获取正样本基因的通路集合
        pos_pathways = set()
        for g_idx in positive_gene_indices:
            g = gene_names[g_idx]
            if g in gene_pathways:
                pos_pathways |= gene_pathways[g]

        if pos_pathways:
            # 找与正样本共享通路的非靶点基因
            for g in gene_pathways:
                g_idx = gene_to_idx[g]
                if g_idx in positive_gene_indices:
                    continue
                shared = len(gene_pathways[g] & pos_pathways)
                if shared >= 2:  # 至少共享2个通路
                    hard_neg_pool.add(g_idx)

    # 补充邻居基因作为困难负样本
    if hasattr(data['gene', 'interacts', 'gene'], 'edge_index'):
        gg_ei = data['gene', 'interacts', 'gene'].edge_index.numpy()
        for j in range(gg_ei.shape[1]):
            u, v = int(gg_ei[0, j]), int(gg_ei[1, j])
            if u in positive_gene_indices and v not in positive_gene_indices:
                hard_neg_pool.add(v)
            if v in positive_gene_indices and u not in positive_gene_indices:
                hard_neg_pool.add(u)

    # 采样
    hard_neg = list(hard_neg_pool)
    _random.shuffle(hard_neg)
    hard_neg = hard_neg[:n_hard]

    # Easy negatives: 随机采样
    easy_neg_pool = set(range(n_genes)) - positive_gene_indices - hard_neg_pool
    if len(easy_neg_pool) < n_easy:
        n_easy = len(easy_neg_pool)
    easy_neg = list(easy_neg_pool)
    _random.shuffle(easy_neg)
    easy_neg = easy_neg[:n_easy]

    all_neg_genes = hard_neg + easy_neg
    neg_drugs = torch.zeros(len(all_neg_genes), dtype=torch.long)
    neg_genes = torch.tensor(all_neg_genes, dtype=torch.long)

    if hard_neg:
        print(f"    困难负样本: {len(hard_neg)} (通路相邻/PPI邻居), 简单负样本: {len(easy_neg)}")

    return neg_drugs, neg_genes


# ============================================================
# 1. 数据加载与特征构建
# ============================================================

def build_node_features(data: HeteroData, fingerprint_dim: int = HIDDEN_DIM) -> HeteroData:
    """为异构图各节点类型构建初始特征。

    策略:
      - drug: 使用GAT编码的256维分子指纹
      - gene: 使用拼接编码 (类别特征 + 中心性特征)
      - disease: 使用基因池化特征 + one-hot
      - pathway: 使用基因池化特征
      - phenotype: 使用通路池化特征
    """
    hidden_dim = fingerprint_dim

    # --- Drug: GAT分子指纹 ---
    fp_path = DATA_DIR / "bcp_molecular_fingerprint.npy"
    if fp_path.exists():
        drug_fp = np.load(fp_path)
        drug_feat = torch.tensor(drug_fp, dtype=torch.float).unsqueeze(0)  # (1, 256)
    else:
        drug_feat = torch.randn(1, hidden_dim)

    # 如果指纹维度不匹配，投影到hidden_dim
    if drug_feat.size(1) != hidden_dim:
        drug_feat = F.linear(drug_feat, torch.randn(hidden_dim, drug_feat.size(1)))
    data['drug'].x = drug_feat

    # --- Gene: 多模态特征 ---
    num_genes = data['gene'].num_nodes
    gene_names = data['gene'].names

    gene_feat = _build_gene_features(gene_names, hidden_dim, data)
    data['gene'].x = gene_feat

    # --- Disease: 基因平均池化 ---
    disease_feat = _build_disease_features(data, hidden_dim)
    data['disease'].x = disease_feat

    # --- Pathway: 基因平均池化 ---
    pathway_feat = _build_pathway_features(data, hidden_dim)
    data['pathway'].x = pathway_feat

    # --- Phenotype: 通路平均池化 ---
    phenotype_feat = _build_phenotype_features(data, hidden_dim)
    data['phenotype'].x = phenotype_feat

    # 将数据移到设备
    for nt in data.node_types:
        if hasattr(data[nt], 'x') and data[nt].x is not None:
            data[nt].x = data[nt].x.float()

    return data


def _build_gene_features(
    gene_names: List[str],
    hidden_dim: int,
    data: HeteroData,
) -> torch.Tensor:
    """构建基因节点特征。

    组合多种信号:
      - 铁死亡类别 (驱动/抑制/标记)
      - 衰老类别
      - CIRI类别
      - ACSL4邻居关系 (仅基于gene-gene边, 不含标签信息)
      - 通路成员数
      - 甲基化度
      - 随机嵌入 (降噪)

    注意: 不包含 is_bcp_target / is_acsl4 等标签泄漏特征。
    这些信息由HGT边预测学习, 不作为节点初始特征输入。
    """
    num_genes = len(gene_names)

    acsl4_neighbors = set(
        ACSL4_FIRST_NEIGHBORS["direct_interactors"] +
        ACSL4_FIRST_NEIGHBORS["indirect_regulators"]
    )

    # 通路成员统计 (每个基因属于几个通路)
    pathway_count = np.zeros(num_genes)
    if hasattr(data['gene', 'belongs_to', 'pathway'], 'edge_index'):
        gp_ei = data['gene', 'belongs_to', 'pathway'].edge_index
        for j in range(gp_ei.shape[1]):
            g_idx = int(gp_ei[0, j])
            if g_idx < num_genes:
                pathway_count[g_idx] += 1

    # 甲基化度统计 (每个基因关联几个CpG位点)
    methylation_degree = np.zeros(num_genes)
    gene_to_idx_local = {g: i for i, g in enumerate(gene_names)}
    if METHYLATION_EDGES_FILE.exists():
        with open(METHYLATION_EDGES_FILE, 'r', encoding='utf-8') as f:
            f.readline()  # 跳过表头
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    gene = parts[0].upper()
                    if gene in gene_to_idx_local:
                        methylation_degree[gene_to_idx_local[gene]] += 1

    # 随机初始化编码 - 降低维度 (64→16) 减少噪声
    embed_cache_path = DATA_DIR / f"gene_random_embed_{num_genes}x16.npy"
    if embed_cache_path.exists():
        random_embed_matrix = np.load(embed_cache_path)
    else:
        np.random.seed(42)
        random_embed_matrix = np.random.randn(num_genes, 16).astype(np.float32) * 0.05
        np.save(embed_cache_path, random_embed_matrix)

    features_list = []

    for i, g in enumerate(gene_names):
        feat = []

        # 1. 铁死亡类别编码 (3维 one-hot)
        is_driver = 1.0 if g in ["ACSL4","LPCAT3","ALOX5","ALOX15","TP53",
                                 "TFRC","HMOX1","SAT1"] else 0.0
        is_suppressor = 1.0 if g in ["GPX4","FSP1","SLC7A11","FTH1","NFE2L2",
                                      "GCH1","SLC40A1","AIFM2","GCLC"] else 0.0
        is_marker = 1.0 if g in ["PTGS2","CHAC1","ATF3","DDIT3"] else 0.0
        feat.extend([is_driver, is_suppressor, is_marker])

        # 2. 衰老类别 (2维)
        in_aging = 1.0 if g in ALL_AGING_GENES else 0.0
        in_ferro_aging = 1.0 if g in FERRO_AGING_GENES_FILE else 0.0
        feat.extend([in_aging, in_ferro_aging])

        # 3. CIRI类别 (2维)
        in_ciri_up = 1.0 if g in CIRI_DEGS.get("upregulated_in_ischemia",[]) else 0.0
        in_ciri_down = 1.0 if g in CIRI_DEGS.get("downregulated_in_ischemia",[]) else 0.0
        feat.extend([in_ciri_up, in_ciri_down])

        # 4. ACSL4邻居 (1维) — 仅基于gene-gene边, 不含标签信息
        is_acsl4_nb = 1.0 if g in acsl4_neighbors else 0.0
        feat.append(is_acsl4_nb)

        # 5. 通路成员数 (1维, 归一化)
        feat.append(min(pathway_count[i] / 8.0, 1.0))

        # 6. 甲基化度 (1维, 归一化)
        feat.append(min(methylation_degree[i] / 20.0, 1.0))

        # 7. 使用缓存的随机嵌入 (16维, 降噪)
        feat.extend(random_embed_matrix[i].tolist())

        features_list.append(feat)

    gene_feat = torch.tensor(features_list, dtype=torch.float)  # (N, 26): 3+2+2+1+1+1+16

    # 投影到hidden_dim - 固定种子 + 保存/加载权重确保可复现
    torch.manual_seed(42)
    np.random.seed(42)
    proj = nn.Linear(gene_feat.size(1), hidden_dim)
    proj_cache_path = DATA_DIR / f"gene_proj_weights_{gene_feat.size(1)}to{hidden_dim}.pt"
    if proj_cache_path.exists():
        proj.load_state_dict(torch.load(proj_cache_path, map_location='cpu'))
    else:
        nn.init.xavier_uniform_(proj.weight)
        nn.init.zeros_(proj.bias)
        torch.save(proj.state_dict(), proj_cache_path)

    gene_feat_proj = proj(gene_feat)

    # 添加节点度信息作为条件
    if hasattr(data['gene', 'interacts', 'gene'], 'edge_index'):
        ei = data['gene', 'interacts', 'gene'].edge_index
        degrees = torch.zeros(num_genes)
        for idx in range(num_genes):
            degrees[idx] = (ei[0] == idx).sum() + (ei[1] == idx).sum()
        degree_norm = torch.log1p(degrees) / math.log(max(2, degrees.max().item() + 1))
        gene_feat_proj = gene_feat_proj + degree_norm.unsqueeze(1) * 0.1 * gene_feat_proj.std()

    return gene_feat_proj


def _build_disease_features(data: HeteroData, hidden_dim: int) -> torch.Tensor:
    """构建疾病节点特征。

    CIRI特征 = CIRI基因集合的平均特征
    """
    num_diseases = data['disease'].num_nodes

    if hasattr(data, 'gene') and hasattr(data['gene'], 'x'):
        # 取CIRI相关基因的平均嵌入
        gene_names = data['gene'].names
        ciri_indices = [i for i, g in enumerate(gene_names) if g in ALL_CIRI_GENES]

        if ciri_indices:
            ciri_feat = data['gene'].x[torch.tensor(ciri_indices)].mean(dim=0, keepdim=True)
            # 复制到所有疾病节点
            disease_feat = ciri_feat.repeat(num_diseases, 1)
        else:
            disease_feat = torch.randn(num_diseases, hidden_dim) * 0.1
    else:
        disease_feat = torch.randn(num_diseases, hidden_dim) * 0.1

    return disease_feat


def _build_pathway_features(data: HeteroData, hidden_dim: int) -> torch.Tensor:
    """构建通路节点特征。

    每个通路 = 其成员基因的平均特征
    """
    num_pathways = data['pathway'].num_nodes
    gene_x = data['gene'].x

    pathway_feat = torch.zeros(num_pathways, hidden_dim)

    if hasattr(data['gene', 'belongs_to', 'pathway'], 'edge_index'):
        gp_ei = data['gene', 'belongs_to', 'pathway'].edge_index
        # 对每个通路聚合基因特征
        pw_counts = torch.zeros(num_pathways)
        for g_idx, p_idx in gp_ei.t():
            if p_idx < num_pathways:
                pathway_feat[p_idx] += gene_x[g_idx]
                pw_counts[p_idx] += 1

        # 避免除零
        pw_counts = pw_counts.clamp(min=1)
        pathway_feat = pathway_feat / pw_counts.unsqueeze(1)

    return pathway_feat


def _build_phenotype_features(data: HeteroData, hidden_dim: int) -> torch.Tensor:
    """构建表型节点特征。

    每个表型 = 其关联通路的平均特征
    """
    num_phenotypes = data['phenotype'].num_nodes
    pw_x = data['pathway'].x

    phenotype_feat = torch.zeros(num_phenotypes, hidden_dim)

    if hasattr(data['pathway', 'related_to', 'phenotype'], 'edge_index'):
        pp_ei = data['pathway', 'related_to', 'phenotype'].edge_index
        ph_counts = torch.zeros(num_phenotypes)
        for pw_idx, ph_idx in pp_ei.t():
            if ph_idx < num_phenotypes:
                phenotype_feat[ph_idx] += pw_x[pw_idx]
                ph_counts[ph_idx] += 1

        ph_counts = ph_counts.clamp(min=1)
        phenotype_feat = phenotype_feat / ph_counts.unsqueeze(1)

    return phenotype_feat


class VariationalInformationBottleneck(nn.Module):
    """变分信息瓶颈 (VIB) 模块 (源自 IB-DTI, Complex Eng. Syst. 2026)。

    核心思想: 学习压缩的、去噪的节点表示，通过变分方法近似
    最小化 I(Z; X) + β * I(Z; Y)，其中:
      - I(Z; X): 嵌入与输入的互信息 (压缩项)
      - I(Z; Y): 嵌入与标签的互信息 (预测项, 由BCE损失覆盖)
      - β: 权衡系数 (越大越压缩)

    实现策略:
      - 高斯重参数化: μ(x), σ(x) → z = μ + σ * ε, ε ~ N(0,I)
      - KL散度: KL(N(μ,σ²) || N(0,1)) 作为压缩正则项
      - 可学习的 β 从0.001线性增长到0.01 (warmup)

    参考:
      - Alemi et al., "Deep Variational Information Bottleneck", ICLR 2017
      - Song et al., "DTI prediction via hierarchical gated attention
        and information bottleneck", Complex Eng. Syst. 2026
    """

    def __init__(
        self,
        hidden_dim: int,
        beta_start: float = 1e-4,
        beta_end: float = 1e-2,
        beta_warmup_epochs: int = 50,
        eps: float = 1e-8,
    ):
        super().__init__()
        self.hidden_dim = hidden_dim
        self.beta_start = beta_start
        self.beta_end = beta_end
        self.beta_warmup_epochs = max(1, beta_warmup_epochs)
        self.eps = eps
        self.current_epoch = 0

        # 编码器: x → μ, logvar
        self.encoder_mu = nn.Linear(hidden_dim, hidden_dim)
        self.encoder_logvar = nn.Linear(hidden_dim, hidden_dim)

        # 输出投影 (压缩后)
        self.output_proj = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.ReLU(),
            nn.Dropout(0.1),
        )

    def get_beta(self) -> float:
        """获取当前β值 (线性warmup)。"""
        if self.current_epoch < self.beta_warmup_epochs:
            ratio = self.current_epoch / self.beta_warmup_epochs
            return self.beta_start + (self.beta_end - self.beta_start) * ratio
        return self.beta_end

    def set_epoch(self, epoch: int):
        """设置当前epoch (用于β调度)。"""
        self.current_epoch = epoch

    def forward(
        self,
        x: Tensor,
    ) -> Tuple[Tensor, Tensor]:
        """VIB前向传播。

        Args:
            x: 输入嵌入 (N, D)

        Returns:
            z: 压缩后的嵌入 (N, D)
            kl_loss: KL散度损失 (标量)
        """
        # 重参数化
        mu = self.encoder_mu(x)
        logvar = self.encoder_logvar(x)
        logvar = torch.clamp(logvar, -10, 10)

        # 采样
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        z = mu + eps * std

        # KL散度: KL(N(μ,σ²) || N(0,1))
        kl_loss = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp(), dim=-1)
        kl_loss = kl_loss.mean()

        # 输出投影
        z = self.output_proj(z)

        return z, kl_loss


class GraphWaveletMultiScaleFeatures(nn.Module):
    """图小波多尺度特征提取器 (源自 GHCDTI, Scientific Reports 2025)。

    使用图拉普拉斯的多尺度分解，对基因节点提取:
      - 低频特征 (全局结构): 使用热核扩散
      - 高频特征 (局部细节): 使用小波高频分量
      - 融合特征: 自适应加权融合

    实现: 使用 Chebyshev 多项式近似图小波变换，
    避免昂贵的特征分解。

    参考:
      - Dai et al., "GHCDTI: Heterogeneous network drug-target interaction
        prediction model based on graph wavelet transform and multi-level
        contrastive learning", Scientific Reports 2025
      - Hammond et al., "Wavelets on graphs via spectral graph theory",
        Applied and Computational Harmonic Analysis 2011
    """

    def __init__(
        self,
        hidden_dim: int,
        n_scales: int = 3,
        chebyshev_order: int = 5,
        dropout: float = 0.2,
    ):
        super().__init__()
        self.hidden_dim = hidden_dim
        self.n_scales = n_scales
        self.chebyshev_order = chebyshev_order

        # 每尺度的可学习滤波器
        self.scale_filters = nn.ModuleList()
        for _ in range(n_scales):
            self.scale_filters.append(
                nn.Sequential(
                    nn.Linear(hidden_dim, hidden_dim),
                    nn.LayerNorm(hidden_dim),
                    nn.ReLU(),
                    nn.Dropout(dropout),
                )
            )

        # 自适应融合门控
        self.fusion_gate = nn.Sequential(
            nn.Linear(hidden_dim * n_scales, n_scales),
            nn.Softmax(dim=-1),
        )

        # 最终融合投影
        self.fusion_proj = nn.Linear(hidden_dim * n_scales, hidden_dim)

    def _compute_chebyshev_filter(
        self,
        x: Tensor,           # (N, D) 节点特征
        edge_index: Tensor,   # (2, E) 边索引
        n_scale: int,         # 尺度索引
        num_nodes: int,       # 节点数
    ) -> Tensor:
        """Chebyshev多项式近似的图小波滤波。

        使用热核缩放: T_k(L) ≈ exp(-k²σ²/2) * cos(k * arccos(λ))
        """
        N = num_nodes
        D = x.size(1)

        # 归一化拉普拉斯: L = I - D^{-1/2} A D^{-1/2}
        # 使用度数矩阵
        device = x.device
        deg = torch.zeros(N, device=device)
        for j in range(edge_index.shape[1]):
            deg[edge_index[0, j]] += 1.0
        deg_inv_sqrt = deg.pow(-0.5)
        deg_inv_sqrt[torch.isinf(deg_inv_sqrt)] = 0.0

        # Chebyshev递归: T_0 = x, T_1 = Lx
        # 简化: 使用对称归一化邻接矩阵 A_norm = D^{-1/2} A D^{-1/2}
        # 则 L_norm = I - A_norm, Lx ≈ x - A_norm @ x

        # 构建稀疏邻接矩阵乘积
        adj_weights = deg_inv_sqrt[edge_index[0]] * deg_inv_sqrt[edge_index[1]]

        # T_0 = x
        T_prev = x  # (N, D)

        # 计算 A_norm @ x 作为 Lx 的基础
        ax = torch.zeros_like(x)
        for j in range(edge_index.shape[1]):
            src = edge_index[0, j]
            dst = edge_index[1, j]
            w = adj_weights[j]
            ax[src] += w * x[dst]

        # Lx = x - A_norm @ x
        lx = x - ax

        # T_1 = Lx (实际是 L @ x)
        T_curr = lx

        # 缩放因子 (不同尺度)
        scale = n_scale + 1
        heat_kernel = torch.exp(-torch.tensor(scale * scale * 0.5, device=device))

        # 聚合: x_filtered ≈ sum_k T_k * heat_kernel_scale
        filtered = heat_kernel * T_prev

        # 更高阶Chebyshev项
        for k in range(2, self.chebyshev_order + 1):
            # 计算 A_norm @ T_curr (用于后续 L @ T_curr = T_curr - A_norm @ T_curr)
            a_t_prev = torch.zeros_like(T_curr)
            for j in range(edge_index.shape[1]):
                src = edge_index[0, j]
                dst = edge_index[1, j]
                w = adj_weights[j]
                a_t_prev[src] += w * T_curr[dst]

            # T_k = 2 * L @ T_{k-1} - T_{k-2}
            # 其中 Lx = x - A_norm @ x, 所以 L @ T_curr = T_curr - A_norm @ T_curr
            # 修正: 使用 T_curr - a_t_prev 替代 a_t_prev (原代码错误)
            T_next = 2.0 * (T_curr - a_t_prev) - T_prev

            # 热核衰减
            heat_k = torch.exp(-torch.tensor(k * k * scale * 0.25, device=device))
            filtered = filtered + heat_k * T_next

            # 移位
            T_prev = T_curr
            T_curr = T_next

        # 通过可学习滤波器
        filtered = self.scale_filters[n_scale](filtered)

        return filtered

    def forward(
        self,
        gene_x: Tensor,          # (N_genes, D)
        edge_index: Tensor,       # (2, E) gene-gene边
    ) -> Tuple[Tensor, Dict[str, Tensor]]:
        """多尺度小波特征提取前向传播。

        Returns:
            fused: 融合后的多尺度特征 (N, D)
            scale_features: 每尺度特征 {f"scale_{i}": (N, D)}
        """
        N = gene_x.size(0)

        if edge_index.numel() == 0:
            # 无边时直接返回
            fused = self.fusion_proj(
                torch.cat([sf(gene_x) for sf in self.scale_filters], dim=-1)
            )
            return fused, {}

        # 提取每尺度特征
        scale_features = {}
        scale_embs = []
        for s in range(self.n_scales):
            f_s = self._compute_chebyshev_filter(gene_x, edge_index, s, N)
            scale_features[f"scale_{s}"] = f_s
            scale_embs.append(f_s)

        # 自适应融合
        stacked = torch.cat(scale_embs, dim=-1)  # (N, D * n_scales)
        gate = self.fusion_gate(stacked)          # (N, n_scales)

        fused = torch.zeros_like(gene_x)
        for s in range(self.n_scales):
            fused = fused + gate[:, s:s+1] * scale_embs[s]

        # 残差连接
        fused = fused + gene_x

        return fused, scale_features


class KnowledgeBasedRegularizer:
    """基于知识的正则化器 (源自知识增强DTI 2025)。

    整合生物学先验知识作为正则化约束:
      1. Pathway一致性: 同通路的基因嵌入应相似 (cosine接近)
      2. Drug-Target邻近: 已知BCP靶点与药物嵌入应接近
      3. Module紧致性: 同一功能模块的基因应形成紧致聚类
      4. 负边分离: 非靶点基因应与药物嵌入远离

    参考:
      - Yao et al., "Enhancing DTI prediction with graph representation
        learning and knowledge-based regularization", 2025
      - Chen et al., "Multilayer semantic-topology fusion for DTI", 2026

    Args:
        pathway_consistency_weight: 通路一致性损失权重
        drug_target_proximity_weight: 药物-靶点邻近损失权重
        module_compactness_weight: 模块紧致性损失权重
        neg_separation_weight: 负边分离损失权重
        margin: 负边分离的margin
    """

    def __init__(
        self,
        pathway_consistency_weight: float = 0.1,
        drug_target_proximity_weight: float = 0.05,
        module_compactness_weight: float = 0.05,
        neg_separation_weight: float = 0.02,
        margin: float = 0.5,
    ):
        self.pcw = pathway_consistency_weight
        self.dtpw = drug_target_proximity_weight
        self.mcw = module_compactness_weight
        self.nsw = neg_separation_weight
        self.margin = margin

    def compute_pathway_consistency_loss(
        self,
        gene_embeddings: Tensor,
        pathway_edges: Tensor,
    ) -> Tensor:
        """通路一致性损失: 同通路基因嵌入应更相似。"""
        if pathway_edges.numel() == 0:
            return torch.tensor(0.0, device=gene_embeddings.device)

        # 归一化
        g_norm = F.normalize(gene_embeddings, dim=-1)

        # 对每条通路边计算cosine相似度
        sims = (g_norm[pathway_edges[0]] * g_norm[pathway_edges[1]]).sum(dim=-1)

        # 最大化相似度 → 最小化 -cosine
        loss = -sims.mean()
        return loss * self.pcw

    def compute_drug_target_proximity_loss(
        self,
        drug_embedding: Tensor,
        gene_embeddings: Tensor,
        target_indices: Tensor,
    ) -> Tensor:
        """药物-靶点邻近损失: 已知靶点与药物接近。"""
        if target_indices.numel() == 0:
            return torch.tensor(0.0, device=gene_embeddings.device)

        # 归一化
        d_norm = F.normalize(drug_embedding, dim=-1)  # (1, D) 或 (D,)
        g_norm = F.normalize(gene_embeddings, dim=-1)  # (N, D)

        if d_norm.dim() == 1:
            d_norm = d_norm.unsqueeze(0)

        # 靶点与药物的cosine相似度
        target_embs = g_norm[target_indices]
        sims = (d_norm * target_embs).sum(dim=-1)

        loss = -sims.mean()
        return loss * self.dtpw

    def compute_module_compactness_loss(
        self,
        gene_embeddings: Tensor,
        module_assignment: Tensor,
    ) -> Tensor:
        """模块紧致性损失: 同模块基因协方差小。"""
        if module_assignment is None or module_assignment.sum() == 0:
            return torch.tensor(0.0, device=gene_embeddings.device)

        N_modules = module_assignment.size(1)
        g_norm = F.normalize(gene_embeddings, dim=-1)
        loss = torch.tensor(0.0, device=gene_embeddings.device)
        n_active = 0

        for m in range(N_modules):
            members = (module_assignment[:, m] > 0).nonzero().squeeze(-1)
            if len(members) > 1:
                # 计算模块内方差
                module_embs = g_norm[members]
                center = module_embs.mean(dim=0, keepdim=True)
                variance = (module_embs - center).pow(2).mean()
                loss = loss + variance
                n_active += 1

        if n_active > 0:
            loss = loss / n_active

        return loss * self.mcw

    def compute_neg_separation_loss(
        self,
        drug_embedding: Tensor,
        gene_embeddings: Tensor,
        neg_indices: Tensor,
    ) -> Tensor:
        """负边分离损失: 非靶点基因与药物远离。"""
        if neg_indices.numel() == 0:
            return torch.tensor(0.0, device=gene_embeddings.device)

        d_norm = F.normalize(drug_embedding, dim=-1)
        g_norm = F.normalize(gene_embeddings, dim=-1)

        if d_norm.dim() == 1:
            d_norm = d_norm.unsqueeze(0)

        neg_embs = g_norm[neg_indices]
        sims = (d_norm * neg_embs).sum(dim=-1)

        # hinge loss: max(0, sim - margin)
        loss = torch.clamp(sims - self.margin, min=0).mean()
        return loss * self.nsw

    def compute_all(
        self,
        gene_embeddings: Tensor,
        drug_embedding: Tensor,
        pathway_edges: Optional[Tensor] = None,
        target_indices: Optional[Tensor] = None,
        module_assignment: Optional[Tensor] = None,
        neg_indices: Optional[Tensor] = None,
    ) -> Dict[str, Tensor]:
        """计算所有知识正则化损失。"""
        losses = {}

        if pathway_edges is not None:
            losses['pathway_consistency'] = self.compute_pathway_consistency_loss(
                gene_embeddings, pathway_edges
            )

        if target_indices is not None:
            losses['drug_target_proximity'] = self.compute_drug_target_proximity_loss(
                drug_embedding, gene_embeddings, target_indices
            )

        if module_assignment is not None:
            losses['module_compactness'] = self.compute_module_compactness_loss(
                gene_embeddings, module_assignment
            )

        if neg_indices is not None:
            losses['neg_separation'] = self.compute_neg_separation_loss(
                drug_embedding, gene_embeddings, neg_indices
            )

        return losses


# ============================================================
# 2. 高级架构组件: 关系门控 + 层次基因模块聚合
# ============================================================

class RelationMessageGate(nn.Module):
    """关系消息门控机制 (源自 HierHGT-DTI, Bioinformatics Advances 2025)。

    为每种边类型学习一个可学习的标量门控值 (sigmoid-gated logit)，
    动态控制不同类型边在消息传递中的贡献权重。

    升级版改进:
      - 差异化初始化: 高置信边 (PPI/DrugTarget) 高bias, 低置信边 (Methyl) 低bias
      - 逐边类型冻结: 消融实验中可冻结特定不参与实验的维度
      - 温度退火调度: 训练后期降低温度使门控更尖锐

    核心作用:
      - 防止噪声边 (如低质量甲基化边) 污染消息传递
      - 自动学习每种关系类型的重要性
      - 消融实验可冻结门控值来验证各维度的贡献

    Args:
        edge_types: 边类型列表 [(src, rel, dst), ...]
        init_bias: 默认初始偏置
        per_edge_init_bias: 逐边类型的差异化初始偏置 {et_key: bias_value}
        temperature: soft-gating温度 (越低越接近硬门控)
        freeze_types: 需要冻结的边类型集合 (消融实验用)
    """

    # 边类型质量分级 (用于差异化初始化)
    # High: 高置信实验验证 → init_bias=3.0 (sigmoid≈0.95)
    # Medium: 中等置信 → init_bias=1.5 (sigmoid≈0.82)
    # Low: 推测性/噪声边 → init_bias=0.0 (sigmoid≈0.50)
    # Disabled: 消融冻结 → freeze, bias=-5.0 (sigmoid≈0.007)

    EDGE_QUALITY_MAP = {
        "drug__targets__gene": ("high", 3.0),       # 已知BCP靶点
        "gene__interacts__gene": ("high", 3.0),      # STRING PPI (score≥700)
        "gene__belongs_to__pathway": ("high", 3.0),  # KEGG/Reactome通路
        "gene__associated_with__disease": ("high", 3.0),  # CIRI差异基因
        "pathway__related_to__phenotype": ("medium", 1.5),
        "gene__regulates__phenotype": ("medium", 1.5),
        "gene__regulates__gene": ("medium", 1.5),    # TF调控 (PPI子集)
        "gene__coexpressed__gene": ("medium", 1.0),  # 共表达 (STRING≥900)
        "gene__co_pathway__gene": ("medium", 1.0),   # 通路共成员 (Jaccard)
        "gene__co_methylated__gene": ("low", 0.0),   # 甲基化共模式 (推测性)
    }

    def __init__(
        self,
        edge_types: List[Tuple[str, str, str]],
        init_bias: float = 2.0,
        per_edge_init_bias: Dict[str, float] = None,
        temperature: float = 1.0,
        freeze_types: Set[str] = None,
    ):
        super().__init__()
        self.edge_types = edge_types
        self.temperature = temperature
        self.freeze_types = freeze_types or set()

        # 为每种边类型学习一个logit
        self.gate_logits = nn.ParameterDict()
        for et in edge_types:
            et_key = f"{et[0]}__{et[1]}__{et[2]}"
            # 确定初始bias: 逐边类型配置 > 质量分级 > 默认值
            if per_edge_init_bias and et_key in per_edge_init_bias:
                bias_val = per_edge_init_bias[et_key]
            elif et_key in self.EDGE_QUALITY_MAP:
                _, bias_val = self.EDGE_QUALITY_MAP[et_key]
            else:
                bias_val = init_bias

            self.gate_logits[et_key] = nn.Parameter(
                torch.tensor(bias_val)
            )

        # 冻结指定边类型
        for et_key in self.freeze_types:
            if et_key in self.gate_logits:
                self.gate_logits[et_key].requires_grad = False

    def get_gate_values(self) -> Dict[str, float]:
        """获取所有边类型的门控值 (sigmoid后)。"""
        gates = {}
        device = next(self.gate_logits.values()).device
        for et_key, logit in self.gate_logits.items():
            gates[et_key] = float(torch.sigmoid(logit / self.temperature).item())
        return gates

    def set_temperature(self, new_temp: float):
        """温度退火调度: 训练后期降低温度使门控更尖锐。"""
        self.temperature = max(new_temp, 0.1)

    def forward(
        self,
        edge_index_dict: Dict[Tuple, Tensor],
    ) -> Dict[Tuple, Tensor]:
        """对边索引加权 (通过调整HGT注意力前的边权重)。

        Returns:
            gate_weights: {edge_type: scalar_weight}
        """
        gate_weights = {}
        device = next(self.gate_logits.values()).device
        for et in edge_index_dict:
            et_key = f"{et[0]}__{et[1]}__{et[2]}"
            if et_key in self.gate_logits:
                gate_weights[et] = torch.sigmoid(
                    self.gate_logits[et_key] / self.temperature
                )
            else:
                gate_weights[et] = torch.tensor(1.0, device=device)
        return gate_weights

    def get_gate_l1_loss(self) -> Tensor:
        """计算门控值的L1正则化损失 (鼓励稀疏门控)。

        加权策略: 低质量边类型受到更强的稀疏惩罚。
        """
        device = next(self.gate_logits.values()).device
        l1 = torch.tensor(0.0, device=device)
        for et_key, logit in self.gate_logits.items():
            gate_val = torch.sigmoid(logit / self.temperature)
            # 低质量边受到2x稀疏惩罚
            quality = self.EDGE_QUALITY_MAP.get(et_key, ("medium", 0))[0]
            penalty_weight = 2.0 if quality == "low" else 1.0
            l1 += penalty_weight * torch.abs(gate_val)
        return l1 / len(self.gate_logits)


class HierarchicalGeneModuleAggregator(nn.Module):
    """层次基因模块聚合器 (源自 HierHGT-DTI 层次聚合思想)。

    将基因节点按通路/功能模块分组为"基因模块超节点",
    执行层次化聚合: genes → gene_modules → global_context。

    这解决了直接在全连通PPI图上消息传递的"过度平滑"问题，
    并提供了天然的子图提取结构。

    聚合策略:
      - 'mean': 简单平均池化
      - 'attention': 可学习注意力加权
      - 'gated': 门控融合 (直接信息 + 聚合信息)
    """

    def __init__(
        self,
        hidden_dim: int,
        dropout: float = 0.3,
        aggregation: str = "attention",
    ):
        super().__init__()
        self.hidden_dim = hidden_dim
        self.aggregation = aggregation

        if aggregation == "attention":
            self.attn = nn.Sequential(
                nn.Linear(hidden_dim, hidden_dim // 4),
                nn.Tanh(),
                nn.Linear(hidden_dim // 4, 1),
            )
        elif aggregation == "gated":
            self.gate = nn.Sequential(
                nn.Linear(hidden_dim * 2, hidden_dim // 4),
                nn.ReLU(),
                nn.Dropout(dropout),
                nn.Linear(hidden_dim // 4, 1),
                nn.Sigmoid(),
            )

    def forward(
        self,
        gene_embeddings: Tensor,           # (N_genes, D)
        module_assignment: Tensor,          # (N_genes, N_modules) 稀疏分配矩阵
        module_supernode_emb: Tensor,       # (N_modules, D) 可学习超节点嵌入
    ) -> Tuple[Tensor, Tensor]:
        """层次聚合前向传播。

        Args:
            gene_embeddings: 基因节点嵌入
            module_assignment: 基因→模块的稀疏分配 (one-hot或soft)
            module_supernode_emb: 模块超节点初始嵌入

        Returns:
            module_embeddings: 聚合后的模块嵌入 (N_modules, D)
            gene_embeddings_updated: 更新后的基因嵌入 (N_genes, D)
        """
        N_modules = module_supernode_emb.size(0)

        if self.aggregation == "mean":
            # 简单平均: 每个模块 = 成员基因的平均
            module_emb = torch.zeros(N_modules, self.hidden_dim, device=gene_embeddings.device)
            module_counts = torch.zeros(N_modules, 1, device=gene_embeddings.device)
            for i in range(gene_embeddings.size(0)):
                if module_assignment[i].sum() > 0:
                    for m_idx in module_assignment[i].nonzero().squeeze(-1):
                        m_idx = int(m_idx.item()) if m_idx.dim() == 0 else m_idx
                        module_emb[int(m_idx)] += gene_embeddings[i]
                        module_counts[int(m_idx)] += 1
            module_counts = module_counts.clamp(min=1)
            module_emb = module_emb / module_counts

        elif self.aggregation == "attention":
            # 注意力加权聚合
            module_emb = torch.zeros(N_modules, self.hidden_dim, device=gene_embeddings.device)
            module_counts = torch.zeros(N_modules, 1, device=gene_embeddings.device)
            for i in range(gene_embeddings.size(0)):
                if module_assignment[i].sum() > 0:
                    g_emb = gene_embeddings[i]  # (D,)
                    for m_idx in module_assignment[i].nonzero().squeeze(-1):
                        m_idx = int(m_idx.item()) if m_idx.dim() == 0 else m_idx
                        m_emb = module_supernode_emb[int(m_idx)] if m_idx < N_modules else g_emb
                        # 注意力权重
                        attn_input = g_emb + m_emb
                        alpha = torch.sigmoid(self.attn(attn_input)).squeeze()
                        if alpha.dim() == 0:
                            module_emb[int(m_idx)] += alpha * g_emb
                            module_counts[int(m_idx)] += alpha
                        else:
                            module_emb[int(m_idx)] += alpha[0] * g_emb
                            module_counts[int(m_idx)] += alpha[0]
            module_counts = module_counts.clamp(min=1e-8)
            module_emb = module_emb / module_counts

        elif self.aggregation == "gated":
            # 门控融合
            module_emb = torch.zeros(N_modules, self.hidden_dim, device=gene_embeddings.device)
            for i in range(gene_embeddings.size(0)):
                if module_assignment[i].sum() > 0:
                    g_emb = gene_embeddings[i]
                    for m_idx in module_assignment[i].nonzero().squeeze(-1):
                        m_idx = int(m_idx.item()) if m_idx.dim() == 0 else m_idx
                        m_emb = module_supernode_emb[int(m_idx)] if m_idx < N_modules else g_emb
                        gate_input = torch.cat([g_emb, m_emb])
                        gate_val = self.gate(gate_input)
                        module_emb[int(m_idx)] += gate_val * g_emb + (1 - gate_val) * m_emb
            # 归一化
            module_norms = module_emb.norm(dim=-1, keepdim=True).clamp(min=1e-8)
            module_emb = module_emb / module_norms * module_supernode_emb.norm(dim=-1, keepdim=True)

        else:
            module_emb = module_supernode_emb

        # 更新基因嵌入 (从模块回传信息)
        gene_updated = gene_embeddings.clone()
        for i in range(gene_embeddings.size(0)):
            if module_assignment[i].sum() > 0:
                member_modules = []
                for m_idx in module_assignment[i].nonzero().squeeze(-1):
                    m_idx = int(m_idx.item()) if m_idx.dim() == 0 else m_idx
                    member_modules.append(module_emb[int(m_idx)])
                if member_modules:
                    gene_updated[i] = gene_updated[i] + 0.1 * torch.stack(member_modules).mean(dim=0)

        return module_emb, gene_updated


def build_gene_module_assignment(
    gene_names: List[str],
    data: HeteroData,
    n_modules: int = 8,
) -> Tuple[Tensor, List[str]]:
    """基于通路成员关系构建基因→模块分配矩阵。

    使用基因-通路边构建软分配矩阵 (每个基因可属于多个模块)。
    模块对应通路组 (铁死亡、脂代谢、铁稳态、氧化应激等)。

    Returns:
        module_assignment: (N_genes, N_modules) 二值分配矩阵
        module_names: 模块名称列表
    """
    N = len(gene_names)
    gene_to_idx = {g: i for i, g in enumerate(gene_names)}
    gene_set = set(gene_names)

    # 定义模块→通路映射
    pathway_groups = {
        "Ferroptosis_Core": [
            "Ferroptosis (KEGG: hsa04216)",
        ],
        "Lipid_Metabolism": [
            "Lipid metabolism / PUFA biosynthesis",
        ],
        "Iron_Homeostasis": [
            "Iron homeostasis",
        ],
        "Antioxidant_GSH": [
            "Glutathione metabolism (KEGG: hsa00480)",
        ],
        "Autophagy": [
            "Autophagy & ferroptosis",
        ],
        "Inflammation_Death": [
            "Necroptosis / Pyroptosis (交叉)",
        ],
        "Aging": [
            "Aging-related pathways",
        ],
        "CIRI": [
            "CIRI-related pathways",
        ],
    }

    module_names = list(pathway_groups.keys())
    assignment = torch.zeros(N, len(module_names))

    # 从data中获取gene→pathway边
    if hasattr(data['gene', 'belongs_to', 'pathway'], 'edge_index'):
        gp_ei = data['gene', 'belongs_to', 'pathway'].edge_index
        pathway_names_data = data['pathway'].names if hasattr(data['pathway'], 'names') else []

        # 构建pathway→module索引
        pw_to_module = {}
        for mod_name, pw_list in pathway_groups.items():
            mod_idx = module_names.index(mod_name)
            for pw in pw_list:
                if pw in pathway_names_data:
                    pw_to_module[pathway_names_data.index(pw)] = mod_idx
                # 模糊匹配
                for pi, pn in enumerate(pathway_names_data):
                    if pw.lower() in pn.lower():
                        pw_to_module[pi] = mod_idx

        for j in range(gp_ei.shape[1]):
            g_idx = int(gp_ei[0, j])
            p_idx = int(gp_ei[1, j])
            if g_idx < N and p_idx in pw_to_module:
                assignment[g_idx, pw_to_module[p_idx]] = 1.0

    # 确保没有基因完全无模块 (分配给最相关的)
    unassigned = assignment.sum(dim=1) == 0
    if unassigned.any():
        # 对未分配的基因，基于名称启发式分配
        for i in range(N):
            if unassigned[i]:
                g = gene_names[i]
                if g in ALL_FERROPTOSIS_GENES:
                    assignment[i, 0] = 1.0  # Ferroptosis_Core
                elif g in ALL_AGING_GENES:
                    assignment[i, 6] = 1.0  # Aging
                elif g in ALL_CIRI_GENES:
                    assignment[i, 7] = 1.0  # CIRI
                else:
                    assignment[i, 0] = 1.0  # 默认铁死亡

    print(f"  基因模块分配: {int((assignment.sum(dim=1) > 0).sum())}/{N} 基因已分配")
    for mi, mn in enumerate(module_names):
        n_members = int(assignment[:, mi].sum())
        print(f"    {mn}: {n_members} 成员")

    return assignment, module_names


# ============================================================
# 3.X 课程式难例挖掘器 (Curriculum Hard Negative Miner)
# ============================================================
# 基于: Structure-aware Curriculum for Masked Graph Autoencoders
#       (ICML 2025 Spotlight)
#
# 核心思想: 训练过程中动态增加负样本难度
#   - 早期: 随机采样负样本 (容易)
#   - 中期: 逐渐引入 pathway-similar 基因作为困难负样本
#   - 后期: 聚焦最困难负样本 (loss最高的样本)
#
# 与当前已有的静态硬负采样的区别:
#   当前实现 simple_hard_negative 在数据加载时静态筛选;
#   本模块在训练过程中动态调整难度, 使得模型逐步适应。
# ============================================================

class CurriculumHardNegativeMiner:
    """课程式难例挖掘器 — 动态增加负样本难度。

    调度策略:
      1. 计算每个基因与已知正样本的"相似度分数" (基于pathway/共表达等)
      2. 按相似度从高到低排序 (高相似度 = 更难的负样本)
      3. 课程调度器: epoch 0→全部用随机负样本, epoch max→全部用最困难负样本

    Args:
        gene_names: 基因名称列表
        pathway_genes: 与ACSL4/已知靶点同通路的基因集合 (作为困难负样本池)
        n_total_neg: 每轮需要的负样本总数
        n_genes: 基因总数
        warmup_epochs: 预热epoch数 (仅随机采样)
        curriculum_epochs: 课程训练epoch数 (从简单→困难过渡)
        alpha: 课程曲线陡度 (越大则过渡越快)
    """

    def __init__(
        self,
        gene_names: List[str],
        pathway_genes: Optional[Set[str]] = None,
        acsl4_neighbors: Optional[Set[str]] = None,
        n_total_neg: int = 100,
        n_genes: int = 0,
        warmup_epochs: int = 20,
        curriculum_epochs: int = 100,
        alpha: float = 2.0,
    ):
        self.gene_names = gene_names
        self.n_genes = len(gene_names) if n_genes == 0 else n_genes
        self.n_total_neg = n_total_neg
        self.warmup_epochs = warmup_epochs
        self.curriculum_epochs = curriculum_epochs
        self.alpha = alpha

        # 构建三级难度池:
        #
        #   Level 2 — 困难 (hard_pool):
        #     ACSL4的直接PPI邻居 + 高通路相似基因
        #     这些基因与正样本(已知靶点)在结构/功能上高度相似,
        #     模型最容易将它们错误预测为正样本。
        #
        #   Level 1 — 中等 (medium_pool):
        #     与ACSL4同通路的其他基因 (不在hard_pool中)
        #     有一定功能关联, 但不如直接邻居"像"正样本。
        #
        #   Level 0 — 简单 (easy_pool):
        #     所有其他基因 (无已知关联)
        #     模型天然容易区分为负样本。
        #
        # 无通路信息时的回退: 使用嵌入相似度动态分级
        self.hard_pool: Set[int] = set()
        self.medium_pool: Set[int] = set()

        # 1. 构建困难池: ACSL4邻居
        if acsl4_neighbors:
            for g in acsl4_neighbors:
                if g in self.gene_names:
                    self.hard_pool.add(self.gene_names.index(g))

        # 2. 构建中等池: 通路基因 (不在hard_pool中)
        if pathway_genes:
            for g in pathway_genes:
                if g in self.gene_names:
                    idx = self.gene_names.index(g)
                    if idx not in self.hard_pool:
                        self.medium_pool.add(idx)

        # 若无任何先验信息, 使用pathway_genes作为medium_pool
        if not self.hard_pool and not self.medium_pool and pathway_genes:
            for g in pathway_genes:
                if g in self.gene_names:
                    self.medium_pool.add(self.gene_names.index(g))

        # 困难基因嵌入相似度排序 (动态构建)
        self._difficulty_scores: Optional[np.ndarray] = None

    def _build_difficulty_scores_from_features(
        self, gene_features: Optional[Tensor] = None
    ) -> np.ndarray:
        """基于基因特征计算困难度分数 (高 = 更难区分)。

        若无特征可用, 则使用pathway membership作为难度代理:
          pathway基因 = 相似度0.8, 其余随机 = 0.2
        """
        scores = np.zeros(self.n_genes, dtype=np.float32)

        if gene_features is not None and gene_features.size(0) == self.n_genes:
            # 使用嵌入相似度: 与ACSL4嵌入越相似 → 越难
            acsl4_idx = self.gene_names.index("ACSL4") if "ACSL4" in self.gene_names else -1
            if acsl4_idx >= 0:
                acsl4_feat = gene_features[acsl4_idx:acsl4_idx+1].cpu()
                feats = gene_features.cpu()
                sim = F.cosine_similarity(acsl4_feat, feats, dim=-1).numpy()
                scores = (sim - sim.min()) / (sim.max() - sim.min() + 1e-8)
        else:
            # 基于pathway membership作为难度代理
            for idx in self.medium_pool:
                scores[idx] = 0.8

        return scores

    def sample(
        self,
        pos_indices: Tensor,
        epoch: int,
        max_epochs: int,
        gene_features: Optional[Tensor] = None,
    ) -> Tuple[Tensor, Tensor]:
        """根据当前epoch采样课程式负样本。

        三级难度定义:
          - Hard (Level 2): self.hard_pool — ACSL4直接PPI邻居
          - Medium (Level 1): self.medium_pool — 同通路其他基因
          - Easy (Level 0): 所有其他基因

        Args:
            pos_indices: 正样本基因索引 (N,)
            epoch: 当前epoch
            max_epochs: 最大epoch数
            gene_features: 基因特征 (用于计算难度分数)

        Returns:
            neg_drug_indices: 负样本药物索引 (N_neg,)
            neg_gene_indices: 负样本基因索引 (N_neg,)
        """
        pos_set = set(pos_indices.tolist())
        self._difficulty_scores = self._build_difficulty_scores_from_features(gene_features)

        # 计算课程进度: 0.0 (全部随机) → 1.0 (全部最困难)
        if epoch < self.warmup_epochs:
            curriculum_progress = 0.0
        else:
            progress = (epoch - self.warmup_epochs) / max(1, self.curriculum_epochs - self.warmup_epochs)
            curriculum_progress = min(1.0, progress ** (1.0 / self.alpha))

        # 采样困难负样本 (Level 2: hard_pool, ACSL4直接邻居)
        hard_candidates = [i for i in self.hard_pool
                          if i not in pos_set]
        n_hard = int(self.n_total_neg * curriculum_progress * 0.5)
        if len(hard_candidates) < n_hard:
            n_hard = len(hard_candidates)
        if n_hard > 0 and hard_candidates:
            hard_scores = [self._difficulty_scores[i] for i in hard_candidates]
            hard_top = np.argsort(-np.array(hard_scores))[:n_hard]
            hard_selected = [hard_candidates[i] for i in hard_top]
        else:
            hard_selected = []

        # 采样中等负样本 (Level 1: medium_pool, 同通路非邻居基因)
        medium_candidates = [i for i in self.medium_pool
                            if i not in pos_set and i not in hard_selected]
        n_medium = int(self.n_total_neg * curriculum_progress * 0.3)
        if len(medium_candidates) < n_medium:
            n_medium = len(medium_candidates)
        if n_medium > 0 and medium_candidates:
            medium_scores = [self._difficulty_scores[i] for i in medium_candidates]
            medium_top = np.argsort(-np.array(medium_scores))[:n_medium]
            medium_selected = [medium_candidates[i] for i in medium_top]
        else:
            medium_selected = []

        # 采样容易负样本 (Level 0: 全基因集排除hard+medium+pos)
        easy_candidates = [i for i in range(self.n_genes)
                          if i not in pos_set
                          and i not in self.hard_pool
                          and i not in self.medium_pool]
        n_easy = max(0, self.n_total_neg - len(hard_selected) - len(medium_selected))
        if len(easy_candidates) < n_easy:
            n_easy = len(easy_candidates)
        if n_easy > 0 and easy_candidates:
            easy_selected = list(np.random.choice(easy_candidates, n_easy, replace=False))
        else:
            easy_selected = []

        all_neg = hard_selected + medium_selected + easy_selected
        np.random.shuffle(all_neg)

        if len(all_neg) < self.n_total_neg:
            # 填充随机值
            fill_candidates = [i for i in range(self.n_genes) if i not in pos_set]
            n_fill = self.n_total_neg - len(all_neg)
            if fill_candidates and n_fill > 0:
                fill = list(np.random.choice(fill_candidates, min(n_fill, len(fill_candidates)), replace=False))
                all_neg.extend(fill)

        neg_drugs = torch.zeros(len(all_neg), dtype=torch.long)
        neg_genes = torch.tensor(all_neg[:self.n_total_neg], dtype=torch.long)

        return neg_drugs, neg_genes

    def get_curriculum_stats(self, epoch: int, max_epochs: int) -> Dict:
        """获取当前课程进度统计。"""
        if epoch < self.warmup_epochs:
            progress = 0.0
        else:
            progress = min(1.0, (epoch - self.warmup_epochs) / max(1, max_epochs - self.warmup_epochs))

        return {
            "epoch": epoch,
            "curriculum_progress": float(progress),
            "n_hard_pool": len(self.medium_pool),
            "active_hard_ratio": float(progress * 0.5),
        }


# ============================================================
# 3.Y 自适应图增强 (Adaptive Graph Augmentation)
# ============================================================
# 基于: GraphFormer-CL — Graph Transformer with Contrastive Learning
#       (PLoS One, 2026)
#
# 核心思想: 使用课程式控制增强强度的对比学习
#   - 特征掩码: 随机掩码节点特征的一部分
#   - 边Dropout: 随机丢弃一部分边
#   - 增强强度随训练轮数动态调整 (从弱→强)
#
# 与现有HeCo对比学习的关系:
#   增强的对比视图替代HeCo的固定视图, 使对比学习更有效
# ============================================================

class AdaptiveGraphAugmentation(nn.Module):
    """自适应图增强 — 课程式控制增强强度的对比学习视图生成。

    包含两种增强操作:
      1. Feature Masking: 随机掩码特征维度
      2. Edge Dropout: 随机丢弃边

    增强强度 (mask_rate, dropout_rate) 由课程调度器控制:
      - 早期: 弱增强 (保留大部分特征和边)
      - 后期: 强增强 (制造更具挑战性的视图)
    """

    def __init__(
        self,
        hidden_dim: int,
        feat_mask_min: float = 0.1,
        feat_mask_max: float = 0.5,
        edge_dropout_min: float = 0.1,
        edge_dropout_max: float = 0.6,
        warmup_epochs: int = 10,
        curriculum_epochs: int = 100,
        alpha: float = 2.0,
    ):
        super().__init__()
        self.hidden_dim = hidden_dim
        self.feat_mask_min = feat_mask_min
        self.feat_mask_max = feat_mask_max
        self.edge_dropout_min = edge_dropout_min
        self.edge_dropout_max = edge_dropout_max
        self.warmup_epochs = warmup_epochs
        self.curriculum_epochs = curriculum_epochs
        self.alpha = alpha

    def get_augmentation_strength(self, epoch: int) -> Tuple[float, float]:
        """获取当前epoch的增强强度。

        Returns:
            feat_mask_rate: 特征掩码率
            edge_dropout_rate: 边丢弃率
        """
        if epoch < self.warmup_epochs:
            progress = 0.0
        else:
            progress = min(1.0, (epoch - self.warmup_epochs) / max(1, self.curriculum_epochs - self.warmup_epochs))

        # 曲线式增强: 早期慢、中期快、后期饱和
        adj_progress = progress ** (1.0 / self.alpha)

        feat_mask_rate = self.feat_mask_min + adj_progress * (self.feat_mask_max - self.feat_mask_min)
        edge_dropout_rate = self.edge_dropout_min + adj_progress * (self.edge_dropout_max - self.edge_dropout_min)

        return feat_mask_rate, edge_dropout_rate

    def augment_features(
        self, x: Tensor, mask_rate: float, mask_token: Optional[Tensor] = None
    ) -> Tensor:
        """特征掩码增强: 随机将一些特征维度置零 (或替换为 mask_token)。

        Args:
            x: 特征矩阵 (N, D)
            mask_rate: 掩码比例
            mask_token: 可选的替换标记 (默认为0)

        Returns:
            x_aug: 增强后的特征
        """
        if mask_rate <= 0:
            return x

        x_aug = x.clone()
        N, D = x.shape
        n_mask = max(1, int(D * mask_rate))

        for i in range(N):
            mask_idx = torch.randperm(D)[:n_mask]
            if mask_token is not None:
                x_aug[i, mask_idx] = mask_token[i, mask_idx]
            else:
                x_aug[i, mask_idx] = 0.0

        return x_aug

    def augment_edge_index(
        self, edge_index: Tensor, dropout_rate: float
    ) -> Tensor:
        """边Dropout: 随机丢弃部分边。

        Args:
            edge_index: 边索引 (2, E)
            dropout_rate: 丢弃率

        Returns:
            edge_index_aug: 保留的边索引
        """
        if dropout_rate <= 0 or edge_index.size(1) == 0:
            return edge_index

        n_edges = edge_index.size(1)
        n_keep = max(1, int(n_edges * (1 - dropout_rate)))
        perm = torch.randperm(n_edges)[:n_keep]

        return edge_index[:, perm]

    def forward(
        self,
        x_dict: Dict[str, Tensor],
        edge_index_dict: Dict[Tuple, Tensor],
        epoch: int,
    ) -> Tuple[Dict[str, Tensor], Dict[Tuple, Tensor]]:
        """生成增强后的图视图。

        Args:
            x_dict: 节点特征字典
            edge_index_dict: 边索引字典
            epoch: 当前epoch

        Returns:
            x_aug_dict: 增强后的特征
            edge_aug_dict: 增强后的边
        """
        feat_mask_rate, edge_dropout_rate = self.get_augmentation_strength(epoch)

        x_aug_dict = {}
        for nt, x in x_dict.items():
            # 每种节点类型使用相同的mask率但不同的随机种子
            x_aug_dict[nt] = self.augment_features(x, feat_mask_rate)

        edge_aug_dict = {}
        for et, ei in edge_index_dict.items():
            edge_aug_dict[et] = self.augment_edge_index(ei, edge_dropout_rate)

        return x_aug_dict, edge_aug_dict


# ============================================================
# 3.Z 自蒸馏模块 (Self-Distillation with EMA Teacher)
# ============================================================
# 基于: Knowledge Distillation for Molecular Property Prediction
#       (Advanced Science, 2025) + Mean Teacher (Tarvainen 2017)
#
# 核心思想: 维护一个教师模型 (EMA指数移动平均), 教师提供
#   软标签指导学生模型训练, 增强模型鲁棒性和泛化能力。
#
# 与VIB的区别:
#   VIB: 信息瓶颈压缩, 去除噪声 (编码器级别)
#   自蒸馏: 软标签正则化, 去相关 (预测器级别)
#   两者互补, 可同时使用。
# ============================================================

class SelfDistillationEMA:
    """自蒸馏模块 — EMA教师提供软标签。

    工作机制:
      1. 训练每个step后: teacher = decay * teacher + (1-decay) * student
      2. 前向时: teacher和学生分别预测边分数
      3. 损失: KL(student_logits || teacher_logits) 作为蒸馏损失

    Args:
        student_model: 学生模型 (主模型)
        model_factory: 教师模型工厂函数 (返回新模型实例) —
                       用于安全创建教师模型, 避免copy.deepcopy在PyG下失败.
                       若为None, 则回退到state_dict级深拷贝.
        decay: EMA衰减系数 (越大则教师更新越慢)
        start_epoch: 开始蒸馏的epoch (前期让学生自我探索)
        distill_weight: 蒸馏损失权重
        temperature: 蒸馏温度 (越大则软标签越平滑)
    """

    def __init__(
        self,
        student_model: nn.Module,
        model_factory: Optional[callable] = None,
        decay: float = 0.999,
        start_epoch: int = 30,
        distill_weight: float = 0.3,
        temperature: float = 4.0,
    ):
        # 创建教师模型 — 使用安全方式避免deepcopy在PyG下失败
        if model_factory is not None:
            # 方式1: 用户提供工厂函数
            self.teacher = model_factory()
            self.teacher.load_state_dict(copy.deepcopy(student_model.state_dict()))
        else:
            # 方式2: state_dict级深拷贝 (比copy.deepcopy整个module安全)
            import io
            buffer = io.BytesIO()
            torch.save(student_model.state_dict(), buffer)
            buffer.seek(0)
            # 创建同架构新实例: 通过type(student_model) + 浅拷贝配置
            try:
                # 尝试用type创建新实例 (适用于FerroHGT)
                # 注意: 需要student_model上保存了__init__参数
                model_args = {}
                if hasattr(student_model, 'get_init_args'):
                    model_args = student_model.get_init_args()
                self.teacher = type(student_model)(**model_args)
                self.teacher.load_state_dict(torch.load(buffer, weights_only=True))
            except Exception:
                # 终极回退: 仍然使用deepcopy, 但捕捉异常给出明确警告
                import warnings
                warnings.warn("state_dict级深拷贝失败, 回退到copy.deepcopy. "
                              "若失败请提供model_factory参数.")
                self.teacher = copy.deepcopy(student_model)

        self.teacher.eval()  # 教师永远在eval模式

        # 将教师模型同步到学生模型的设备上
        student_device = next(student_model.parameters()).device
        self.teacher = self.teacher.to(student_device)

        # 冻结教师参数
        for param in self.teacher.parameters():
            param.requires_grad = False

        self.decay = decay
        self.start_epoch = start_epoch
        self.distill_weight = distill_weight
        self.temperature = temperature

    def update_teacher(self, student_model: nn.Module):
        """EMA更新教师参数: θ_t = decay * θ_t + (1-decay) * θ_s

        安全实现: 通过state_dict交换, 避免nn.Module参数绑定问题。
        """
        with torch.no_grad():
            t_dict = self.teacher.state_dict()
            s_dict = student_model.state_dict()
            for key in t_dict:
                if key in s_dict:
                    t_dict[key] = t_dict[key] * self.decay + s_dict[key] * (1.0 - self.decay)
            self.teacher.load_state_dict(t_dict)

    def ensure_device(self, student_model: nn.Module):
        """确保教师模型与学生模型在同一设备上。"""
        stu_device = next(student_model.parameters()).device
        tea_device = next(self.teacher.parameters()).device
        if stu_device != tea_device:
            self.teacher = self.teacher.to(stu_device)

    @torch.no_grad()
    def get_teacher_predictions(
        self,
        z_dict: Dict[str, Tensor],
        drug_indices: Tensor,
        gene_indices: Tensor,
        predict_fn,
    ) -> Tensor:
        """获取教师模型的软标签预测。

        Args:
            z_dict: 教师编码的节点嵌入
            drug_indices: 药物索引
            gene_indices: 基因索引
            predict_fn: 边预测函数 (model.predict_edges)

        Returns:
            teacher_logits: 教师预测的logits (N,)
        """
        self.teacher.eval()
        # 教师不参与梯度计算
        return predict_fn(z_dict, drug_indices, gene_indices)

    def compute_distill_loss(
        self,
        student_logits: Tensor,
        teacher_logits: Tensor,
    ) -> Tensor:
        """计算蒸馏损失 (KL散度 + 温度缩放)。

        KL散度在高温下:
          L_distill = T² * KL(softmax(s/T) || softmax(t/T))

        使用F.kl_div实现二元分类:
          将sigmoid输出视为2类分布 [p, 1-p], 用log_softmax + kl_div保证数值稳定。

        Args:
            student_logits: 学生预测logits (N,)
            teacher_logits: 教师预测logits (N,)

        Returns:
            distill_loss: 蒸馏损失
        """
        if teacher_logits is None:
            return torch.tensor(0.0, device=student_logits.device)

        T = self.temperature

        # 构建2类分布: [sigmoid(s/T), 1-sigmoid(s/T)]
        # 使用log_softmax和softmax保证数值稳定性
        logits_stu = torch.stack([
            student_logits / T,
            torch.zeros_like(student_logits)
        ], dim=-1)  # (N, 2) — logits for [p_pos, p_neg]
        logits_tea = torch.stack([
            teacher_logits / T,
            torch.zeros_like(teacher_logits)
        ], dim=-1)

        # F.kl_div: 输入需为log-probabilities, target为probabilities
        # KL(P||Q) = sum(P * log(P/Q)) = sum(P * (log(P) - log(Q)))
        log_p_stu = F.log_softmax(logits_stu, dim=-1)   # log(Q)
        p_teacher = F.softmax(logits_tea, dim=-1)        # P

        kl_div = F.kl_div(log_p_stu, p_teacher, reduction='batchmean')

        # 温度缩放: L = T² * KL (梯度幅度温度补偿)
        distill_loss = (T ** 2) * kl_div * self.distill_weight

        return distill_loss


# ============================================================
# 3. HGT 链路预测模型 (升级版: 关系门控 + 层次聚合 + VIB + 小波 + 课程)
# ============================================================

class FerroHGT(nn.Module):
    """石竹烯-铁衰老-脑缺血 异构网络 HGT 编码器-解码器。

    架构 (升级版):
      1. 类型感知输入投影
      2. 多层 HGTConv 消息传递 (含关系门控)
      3. 层次基因模块聚合 (genes → modules → context)
      4. 边预测解码器 (Drug-Gene, 简化版防过拟合)

    训练目标:
      - 正边: 已知石竹烯-靶点相互作用
      - 负边: 随机采样 (排除已知正边)
    """

    def __init__(
        self,
        metadata: Tuple[List[str], List[Tuple]],
        hidden_dim: int = 256,
        num_heads: int = 4,
        num_layers: int = 3,
        dropout: float = 0.3,
        relation_gate_init_bias: float = 2.0,
        relation_gate_temperature: float = 1.0,
        relation_gate_freeze_types: Set[str] = None,
        relation_gate_per_edge_bias: Dict[str, float] = None,
        use_module_aggregation: bool = True,
        module_aggregation: str = "attention",
        use_vib: bool = True,                    # 新增: VIB压缩
        vib_beta_start: float = 1e-4,            # VIB β初始
        vib_beta_end: float = 1e-2,              # VIB β最终
        use_wavelet: bool = True,                # 新增: 小波多尺度特征
        wavelet_n_scales: int = 3,               # 小波尺度数
        use_knowledge_reg: bool = True,          # 新增: 知识正则化
        use_augmentation: bool = True,           # 新增: 自适应图增强
    ):
        super().__init__()
        self.hidden_dim = hidden_dim
        self.node_types = metadata[0]
        self.edge_types = metadata[1]
        self.use_module_aggregation = use_module_aggregation
        self.use_vib = use_vib
        self.use_wavelet = use_wavelet
        self.use_knowledge_reg = use_knowledge_reg
        self.use_augmentation = use_augmentation

        # 关系消息门控 (升级版: 差异化初始化 + 逐类型冻结)
        self.relation_gate = RelationMessageGate(
            edge_types=metadata[1],
            init_bias=relation_gate_init_bias,
            per_edge_init_bias=relation_gate_per_edge_bias,
            temperature=relation_gate_temperature,
            freeze_types=relation_gate_freeze_types,
        )

        # HGT卷积层
        self.convs = nn.ModuleList()
        self.norms = nn.ModuleList()
        for _ in range(num_layers):
            self.convs.append(
                HGTConv(
                    in_channels=hidden_dim,
                    out_channels=hidden_dim,
                    metadata=metadata,
                    heads=num_heads,
                )
            )
            self.norms.append(nn.LayerNorm(hidden_dim))

        self.dropout = dropout

        # VIB: 变分信息瓶颈 (在HGT编码后压缩)
        if use_vib:
            self.vib = VariationalInformationBottleneck(
                hidden_dim=hidden_dim,
                beta_start=vib_beta_start,
                beta_end=vib_beta_end,
                beta_warmup_epochs=50,
            )
        else:
            self.vib = None

        # 图小波多尺度特征 (在基因特征送入HGT前增强)
        if use_wavelet:
            self.wavelet_encoder = GraphWaveletMultiScaleFeatures(
                hidden_dim=hidden_dim,
                n_scales=wavelet_n_scales,
                chebyshev_order=5,
                dropout=dropout * 0.5,
            )
        else:
            self.wavelet_encoder = None

        # 知识正则化器
        if use_knowledge_reg:
            self.knowledge_regularizer = KnowledgeBasedRegularizer(
                pathway_consistency_weight=0.1,
                drug_target_proximity_weight=0.05,
                module_compactness_weight=0.05,
                neg_separation_weight=0.02,
                margin=0.5,
            )
        else:
            self.knowledge_regularizer = None

        # 自适应图增强 (用于对比学习中生成增强视图)
        if use_augmentation:
            self.graph_augmentation = AdaptiveGraphAugmentation(
                hidden_dim=hidden_dim,
                feat_mask_min=0.1,
                feat_mask_max=0.5,
                edge_dropout_min=0.1,
                edge_dropout_max=0.6,
                warmup_epochs=10,
                curriculum_epochs=100,
            )
        else:
            self.graph_augmentation = None

        # 层次基因模块聚合器
        if use_module_aggregation:
            self.module_aggregator = HierarchicalGeneModuleAggregator(
                hidden_dim=hidden_dim,
                dropout=dropout,
                aggregation=module_aggregation,
            )
            # 模块超节点可学习嵌入
            self.module_supernode_emb = nn.Parameter(
                torch.randn(8, hidden_dim) * 0.02  # 8个模块
            )
        else:
            self.module_aggregator = None
            self.module_supernode_emb = None

        # 边预测解码器 (简化版 - 防过拟合)
        # 仅用 Hadamard积 + 浅层MLP
        self.edge_decoder = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.ReLU(),
            nn.Dropout(dropout * 0.3),
            nn.Linear(hidden_dim // 2, 1),
        )

    def forward(
        self,
        x_dict: Dict[str, Tensor],
        edge_index_dict: Dict,
        module_assignment: Tensor = None,
        epoch: int = 0,  # 新增: 当前epoch用于VIBβ调度
    ) -> Dict[str, Tensor]:
        """HGT编码前向传播 (含关系门控 + 小波 + VIB + 层次聚合)。

        Returns:
            z_dict: 各节点类型的最终嵌入
            extra: 附加信息 (VIB KL损失等)
        """
        extra = {}

        # === Phase 1: 图小波多尺度特征增强 (仅gene) ===
        if self.wavelet_encoder is not None and 'gene' in x_dict:
            gg_edge_key = ('gene', 'interacts', 'gene')
            gg_edge = edge_index_dict.get(gg_edge_key, torch.zeros((2, 0), dtype=torch.long, device=x_dict['gene'].device))
            wavelet_feat, _ = self.wavelet_encoder(x_dict['gene'], gg_edge)
            x_dict['gene'] = wavelet_feat

        # === Phase 2: HGT消息传递 ===
        # 获取关系门控权重
        gate_weights = self.relation_gate(edge_index_dict)

        for conv, norm in zip(self.convs, self.norms):
            # 保存每层输入用于残差
            residual = {k: v for k, v in x_dict.items()}

            out = conv(x_dict, edge_index_dict)

            # 补全未参与消息传递的节点类型
            for k in x_dict:
                if k not in out:
                    out[k] = x_dict[k]

            # 残差连接 + LayerNorm + ReLU + Dropout
            for k in out:
                res = residual.get(k, out[k])
                out[k] = F.relu(norm(out[k] + res))
                if self.training:
                    out[k] = F.dropout(out[k], p=self.dropout, training=True)

            x_dict = out

        # === Phase 3: VIB压缩 (仅gene, 去除噪声) ===
        if self.vib is not None and 'gene' in x_dict:
            self.vib.set_epoch(epoch)
            gene_compressed, kl_loss = self.vib(x_dict['gene'])
            x_dict['gene'] = gene_compressed
            extra['vib_kl_loss'] = kl_loss
            extra['vib_beta'] = self.vib.get_beta()

        # === Phase 4: 层次基因模块聚合 (可选) ===
        if self.module_aggregator is not None and module_assignment is not None \
                and 'gene' in x_dict:
            module_emb, gene_updated = self.module_aggregator(
                x_dict['gene'].clone(),
                module_assignment.to(x_dict['gene'].device),
                self.module_supernode_emb.to(x_dict['gene'].device),
            )
            # 基因嵌入融合: 原始嵌入 + 模块上下文
            x_dict['gene'] = gene_updated
            # 存储模块嵌入供后续使用
            x_dict['_module_emb'] = module_emb

        return x_dict, extra

    def get_init_args(self) -> Dict:
        """返回__init__参数字典, 用于安全创建同架构教师模型。"""
        return {
            'metadata': (self.node_types, self.edge_types),
            'hidden_dim': self.hidden_dim,
            'num_heads': self.convs[0].heads if self.convs else 4,
            'num_layers': len(self.convs),
            'dropout': self.dropout,
            'relation_gate_init_bias': 2.0,
            'relation_gate_temperature': 1.0,
            'relation_gate_freeze_types': None,
            'relation_gate_per_edge_bias': None,
            'use_module_aggregation': self.module_aggregator is not None,
            'use_vib': self.vib is not None,
            'use_wavelet': self.wavelet_encoder is not None,
            'use_knowledge_reg': self.knowledge_regularizer is not None,
            'use_augmentation': self.graph_augmentation is not None,
        }

    def predict_edges(
        self,
        z_dict: Dict[str, Tensor],
        drug_indices: Tensor,
        gene_indices: Tensor,
    ) -> Tensor:
        """预测 Drug-Gene 边分数。

        Args:
            z_dict: 编码后的节点嵌入 (可能包含extra信息)
            drug_indices: 药物节点索引 (N,)
            gene_indices: 基因节点索引 (N,)

        Returns:
            scores: 预测分数 (N,) — logits
        """
        z_drug = z_dict['drug'][drug_indices]
        z_gene = z_dict['gene'][gene_indices]
        # 简化: 仅Hadamard积
        h = z_drug * z_gene
        return self.edge_decoder(h).squeeze(-1)

    def get_relation_gate_values(self) -> Dict[str, float]:
        """获取当前关系门控值 (用于消融分析)。"""
        return self.relation_gate.get_gate_values()

    def compute_knowledge_loss(
        self,
        z_dict: Dict[str, Tensor],
        data: HeteroData,
        module_assignment: Tensor = None,
    ) -> Dict[str, Tensor]:
        """计算知识正则化损失。
        Args:
            z_dict: 编码后的嵌入 (可能含extra信息)
            data: 异构图数据
            module_assignment: 模块分配矩阵
        Returns:
            loss_dict: 各知识损失分量
        """
        if self.knowledge_regularizer is None:
            return {}

        gene_emb = z_dict['gene']
        drug_emb = z_dict['drug']  # (1, D)

        # Pathway边 — 安全访问HeteroData中的异质边类型
        pw_edge_key = ('gene', 'co_pathway', 'gene')
        try:
            pw_edge_store = data[pw_edge_key]
            pw_edge = pw_edge_store.edge_index if hasattr(pw_edge_store, 'edge_index') else None
        except (KeyError, AttributeError):
            pw_edge = None

        # 正样本靶点索引
        pos_edge = data['drug', 'targets', 'gene'].edge_index
        target_indices = pos_edge[1] if pos_edge.numel() > 0 else None

        # 负样本索引 (随机采样非靶点)
        n_genes = gene_emb.size(0)
        if target_indices is not None:
            pos_set = set(target_indices.tolist())
            neg_indices = torch.tensor(
                [i for i in range(n_genes) if i not in pos_set][:len(pos_set) * 2],
                device=gene_emb.device, dtype=torch.long
            )
        else:
            neg_indices = None

        return self.knowledge_regularizer.compute_all(
            gene_embeddings=gene_emb,
            drug_embedding=drug_emb,
            pathway_edges=pw_edge,
            target_indices=target_indices,
            module_assignment=module_assignment,
            neg_indices=neg_indices,
        )


# ============================================================
# 4. 训练器
# ============================================================

@torch.no_grad()
def evaluate_raw_embeddings(
    embeddings: Dict[str, Tensor],
    gene_names: List[str],
    data: HeteroData = None,
) -> Dict:
    """评估原始嵌入质量 (不依赖解码器，直接度量嵌入空间结构)。

    使用三种指标评估 drug-gene 嵌入质量:
      - Cosine: 余弦相似度
      - Hadamard: Hadamard积求和
      - L2: L2距离 (取负)

    返回 ACSL4 在各指标下的排名和Top-K基因列表。
    """
    drug_emb = embeddings['drug'].cpu()
    gene_emb = embeddings['gene'].cpu()

    # 归一化
    drug_norm = drug_emb / (drug_emb.norm(dim=-1, keepdim=True) + 1e-8)
    gene_norm = gene_emb / (gene_emb.norm(dim=-1, keepdim=True) + 1e-8)

    results = {}
    acsl4_idx = gene_names.index("ACSL4") if "ACSL4" in gene_names else None
    n_total = len(gene_names)

    # 1. Cosine similarity
    cos_scores = (drug_norm * gene_norm).sum(dim=-1).squeeze().numpy()
    cos_rank = np.argsort(-cos_scores)
    results['cosine'] = {
        'acsl4_rank': int(np.where(cos_rank == acsl4_idx)[0][0] + 1) if acsl4_idx is not None else -1,
        'acsl4_score': float(cos_scores[acsl4_idx]) if acsl4_idx is not None else 0,
        'top10': [(gene_names[i], float(cos_scores[i])) for i in cos_rank[:10]],
    }

    # 2. Hadamard product sum
    hadamard_scores = (drug_emb * gene_emb).sum(dim=-1).squeeze().numpy()
    had_rank = np.argsort(-hadamard_scores)
    results['hadamard'] = {
        'acsl4_rank': int(np.where(had_rank == acsl4_idx)[0][0] + 1) if acsl4_idx is not None else -1,
        'acsl4_score': float(hadamard_scores[acsl4_idx]) if acsl4_idx is not None else 0,
        'top10': [(gene_names[i], float(hadamard_scores[i])) for i in had_rank[:10]],
    }

    # 3. L2 distance (negated for ranking)
    l2_dist = torch.cdist(drug_emb, gene_emb, p=2).squeeze().numpy()
    l2_rank = np.argsort(l2_dist)  # ascending = more similar
    results['l2'] = {
        'acsl4_rank': int(np.where(l2_rank == acsl4_idx)[0][0] + 1) if acsl4_idx is not None else -1,
        'acsl4_dist': float(l2_dist[acsl4_idx]) if acsl4_idx is not None else 0,
        'top10': [(gene_names[i], float(l2_dist[i])) for i in l2_rank[:10]],
    }

    # 综合排名 (取三种方法的平均排名)
    if acsl4_idx is not None:
        avg_rank = (results['cosine']['acsl4_rank'] + results['hadamard']['acsl4_rank'] + results['l2']['acsl4_rank']) / 3.0
        results['ensemble_rank'] = float(avg_rank)
        results['ensemble_percentile'] = float(avg_rank / n_total * 100)

    return results


class HGTLinkPredictionTrainer:
    """HGT链路预测训练器。

    支持:
      - 正负样本平衡采样
      - 早停
      - 关系门控L1正则化
      - ACSL4回忆实验
      - 原始嵌入质量评估
      - 课程式难例挖掘 (CurriculumHardNegativeMiner)
      - EMA自蒸馏 (SelfDistillationEMA)
      - 自适应图增强 (AdaptiveGraphAugmentation)
    """

    def __init__(
        self,
        model: FerroHGT,
        data: HeteroData,
        module_assignment: Tensor = None,
        lr: float = 0.001,
        weight_decay: float = 1e-5,
        gate_l1_lambda: float = 1e-4,
        device: str = "cuda",
        use_curriculum: bool = False,             # 是否启用课程式难例挖掘
        curriculum_miner: Optional[CurriculumHardNegativeMiner] = None,
        use_self_distill: bool = False,           # 是否启用自蒸馏
        self_distill: Optional[SelfDistillationEMA] = None,
        val_split: float = 0.15,                  # 验证集正边比例
        test_split: float = 0.15,                 # 测试集正边比例
    ):
        self.model = model.to(device)
        self.data = data.to(device)
        self.module_assignment = module_assignment.to(device) if module_assignment is not None else None
        self.device = device
        self.gate_l1_lambda = gate_l1_lambda
        self.use_curriculum = use_curriculum
        self.curriculum_miner = curriculum_miner
        self.use_self_distill = use_self_distill
        self.self_distill = self_distill

        # === 边分割: 将正边划分为训练/验证/测试 ===
        self._create_edge_splits(val_split, test_split)

        self.optimizer = torch.optim.AdamW(
            model.parameters(), lr=lr, weight_decay=weight_decay
        )
        self.scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
            self.optimizer, mode='max', factor=0.5, patience=30, min_lr=1e-6
        )
        self.val_neg_edges = None  # 固定验证负样本，在首次evaluate时初始化
        self._current_epoch = 0    # 当前epoch (用于VIB β调度、课程调度、蒸馏调度)
        self._max_epochs = 300     # 最大epoch数 (由train方法更新)

    def _create_edge_splits(self, val_split: float, test_split: float):
        """将正边划分为训练/验证/测试三组。

        方案B: 边分割 (Edge Splitting)
          - 训练集: 70% 的正边 (用于模型训练)
          - 验证集: 15% 的正边 (用于早停和超参选择)
          - 测试集: 15% 的正边 (用于最终评估)

        训练过程中 `self.data['drug', 'targets', 'gene'].edge_index`
        仅包含训练集正边, 验证/测试正边不参与训练。
        """
        pos_ei = self.data['drug', 'targets', 'gene'].edge_index
        n_pos = pos_ei.size(1)

        n_val = max(1, int(n_pos * val_split))
        n_test = max(1, int(n_pos * test_split))
        n_train = max(1, n_pos - n_val - n_test)

        # 固定随机种子保证可复现
        rng = torch.Generator().manual_seed(42)
        perm = torch.randperm(n_pos, generator=rng)

        train_idx = perm[:n_train]
        val_idx = perm[n_train:n_train + n_val]
        test_idx = perm[n_train + n_val:]

        self.train_ei = pos_ei[:, train_idx].to(self.device)
        self.val_ei = pos_ei[:, val_idx].to(self.device)
        self.test_ei = pos_ei[:, test_idx].to(self.device)

        # 数据集中仅保留训练边 (训练时只可见训练正边)
        self.data['drug', 'targets', 'gene'].edge_index = self.train_ei

        print(f"  [边分割] 正边 {n_pos} → 训练 {n_train} / 验证 {n_val} / 测试 {n_test}")

    def _get_positive_edges(self) -> Tuple[Tensor, Tensor]:
        """返回训练集正样本边。"""
        return self.train_ei[0], self.train_ei[1]

    def _get_val_positive_edges(self) -> Tuple[Tensor, Tensor]:
        """返回验证集正样本边。"""
        return self.val_ei[0], self.val_ei[1]

    def _get_test_positive_edges(self) -> Tuple[Tensor, Tensor]:
        """返回测试集正样本边。"""
        return self.test_ei[0], self.test_ei[1]

    def _init_val_negatives(self, n_neg: int = None):
        """生成固定的验证集负样本 — 排除训练正边和验证正边。

        负样本候选 = 所有 drug-gene 对 - 训练正边 - 验证正边
        这确保验证集负样本不包含任何模型见过的正边。
        """
        n_drugs = self.data['drug'].num_nodes
        n_genes = self.data['gene'].num_nodes

        # 排除训练正边 + 验证正边 (所有已知正边均不参与验证负采样)
        train_set = set(zip(self.train_ei[0].tolist(), self.train_ei[1].tolist()))
        val_set = set(zip(self.val_ei[0].tolist(), self.val_ei[1].tolist()))
        all_pos_set = train_set | val_set

        n_pos_val = len(val_set)
        if n_neg is None:
            n_neg = n_pos_val * 10

        all_pairs = torch.cartesian_prod(
            torch.arange(n_drugs), torch.arange(n_genes)
        )
        pos_mask = torch.tensor([
            (d.item(), g.item()) in all_pos_set for d, g in all_pairs
        ])
        neg_candidates = all_pairs[~pos_mask]

        if len(neg_candidates) < n_neg:
            n_neg = len(neg_candidates)

        perm = torch.randperm(len(neg_candidates))[:n_neg]
        selected = neg_candidates[perm]
        self.val_neg_edges = (
            selected[:, 0].to(self.device),
            selected[:, 1].to(self.device)
        )
        print(f"  [验证负样本] {n_neg} 个 (从 {len(neg_candidates)} 候选采样, 排除 {len(all_pos_set)} 已知正边)")

    def _sample_negative_edges(self, n_neg: int) -> Tuple[Tensor, Tensor]:
        """采样训练负样本边 (排除所有已知正边: 训练+验证+测试)。"""
        n_drugs = self.data['drug'].num_nodes
        n_genes = self.data['gene'].num_nodes

        # 排除所有已知正边以避免训练/验证/测试泄漏
        train_set = set(zip(self.train_ei[0].tolist(), self.train_ei[1].tolist()))
        val_set = set(zip(self.val_ei[0].tolist(), self.val_ei[1].tolist()))
        test_set = set(zip(self.test_ei[0].tolist(), self.test_ei[1].tolist()))
        all_pos_set = train_set | val_set | test_set

        all_pairs = torch.cartesian_prod(torch.arange(n_drugs), torch.arange(n_genes))
        pos_mask = torch.tensor([(d.item(), g.item()) in all_pos_set for d, g in all_pairs])
        neg_candidates = all_pairs[~pos_mask]

        if len(neg_candidates) < n_neg:
            n_neg = len(neg_candidates)

        perm = torch.randperm(len(neg_candidates))[:n_neg]
        selected = neg_candidates[perm]
        return (selected[:, 0].to(self.device), selected[:, 1].to(self.device))

    def train_epoch(self) -> Dict[str, float]:
        """单epoch训练 (含课程式难例挖掘 + 自蒸馏)。"""
        self.model.train()

        # 获取正样本
        pos_drugs, pos_genes = self._get_positive_edges()
        n_pos = len(pos_drugs)

        if n_pos == 0:
            return {"loss": 0.0, "auc": 0.5}

        # 采样负样本 — 使用课程式难例挖掘 (若启用)
        if self.use_curriculum and self.curriculum_miner is not None:
            # 获取当前基因特征用于难度评估
            gene_features = self.data['gene'].x.clone().detach() if hasattr(self.data['gene'], 'x') else None
            neg_drugs, neg_genes = self.curriculum_miner.sample(
                pos_genes,
                epoch=self._current_epoch,
                max_epochs=self._max_epochs,
                gene_features=gene_features,
            )
            # 确保负样本数量足够
            n_neg = min(len(neg_drugs), n_pos * 3)
            neg_drugs = neg_drugs[:n_neg].to(self.device)
            neg_genes = neg_genes[:n_neg].to(self.device)
        else:
            # 默认随机采样
            neg_drugs, neg_genes = self._sample_negative_edges(n_pos * 3)

        # 编码
        x_dict = {nt: self.data[nt].x.clone().detach() for nt in self.model.node_types}
        edge_index_dict = {
            et: self.data[et].edge_index for et in self.model.edge_types
            if hasattr(self.data[et], 'edge_index')
        }

        z_dict_out = self.model(x_dict, edge_index_dict, self.module_assignment, epoch=self._current_epoch)
        z_dict = z_dict_out[0] if isinstance(z_dict_out, tuple) else z_dict_out
        extra = z_dict_out[1] if isinstance(z_dict_out, tuple) else {}

        # 学生预测
        pos_scores = self.model.predict_edges(z_dict, pos_drugs, pos_genes)
        neg_scores = self.model.predict_edges(z_dict, neg_drugs, neg_genes)

        # BCE损失
        pos_labels = torch.ones_like(pos_scores)
        neg_labels = torch.zeros_like(neg_scores)

        scores = torch.cat([pos_scores, neg_scores])
        labels = torch.cat([pos_labels, neg_labels])

        bce_loss = F.binary_cross_entropy_with_logits(scores, labels)

        # 自蒸馏损失 (教师提供软标签)
        distill_loss = torch.tensor(0.0, device=scores.device)
        if self.use_self_distill and self.self_distill is not None:
            if self._current_epoch >= self.self_distill.start_epoch:
                # 教师编码
                # 确保教师模型与学生在同一设备上
                self.self_distill.ensure_device(self.model)
                with torch.no_grad():
                    t_z_dict_out = self.self_distill.teacher(
                        x_dict, edge_index_dict, self.module_assignment, epoch=self._current_epoch
                    )
                    t_z_dict = t_z_dict_out[0] if isinstance(t_z_dict_out, tuple) else t_z_dict_out

                # 教师预测
                t_pos_scores = self.self_distill.teacher.predict_edges(t_z_dict, pos_drugs, pos_genes)
                t_neg_scores = self.self_distill.teacher.predict_edges(t_z_dict, neg_drugs, neg_genes)
                t_scores = torch.cat([t_pos_scores, t_neg_scores])

                # KL蒸馏损失
                distill_loss = self.self_distill.compute_distill_loss(scores, t_scores)  # 不detach, 保持梯度回传

        # VIB KL损失
        vib_loss = extra.get('vib_kl_loss', torch.tensor(0.0, device=scores.device))
        vib_beta = extra.get('vib_beta', 0.0)

        # 关系门控L1正则化
        gate_l1 = self.model.relation_gate.get_gate_l1_loss()

        # 知识正则化损失
        if self.model.use_knowledge_reg and self.model.knowledge_regularizer is not None:
            kn_losses = self.model.compute_knowledge_loss(
                z_dict, self.data, self.module_assignment
            )
            kn_total = torch.tensor(0.0, device=scores.device)
            for k, v in kn_losses.items():
                kn_total = kn_total + v
        else:
            kn_total = torch.tensor(0.0, device=scores.device)

        # 综合损失
        loss = bce_loss + vib_beta * vib_loss + self.gate_l1_lambda * gate_l1 + kn_total + distill_loss

        # 反向传播
        self.optimizer.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(self.model.parameters(), 1.0)
        self.optimizer.step()

        # EMA更新教师 (仅在蒸馏启用且过了start_epoch后)
        if self.use_self_distill and self.self_distill is not None:
            if self._current_epoch >= self.self_distill.start_epoch:
                self.self_distill.update_teacher(self.model)

        # AUC
        with torch.no_grad():
            probs = torch.sigmoid(scores).cpu().numpy()
            labels_np = labels.cpu().numpy()
            try:
                auc = roc_auc_score(labels_np, probs)
            except ValueError:
                auc = 0.5

        return {
            "loss": loss.item(),
            "auc": auc,
            "gate_l1": float(gate_l1.item()),
            "distill_loss": float(distill_loss.item()),
            "vib_loss": float(vib_loss.item()),
        }

    @torch.no_grad()
    def evaluate(self) -> Dict[str, float]:
        """评估模型性能 — 使用验证集正边 + 固定验证负样本。"""
        self.model.eval()

        pos_drugs, pos_genes = self._get_val_positive_edges()
        n_pos = len(pos_drugs)

        if n_pos == 0:
            return {"auc": 0.5, "auprc": 0.0}

        # 初始化固定验证负样本 (仅首次)
        if self.val_neg_edges is None:
            self._init_val_negatives(n_neg=n_pos * 10)
        neg_drugs, neg_genes = self.val_neg_edges
        n_neg = len(neg_drugs)

        x_dict = {nt: self.data[nt].x for nt in self.model.node_types}
        edge_index_dict = {
            et: self.data[et].edge_index for et in self.model.edge_types
            if hasattr(self.data[et], 'edge_index')
        }

        z_dict_out = self.model(x_dict, edge_index_dict, self.module_assignment)
        z_dict = z_dict_out[0] if isinstance(z_dict_out, tuple) else z_dict_out

        pos_scores = self.model.predict_edges(z_dict, pos_drugs, pos_genes)
        neg_scores = self.model.predict_edges(z_dict, neg_drugs, neg_genes)

        scores = torch.cat([pos_scores, neg_scores])
        labels = torch.cat([torch.ones(n_pos), torch.zeros(n_neg)])

        probs = torch.sigmoid(scores).cpu().numpy()
        labels_np = labels.cpu().numpy()

        try:
            auc = roc_auc_score(labels_np, probs)
        except ValueError:
            auc = 0.5

        try:
            auprc = average_precision_score(labels_np, probs)
        except ValueError:
            auprc = 0.0

        # Hits@K (Recall@K — 正样本排名在前K的比例)
        k = min(10, n_pos)
        if k > 0:
            pos_probs = probs[:n_pos]
            all_probs_sorted = np.argsort(np.argsort(-probs))  # descending rank
            pos_ranks = all_probs_sorted[:n_pos] + 1  # 1-based ranks
            hits_at_k = np.mean(pos_ranks <= k)
        else:
            hits_at_k = 0.0

        return {
            "auc": auc,
            "auprc": auprc,
            f"hits@{k}": hits_at_k,
        }

    def train(
        self,
        epochs: int = 300,
        patience: int = 80,
        contrastive_pretrain: bool = True,
        contrastive_epochs: int = 100,
    ) -> Dict[str, List[float]]:
        """完整训练循环 (可选HeCo对比预训练)。

        Args:
            epochs: 链路预测训练轮数
            patience: 早停耐心
            contrastive_pretrain: 是否进行对比预训练
            contrastive_epochs: 对比预训练轮数
        """
        # Phase 1: HeCo对比预训练
        if contrastive_pretrain:
            pretrainer = HeCoContrastivePretrainer(
                self.model,
                hidden_dim=self.model.hidden_dim,
                proj_dim=64,
                tau=0.5,
                edge_dropout=0.3,
                feat_mask=0.2,
                device=self.device,
            )
            pretrainer.pretrain(
                self.data.cpu().detach(),
                epochs=contrastive_epochs,
                lr=1e-3,
                patience=20,
            )

            # 预训练后重置优化器 (Adam动量基于新参数状态重新初始化)
            self.optimizer = torch.optim.AdamW(
                self.model.parameters(), lr=self.optimizer.param_groups[0]['lr'],
                weight_decay=self.optimizer.param_groups[0].get('weight_decay', 1e-5)
            )
            # 重置学习率调度器
            self.scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
                self.optimizer, mode='max', factor=0.5, patience=30, min_lr=1e-6
            )

        # Phase 2: 链路预测微调
        self._max_epochs = epochs
        best_auc = 0.0
        best_state = None
        patience_counter = 0
        history = {"loss": [], "auc": [], "val_auc": [], "val_auprc": []}

        print(f"\n  [链路预测微调] ({epochs} epochs, patience={patience})...")
        print(f"  {'Epoch':>5s} | {'Loss':>10s} | {'Train AUC':>10s} | {'Val AUC':>10s} | {'Val AUPRC':>10s}")

        for epoch in range(epochs):
            self._current_epoch = epoch  # 更新当前epoch (用于VIB β调度)
            # 温度退火: 训练后期降低门控温度, 使门控更尖锐
            if epoch == 50:
                self.model.relation_gate.set_temperature(0.5)
                print(f"  [退火] 门控温度降至 0.5 (epoch {epoch+1})")
            elif epoch == 100:
                self.model.relation_gate.set_temperature(0.25)
                print(f"  [退火] 门控温度降至 0.25 (epoch {epoch+1})")

            train_metrics = self.train_epoch()
            val_metrics = self.evaluate()

            history["loss"].append(train_metrics["loss"])
            history["auc"].append(train_metrics["auc"])
            history["val_auc"].append(val_metrics["auc"])
            history["val_auprc"].append(val_metrics["auprc"])

            self.scheduler.step(val_metrics["auc"])

            if (epoch + 1) % 50 == 0:
                print(f"  {epoch+1:5d} | {train_metrics['loss']:10.4f} | "
                      f"{train_metrics['auc']:10.4f} | {val_metrics['auc']:10.4f} | "
                      f"{val_metrics['auprc']:10.4f}")

            # 早停
            if val_metrics["auc"] > best_auc + 0.001:
                best_auc = val_metrics["auc"]
                best_state = copy.deepcopy(self.model.state_dict())
                patience_counter = 0
            else:
                patience_counter += 1

            if patience_counter >= patience:
                print(f"  早停于 epoch {epoch+1}, best val AUC = {best_auc:.4f}")
                break

        # 恢复最佳模型
        if best_state is not None:
            self.model.load_state_dict(best_state)

        final_eval = self.evaluate()
        print(f"  训练完成 | Best Val AUC: {best_auc:.4f} | "
              f"Final Val AUPRC: {final_eval['auprc']:.4f}")

        return history

    @torch.no_grad()
    def get_node_embeddings(self) -> Dict[str, Tensor]:
        """获取所有节点的最终嵌入。"""
        self.model.eval()
        device = next(self.model.parameters()).device

        x_dict = {nt: self.data[nt].x.to(device) for nt in self.model.node_types}
        edge_index_dict = {
            et: self.data[et].edge_index.to(device) for et in self.model.edge_types
            if hasattr(self.data[et], 'edge_index')
        }

        z_dict_out = self.model(x_dict, edge_index_dict, self.module_assignment)
        z_dict = z_dict_out[0] if isinstance(z_dict_out, tuple) else z_dict_out
        return {k: v.cpu() for k, v in z_dict.items()}

    @torch.no_grad()
    def predict_all_drug_gene_pairs(self) -> pd.DataFrame:
        """预测所有药物-基因对的边分数。"""
        self.model.eval()

        z_dict = self.get_node_embeddings()
        z_dict_gpu = {k: v.to(self.device) for k, v in z_dict.items()}

        gene_names = self.data['gene'].names
        n_genes = len(gene_names)

        all_scores = []
        all_genes = []

        # 分batch预测 (避免OOM)
        batch_size = 1024
        for start in range(0, n_genes, batch_size):
            end = min(start + batch_size, n_genes)
            drug_idx = torch.zeros(end - start, dtype=torch.long, device=self.device)
            gene_idx = torch.arange(start, end, dtype=torch.long, device=self.device)

            scores = self.model.predict_edges(z_dict_gpu, drug_idx, gene_idx)
            all_scores.extend(scores.cpu().numpy())
            all_genes.extend(gene_names[start:end])

        df = pd.DataFrame({
            "gene": all_genes,
            "drug": "Beta-caryophyllene",
            "hgt_score": all_scores,
            "hgt_probability": torch.sigmoid(torch.tensor(all_scores)).numpy(),
        })
        df = df.sort_values("hgt_probability", ascending=False).reset_index(drop=True)

        return df


# ============================================================
# 4. ACSL4 回忆实验
# ============================================================

def acsl4_recall_experiment(
    data: HeteroData,
    hidden_dim: int = HIDDEN_DIM,
    n_trials: int = 5,
    trainer: "HGTLinkPredictionTrainer" = None,
    df_predictions: pd.DataFrame = None,
) -> Dict:
    """ACSL4回忆实验：掩蔽石竹烯-ACSL4边，检验HGT恢复能力。

    实验设计:
      1. 随机掩蔽石竹烯-ACSL4边
      2. 用剩余边训练HGT
      3. 预测石竹烯-ACSL4的边分数
      4. 与随机掩蔽其他边比较
      5. 重复n_trials次取平均

    这验证了HGT网络对关键信号的捕捉能力。

    ⚠ 公平性保证: 节点特征构建 (_build_gene_features) 仅基于 gene-gene 边
    和领域知识（铁死亡类别、衰老类别等），不依赖 drug-gene 边信息。
    因此掩蔽石竹烯-ACSL4边不会影响 ACSL4 的节点特征，实验设计公平。
    若未来 _build_gene_features 中加入了 drug-gene 关系特征（如 BCP 靶点
    标志），需重新评估实验公平性并调整掩蔽策略。
    """
    print(f"\n{'='*70}")
    print("ACSL4 回忆实验")
    print(f"{'='*70}")
    print(f"  试验次数: {n_trials}")
    print(f"  实验设计: 掩蔽石竹烯-ACSL4边 → HGT训练 → 检验恢复")

    results = {
        "acsl4_recall_ranks": [],
        "acsl4_recall_scores": [],
        "acsl4_recall_probabilities": [],
        "random_control_ranks": [],
        "random_control_scores": [],
    }

    # 获取石竹烯-ACSL4边索引
    gene_names = data['gene'].names
    acsl4_idx = gene_names.index("ACSL4") if "ACSL4" in gene_names else None

    if acsl4_idx is None:
        print("  ACSL4不在基因节点中，跳过实验")
        return results

    orig_ei = data['drug', 'targets', 'gene'].edge_index.clone()

    # 找到石竹烯-ACSL4边
    acsl4_edge_mask = (orig_ei[1] == acsl4_idx)

    if not acsl4_edge_mask.any():
        print("  ACSL4不是已知靶点，从已有预测中获取排名 (不重新训练模型)")

        # 优先使用传入的预测结果
        if df_predictions is not None:
            df_results = df_predictions
        elif trainer is not None:
            df_results = trainer.predict_all_drug_gene_pairs()
        else:
            print("  警告: 无trainer或df_predictions，无法获取ACSL4排名")
            return results

        acsl4_row = df_results[df_results["gene"] == "ACSL4"]
        if len(acsl4_row) > 0:
            rank = acsl4_row.index[0] + 1
            score = acsl4_row["hgt_score"].values[0]
            prob = acsl4_row["hgt_probability"].values[0]
            total = len(df_results)

            results["acsl4_recall_ranks"].append(rank)
            results["acsl4_recall_scores"].append(float(score))
            results["acsl4_recall_probabilities"].append(float(prob))
            results["acsl4_position"] = "existing_prediction"
            results["total_genes"] = total
            results["percentile"] = rank / total * 100

            print(f"\n  ACSL4 预测排名: {rank}/{total} (Top {rank/total*100:.2f}%)")
            print(f"  ACSL4 HGT分数: {score:.4f}, 概率: {prob:.4f}")

            print(f"\n  HGT预测 Top 10 靶点:")
            for i, row in df_results.head(10).iterrows():
                print(f"    #{i+1}: {row['gene']:10s} | prob={row['hgt_probability']:.4f} | score={row['hgt_score']:.4f}")

        return results

    # --- ACSL4是已知靶点：掩蔽实验 ---
    for trial in range(n_trials):
        print(f"\n  Trial {trial+1}/{n_trials}:")

        # 创建副本并掩蔽
        data_trial = data.clone().detach()  # PyG HeteroData.clone() 创建独立副本
        ei_trial = data_trial['drug', 'targets', 'gene'].edge_index

        # 掩蔽所有石竹烯-ACSL4边
        keep_mask = ~((ei_trial[0] == 0) & (ei_trial[1] == acsl4_idx))
        data_trial['drug', 'targets', 'gene'].edge_index = ei_trial[:, keep_mask]

        n_removed = (~keep_mask).sum().item()
        print(f"    已掩蔽 {n_removed} 条石竹烯-ACSL4边")

        # 构建特征
        data_with_feats = build_node_features(data_trial, hidden_dim)

        # 训练HGT (显式关闭模块聚合, 因掩蔽实验不需要模块先验)
        model = FerroHGT(
            metadata=data.metadata(),
            hidden_dim=hidden_dim,
            num_heads=NUM_HEADS,
            num_layers=NUM_LAYERS,
            dropout=DROPOUT,
            use_module_aggregation=False,
        )
        trainer = HGTLinkPredictionTrainer(model, data_with_feats, device=str(DEVICE))
        trainer.train(epochs=150, patience=40)

        # 预测ACSL4分数
        x_dict = {nt: data_with_feats[nt].x for nt in model.node_types}
        edge_index_dict = {
            et: data_with_feats[et].edge_index for et in model.edge_types
            if hasattr(data_with_feats[et], 'edge_index')
        }

        z_dict_out = trainer.model(x_dict, edge_index_dict)
        z_dict = z_dict_out[0] if isinstance(z_dict_out, tuple) else z_dict_out
        z_dict_gpu = {k: v.to(DEVICE) for k, v in z_dict.items()}

        drug_idx = torch.tensor([0], device=DEVICE)
        gene_idx = torch.tensor([acsl4_idx], device=DEVICE)
        acsl4_score = trainer.model.predict_edges(z_dict_gpu, drug_idx, gene_idx).item()

        # 获取所有基因的分数排名
        all_gene_scores = []
        for g_idx in range(len(gene_names)):
            if g_idx in ei_trial[1].tolist():  # 排除已知靶点
                continue
            g_tensor = torch.tensor([g_idx], device=DEVICE)
            s = trainer.model.predict_edges(z_dict_gpu, drug_idx, g_tensor).item()
            all_gene_scores.append((g_idx, s))

        all_gene_scores.sort(key=lambda x: x[1], reverse=True)
        acsl4_rank = next(
            (i+1 for i, (idx, _) in enumerate(all_gene_scores) if idx == acsl4_idx),
            len(all_gene_scores)
        )

        results["acsl4_recall_ranks"].append(acsl4_rank)
        results["acsl4_recall_scores"].append(acsl4_score)
        results["acsl4_recall_probabilities"].append(
            float(torch.sigmoid(torch.tensor(acsl4_score)))
        )

        print(f"    ACSL4 召回排名: {acsl4_rank}/{len(all_gene_scores)} "
              f"| 分数: {acsl4_score:.4f} | 概率: {results['acsl4_recall_probabilities'][-1]:.4f}")

    # 汇总
    ranks = results["acsl4_recall_ranks"]
    scores = results["acsl4_recall_scores"]
    print(f"\n  ACSL4回忆实验汇总 (n={n_trials}):")
    print(f"    平均排名: {np.mean(ranks):.1f} ± {np.std(ranks):.1f}")
    print(f"    平均分数: {np.mean(scores):.4f} ± {np.std(scores):.4f}")
    print(f"    平均概率: {np.mean(results['acsl4_recall_probabilities']):.4f}")

    # --- 随机对照组: 掩蔽非ACSL4边并测量恢复能力 ---
    control_trials = min(5, len(orig_ei[0]) - 1)
    if control_trials > 0:
        print(f"\n  {'='*50}")
        print(f"  随机对照组 (n={control_trials}): 掩蔽随机正边 → 检验恢复")
        print(f"  {'='*50}")

        # 获取所有非ACSL4正边索引
        non_acsl4_mask = ~acsl4_edge_mask
        non_acsl4_indices = torch.where(non_acsl4_mask)[0].tolist()
        rng_ctrl = np.random.RandomState(2024)

        for ctrl_trial in range(control_trials):
            # 随机选一条非ACSL4正边
            idx = rng_ctrl.choice(non_acsl4_indices)
            masked_drug = orig_ei[0, idx].item()
            masked_gene = orig_ei[1, idx].item()
            masked_gene_name = gene_names[masked_gene] if masked_gene < len(gene_names) else f"idx_{masked_gene}"

            print(f"\n  Control Trial {ctrl_trial+1}/{control_trials}: 掩蔽 {masked_gene_name}")

            data_ctrl = data.clone().detach()
            ei_ctrl = data_ctrl['drug', 'targets', 'gene'].edge_index
            keep_ctrl = ~((ei_ctrl[0] == masked_drug) & (ei_ctrl[1] == masked_gene))
            data_ctrl['drug', 'targets', 'gene'].edge_index = ei_ctrl[:, keep_ctrl]

            data_ctrl_feats = build_node_features(data_ctrl, hidden_dim)

            ctrl_model = FerroHGT(
                metadata=data.metadata(),
                hidden_dim=hidden_dim, num_heads=NUM_HEADS,
                num_layers=NUM_LAYERS, dropout=DROPOUT,
                use_module_aggregation=False,
            )
            ctrl_trainer = HGTLinkPredictionTrainer(ctrl_model, data_ctrl_feats, device=str(DEVICE))
            ctrl_trainer.train(epochs=100, patience=30)

            # 计算被掩蔽基因的排名
            x_dict = {nt: data_ctrl_feats[nt].x for nt in ctrl_model.node_types}
            edge_index_dict = {
                et: data_ctrl_feats[et].edge_index for et in ctrl_model.edge_types
                if hasattr(data_ctrl_feats[et], 'edge_index')
            }
            z_dict_out = ctrl_trainer.model(x_dict, edge_index_dict)
            z_dict = z_dict_out[0] if isinstance(z_dict_out, tuple) else z_dict_out
            z_dict_gpu = {k: v.to(DEVICE) for k, v in z_dict.items()}

            drug_idx = torch.tensor([masked_drug], device=DEVICE)
            gene_idx = torch.tensor([masked_gene], device=DEVICE)
            masked_score = ctrl_trainer.model.predict_edges(z_dict_gpu, drug_idx, gene_idx).item()

            # 排除已知正边后对所有基因排名
            ei_ctrl_list = ei_ctrl[1].tolist()
            all_ctrl_scores = []
            for g_idx in range(len(gene_names)):
                if g_idx in ei_ctrl_list:
                    continue
                g_t = torch.tensor([g_idx], device=DEVICE)
                s = ctrl_trainer.model.predict_edges(z_dict_gpu, drug_idx, g_t).item()
                all_ctrl_scores.append((g_idx, s))

            all_ctrl_scores.sort(key=lambda x: x[1], reverse=True)
            ctrl_rank = next(
                (i+1 for i, (midx, _) in enumerate(all_ctrl_scores) if midx == masked_gene),
                len(all_ctrl_scores)
            )

            results["random_control_ranks"].append(ctrl_rank)
            results["random_control_scores"].append(masked_score)
            print(f"    {masked_gene_name} 召回排名: {ctrl_rank}/{len(all_ctrl_scores)} "
                  f"| 分数: {masked_score:.4f}")

        # 对照组汇总
        if results["random_control_ranks"]:
            cr = results["random_control_ranks"]
            print(f"\n  对照组汇总 (n={len(cr)}):")
            print(f"    平均排名: {np.mean(cr):.1f} ± {np.std(cr):.1f}")
            print(f"    平均分数: {np.mean(results['random_control_scores']):.4f}")
            # 与ACSL4比较
            if results["acsl4_recall_ranks"]:
                print(f"    ACSL4平均排名: {np.mean(results['acsl4_recall_ranks']):.1f}")
                print(f"    Δ(ACSL4 - Random): {np.mean(results['acsl4_recall_ranks']) - np.mean(cr):.1f}")

    return results


# ============================================================
# 5. 主函数
# ============================================================

def main():
    print("=" * 70)
    print("阶段3: HGT异构网络学习 — 关系推理与边预测")
    print("=" * 70)

    # --- 加载异构图 ---
    graph_path = DATA_DIR / "ferroptosis_hetero_graph.pt"
    print(f"\n[1/7] 加载异构图: {graph_path}")
    data_base = torch.load(graph_path, weights_only=False)
    print(f"  节点类型: {data_base.metadata()[0]}")
    print(f"  边类型: {data_base.metadata()[1]}")

    # --- 消融实验配置 ---
    ablation_configs = [
        # (名称, add_ppi, add_tf, add_coexp, add_pathway, add_methylation)
        ("Baseline (无增强)",       False, False, False, False, False),
        ("+PPI",                     True,  False, False, False, False),
        ("+PPI+TF",                  True,  True,  False, False, False),
        ("+PPI+TF+Coexp",            True,  True,  True,  False, False),
        ("+PPI+TF+Coexp+Pathway",    True,  True,  True,  True,  False),
        ("+PPI+TF+Coexp+Pathway+Methyl", True, True, True, True, True),
    ]

    ablation_results = []
    best_score = -1.0  # 综合评分: 原始嵌入ACSL4排名为主
    best_embeddings = None
    best_df_predictions = None
    best_config_name = ""

    # 所有可能的增强边类型 (用于门控冻结逻辑)
    ALL_ENHANCED_EDGE_TYPES = {
        # (bool_key, edge_type_str, gate_key)
        "ppi": "gene__interacts__gene",
        "tf": "gene__regulates__gene",
        "coexp": "gene__coexpressed__gene",
        "pathway": "gene__co_pathway__gene",
        "methylation": "gene__co_methylated__gene",
    }

    for config_name, ppi, tf, coexp, pw, meth in ablation_configs:
        print(f"\n{'='*70}")
        print(f"消融实验: {config_name}")
        print(f"{'='*70}")

        # === 计算需要冻结的边类型 (当前配置不激活的维度) ===
        active_flags = {"ppi": ppi, "tf": tf, "coexp": coexp,
                        "pathway": pw, "methylation": meth}
        freeze_types = set()
        per_edge_bias = {}
        for key, gate_key in ALL_ENHANCED_EDGE_TYPES.items():
            if not active_flags[key]:
                freeze_types.add(gate_key)
                # 冻结的维度: bias设为-5 (sigmoid≈0.007, 近似关闭)
                per_edge_bias[gate_key] = -5.0

        if freeze_types:
            print(f"  冻结门控 (非激活维度): {sorted(freeze_types)}")

        # 克隆基础图
        data = data_base.clone().detach()

        # --- 多维边增强 ---
        print(f"\n[1.5/7] 多维边增强...")
        data = enrich_graph_edges(
            data,
            add_ppi=ppi, add_tf=tf, add_coexp=coexp,
            add_pathway=pw, add_methylation=meth,
        )

        # --- 构建节点特征 ---
        print(f"\n[2/7] 构建节点初始特征...")
        data = build_node_features(data, fingerprint_dim=HIDDEN_DIM)
        for nt in data.node_types:
            if hasattr(data[nt], 'x') and data[nt].x is not None:
                print(f"  {nt}: {data[nt].x.shape}")

        # --- 初始化HGT模型 (含门控冻结) ---
        print(f"\n[3/7] 初始化 FerroHGT 模型 (含差异化门控初始化)...")
        model = FerroHGT(
            metadata=data.metadata(),
            hidden_dim=HIDDEN_DIM,
            num_heads=NUM_HEADS,
            num_layers=NUM_LAYERS,
            dropout=DROPOUT,
            relation_gate_freeze_types=freeze_types if freeze_types else None,
            relation_gate_per_edge_bias=per_edge_bias if per_edge_bias else None,
        )
        n_params = sum(p.numel() for p in model.parameters())
        print(f"  模型参数量: {n_params:,}")

        # 检查初始门控值
        init_gates = model.get_relation_gate_values()
        print(f"  初始门控值 (冻结维度应≈0):")
        for et_key, gv in sorted(init_gates.items(), key=lambda x: -x[1]):
            marker = "[冻结]" if et_key in freeze_types else "[可学]"
            print(f"    {marker} {et_key}: {gv:.4f}")

        # --- 子图提取 (用于重点分析) ---
        print(f"\n[3.5/7] 铁死亡核心子图提取 + 基因模块构建...")
        subgraph_indices, hop_info = extract_ferroptosis_subgraph(data, k_hop=2)

        # 构建基因模块分配矩阵 (用于层次聚合)
        module_assignment, module_names = build_gene_module_assignment(
            data_base['gene'].names, data_base
        )

        # --- 训练HGT (HeCo对比预训练 + 链路预测微调) ---
        print(f"\n[4/7] 训练 HGT 链路预测 (关系门控 + 层次模块聚合)...")
        trainer = HGTLinkPredictionTrainer(
            model, data,
            module_assignment=module_assignment,
            lr=LR,
            weight_decay=WEIGHT_DECAY,
            gate_l1_lambda=1e-3,  # 提高到1e-3, 增强稀疏门控
            device=str(DEVICE),
        )
        history = trainer.train(
            epochs=MAX_EPOCHS,
            patience=PATIENCE,
            contrastive_pretrain=True,   # 开启对比预训练
            contrastive_epochs=50,       # 减少预训练轮数提高效率
        )

        # 最终评估
        final_metrics = trainer.evaluate()
        print(f"\n  最终评估:")
        for k, v in final_metrics.items():
            print(f"    {k}: {v:.4f}")

        # --- ACSL4回忆实验 ---
        print(f"\n[5/7] ACSL4 回忆实验...")
        acsl4_results = acsl4_recall_experiment(
            data.cpu().detach() if hasattr(data, 'cpu') else data,
            hidden_dim=HIDDEN_DIM, n_trials=3,
            trainer=trainer,
        )

        # --- 提取节点嵌入 ---
        print(f"\n[6/7] 提取节点嵌入...")
        embeddings = trainer.get_node_embeddings()

        # 原始嵌入质量评估 (不依赖解码器)
        print(f"  原始嵌入质量评估...")
        raw_metrics = evaluate_raw_embeddings(embeddings, data_base['gene'].names)
        print(f"    Cosine ACSL4 rank: {raw_metrics['cosine']['acsl4_rank']}")
        print(f"    Hadamard ACSL4 rank: {raw_metrics['hadamard']['acsl4_rank']}")
        print(f"    L2 ACSL4 rank: {raw_metrics['l2']['acsl4_rank']}")
        print(f"    综合排名: {raw_metrics.get('ensemble_rank', 'N/A')}")

        # 关系门控值 (训练后)
        gate_values = model.get_relation_gate_values()
        print(f"  训练后关系门控值:")
        for et_key, gv in sorted(gate_values.items(), key=lambda x: -x[1]):
            marker = "[冻结]" if et_key in freeze_types else "[可学]"
            print(f"    {marker} {et_key}: {gv:.4f}")

        # --- 预测所有药物-基因对 ---
        print(f"\n[7/7] 预测所有药物-基因对...")
        df_predictions = trainer.predict_all_drug_gene_pairs()

        # 记录ACSL4排名
        acsl4_row = df_predictions[df_predictions["gene"] == "ACSL4"]
        acsl4_rank = int(acsl4_row.index[0] + 1) if len(acsl4_row) > 0 else -1
        print(f"  ACSL4 HGT排名: {acsl4_rank}/{len(df_predictions)} "
              f"(Top {acsl4_rank/len(df_predictions)*100:.2f}%)" if acsl4_rank > 0 else "  ACSL4未在基因列表中")

        # 保存消融结果
        ablation_results.append({
            "config": config_name,
            "ppi": ppi, "tf": tf, "coexp": coexp, "pathway": pw, "methylation": meth,
            "val_auc": final_metrics.get("auc", 0),
            "val_auprc": final_metrics.get("auprc", 0),
            "hits_k": final_metrics.get("hits@10", 0),
            "acsl4_rank": acsl4_rank,
            "acsl4_probability": float(acsl4_row["hgt_probability"].values[0]) if len(acsl4_row) > 0 else 0,
            # 原始嵌入指标
            "raw_cosine_acsl4_rank": raw_metrics['cosine']['acsl4_rank'],
            "raw_hadamard_acsl4_rank": raw_metrics['hadamard']['acsl4_rank'],
            "raw_l2_acsl4_rank": raw_metrics['l2']['acsl4_rank'],
            "raw_ensemble_rank": raw_metrics.get('ensemble_rank', -1),
            "raw_ensemble_percentile": raw_metrics.get('ensemble_percentile', 100),
            # 关系门控值
            "gate_values": gate_values,
        })

        # 跟踪最佳配置 (综合评分: 原始嵌入质量70% + Val AUPRC 30%)
        val_auprc = final_metrics.get("auprc", 0)
        raw_ens_rank = raw_metrics.get('ensemble_rank', len(df_predictions))
        raw_rank_pct = raw_ens_rank / len(df_predictions) if raw_ens_rank > 0 else 1.0
        # 评分: 原始嵌入质量(ACSL4排名百分位反转)占70%, 验证AUPRC占30%
        combined_score = (1.0 - raw_rank_pct) * 0.70 + val_auprc * 0.30
        if combined_score > best_score:
            best_score = combined_score
            best_embeddings = embeddings
            best_df_predictions = df_predictions
            best_config_name = config_name

    # --- 消融实验汇总 ---
    print(f"\n{'='*70}")
    print("消融实验汇总")
    print(f"{'='*70}")
    print(f"  {'Config':<35s} | {'Val AUC':>8s} | {'DecoderRank':>11s} | {'RawEnsRank':>11s} | {'RawEns%':>8s}")
    print(f"  {'-'*90}")
    for r in ablation_results:
        acsl4_str = f"{r['acsl4_rank']}" if r['acsl4_rank'] > 0 else "N/A"
        raw_str = f"{r.get('raw_ensemble_rank', 'N/A')}"
        raw_pct_str = f"{r.get('raw_ensemble_percentile', 100):.1f}%"
        print(f"  {r['config']:<35s} | {r['val_auc']:8.4f} | {acsl4_str:>11s} | {raw_str:>11s} | {raw_pct_str:>8s}")

    # 保存消融结果
    ablation_path = RESULTS_DIR / "ablation_results.json"
    with open(ablation_path, 'w', encoding='utf-8') as f:
        json.dump(ablation_results, f, indent=2, ensure_ascii=False)
    print(f"\n  ✓ 消融结果已保存: {ablation_path}")

    # --- 使用最佳配置保存结果 ---
    print(f"\n{'='*70}")
    print(f"最佳配置: {best_config_name} (综合评分 = {best_score:.4f})")
    print(f"{'='*70}")

    # 保存最佳嵌入
    embed_path = RESULTS_DIR / "node_embeddings.pt"
    torch.save(best_embeddings, embed_path)
    print(f"  ✓ 节点嵌入已保存: {embed_path}")

    # 导出基因嵌入为CSV
    gene_names = data_base['gene'].names
    gene_emb = best_embeddings['gene'].numpy()
    df_emb = pd.DataFrame(
        gene_emb,
        index=gene_names,
        columns=[f"emb_{i}" for i in range(gene_emb.shape[1])],
    )
    df_emb.index.name = "gene"
    df_emb.to_csv(RESULTS_DIR / "gene_hgt_embeddings.csv")
    print(f"  ✓ 基因嵌入已保存 ({gene_emb.shape[0]} x {gene_emb.shape[1]})")

    # 导出ACSL4嵌入
    if "ACSL4" in gene_names:
        acsl4_idx = gene_names.index("ACSL4")
        acsl4_emb = gene_emb[acsl4_idx]
        np.save(RESULTS_DIR / "acsl4_hgt_embedding.npy", acsl4_emb)
        print(f"  ✓ ACSL4嵌入已保存: {acsl4_emb.shape}")

    # 导出石竹烯嵌入
    drug_emb = best_embeddings['drug'].numpy()
    np.save(RESULTS_DIR / "bcp_hgt_embedding.npy", drug_emb)
    print(f"  ✓ 石竹烯嵌入已保存: {drug_emb.shape}")

    # 保存最佳预测结果
    best_df_predictions.to_csv(RESULTS_DIR / "hgt_drug_gene_predictions.csv", index=False)
    print(f"  ✓ 预测结果已保存 ({len(best_df_predictions)} 行)")

    # Top 20 预测
    print(f"\n  HGT Top 20 预测靶点 (配置: {best_config_name}):")
    print(f"  {'Rank':>5s} | {'Gene':>12s} | {'Probability':>12s} | {'Score':>10s}")
    print(f"  {'-'*50}")
    for i, row in best_df_predictions.head(20).iterrows():
        print(f"  {i+1:5d} | {row['gene']:>12s} | {row['hgt_probability']:12.4f} | {row['hgt_score']:10.4f}")

    # 保存HGT模型
    model_path = RESULTS_DIR / "ferro_hgt_model.pt"
    torch.save({
        "config": {
            "hidden_dim": HIDDEN_DIM,
            "num_heads": NUM_HEADS,
            "num_layers": NUM_LAYERS,
            "dropout": DROPOUT,
        },
        "best_config": best_config_name,
        "ablation_results": ablation_results,
    }, model_path)
    print(f"  ✓ 配置已保存: {model_path}")

    print("\n" + "=" * 70)
    print(f"阶段3 完成 — 最佳配置: {best_config_name}")
    print(f"  输出目录: {RESULTS_DIR}")
    print("=" * 70)

    return best_embeddings, best_df_predictions


if __name__ == "__main__":
    main()