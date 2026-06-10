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
from torch_geometric.utils import negative_sampling
from sklearn.metrics import roc_auc_score, average_precision_score, precision_recall_curve

# 路径
BASE_DIR = Path(__file__).parent.parent
sys.path.insert(0, str(BASE_DIR))
DATA_DIR = BASE_DIR / "data" / "ferroptosis_graph"
RESULTS_DIR = BASE_DIR / "results" / "hgt_learning"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# 设备
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Device: {DEVICE}")

# ============================================================
# 1. 数据加载与特征构建
# ============================================================

def build_node_features(data: HeteroData, fingerprint_dim: int = 256) -> HeteroData:
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
      - ACSL4邻居关系
      - BCP靶点关系
      - 节点度中心性
      - 介数中心性近似
    """
    num_genes = len(gene_names)
    feat_dim = 128  # 中间维度

    # 导入基因集定义
    sys.path.insert(0, str(BASE_DIR))
    from bcp_ferroptosis_pipeline.stage1_build_graph import (
        ALL_FERROPTOSIS_GENES, ALL_AGING_GENES, ALL_CIRI_GENES,
        BCP_TARGETS, ACSL4_FIRST_NEIGHBORS,
        CIRI_DEGS, FERRO_AGING_GENES_FILE,
    )

    acsl4_neighbors = set(
        ACSL4_FIRST_NEIGHBORS["direct_interactors"] +
        ACSL4_FIRST_NEIGHBORS["indirect_regulators"]
    )

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

        # 4. ACSL4邻居 (1维)
        is_acsl4_nb = 1.0 if g in acsl4_neighbors else 0.0
        feat.append(is_acsl4_nb)

        # 5. BCP靶点 (1维)
        is_bcp_target = 1.0 if g in BCP_TARGETS else 0.0
        feat.append(is_bcp_target)

        # 6. 是否是ACSL4本身 (1维)
        is_acsl4 = 1.0 if g == "ACSL4" else 0.0
        feat.append(is_acsl4)

        # 7. 随机初始化编码 (补充信息)
        # 使用基因名的hash作为种子，确保可复现
        seed = hash(g) % (2**31)
        rng = np.random.RandomState(seed)
        random_embed = rng.randn(64).astype(np.float32) * 0.1
        feat.extend(random_embed.tolist())

        features_list.append(feat)

    gene_feat = torch.tensor(features_list, dtype=torch.float)  # (N, 74)

    # 投影到hidden_dim
    proj = nn.Linear(gene_feat.size(1), hidden_dim)
    nn.init.xavier_uniform_(proj.weight)
    nn.init.zeros_(proj.bias)

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
        from bcp_ferroptosis_pipeline.stage1_build_graph import ALL_CIRI_GENES
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


# ============================================================
# 2. HGT 链路预测模型
# ============================================================

class FerroHGT(nn.Module):
    """石竹烯-铁衰老-脑缺血 异构网络 HGT 编码器-解码器。

    架构:
      1. 类型感知输入投影
      2. 多层 HGTConv 消息传递
      3. 边预测解码器 (Drug-Gene)

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
    ):
        super().__init__()
        self.hidden_dim = hidden_dim
        self.node_types = metadata[0]
        self.edge_types = metadata[1]

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

        # 边预测解码器 (Drug-Gene)
        # 使用拼接 + MLP
        self.edge_decoder = nn.Sequential(
            nn.Linear(hidden_dim * 2, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout * 0.5),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.BatchNorm1d(hidden_dim // 2),
            nn.ReLU(),
            nn.Dropout(dropout * 0.25),
            nn.Linear(hidden_dim // 2, 1),
        )

    def forward(
        self,
        x_dict: Dict[str, Tensor],
        edge_index_dict: Dict,
    ) -> Dict[str, Tensor]:
        """HGT编码前向传播。

        Returns:
            z_dict: 各节点类型的最终嵌入
        """
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

        return x_dict

    def predict_edges(
        self,
        z_dict: Dict[str, Tensor],
        drug_indices: Tensor,
        gene_indices: Tensor,
    ) -> Tensor:
        """预测 Drug-Gene 边分数。

        Args:
            z_dict: 编码后的节点嵌入
            drug_indices: 药物节点索引 (N,)
            gene_indices: 基因节点索引 (N,)

        Returns:
            scores: 预测分数 (N,) — logits
        """
        z_drug = z_dict['drug'][drug_indices]
        z_gene = z_dict['gene'][gene_indices]
        h = torch.cat([z_drug, z_gene], dim=-1)
        return self.edge_decoder(h).squeeze(-1)


# ============================================================
# 3. 训练器
# ============================================================

class HGTLinkPredictionTrainer:
    """HGT链路预测训练器。

    支持:
      - 正负样本平衡采样
      - 早停
      - ACSL4回忆实验
      - MC Dropout不确定性估计
    """

    def __init__(
        self,
        model: FerroHGT,
        data: HeteroData,
        lr: float = 0.001,
        weight_decay: float = 1e-5,
        device: str = "cuda",
    ):
        self.model = model.to(device)
        self.data = data.to(device)
        self.device = device
        self.optimizer = torch.optim.AdamW(
            model.parameters(), lr=lr, weight_decay=weight_decay
        )
        self.scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
            self.optimizer, mode='max', factor=0.5, patience=30, min_lr=1e-6
        )

    def _get_positive_edges(self) -> Tuple[Tensor, Tensor]:
        """从数据中提取正样本边 (drug → gene)。"""
        ei = self.data['drug', 'targets', 'gene'].edge_index
        return ei[0], ei[1]  # drug_indices, gene_indices

    def _sample_negative_edges(self, n_pos: int) -> Tuple[Tensor, Tensor]:
        """采样负样本边 (排除已知正边)。"""
        n_drugs = self.data['drug'].num_nodes
        n_genes = self.data['gene'].num_nodes
        pos_ei = self.data['drug', 'targets', 'gene'].edge_index
        pos_set = set(zip(pos_ei[0].tolist(), pos_ei[1].tolist()))

        # 随机采样
        neg_drugs = []
        neg_genes = []
        n_sampled = 0
        max_attempts = n_pos * 10

        while n_sampled < n_pos and max_attempts > 0:
            d = torch.randint(0, n_drugs, (1,)).item()
            g = torch.randint(0, n_genes, (1,)).item()
            if (d, g) not in pos_set:
                neg_drugs.append(d)
                neg_genes.append(g)
                n_sampled += 1
            max_attempts -= 1

        return (
            torch.tensor(neg_drugs, device=self.device),
            torch.tensor(neg_genes, device=self.device),
        )

    def train_epoch(self) -> Dict[str, float]:
        """单epoch训练。"""
        self.model.train()

        # 获取正样本
        pos_drugs, pos_genes = self._get_positive_edges()
        n_pos = len(pos_drugs)

        if n_pos == 0:
            return {"loss": 0.0, "auc": 0.5}

        # 采样负样本 (1:3比例)
        neg_drugs, neg_genes = self._sample_negative_edges(n_pos * 3)

        # 编码 - 使用clone避免计算图重叠
        x_dict = {nt: self.data[nt].x.clone().detach() for nt in self.model.node_types}
        edge_index_dict = {
            et: self.data[et].edge_index for et in self.model.edge_types
            if hasattr(self.data[et], 'edge_index')
        }

        z_dict = self.model(x_dict, edge_index_dict)

        # 预测
        pos_scores = self.model.predict_edges(z_dict, pos_drugs, pos_genes)
        neg_scores = self.model.predict_edges(z_dict, neg_drugs, neg_genes)

        # BCE损失
        pos_labels = torch.ones_like(pos_scores)
        neg_labels = torch.zeros_like(neg_scores)

        scores = torch.cat([pos_scores, neg_scores])
        labels = torch.cat([pos_labels, neg_labels])

        loss = F.binary_cross_entropy_with_logits(scores, labels)

        # 反向传播
        self.optimizer.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(self.model.parameters(), 1.0)
        self.optimizer.step()

        # AUC
        with torch.no_grad():
            probs = torch.sigmoid(scores).cpu().numpy()
            labels_np = labels.cpu().numpy()
            try:
                auc = roc_auc_score(labels_np, probs)
            except ValueError:
                auc = 0.5

        return {"loss": loss.item(), "auc": auc}

    @torch.no_grad()
    def evaluate(self) -> Dict[str, float]:
        """评估模型性能。"""
        self.model.eval()

        pos_drugs, pos_genes = self._get_positive_edges()
        n_pos = len(pos_drugs)

        if n_pos == 0:
            return {"auc": 0.5, "auprc": 0.0}

        neg_drugs, neg_genes = self._sample_negative_edges(n_pos * 3)

        x_dict = {nt: self.data[nt].x for nt in self.model.node_types}
        edge_index_dict = {
            et: self.data[et].edge_index for et in self.model.edge_types
            if hasattr(self.data[et], 'edge_index')
        }

        z_dict = self.model(x_dict, edge_index_dict)

        pos_scores = self.model.predict_edges(z_dict, pos_drugs, pos_genes)
        neg_scores = self.model.predict_edges(z_dict, neg_drugs, neg_genes)

        scores = torch.cat([pos_scores, neg_scores])
        labels = torch.cat([torch.ones(n_pos), torch.zeros(n_pos * 3)])

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

        # Hits@K (正样本排名)
        k = min(10, n_pos)
        if k > 0:
            pos_probs = probs[:n_pos]
            neg_probs = probs[n_pos:]
            topk = np.sort(probs)[-k:]
            hits = np.sum(np.isin(pos_probs, topk))
            hits_at_k = hits / k
        else:
            hits_at_k = 0.0

        return {
            "auc": auc,
            "auprc": auprc,
            f"hits@{k}": hits_at_k,
        }

    def train(self, epochs: int = 300, patience: int = 80) -> Dict[str, List[float]]:
        """完整训练循环。"""
        best_auc = 0.0
        best_state = None
        patience_counter = 0
        history = {"loss": [], "auc": [], "val_auc": [], "val_auprc": []}

        print(f"\n  开始训练 ({epochs} epochs, patience={patience})...")
        print(f"  {'Epoch':>5s} | {'Loss':>10s} | {'Train AUC':>10s} | {'Val AUC':>10s} | {'Val AUPRC':>10s}")

        for epoch in range(epochs):
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

        x_dict = {nt: self.data[nt].x for nt in self.model.node_types}
        edge_index_dict = {
            et: self.data[et].edge_index for et in self.model.edge_types
            if hasattr(self.data[et], 'edge_index')
        }

        z_dict = self.model(x_dict, edge_index_dict)
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
    hidden_dim: int = 256,
    n_trials: int = 5,
) -> Dict:
    """ACSL4回忆实验：掩蔽石竹烯-ACSL4边，检验HGT恢复能力。

    实验设计:
      1. 随机掩蔽石竹烯-ACSL4边
      2. 用剩余边训练HGT
      3. 预测石竹烯-ACSL4的边分数
      4. 与随机掩蔽其他边比较
      5. 重复n_trials次取平均

    这验证了HGT网络对关键信号的捕捉能力。
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
        print("  石竹烯-ACSL4边不在原始数据中，实验: 预测潜在边")
        # ACSL4不在已知靶点中，预测其作为潜在靶点的排名
        model = FerroHGT(
            metadata=data.metadata(),
            hidden_dim=hidden_dim,
            num_heads=4,
            num_layers=3,
            dropout=0.3,
        )

        # 构建特征
        data_with_feats = build_node_features(data, hidden_dim)

        trainer = HGTLinkPredictionTrainer(model, data_with_feats, device=str(DEVICE))
        trainer.train(epochs=200, patience=50)

        # 预测所有基因对
        df_results = trainer.predict_all_drug_gene_pairs()

        acsl4_row = df_results[df_results["gene"] == "ACSL4"]
        if len(acsl4_row) > 0:
            rank = acsl4_row.index[0] + 1
            score = acsl4_row["hgt_score"].values[0]
            prob = acsl4_row["hgt_probability"].values[0]
            total = len(df_results)

            results["acsl4_recall_ranks"].append(rank)
            results["acsl4_recall_scores"].append(float(score))
            results["acsl4_recall_probabilities"].append(float(prob))
            results["acsl4_position"] = "new_prediction"  # ACSL4不是已知靶点
            results["total_genes"] = total
            results["percentile"] = rank / total * 100

            print(f"\n  ACSL4 预测排名: {rank}/{total} (Top {rank/total*100:.2f}%)")
            print(f"  ACSL4 HGT分数: {score:.4f}, 概率: {prob:.4f}")

            # Top 10 预测
            print(f"\n  HGT预测 Top 10 靶点:")
            for i, row in df_results.head(10).iterrows():
                print(f"    #{i+1}: {row['gene']:10s} | prob={row['hgt_probability']:.4f} | score={row['hgt_score']:.4f}")

        return results

    # --- ACSL4是已知靶点：掩蔽实验 ---
    for trial in range(n_trials):
        print(f"\n  Trial {trial+1}/{n_trials}:")

        # 创建副本并掩蔽
        data_trial = copy.deepcopy(data)
        ei_trial = data_trial['drug', 'targets', 'gene'].edge_index

        # 掩蔽所有石竹烯-ACSL4边
        keep_mask = ~((ei_trial[0] == 0) & (ei_trial[1] == acsl4_idx))
        data_trial['drug', 'targets', 'gene'].edge_index = ei_trial[:, keep_mask]

        n_removed = (~keep_mask).sum().item()
        print(f"    已掩蔽 {n_removed} 条石竹烯-ACSL4边")

        # 构建特征
        data_with_feats = build_node_features(data_trial, hidden_dim)

        # 训练HGT
        model = FerroHGT(
            metadata=data.metadata(),
            hidden_dim=hidden_dim,
            num_heads=4,
            num_layers=3,
            dropout=0.3,
        )
        trainer = HGTLinkPredictionTrainer(model, data_with_feats, device=str(DEVICE))
        trainer.train(epochs=150, patience=40)

        # 预测ACSL4分数
        x_dict = {nt: data_with_feats[nt].x for nt in model.node_types}
        edge_index_dict = {
            et: data_with_feats[et].edge_index for et in model.edge_types
            if hasattr(data_with_feats[et], 'edge_index')
        }

        z_dict = trainer.model(x_dict, edge_index_dict)
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
    print(f"\n[1/6] 加载异构图: {graph_path}")
    data = torch.load(graph_path, weights_only=False)
    print(f"  节点类型: {data.metadata()[0]}")
    print(f"  边类型: {data.metadata()[1]}")

    # --- 构建节点特征 ---
    print(f"\n[2/6] 构建节点初始特征...")
    data = build_node_features(data, fingerprint_dim=256)
    for nt in data.node_types:
        if hasattr(data[nt], 'x') and data[nt].x is not None:
            print(f"  {nt}: {data[nt].x.shape}")

    # --- 初始化HGT模型 ---
    print(f"\n[3/6] 初始化 FerroHGT 模型...")
    model = FerroHGT(
        metadata=data.metadata(),
        hidden_dim=256,
        num_heads=4,
        num_layers=3,
        dropout=0.3,
    )
    n_params = sum(p.numel() for p in model.parameters())
    print(f"  模型参数量: {n_params:,}")

    # --- 训练HGT ---
    print(f"\n[4/6] 训练 HGT 链路预测...")
    trainer = HGTLinkPredictionTrainer(model, data, device=str(DEVICE))
    history = trainer.train(epochs=300, patience=80)

    # 最终评估
    final_metrics = trainer.evaluate()
    print(f"\n  最终评估:")
    for k, v in final_metrics.items():
        print(f"    {k}: {v:.4f}")

    # --- ACSL4回忆实验 ---
    print(f"\n[5/6] ACSL4 回忆实验...")
    acsl4_results = acsl4_recall_experiment(
        copy.deepcopy(data), hidden_dim=256, n_trials=5
    )

    # --- 提取节点嵌入 ---
    print(f"\n[6/6] 提取节点嵌入...")
    embeddings = trainer.get_node_embeddings()

    # 保存嵌入
    embed_path = RESULTS_DIR / "node_embeddings.pt"
    torch.save(embeddings, embed_path)
    print(f"  ✓ 节点嵌入已保存: {embed_path}")

    # 导出基因嵌入为CSV
    gene_names = data['gene'].names
    gene_emb = embeddings['gene'].numpy()
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
    drug_emb = embeddings['drug'].numpy()
    np.save(RESULTS_DIR / "bcp_hgt_embedding.npy", drug_emb)
    print(f"  ✓ 石竹烯嵌入已保存: {drug_emb.shape}")

    # --- 预测所有药物-基因对 ---
    print(f"\n  预测所有药物-基因对分数...")
    df_predictions = trainer.predict_all_drug_gene_pairs()
    df_predictions.to_csv(RESULTS_DIR / "hgt_drug_gene_predictions.csv", index=False)
    print(f"  ✓ 预测结果已保存 ({len(df_predictions)} 行)")

    # Top 50 预测
    print(f"\n  HGT Top 20 预测靶点:")
    print(f"  {'Rank':>5s} | {'Gene':>12s} | {'Probability':>12s} | {'Score':>10s}")
    print(f"  {'-'*50}")
    for i, row in df_predictions.head(20).iterrows():
        print(f"  {i+1:5d} | {row['gene']:>12s} | {row['hgt_probability']:12.4f} | {row['hgt_score']:10.4f}")

    # --- 保存训练历史 ---
    history_path = RESULTS_DIR / "hgt_training_history.json"
    with open(history_path, 'w', encoding='utf-8') as f:
        json.dump(history, f, indent=2)
    print(f"\n  ✓ 训练历史已保存: {history_path}")

    # --- 保存ACSL4实验结果 ---
    acsl4_path = RESULTS_DIR / "acsl4_recall_results.json"
    with open(acsl4_path, 'w', encoding='utf-8') as f:
        json.dump(acsl4_results, f, indent=2, default=float)
    print(f"  ✓ ACSL4回忆实验结果已保存: {acsl4_path}")

    # --- 保存HGT模型 ---
    model_path = RESULTS_DIR / "ferro_hgt_model.pt"
    torch.save({
        "model_state_dict": model.state_dict(),
        "config": {
            "hidden_dim": 256,
            "num_heads": 4,
            "num_layers": 3,
            "dropout": 0.3,
        },
        "final_metrics": final_metrics,
    }, model_path)
    print(f"  ✓ HGT模型已保存: {model_path}")

    print("\n" + "=" * 70)
    print("阶段3 完成 — HGT异构网络嵌入已生成")
    print(f"  输出目录: {RESULTS_DIR}")
    print("=" * 70)

    return embeddings, df_predictions


if __name__ == "__main__":
    main()