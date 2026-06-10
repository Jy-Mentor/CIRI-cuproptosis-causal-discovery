# -*- coding: utf-8 -*-
"""
注意力权重可解释性分析 — 独立脚本

从 gat_hetero_link_prediction.py 导入模型和数据管线，
训练一个单模型（非集成），提取每条边的注意力系数 alpha，
输出:
  attention_analysis/attention_summary.csv     — 边类型/层粒度汇总
  attention_analysis/attention_per_gene.csv    — 桥梁基因 Top-K 注意力邻居
  attention_analysis/attention_heatmap.png     — 注意力热度图（PDF出版级）

依赖:
  - 运行前确保已有 all_bridge_genes.csv（由 gat_hetero_link_prediction.py 生成）
  - 如需从零训练，直接运行本脚本即可
"""

import os, sys, gc, math, copy, warnings
import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F

# 导入 GAT 管线中的函数和类
from gat_hetero_link_prediction import (
    # Config
    GAT_DATA_DIR, SCRIPT_DIR, OUTPUT, device,
    HIDDEN_DIM, OUT_DIM, GAT_HEADS, DROPOUT, LR, WEIGHT_DECAY,
    EPOCHS, PATIENCE, N_FOLDS, N_LAYERS, VAL_RATIO, RANDOM_SEED,
    SWA_START, DROPEDGE_P, SIGMOID_TEMP, LOGIT_CLAMP,
    # Data loading
    load_drug_targets, load_disease_genes, load_ppi,
    load_gene_features, load_drug_fingerprint, load_all_genes,
    load_toxirna_features, build_hetero_data,
    # Model
    HeteroGAT, SigmoidAttnConv,
    # Training
    train_epoch, evaluate, sample_negative_edges_via_set,
    NegativeSampler, build_ppi_neighbor_set,
    update_swa_model, finalize_swa, get_dropedge_p,
    # Utils
    set_seed, logger, compute_ensemble_jaccard,
    analyze_attention_weights,
)

warnings.filterwarnings("ignore")
set_seed(RANDOM_SEED)

# ============================================================================
# 管道：数据加载 + 单模型训练 + 注意力分析
# ============================================================================

def main():
    logger.info("=" * 60)
    logger.info("  Attention Interpretability Analysis")
    logger.info("  Single-run model for attention extraction")
    logger.info("=" * 60)

    # ---- 加载数据 ----
    logger.info("[1/5] Loading data...")
    drug_target_genes = load_drug_targets(os.path.join(GAT_DATA_DIR, "drug_targets.txt"))
    disease_genes_raw = load_disease_genes(os.path.join(GAT_DATA_DIR, "disease_genes.txt"))
    ppi_edges = load_ppi(os.path.join(GAT_DATA_DIR, "ppi_subgraph.csv"))
    gene_features_arr, _ = load_gene_features(os.path.join(GAT_DATA_DIR, "subgraph_embeddings.csv"))
    drug_fingerprint_arr = load_drug_fingerprint(os.path.join(GAT_DATA_DIR, "drug_fingerprint.csv"))
    all_genes = load_all_genes(os.path.join(GAT_DATA_DIR, "subgraph_genes.txt"))
    toxirna_df = load_toxirna_features(os.path.join(SCRIPT_DIR, "toxirna_enhanced_features.csv"))

    logger.info("[2/5] Building heterogeneous graph...")
    data, gene_to_idx = build_hetero_data(
        drug_target_genes, disease_genes_raw, ppi_edges,
        gene_features_arr, _, drug_fingerprint_arr,
        all_genes, toxirna_df,
    )
    n_genes = data["gene"].x.size(0)
    gene_list_sorted = list(gene_to_idx.keys())
    data = data.to(device)

    # 提取边索引
    full_dt_ei = data["drug", "targets", "gene"].edge_index
    full_td_ei = data["gene", "associated_with", "disease"].edge_index
    full_ppi_ei = data["gene", "interacts", "gene"].edge_index

    # ---- 训练集/验证集划分 ----
    dt_all = torch.arange(full_dt_ei.size(1))
    td_all = torch.arange(full_td_ei.size(1))
    n_dt_val = max(1, int(len(dt_all) * VAL_RATIO))
    n_td_val = max(1, int(len(td_all) * VAL_RATIO))
    dt_perm = torch.randperm(len(dt_all), generator=torch.Generator().manual_seed(RANDOM_SEED))
    td_perm = torch.randperm(len(td_all), generator=torch.Generator().manual_seed(RANDOM_SEED))
    dt_train_idx = dt_all[dt_perm[n_dt_val:]]
    dt_val_idx = dt_all[dt_perm[:n_dt_val]]
    td_train_idx = td_all[td_perm[n_td_val:]]
    td_val_idx = td_all[td_perm[:n_td_val]]

    # 训练子图（隔离验证边）
    td_train_edges = full_td_ei[:, td_train_idx]
    train_disease_nodes = set(td_train_edges[0].tolist())
    gene_feat_full = data["gene"].x.cpu().numpy()
    disease_feat = np.mean(gene_feat_full[list(train_disease_nodes)], axis=0, keepdims=True) \
        if len(train_disease_nodes) > 0 else np.mean(gene_feat_full, axis=0, keepdims=True)

    train_data = type(data)()
    train_data["drug"].x = data["drug"].x.clone()
    train_data["gene"].x = data["gene"].x.clone()
    train_data["disease"].x = torch.from_numpy(disease_feat).to(device)
    train_data["drug", "targets", "gene"].edge_index = full_dt_ei[:, dt_train_idx].clone()
    train_data["gene", "associated_with", "disease"].edge_index = full_td_ei[:, td_train_idx].clone()
    train_data["gene", "interacts", "gene"].edge_index = full_ppi_ei.clone()

    # 候选池
    ppi_neighbors = build_ppi_neighbor_set(full_ppi_ei, n_genes)
    dt_target_indices = [gene_to_idx[g] for g in drug_target_genes if g in gene_to_idx]
    td_disease_indices = [gene_to_idx[g] for g in disease_genes_raw if g in gene_to_idx]

    dt_candidates = set()
    for idx in dt_target_indices: dt_candidates.update(ppi_neighbors.get(idx, set()))
    dt_candidates -= set(dt_target_indices)
    td_candidates = set()
    for idx in td_disease_indices: td_candidates.update(ppi_neighbors.get(idx, set()))
    td_candidates -= set(td_disease_indices)

    dt_sampler = NegativeSampler(np.array(sorted(dt_candidates), dtype=np.int64), n_genes)
    td_sampler = NegativeSampler(np.array(sorted(td_candidates), dtype=np.int64), n_genes)

    # 预缓存验证负边
    cached_neg_dt = sample_negative_edges_via_set(
        full_dt_ei[:, dt_val_idx], n_genes, dt_sampler, is_dt=True).to(device)
    cached_neg_td = sample_negative_edges_via_set(
        full_td_ei[:, td_val_idx], n_genes, td_sampler, is_dt=False).to(device)

    # ---- 训练单模型 ----
    logger.info("[3/5] Training single model (no ensemble)...")
    model = HeteroGAT(
        metadata=data.metadata(),
        drug_dim=data["drug"].x.size(1),
        gene_dim=data["gene"].x.size(1),
        disease_dim=data["disease"].x.size(1),
        hidden_dim=HIDDEN_DIM, out_dim=OUT_DIM,
        heads=GAT_HEADS, dropout=DROPOUT, n_layers=N_LAYERS,
    ).to(device)

    optimizer = torch.optim.AdamW(model.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)
    scaler = torch.cuda.amp.GradScaler()

    best_val = float("inf")
    patience_cnt = 0
    swa_model, swa_n = None, 0
    fold_best = {}

    for epoch in range(EPOCHS):
        de_p = get_dropedge_p(epoch, EPOCHS)
        loss = train_epoch(model, train_data, full_dt_ei, full_td_ei,
                           dt_train_idx, td_train_idx, n_genes,
                           optimizer, dt_sampler, td_sampler, scaler,
                           dropedge_p=de_p)

        if epoch >= SWA_START:
            swa_model, swa_n = update_swa_model(swa_model, model, swa_n)

        if (epoch + 1) % 50 == 0 or epoch == 0:
            dt_auroc, dt_auprc, td_auroc, td_auprc, td_p20, td_mrr, hits_k = evaluate(
                model, train_data, full_dt_ei, full_td_ei,
                dt_val_idx, td_val_idx, n_genes, dt_sampler, td_sampler,
                cached_neg_dt=cached_neg_dt, cached_neg_td=cached_neg_td)
            hits_at_k = hits_k.get(100, 0.0) if hits_k else 0.0
            val_metric = -((dt_auroc + td_auroc) / 2.0 + hits_at_k)
            logger.info(f"  Epoch {epoch+1:3d}/{EPOCHS} | Loss: {loss:.4f} | "
                  f"DT AUC: {dt_auroc:.4f} | TD AUC: {td_auroc:.4f} | "
                  f"H@100: {hits_at_k:.3f}")

            if val_metric < best_val:
                best_val = val_metric
                patience_cnt = 0
                fold_best = dict(epoch=epoch+1, dt_auroc=dt_auroc, td_auroc=td_auroc)
                # 保存最佳模型
                best_state = copy.deepcopy(model.state_dict())
            else:
                patience_cnt += 1
                if patience_cnt >= PATIENCE:
                    logger.info(f"  [Early stop] epoch {epoch+1}")
                    break

    # 恢复最佳模型（早停）
    if best_state:
        model.load_state_dict(best_state)
    logger.info(f"  [Best] Epoch {fold_best.get('epoch','?')}: "
          f"DT AUC={fold_best.get('dt_auroc',0):.4f}, "
          f"TD AUC={fold_best.get('td_auroc',0):.4f}")

    # ---- 全图推理（启用注意力存储）----
    logger.info("[4/5] Extracting attention weights from full graph...")
    data = data.to(device)
    analyze_attention_weights(model, data, gene_to_idx, gene_list_sorted,
                               top_n_bridge=20, top_k_neighbors=10)

    # ---- 生成注意力热度图 ----
    logger.info("[5/5] Generating attention visualizations...")
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import seaborn as sns

        out_dir = os.path.join(SCRIPT_DIR, "attention_analysis")
        os.makedirs(out_dir, exist_ok=True)

        # 加载注意力汇总
        summary_path = os.path.join(out_dir, "attention_summary.csv")
        if os.path.exists(summary_path):
            summary_df = pd.read_csv(summary_path)

            # 热度图: 边类型 × 层
            fig, axes = plt.subplots(1, 2, figsize=(14, 5))

            # 左图: 平均 α × 边类型
            ax1 = axes[0]
            pivot1 = summary_df.pivot_table(
                index="edge_type", columns="layer", values="mean_alpha")
            sns.heatmap(pivot1, annot=True, fmt=".4f", cmap="YlOrRd",
                       ax=ax1, cbar_kws={"label": "Mean Attention α"})
            ax1.set_title("Mean Attention Coefficient per Edge Type × Layer")
            ax1.set_ylabel("Edge Type")
            ax1.set_xlabel("GNN Layer")

            # 右图: α 分布 (mean ± std)
            ax2 = axes[1]
            x_pos = np.arange(len(summary_df))
            ax2.bar(x_pos - 0.2, summary_df["mean_alpha"], 0.4,
                   label="Mean α", color="steelblue")
            ax2.bar(x_pos + 0.2, summary_df["std_alpha"], 0.4,
                   label="Std α", color="tomato")
            ax2.set_xticks(x_pos)
            ax2.set_xticklabels(
                [f"{r['edge_type'][:20]}·{r['layer']}" for _, r in summary_df.iterrows()],
                rotation=45, ha="right", fontsize=8)
            ax2.set_ylabel("Attention α")
            ax2.set_title("Attention Distribution (Mean ± Std)")
            ax2.legend()
            ax2.grid(axis="y", alpha=0.3)

            plt.tight_layout()
            heatmap_path = os.path.join(out_dir, "attention_heatmap.png")
            plt.savefig(heatmap_path, dpi=150, bbox_inches="tight")
            plt.close()
            logger.info(f"  [Output] 注意力热度图 → {heatmap_path}")

        # 加载基因级注意力 → 打印 Top-5 桥梁基因的驱动因子
        per_gene_path = os.path.join(out_dir, "attention_per_gene.csv")
        if os.path.exists(per_gene_path):
            attn_df = pd.read_csv(per_gene_path)
            top_bridge = attn_df["gene_symbol"].unique()[:5]
            for g in top_bridge:
                gdf = attn_df[attn_df["gene_symbol"] == g].head(10)
                logger.info(f"\n  Top-10 attention edges for {g}:")
                for _, row in gdf.iterrows():
                    logger.info(f"    rank={row['rank']}: "
                          f"{row['neighbor']:>20s} "
                          f"({row['edge_type']:30s} "
                          f"layer={row['layer']:5s}) "
                          f"α={row['attention_score']:.4f}")

    except ImportError:
        logger.warning("[Attention] matplotlib/seaborn not installed, skipping heatmap")

    logger.info(f"\n{'='*60}")
    logger.info(f"  Attention analysis complete!")
    logger.info(f"  Results in: {os.path.join(SCRIPT_DIR, 'attention_analysis')}")
    logger.info(f"{'='*60}")


if __name__ == "__main__":
    main()