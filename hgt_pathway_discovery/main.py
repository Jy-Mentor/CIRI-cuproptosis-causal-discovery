# -*- coding: utf-8 -*-
"""
HGT (Heterogeneous Graph Transformer) 基因‑通路关联预测模型 v2.0
用于网络药理学第二阶段通路发现任务。

============================================================================
参考文献:
  - HGT: Heterogeneous Graph Transformer, Hu et al., WWW 2020
  - GCNII: Simple and Deep GCN, Chen et al., ICML 2020
  - MHGTMDA: Molecular Heterogeneous Graph Transformer for miRNA-disease, Zou et al., 2024
  - HGTDR: Advancing Drug Repurposing with Heterogeneous Graph Transformers, Gharizadeh et al., 2024
  - PyG official hgt_dblp example (https://github.com/pyg-team/pytorch_geometric)
  - OGB Link Prediction Benchmark, Hu et al., NeurIPS 2020
  - Feature Propagation for Missing Node Features, Rossi et al., NeurIPS 2021
  - Disentangling Node Attributes for Link Prediction, Chatterjee et al., arXiv:2307.08877
  - Implicit degree bias in the link prediction task, Aiyappa et al., WWW 2024
  - HGNN-IMA: Multi-modal Heterogeneous Networks with Missing Modalities, Li et al., 2025
  - Deep Ensembles, Lakshminarayanan et al., NeurIPS 2017
  - Dropout as Bayesian Approximation, Gal & Ghahramani, ICML 2016
  - Platt Scaling, Platt, 1999

v2.0 改进:
  1. 模块化重构: config / data_loader / build_graph / model / train / inference
  2. YAML 配置管理: 消除硬编码路径
  3. 特征 Z-score 标准化: 统一多模态特征分布
  4. PPI 权重拆分: strong_ppi / weak_ppi 边类型
  5. 有向甲基化边: gene → CpG (单向)
  6. CpG 可学习参数: 特征传播初值 + 梯度优化
  7. 3 层 HGT + GCNII 初始残差: 缓解过平滑
  8. 解耦解码器: gene_bias + pathway_bias + DistMult 因子分解
  9. 复合早停: 0.7*AUROC + 0.3*AUPRC
  10. 自适应负采样比: 余弦退火 1:1 → 10:1
  11. Platt Scaling 校准
  12. 验证负边定时刷新
  13. CV 模型保存与集成预测
  14. MC Dropout 不确定性量化
  15. logging + TensorBoard
============================================================================
"""

import sys
import warnings
from pathlib import Path

import numpy as np
import torch

# 将项目根目录和 hgt_pathway_discovery 加入 sys.path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "hgt_pathway_discovery"))

from hgt_pathway_discovery.config import load_config
from hgt_pathway_discovery.utils import (
    set_seed, setup_logging, setup_tensorboard, try_compile,
    pca_reduce, pca_then_umap, normalize_features, setup_wandb,
)
from hgt_pathway_discovery.data_loader import load_all_data
from hgt_pathway_discovery.build_graph import build_hetero_graph
from hgt_pathway_discovery.model import HGTModel, NegEdgeSampler
from hgt_pathway_discovery.pretrain import pretrain_gae
from hgt_pathway_discovery.train import cross_validate, train_final
from hgt_pathway_discovery.inference import (
    predict_bridge_pathways,
    ensemble_predict_bridge_pathways,
    evaluate_inductive_subset,
)

warnings.filterwarnings("ignore")


def main() -> None:
    # --- 加载配置 ---
    config_path = Path(__file__).resolve().parent / "config.yaml"
    config = load_config(str(config_path))

    # --- 日志 ---
    log_cfg = config.logging
    log_file = log_cfg.log_file if log_cfg.log_file else None
    logger = setup_logging(
        "hgt_pipeline", log_cfg.level, log_file,
        enable_tensorboard=log_cfg.tensorboard,
    )

    # --- TensorBoard ---
    tb_writer = None
    if log_cfg.tensorboard:
        log_dir = Path(config.paths.log_dir)
        tb_writer = setup_tensorboard(str(log_dir))

    # --- WandB ---
    wandb_run = setup_wandb(config, logger)

    # --- 固定种子 ---
    set_seed(config.seed)
    logger.info(f"Device: {config.device}, Seed: {config.seed}")
    logger.info(f"Config loaded from {config_path}")

    # --- 打印关键配置 ---
    m_cfg = config.model
    t_cfg = config.training
    logger.info(f"Model: layers={m_cfg.num_layers}, hidden={m_cfg.hidden_dim}, "
                f"heads={m_cfg.num_heads}, dropout={m_cfg.dropout}")
    logger.info(f"Model: initial_residual={m_cfg.initial_residual}, "
                f"decoder_bias={m_cfg.decoder_bias}, factorization={m_cfg.decoder_factorization}")
    logger.info(f"Training: epochs={t_cfg.epochs}, patience={t_cfg.patience}, "
                f"lr={t_cfg.lr}, wd={t_cfg.weight_decay}")
    logger.info(f"Training: neg_mode={t_cfg.neg_sampling_mode}, "
                f"adaptive_neg={t_cfg.adaptive_neg_ratio}, "
                f"early_stop={t_cfg.early_stop_metric}")
    logger.info(f"Graph: ppi_split={config.graph.ppi_split_by_score}, "
                f"methyl_directed={config.graph.methylation_directed}")

    # ========================================================================
    # 1. 数据加载
    # ========================================================================
    logger.info("=" * 60)
    logger.info("[1/6] Loading data...")

    data_dict = load_all_data(
        Path(config.paths.data_dir),
        Path(config.paths.gat_data_dir),
        Path(config.paths.cache_dir),
        Path(config.paths.bridge_genes),
        config,
    )

    gene_feat_arr = data_dict["gene_feat_arr"]
    gene_feat_names = data_dict["gene_feat_names"]
    pathway_feat_arr = data_dict["pathway_feat_arr"]
    drug_fp_arr = data_dict["drug_fp_arr"]
    disease_feat_arr = data_dict["disease_feat_arr"]
    pathway_names = data_dict["pathway_names"]
    bridge_genes = data_dict["bridge_genes"]
    ppi_edges = data_dict["ppi_edges"]
    coexp_edges = data_dict["coexp_edges"]
    tf_edges = data_dict["tf_edges"]
    gene_pathway_edges = data_dict["gene_pathway_edges"]
    all_genes_list = data_dict["all_genes_list"]
    methyl_edges = data_dict.get("methyl_edges")
    mirna_edges = data_dict.get("mirna_edges")
    pathway_hierarchy = data_dict.get("pathway_hierarchy", [])
    disease_pathway_edges = data_dict.get("disease_pathway_edges", [])

    logger.info(f"  Gene features: {gene_feat_arr.shape}")
    logger.info(f"  Pathway features: {pathway_feat_arr.shape}")
    logger.info(f"  Bridge genes: {len(bridge_genes)}")
    logger.info(f"  PPI edges: {len(ppi_edges)}")
    logger.info(f"  Gene-pathway edges: {len(gene_pathway_edges)}")
    logger.info(f"  Pathway hierarchy: {len(pathway_hierarchy)}")

    # ========================================================================
    # 2. 特征预处理 (降维 + 标准化)
    # ========================================================================
    logger.info("=" * 60)
    logger.info("[2/6] Preprocessing features...")

    preproc_cfg = config.preprocessing
    feature_dim = preproc_cfg.feature_dim

    # 基因特征降维
    if preproc_cfg.use_pca and gene_feat_arr.shape[1] > feature_dim:
        if preproc_cfg.use_umap and gene_feat_arr.shape[1] > preproc_cfg.pca_intermediate_dim:
            logger.info(f"  Gene features: PCA({preproc_cfg.pca_intermediate_dim}) → UMAP({feature_dim})")
            gene_feat_arr = pca_then_umap(
                gene_feat_arr,
                preproc_cfg.pca_intermediate_dim,
                feature_dim,
                seed=config.seed,
            )
        else:
            gene_feat_arr, _ = pca_reduce(gene_feat_arr, feature_dim, config.seed)
            logger.info(f"  Gene features PCA: {gene_feat_arr.shape}")

    # 通路特征降维 (使用缓存)
    pathway_pca_cache = Path(config.paths.pathway_pca_cache)
    if preproc_cfg.use_pca and pathway_feat_arr.shape[1] > feature_dim:
        if pathway_pca_cache.exists():
            pathway_feat_arr = np.load(str(pathway_pca_cache)).astype(np.float32)
            logger.info(f"  Pathway PCA cache loaded: {pathway_feat_arr.shape}")
        else:
            pathway_feat_arr, _ = pca_reduce(pathway_feat_arr, feature_dim, config.seed)
            np.save(str(pathway_pca_cache), pathway_feat_arr)
            logger.info(f"  Pathway PCA: {pathway_feat_arr.shape} (cached)")

    # 通路名数量对齐
    if len(pathway_names) != pathway_feat_arr.shape[0]:
        logger.warning(f"  pathway_names ({len(pathway_names)}) != features ({pathway_feat_arr.shape[0]})")
        pathway_names = [f"pathway_{i}" for i in range(pathway_feat_arr.shape[0])]

    # ========================================================================
    # 3. 构建异构图
    # ========================================================================
    logger.info("=" * 60)
    logger.info("[3/6] Building heterogeneous graph...")

    data, gene_to_idx, gene_list, pathway_name_to_idx = build_hetero_graph(
        gene_feat_arr, gene_feat_names,
        drug_fp_arr, disease_feat_arr,
        pathway_feat_arr, pathway_names,
        ppi_edges, coexp_edges, tf_edges,
        gene_pathway_edges,
        all_genes_list=all_genes_list,
        bridge_genes=bridge_genes,
        methyl_edges=methyl_edges,
        mirna_edges=mirna_edges,
        pathway_hierarchy=pathway_hierarchy,
        disease_pathway_edges=disease_pathway_edges,
        config=config,
    )

    gp_edge_index = data["gene", "involved_in", "pathway"].edge_index
    n_genes = data["gene"].x.size(0)
    n_pathways = data["pathway"].x.size(0)

    if gp_edge_index.size(1) == 0:
        logger.error("No gene-pathway edges found! Aborting.")
        return

    logger.info(f"  Gene nodes: {n_genes}, Pathway nodes: {n_pathways}")
    logger.info(f"  Gene-pathway edges (positive): {gp_edge_index.size(1)}")
    logger.info(f"  Edge types: {list(data.edge_types)}")

    # ========================================================================
    # 4. 节点度计算 (度感知负采样)
    # ========================================================================
    logger.info("[4/6] Computing node degrees...")

    gp_ei_np = gp_edge_index.cpu().numpy()
    gene_degrees = np.ones(n_genes, dtype=np.float64)
    pathway_degrees = np.ones(n_pathways, dtype=np.float64)
    for i in range(gp_ei_np.shape[1]):
        gene_degrees[gp_ei_np[0, i]] += 1
        pathway_degrees[gp_ei_np[1, i]] += 1
    logger.info(f"  Gene degree: max={gene_degrees.max():.0f}, mean={gene_degrees.mean():.1f}")
    logger.info(f"  Pathway degree: max={pathway_degrees.max():.0f}, mean={pathway_degrees.mean():.1f}")

    # ========================================================================
    # 4.5. GAE 预训练 (可选)
    # ========================================================================
    pretrained_gae = None
    if config.pretraining.enabled:
        pretrained_gae = pretrain_gae(data, config, torch.device(config.device))

    # ========================================================================
    # 5. 交叉验证
    # ========================================================================
    logger.info("=" * 60)
    logger.info(f"[5/6] {config.cv.n_folds}-fold Cross Validation...")

    cv_results, cv_models, cv_scalers, cv_heldout_edges = cross_validate(
        data, gp_edge_index, n_genes, n_pathways,
        gene_degrees, pathway_degrees, config, tb_writer, wandb_run,
    )

    if cv_results:
        aurocs = [m["auroc"] for m in cv_results]
        auprcs = [m["auprc"] for m in cv_results]
        logger.info(f"  CV AUROC: {np.mean(aurocs):.4f} +/- {np.std(aurocs):.4f}")
        logger.info(f"  CV AUPRC: {np.mean(auprcs):.4f} +/- {np.std(auprcs):.4f}")

    # ========================================================================
    # 6. 最终训练 + 推理
    # ========================================================================
    logger.info("=" * 60)
    logger.info("[6/6] Final training & inference...")

    # 最终负采样器
    final_neg_sampler = NegEdgeSampler(
        pos_edges=gp_edge_index.t().tolist(),
        n_src=n_genes, n_dst=n_pathways,
        seed=config.seed + 999,
        mode=config.training.neg_sampling_mode,
        src_degrees=gene_degrees,
        dst_degrees=pathway_degrees,
        degree_power=config.training.neg_degree_power,
    )

    # 最终模型
    final_model = HGTModel(
        metadata=data.metadata(),
        dim_dict={nt: data[nt].x.size(-1) for nt in data.node_types},
        hidden_dim=m_cfg.hidden_dim, num_heads=m_cfg.num_heads,
        num_layers=m_cfg.num_layers, dropout=m_cfg.dropout,
        initial_residual=m_cfg.initial_residual,
        drop_edge_p=getattr(m_cfg, "drop_edge_p", 0.0),
        decoder_bias=m_cfg.decoder_bias,
        decoder_factorization=m_cfg.decoder_factorization,
        use_input_bn=getattr(m_cfg, "use_input_bn", True),
    ).to(config.device)

    # CpG 可学习参数注册
    if "cpg" in data.node_types:
        final_model.to_cpg_learnable(
            data["cpg"].x,
            quality_mask=data["cpg"].propagation_mask if hasattr(data["cpg"], "propagation_mask") else None,
        )
        logger.info("  CpG learnable parameters registered")

    final_model = try_compile(final_model)

    # GAE 预训练权重迁移
    if pretrained_gae is not None:
        logger.info("  Transferring pretrained GAE weights to HGTModel...")
        pretrained_gae.transfer_weights_to(final_model)
        logger.info("  Pretrained weights transferred successfully")

    final_auroc, final_auprc, final_scaler = train_final(
        final_model, data, gp_edge_index.to(config.device),
        final_neg_sampler, config, tb_writer,
    )
    logger.info(f"  Final: AUROC={final_auroc:.4f}, AUPRC={final_auprc:.4f}")

    # 保存最终模型
    model_path = Path(config.paths.model_save)
    torch.save(final_model.state_dict(), str(model_path))
    logger.info(f"  Model saved to {model_path}")

    # 基因索引映射
    gene_idx_to_name = {v: k for k, v in gene_to_idx.items()}

    # 最终推理 (转导 + Platt 校准 + MC Dropout)
    if bridge_genes:
        logger.info("  Running final inference...")
        mc_samples = config.inference.mc_dropout_samples
        if mc_samples > 0:
            logger.info(f"  MC Dropout enabled: {mc_samples} samples per gene")

        predict_bridge_pathways(
            final_model, data.to(config.device),
            bridge_genes, gene_to_idx, pathway_names,
            gene_idx_to_name=gene_idx_to_name,
            gp_edge_index=gp_edge_index,
            config=config,
            platt_scaler=final_scaler,
        )

        # 归纳子集评估
        if config.inference.eval_inductive:
            evaluate_inductive_subset(
                final_model, data.to(config.device),
                gp_edge_index, bridge_genes,
                gene_to_idx, pathway_names,
                gene_idx_to_name,
                final_neg_sampler, config,
                heldout_gp_edges=cv_heldout_edges,
            )

    # 集成推理
    if config.cv.ensemble_inference and cv_models:
        logger.info("  Running ensemble inference...")
        ensemble_predict_bridge_pathways(
            cv_models, data.to(config.device),
            bridge_genes, gene_to_idx, pathway_names,
            gp_edge_index=gp_edge_index,
            config=config,
            platt_scalers=cv_scalers,
        )

    # --- 清理 ---
    if tb_writer is not None:
        tb_writer.close()
    if wandb_run is not None:
        wandb_run.finish()
    if config.device == "cuda":
        torch.cuda.empty_cache()

    logger.info("=" * 60)
    logger.info("  Pipeline Complete!")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()