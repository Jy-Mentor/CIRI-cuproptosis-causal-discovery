# -*- coding: utf-8 -*-
"""快速诊断脚本：测试模型前向传播是否正常。"""
import sys
from pathlib import Path

_PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_PROJECT_ROOT))

from hgt_pathway_discovery.config import load_config
from hgt_pathway_discovery.utils import set_seed, pca_reduce
from hgt_pathway_discovery.data_loader import load_all_data
from hgt_pathway_discovery.build_graph import build_hetero_graph
from hgt_pathway_discovery.model import HGTModel
import numpy as np
import torch
import time

cfg = load_config(str(_PROJECT_ROOT / "hgt_pathway_discovery" / "config.yaml"))
set_seed(cfg.seed)

print("Loading data...")
data_dict = load_all_data(
    Path(cfg.paths.data_dir), Path(cfg.paths.gat_data_dir),
    Path(cfg.paths.cache_dir), Path(cfg.paths.bridge_genes), cfg,
)

gene_feat_arr = data_dict["gene_feat_arr"]
pathway_feat_arr = data_dict["pathway_feat_arr"]
gene_feat_arr, _ = pca_reduce(gene_feat_arr, 256, cfg.seed)

# 通路缓存
pathway_pca_cache = Path(cfg.paths.pathway_pca_cache)
if pathway_pca_cache.exists():
    pathway_feat_arr = np.load(str(pathway_pca_cache)).astype(np.float32)
    print(f"  Pathway PCA cache loaded: {pathway_feat_arr.shape}")

print("Building graph...")
data, gene_to_idx, gene_list, pathway_name_to_idx = build_hetero_graph(
    gene_feat_arr, data_dict["gene_feat_names"],
    data_dict["drug_fp_arr"], data_dict["disease_feat_arr"],
    pathway_feat_arr, data_dict["pathway_names"],
    data_dict["ppi_edges"], data_dict["coexp_edges"],
    data_dict["tf_edges"], data_dict["gene_pathway_edges"],
    all_genes_list=data_dict["all_genes_list"],
    bridge_genes=data_dict["bridge_genes"],
    methyl_edges=data_dict.get("methyl_edges"),
    mirna_edges=data_dict.get("mirna_edges"),
    pathway_hierarchy=data_dict.get("pathway_hierarchy", []),
    config=cfg,
)

print(f"Nodes: {data.node_types}")
for et in data.edge_types:
    ei = data[et].edge_index
    print(f"  {et}: {ei.size(1)} edges")

for nt in data.node_types:
    print(f"  {nt}: x.shape={data[nt].x.shape}")

if torch.cuda.is_available():
    print(f"GPU: allocated={torch.cuda.memory_allocated()/1e9:.2f} GB, "
          f"reserved={torch.cuda.memory_reserved()/1e9:.2f} GB")

device = cfg.device
print(f"\nDevice: {device}")

model = HGTModel(
    metadata=data.metadata(),
    dim_dict={nt: data[nt].x.size(-1) for nt in data.node_types},
    hidden_dim=cfg.model.hidden_dim, num_heads=cfg.model.num_heads,
    num_layers=cfg.model.num_layers, dropout=cfg.model.dropout,
    initial_residual=cfg.model.initial_residual,
    decoder_bias=cfg.model.decoder_bias,
    decoder_factorization=cfg.model.decoder_factorization,
    use_input_bn=getattr(cfg.model, "use_input_bn", True),
)
total_params = sum(p.numel() for p in model.parameters())
print(f"Model params: {total_params:,}")

if "cpg" in data.node_types:
    model.to_cpg_learnable(data["cpg"].x)
    print(f"CpG learnable params registered")

model = model.to(device)
print(f"Moving data to {device}...")
data = data.to(device)

print(f"GPU after model+data: allocated={torch.cuda.memory_allocated()/1e9:.2f} GB, "
      f"reserved={torch.cuda.memory_reserved()/1e9:.2f} GB")

print("\nRunning forward pass...")
model.eval()
with torch.inference_mode():
    t0 = time.time()
    z_dict = model(data.x_dict, data.edge_index_dict)
    t1 = time.time()

print(f"Forward pass: {t1-t0:.2f}s")
print(f"Gene embeddings: {z_dict['gene'].shape}")
if "pathway" in z_dict:
    print(f"Pathway embeddings: {z_dict['pathway'].shape}")
print("\nSUCCESS: Model forward pass works correctly!")