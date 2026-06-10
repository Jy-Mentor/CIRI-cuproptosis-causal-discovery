# -*- coding: utf-8 -*-
"""
GAT Embedding Extractor v4 — GATv2 + Residual + DropEdge + Focal Loss (synced w/ main)
  Model_A (no TD edges): for DG classifier  — blind to disease-gene links
  Model_B (no DT edges): for DT classifier  — blind to drug-target links
  v4 upgrades from gat_hetero_link_prediction.py:
    - GATv2Conv (dynamic attention) instead of GATConv
    - Residual connections (gene, disease)
    - DropEdge on PPI during training
    - Disease features = mean of disease-gene embeddings (not ones)
    - Focal BCE Loss (alpha=0.5, gamma=2.0)
    - share_weights=False
    - AMP mixed precision (float16)
"""
import os, sys, warnings
import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from torch.cuda.amp import GradScaler, autocast
from torch_geometric.data import HeteroData
from sklearn.decomposition import PCA

warnings.filterwarnings("ignore")

GAT_DIR = r"C:\Users\Jy-Mentor-7\Desktop\GAT"
RF_DIR = r"C:\Users\Jy-Mentor-7\Desktop\随机森林"
FEATURE_TABLE = os.path.join(RF_DIR, "gene_features_table_with_gat.csv")

HIDDEN_DIM = 128
OUT_DIM = 64
GAT_HEADS = 4
DROPOUT = 0.3
CONV_DROPOUT = 0.18
PPI_DROPEDGE = 0.15
LR = 5e-4
WEIGHT_DECAY = 5e-5
EPOCHS_A = 400
EPOCHS_B = 600
PATIENCE = 60
RANDOM_SEED = 42
N_PCA = 12
NEG_RATIO = 3
FOCAL_ALPHA = 0.5
FOCAL_GAMMA = 2.0

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"[Init] Device: {device}")

import random
def set_seed(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)

set_seed(RANDOM_SEED)

def _read_lines(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Not found: {path}")
    for enc in ("utf-8", "gbk", "latin1"):
        try:
            with open(path, "r", encoding=enc) as fh:
                lines = [l.strip() for l in fh if l.strip()]
            return lines
        except UnicodeDecodeError:
            continue
    return []

def load_gene_list(path):
    lines = _read_lines(path)
    genes = []
    if lines and ("," in lines[0] or "\t" in lines[0]):
        for ln in lines:
            parts = ln.replace("\t", ",").split(",")
            g = parts[0].strip().upper()
            if g:
                genes.append(g)
    else:
        genes = [g.upper() for g in lines if g]
    return list(dict.fromkeys(genes))

def load_ppi(path):
    lines = _read_lines(path)
    ppi_edges = []
    for line in lines:
        if "protein" in line.lower():
            continue
        parts = line.replace("\t", ",").split(",")
        if len(parts) >= 2:
            a, b = parts[0].strip().upper(), parts[1].strip().upper()
            if a and b:
                ppi_edges.append((a, b))
    return ppi_edges

def load_gene_features(path):
    ft = pd.read_csv(path)
    if "gene_symbol" in ft.columns:
        ft = ft.set_index("gene_symbol")
    arr = ft.values.astype(np.float32)
    names = ft.index.tolist()
    arr = np.nan_to_num(arr, 0.0)
    return arr, names

def load_drug_fingerprint(path):
    df = pd.read_csv(path)
    fp = df.drop(columns=df.columns[:1]).values.astype(np.float32).flatten()
    return fp

def load_all_genes(path):
    with open(path) as f:
        return [l.strip().upper() for l in f if l.strip()]

def build_graph(dt_list, dg_list, ppi_edges, gene_feat, gene_names, drug_fp, all_genes):
    data = HeteroData()

    all_genes_set = set(g.upper() for g in all_genes)
    gene_to_idx = {g: i for i, g in enumerate(all_genes)}
    n_genes = len(all_genes)
    gene_dim = gene_feat.shape[1]

    gene_feat_dict = {g.upper(): gene_feat[i] for i, g in enumerate(gene_names) if i < len(gene_feat)}
    gx = np.zeros((n_genes, gene_dim), dtype=np.float32)
    for g in all_genes:
        if g in gene_feat_dict:
            gx[gene_to_idx[g]] = gene_feat_dict[g]
    data["gene"].x = torch.from_numpy(gx).float()
    data["gene"].num_nodes = n_genes

    data["drug"].x = torch.from_numpy(drug_fp.reshape(1, -1)).float()
    data["drug"].num_nodes = 1
    drug_dim = drug_fp.shape[0]

    dg_indices = np.array([gene_to_idx[g] for g in dg_list if g in gene_to_idx], dtype=np.int64)
    if len(dg_indices) > 0:
        disease_feat = np.mean(gx[dg_indices], axis=0, keepdims=True)
    else:
        disease_feat = np.zeros((1, gene_dim), dtype=np.float32)
    data["disease"].x = torch.from_numpy(disease_feat).float()
    data["disease"].num_nodes = 1

    dt_src, dt_dst = [], []
    for g in dt_list:
        if g in gene_to_idx:
            dt_src.append(0)
            dt_dst.append(gene_to_idx[g])
    data["drug", "targets", "gene"].edge_index = torch.tensor([dt_src, dt_dst], dtype=torch.long)

    td_src, td_dst = [], []
    for g in dg_list:
        if g in gene_to_idx:
            td_src.append(gene_to_idx[g])
            td_dst.append(0)
    data["gene", "associated_with", "disease"].edge_index = torch.tensor([td_src, td_dst], dtype=torch.long)

    ppi_src, ppi_dst = [], []
    for a, b in ppi_edges:
        if a in gene_to_idx and b in gene_to_idx:
            ppi_src.append(gene_to_idx[a])
            ppi_dst.append(gene_to_idx[b])
    data["gene", "interacts", "gene"].edge_index = torch.tensor([ppi_src, ppi_dst], dtype=torch.long)

    print(f"[Build] DT:{len(dt_src)} TD:{len(td_src)} PPI:{len(ppi_src)} Genes:{n_genes}")

    meta = {
        "gene_dim": gene_dim,
        "drug_dim": drug_dim,
        "disease_dim": gene_dim,
    }
    return data, gene_to_idx, all_genes, meta

from torch_geometric.nn import GATv2Conv, HeteroConv
from torch_geometric.utils import dropout_edge

class HeteroGATv2(torch.nn.Module):
    def __init__(self, metadata, drug_dim, gene_dim, disease_dim,
                 hidden_dim, out_dim, heads, dropout, conv_dropout):
        super().__init__()
        self.dropout_rate = dropout
        self.conv_dropout = conv_dropout

        self.drug_proj = torch.nn.Linear(drug_dim, hidden_dim)
        self.gene_proj = torch.nn.Linear(gene_dim, hidden_dim)
        self.disease_proj = torch.nn.Linear(disease_dim, hidden_dim)

        self.gene_residual = torch.nn.Linear(gene_dim, hidden_dim)
        self.disease_residual = torch.nn.Linear(disease_dim, hidden_dim)

        convs1 = {}
        for etype in metadata[1]:
            convs1[etype] = GATv2Conv(
                (-1, -1), hidden_dim // heads, heads=heads, concat=True,
                dropout=dropout, add_self_loops=(etype[0] == "gene" and etype[2] == "gene"),
                share_weights=False)
        self.conv1 = HeteroConv(convs1, aggr="mean")

        convs2 = {}
        for etype in metadata[1]:
            convs2[etype] = GATv2Conv(
                (-1, -1), out_dim, heads=1, concat=False,
                dropout=dropout, add_self_loops=(etype[0] == "gene" and etype[2] == "gene"),
                share_weights=False)
        self.conv2 = HeteroConv(convs2, aggr="mean")

        self.drug_out = torch.nn.Linear(hidden_dim, out_dim)
        self.disease_out = torch.nn.Linear(hidden_dim, out_dim)

        self.dt_decoder = torch.nn.Bilinear(out_dim, out_dim, 1)
        self.td_decoder = torch.nn.Bilinear(out_dim, out_dim, 1)

    def forward(self, x_dict, edge_index_dict, ppi_dropedge_prob=0.0):
        eid = edge_index_dict
        if ppi_dropedge_prob > 0 and self.training:
            ppi_key = ("gene", "interacts", "gene")
            if ppi_key in eid:
                ei, mask = dropout_edge(eid[ppi_key], p=ppi_dropedge_prob,
                                        force_undirected=True, training=True)
                eid = dict(eid)
                eid[ppi_key] = ei
                edge_index_dict = eid

        x_proj = {
            "drug": self.drug_proj(x_dict["drug"]),
            "gene": self.gene_proj(x_dict["gene"]),
            "disease": self.disease_proj(x_dict["disease"]),
        }

        gene_res = self.gene_residual(x_dict["gene"])
        disease_res = self.disease_residual(x_dict["disease"])

        out1 = self.conv1(x_proj, edge_index_dict)
        for k in ["drug", "gene", "disease"]:
            if k not in out1:
                out1[k] = x_proj[k]

        out1["gene"] = out1["gene"] + gene_res
        out1["disease"] = out1["disease"] + disease_res

        out1["drug"] = self.drug_out(out1["drug"])
        out1["disease"] = self.disease_out(out1["disease"])

        out1 = {k: F.relu(v) for k, v in out1.items()}
        out1 = {k: F.dropout(v, p=self.dropout_rate, training=self.training)
                for k, v in out1.items()}

        out2 = self.conv2(out1, edge_index_dict)
        for k in ["drug", "gene", "disease"]:
            if k not in out2:
                out2[k] = out1[k]

        out2["gene"] = F.dropout(out2["gene"], p=self.conv_dropout, training=self.training)

        return out2

class NegSampler:
    def __init__(self, candidates, n_genes, seed=42):
        self.candidates = torch.from_numpy(candidates)
        self.rng = torch.Generator()
        self.rng.manual_seed(seed)

    def sample(self, pos_edge, is_dt, neg_ratio=1):
        n_pos = pos_edge.size(1)
        n_neg = n_pos * neg_ratio
        idx = torch.randint(0, len(self.candidates), (n_neg,), generator=self.rng)
        neg = self.candidates[idx].to(pos_edge.device)
        if is_dt:
            src = pos_edge[0].repeat_interleave(neg_ratio)
            return torch.stack([src, neg])
        else:
            dst = pos_edge[1].repeat_interleave(neg_ratio)
            return torch.stack([neg, dst])

def build_neg_candidates(ppi_ei, pos_indices, n_genes):
    ppi_src = ppi_ei[0].cpu().numpy()
    ppi_dst = ppi_ei[1].cpu().numpy()
    nb = {i: set() for i in range(n_genes)}
    for s, d in zip(ppi_src, ppi_dst):
        nb[s].add(d)
        nb[d].add(s)
    cand = set()
    for idx in pos_indices:
        cand.update(nb.get(idx, set()))
    cand -= set(pos_indices)
    return np.array(sorted(cand), dtype=np.int64)

def focal_bce_loss(logits, targets, alpha=FOCAL_ALPHA, gamma=FOCAL_GAMMA):
    bce = F.binary_cross_entropy_with_logits(logits, targets, reduction="none")
    pt = torch.exp(-bce)
    alpha_t = targets * alpha + (1 - targets) * (1 - alpha)
    focal = alpha_t * ((1 - pt) ** gamma) * bce
    return focal.mean()

def train_model(model, data, dt_edges, td_edges, dt_sampler, td_sampler,
                n_genes, neg_ratio, epochs, patience, model_name):
    opt = torch.optim.AdamW(model.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)
    scaler = GradScaler(enabled=(device.type == "cuda"))
    best_loss, wait = float("inf"), 0
    xd = data.x_dict
    eid = data.edge_index_dict

    for ep in range(epochs):
        model.train()
        opt.zero_grad()

        with autocast(enabled=(device.type == "cuda")):
            z = model(xd, eid, ppi_dropedge_prob=PPI_DROPEDGE)

            loss_dt = torch.tensor(0.0, device=device)
            if dt_edges is not None:
                pos = dt_edges.to(device)
                neg = dt_sampler.sample(pos, is_dt=True, neg_ratio=neg_ratio)
                e = torch.cat([pos, neg], dim=1)
                lbl = torch.cat([torch.ones(pos.size(1)), torch.zeros(neg.size(1))]).to(device)
                logits = model.dt_decoder(z["drug"][e[0]], z["gene"][e[1]]).squeeze(-1)
                loss_dt = focal_bce_loss(logits, lbl)

            loss_td = torch.tensor(0.0, device=device)
            if td_edges is not None:
                pos = td_edges.to(device)
                neg = td_sampler.sample(pos, is_dt=False, neg_ratio=neg_ratio)
                e = torch.cat([pos, neg], dim=1)
                lbl = torch.cat([torch.ones(pos.size(1)), torch.zeros(neg.size(1))]).to(device)
                logits = model.td_decoder(z["gene"][e[0]], z["disease"][e[1]]).squeeze(-1)
                loss_td = focal_bce_loss(logits, lbl)

            loss = loss_dt + loss_td

        scaler.scale(loss).backward()
        scaler.unscale_(opt)
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        scaler.step(opt)
        scaler.update()

        loss_val = loss.item()
        if ep % 100 == 0:
            print(f"  E{ep:4d}: loss={loss_val:.4f}  DT={loss_dt.item():.4f}  TD={loss_td.item():.4f}")

        if loss_val < best_loss - 1e-4:
            best_loss = loss_val
            wait = 0
        else:
            wait += 1
        if wait >= patience:
            print(f"  Early stop @ epoch {ep}")
            break

    print(f"  [{model_name}] Best loss: {best_loss:.4f}, epochs: {ep+1}")
    return best_loss

def extract_embeddings(model, data):
    model.eval()
    with torch.no_grad():
        with autocast(enabled=(device.type == "cuda")):
            z = model(data.x_dict, data.edge_index_dict, ppi_dropedge_prob=0.0)
        emb = z["gene"].cpu().float().numpy()
    return emb

def main():
    print("=" * 60)
    print("  GAT Embedding v4 — GATv2 + Residual + DropEdge + Focal (synced)")
    print("  Model_A: no-TD (for DG classifier)")
    print("  Model_B: no-DT (for DT/DG classifiers)")
    print("=" * 60)

    dt_genes = load_gene_list(os.path.join(GAT_DIR, "drug_targets.txt"))
    dg_genes = load_gene_list(os.path.join(GAT_DIR, "disease_genes.txt"))
    ppi_edges = load_ppi(os.path.join(GAT_DIR, "ppi_subgraph.csv"))
    gene_feat, gene_names = load_gene_features(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "subgraph_embeddings.csv"))
    drug_fp = load_drug_fingerprint(os.path.join(GAT_DIR, "drug_fingerprint.csv"))
    all_genes = load_all_genes(os.path.join(GAT_DIR, "subgraph_genes.txt"))

    print(f"[Data] DT genes:{len(dt_genes)} DG genes:{len(dg_genes)} PPI edges:{len(ppi_edges)} All genes:{len(all_genes)}")

    data, gene_to_idx, gene_list, meta = build_graph(
        dt_genes, dg_genes, ppi_edges, gene_feat, gene_names, drug_fp, all_genes)
    n_genes = len(gene_list)

    full_dt = data["drug", "targets", "gene"].edge_index
    full_td = data["gene", "associated_with", "disease"].edge_index
    full_ppi = data["gene", "interacts", "gene"].edge_index

    dt_idx = [gene_to_idx[g] for g in dt_genes if g in gene_to_idx]
    dg_idx = [gene_to_idx[g] for g in dg_genes if g in gene_to_idx]
    dt_cand = build_neg_candidates(full_ppi, dt_idx, n_genes)
    dg_cand = build_neg_candidates(full_ppi, dg_idx, n_genes)
    print(f"[NegCand] DT candidates:{len(dt_cand)} DG candidates:{len(dg_cand)}")

    g_dim = meta["gene_dim"]
    d_dim = meta["drug_dim"]
    dd_dim = meta["disease_dim"]

    # ════════════════════════════════════════
    #  MODEL A: no TD edges (for DG classifier)
    # ════════════════════════════════════════
    print("\n" + "=" * 50)
    print("  MODEL A: DT + PPI only (no TD) — for DG classifier")
    print("=" * 50)

    data_a = HeteroData()
    data_a["drug"].x = data["drug"].x
    data_a["gene"].x = data["gene"].x
    data_a["disease"].x = data["disease"].x
    data_a["drug", "targets", "gene"].edge_index = full_dt
    data_a["gene", "interacts", "gene"].edge_index = full_ppi
    data_a = data_a.to(device)

    dt_sampler_a = NegSampler(dt_cand, n_genes, RANDOM_SEED)

    set_seed(RANDOM_SEED)
    model_a = HeteroGATv2(
        metadata=data_a.metadata(),
        drug_dim=d_dim, gene_dim=g_dim, disease_dim=dd_dim,
        hidden_dim=HIDDEN_DIM, out_dim=OUT_DIM, heads=GAT_HEADS,
        dropout=DROPOUT, conv_dropout=CONV_DROPOUT,
    ).to(device)

    train_model(model_a, data_a, full_dt, None, dt_sampler_a, None,
                n_genes, NEG_RATIO, EPOCHS_A, PATIENCE, "Model-A")
    emb_a = extract_embeddings(model_a, data_a)

    # ════════════════════════════════════════
    #  MODEL B: no DT edges (for DT/DG classifiers)
    # ════════════════════════════════════════
    print("\n" + "=" * 50)
    print("  MODEL B: TD + PPI only (no DT) — for DT/DG classifiers")
    print("=" * 50)

    data_b = HeteroData()
    data_b["drug"].x = data["drug"].x
    data_b["gene"].x = data["gene"].x
    data_b["disease"].x = data["disease"].x
    data_b["gene", "associated_with", "disease"].edge_index = full_td
    data_b["gene", "interacts", "gene"].edge_index = full_ppi
    data_b = data_b.to(device)

    td_sampler_b = NegSampler(dg_cand, n_genes, RANDOM_SEED)

    set_seed(RANDOM_SEED)
    model_b = HeteroGATv2(
        metadata=data_b.metadata(),
        drug_dim=d_dim, gene_dim=g_dim, disease_dim=dd_dim,
        hidden_dim=HIDDEN_DIM, out_dim=OUT_DIM, heads=GAT_HEADS,
        dropout=DROPOUT, conv_dropout=CONV_DROPOUT,
    ).to(device)

    train_model(model_b, data_b, None, full_td, None, td_sampler_b,
                n_genes, NEG_RATIO, EPOCHS_B, PATIENCE, "Model-B")
    emb_b = extract_embeddings(model_b, data_b)

    # ════════════════════════════════════════
    #  Embedding Diagnostics
    # ════════════════════════════════════════
    print("\n" + "=" * 50)
    print("  EMBEDDING DIAGNOSTICS")
    print("=" * 50)

    var_raw_a = emb_a.var(axis=0).sum()
    var_raw_b = emb_b.var(axis=0).sum()
    mean_norm_a = np.linalg.norm(emb_a, axis=1).mean()
    mean_norm_b = np.linalg.norm(emb_b, axis=1).mean()
    zero_var_a = (emb_a.var(axis=0) < 1e-8).sum()
    zero_var_b = (emb_b.var(axis=0) < 1e-8).sum()
    print(f"[Raw-Emb-A] no-TD: total_var={var_raw_a:.4f}  mean_norm={mean_norm_a:.4f}  zero_var_cols={zero_var_a}")
    print(f"[Raw-Emb-B] no-DT: total_var={var_raw_b:.4f}  mean_norm={mean_norm_b:.4f}  zero_var_cols={zero_var_b}")

    pca_a = PCA(n_components=N_PCA, random_state=RANDOM_SEED)
    emb_a_pca = pca_a.fit_transform(emb_a)
    var_a = pca_a.explained_variance_ratio_.sum()
    print(f"\n[PCA-A] no-TD: 64d → {N_PCA}d, cumulative var={var_a:.4f}")
    for i, v in enumerate(pca_a.explained_variance_ratio_):
        print(f"  PC{i:2d}: {v:.4f}")

    pca_b = PCA(n_components=N_PCA, random_state=RANDOM_SEED)
    emb_b_pca = pca_b.fit_transform(emb_b)
    var_b = pca_b.explained_variance_ratio_.sum()
    print(f"\n[PCA-B] no-DT: 64d → {N_PCA}d, cumulative var={var_b:.4f}")
    for i, v in enumerate(pca_b.explained_variance_ratio_):
        print(f"  PC{i:2d}: {v:.4f}")

    var_pca_a = emb_a_pca.var(axis=0).sum()
    var_pca_b = emb_b_pca.var(axis=0).sum()
    print(f"\n[PCA-Var-A] no-TD: total_var={var_pca_a:.4f}")
    print(f"[PCA-Var-B] no-DT: total_var={var_pca_b:.4f}")

    cols_a = [f"gat_nodg_{i}" for i in range(N_PCA)]
    cols_b = [f"gat_nodt_{i}" for i in range(N_PCA)]

    df_a = pd.DataFrame(emb_a_pca, columns=cols_a, index=gene_list).reset_index()
    df_a.columns = ["gene_symbol"] + cols_a
    df_a["gene_symbol"] = df_a["gene_symbol"].str.upper()

    df_b = pd.DataFrame(emb_b_pca, columns=cols_b, index=gene_list).reset_index()
    df_b.columns = ["gene_symbol"] + cols_b
    df_b["gene_symbol"] = df_b["gene_symbol"].str.upper()

    ft = pd.read_csv(FEATURE_TABLE)
    ft["gene_symbol"] = ft["gene_symbol"].str.upper()
    for c in cols_a + cols_b:
        if c in ft.columns:
            ft = ft.drop(columns=[c])

    merged = ft.merge(df_a, on="gene_symbol", how="left").merge(df_b, on="gene_symbol", how="left")
    for c in cols_a + cols_b:
        merged[c] = merged[c].fillna(0.0)

    out_path = os.path.join(RF_DIR, "gene_features_table_with_gat_emb.csv")
    merged.to_csv(out_path, index=False, encoding="utf-8-sig")

    n_matched_a = (merged[cols_a].var(axis=1) > 1e-8).sum()
    n_matched_b = (merged[cols_b].var(axis=1) > 1e-8).sum()
    print(f"\n[Save] → {out_path}  ({merged.shape[0]} x {merged.shape[1]})")
    print(f"[Match] no-TD: {n_matched_a}/{n_genes}  no-DT: {n_matched_b}/{n_genes}")

    vt = var_raw_a, var_raw_b, var_pca_a, var_pca_b
    print(f"\n  Variance comparison | v3 → v4")
    print(f"  Model-A (no-TD):  raw={var_raw_a:.4f}  pca={var_pca_a:.4f}")
    print(f"  Model-B (no-DT):  raw={var_raw_b:.4f}  pca={var_pca_b:.4f}")

    print("\n[Done] GAT v4 embeddings ready with GATv2 + residual + focal.")

if __name__ == "__main__":
    main()