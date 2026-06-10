# -*- coding: utf-8 -*-
"""
GAT Embedding Extractor v6 — All-OOF (No Full Model, Uniform Scale)
  v5 problem: positives got OOF embeddings (norm ~1.93), negatives got
  full-model embeddings (norm ~0.15). PCA captured this 12x scale gap
  as PC1, creating pseudo-leakage (gat_nodg_0 AUC vs DT = 1.0).

  v6 fix:
    - Positives (DT/DG genes): OOF from their validation fold only (unchanged)
    - Negatives (non-DT/non-DG genes): MEAN across all 5 OOF folds
    - NO full model used anywhere → all embeddings from same-scale models
    - Result: uniform embedding norm, PCA captures genuine biology

  Model_A (no-TD): DT + PPI — OOF on DT edges → gat_nodg_* for DT classifier
  Model_B (no-DT): TD + PPI — OOF on TD edges → gat_nodt_* for DG classifier
"""
import os, sys, warnings
import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from torch.cuda.amp import GradScaler, autocast
from torch_geometric.data import HeteroData
from sklearn.decomposition import PCA
from sklearn.model_selection import KFold

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
N_FOLDS = 5

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"[Init] Device: {device}")

import random
def set_seed(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)

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
            if dt_edges is not None and dt_edges.size(1) > 0:
                pos = dt_edges.to(device)
                neg = dt_sampler.sample(pos, is_dt=True, neg_ratio=neg_ratio)
                e = torch.cat([pos, neg], dim=1)
                lbl = torch.cat([torch.ones(pos.size(1)), torch.zeros(neg.size(1))]).to(device)
                logits = model.dt_decoder(z["drug"][e[0]], z["gene"][e[1]]).squeeze(-1)
                loss_dt = focal_bce_loss(logits, lbl)

            loss_td = torch.tensor(0.0, device=device)
            if td_edges is not None and td_edges.size(1) > 0:
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

def make_sub_data(data, dt_ei_to_keep, td_ei_to_keep):
    sub = HeteroData()
    sub["drug"].x = data["drug"].x
    sub["gene"].x = data["gene"].x
    sub["disease"].x = data["disease"].x
    sub["gene", "interacts", "gene"].edge_index = data["gene", "interacts", "gene"].edge_index
    if dt_ei_to_keep is not None and dt_ei_to_keep.size(1) > 0:
        sub["drug", "targets", "gene"].edge_index = dt_ei_to_keep
    if td_ei_to_keep is not None and td_ei_to_keep.size(1) > 0:
        sub["gene", "associated_with", "disease"].edge_index = td_ei_to_keep
    return sub.to(device)

def main():
    print("=" * 60)
    print("  GAT Embedding v6 — All-OOF (No Full Model)")
    print("  Positives: OOF from validation fold only")
    print("  Negatives: MEAN across all 5 OOF folds")
    print("  → UNIFORM embedding scale → no PCA magnitude artifact")
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

    dt_pos_indices = full_dt[1].cpu().numpy()
    td_pos_indices = full_td[0].cpu().numpy()

    dt_idx_all = [gene_to_idx[g] for g in dt_genes if g in gene_to_idx]
    dg_idx_all = [gene_to_idx[g] for g in dg_genes if g in gene_to_idx]
    dt_cand = build_neg_candidates(full_ppi, dt_idx_all, n_genes)
    dg_cand = build_neg_candidates(full_ppi, dg_idx_all, n_genes)

    g_dim = meta["gene_dim"]
    d_dim = meta["drug_dim"]
    dd_dim = meta["disease_dim"]

    # ════════════════════════════════════════════════════
    #  STAGE 1: Model A (no-TD) — All-OOF
    #  Positives: OOF from validation fold
    #  Negatives: MEAN across 5 OOF folds
    # ════════════════════════════════════════════════════
    print("\n" + "=" * 60)
    print("  STAGE 1: Model A (no-TD) — All-OOF (no full model)")
    print("  Output: gat_nodg_* → for DG classifier")
    print("=" * 60)

    emb_a_oof = np.zeros((n_genes, OUT_DIM), dtype=np.float32)
    emb_a_count = np.zeros(n_genes, dtype=np.int32)
    dt_pos_set = set(dt_pos_indices)

    kf = KFold(n_splits=N_FOLDS, shuffle=True, random_state=RANDOM_SEED)
    dt_indices = np.arange(full_dt.size(1))

    for fold, (train_idx, val_idx) in enumerate(kf.split(dt_indices)):
        print(f"\n  ── Fold {fold+1}/{N_FOLDS} ──")
        val_edges = full_dt[:, val_idx]
        train_edges = full_dt[:, train_idx]
        val_gene_indices = val_edges[1].cpu().numpy()
        unique_val_genes = np.unique(val_gene_indices)
        print(f"  Train edges: {train_edges.size(1)}  Val edges: {val_edges.size(1)}  Val genes: {len(unique_val_genes)}")

        data_a_fold = make_sub_data(data, train_edges, None)
        dt_sampler = NegSampler(dt_cand, n_genes, RANDOM_SEED + fold)

        set_seed(RANDOM_SEED + fold)
        model = HeteroGATv2(
            metadata=data_a_fold.metadata(),
            drug_dim=d_dim, gene_dim=g_dim, disease_dim=dd_dim,
            hidden_dim=HIDDEN_DIM, out_dim=OUT_DIM, heads=GAT_HEADS,
            dropout=DROPOUT, conv_dropout=CONV_DROPOUT,
        ).to(device)

        train_model(model, data_a_fold, train_edges, None, dt_sampler, None,
                    n_genes, NEG_RATIO, EPOCHS_A, PATIENCE, f"Model-A fold{fold+1}")

        emb_fold = extract_embeddings(model, data_a_fold)

        for g_idx in unique_val_genes:
            if g_idx in dt_pos_set:
                emb_a_oof[g_idx] = emb_fold[g_idx]
                emb_a_count[g_idx] = 1

        neg_mask = np.array([i for i in range(n_genes) if i not in dt_pos_set])
        emb_a_oof[neg_mask] += emb_fold[neg_mask]
        emb_a_count[neg_mask] += 1

        print(f"  Positives (OOF): {len(unique_val_genes)} genes, Negatives (accum): {len(neg_mask)}")

    neg_mean_mask = emb_a_count > 0
    for i in range(n_genes):
        if emb_a_count[i] > 0 and i not in dt_pos_set:
            emb_a_oof[i] /= emb_a_count[i]

    pos_covered = sum(1 for i in dt_pos_set if emb_a_count[i] > 0)
    print(f"\n  [Model-A All-OOF] Positives covered: {pos_covered}/{len(dt_pos_set)}")
    print(f"  Negatives: mean of {emb_a_count[list(set(range(n_genes)) - dt_pos_set).__iter__().__next__() if len(set(range(n_genes)) - dt_pos_set) > 0 else 0]} folds")
    print(f"  Actual negative mean folds: {emb_a_count[list(set(range(n_genes)) - dt_pos_set)][0]}")
    dead = (emb_a_count == 0).sum()
    if dead > 0:
        print(f"  Dead genes (no embedding): {dead}")

    # ════════════════════════════════════════════════════
    #  STAGE 2: Model B (no-DT) — All-OOF
    #  Positives: OOF from validation fold
    #  Negatives: MEAN across 5 OOF folds
    # ════════════════════════════════════════════════════
    print("\n" + "=" * 60)
    print("  STAGE 2: Model B (no-DT) — All-OOF (no full model)")
    print("  Output: gat_nodt_* → for DT classifier")
    print("=" * 60)

    emb_b_oof = np.zeros((n_genes, OUT_DIM), dtype=np.float32)
    emb_b_count = np.zeros(n_genes, dtype=np.int32)
    td_pos_set = set(td_pos_indices)

    td_indices = np.arange(full_td.size(1))

    for fold, (train_idx, val_idx) in enumerate(kf.split(td_indices)):
        print(f"\n  ── Fold {fold+1}/{N_FOLDS} ──")
        val_edges = full_td[:, val_idx]
        train_edges = full_td[:, train_idx]
        val_gene_indices = val_edges[0].cpu().numpy()
        unique_val_genes = np.unique(val_gene_indices)
        print(f"  Train edges: {train_edges.size(1)}  Val edges: {val_edges.size(1)}  Val genes: {len(unique_val_genes)}")

        data_b_fold = make_sub_data(data, None, train_edges)
        td_sampler = NegSampler(dg_cand, n_genes, RANDOM_SEED + fold)

        set_seed(RANDOM_SEED + fold)
        model = HeteroGATv2(
            metadata=data_b_fold.metadata(),
            drug_dim=d_dim, gene_dim=g_dim, disease_dim=dd_dim,
            hidden_dim=HIDDEN_DIM, out_dim=OUT_DIM, heads=GAT_HEADS,
            dropout=DROPOUT, conv_dropout=CONV_DROPOUT,
        ).to(device)

        train_model(model, data_b_fold, None, train_edges, None, td_sampler,
                    n_genes, NEG_RATIO, EPOCHS_B, PATIENCE, f"Model-B fold{fold+1}")

        emb_fold = extract_embeddings(model, data_b_fold)

        for g_idx in unique_val_genes:
            if g_idx in td_pos_set:
                emb_b_oof[g_idx] = emb_fold[g_idx]
                emb_b_count[g_idx] = 1

        neg_mask = np.array([i for i in range(n_genes) if i not in td_pos_set])
        emb_b_oof[neg_mask] += emb_fold[neg_mask]
        emb_b_count[neg_mask] += 1

        print(f"  Positives (OOF): {len(unique_val_genes)} genes, Negatives (accum): {len(neg_mask)}")

    for i in range(n_genes):
        if emb_b_count[i] > 0 and i not in td_pos_set:
            emb_b_oof[i] /= emb_b_count[i]

    pos_covered_b = sum(1 for i in td_pos_set if emb_b_count[i] > 0)
    print(f"\n  [Model-B All-OOF] Positives covered: {pos_covered_b}/{len(td_pos_set)}")
    dead_b = (emb_b_count == 0).sum()
    if dead_b > 0:
        print(f"  Dead genes (no embedding): {dead_b}")

    # ════════════════════════════════════════════════════
    #  Norm Diagnostics
    # ════════════════════════════════════════════════════
    print("\n" + "=" * 60)
    print("  NORM DIAGNOSTICS (should be uniform pos vs neg)")
    print("=" * 60)

    norm_a = np.linalg.norm(emb_a_oof, axis=1)
    norm_b = np.linalg.norm(emb_b_oof, axis=1)

    dt_set_all = {gene_to_idx[g] for g in dt_genes if g in gene_to_idx}
    dg_set_all = {gene_to_idx[g] for g in dg_genes if g in gene_to_idx}

    for name, norm_vec, pos_set in [
        ("gat_nodg_* (Model A)", norm_a, dt_set_all),
        ("gat_nodt_* (Model B)", norm_b, dg_set_all),
    ]:
        pos_mask_arr = np.array([i in pos_set for i in range(n_genes)])
        neg_mask_arr = ~pos_mask_arr
        pn = norm_vec[pos_mask_arr]
        nn = norm_vec[neg_mask_arr]
        pn_valid = pn[pn > 1e-6]
        nn_valid = nn[nn > 1e-6]
        ratio = pn_valid.mean() / max(nn_valid.mean(), 1e-6)
        print(f"\n  {name}:")
        print(f"    Pos norm: mean={pn_valid.mean():.4f} std={pn_valid.std():.4f} (n={len(pn_valid)})")
        print(f"    Neg norm: mean={nn_valid.mean():.4f} std={nn_valid.std():.4f} (n={len(nn_valid)})")
        print(f"    Pos/Neg ratio: {ratio:.2f}x  {'⚠️  LARGE GAP' if ratio > 5 else '✓  UNIFORM'}")

        dead_p = (pn < 1e-6).sum()
        dead_n = (nn < 1e-6).sum()
        if dead_p > 0:
            print(f"    Dead positives: {dead_p}")
        if dead_n > 0:
            print(f"    Dead negatives: {dead_n}")

    # ════════════════════════════════════════════════════
    #  PCA + Merge
    # ════════════════════════════════════════════════════
    print("\n" + "=" * 60)
    print("  PCA + MERGE")
    print("=" * 60)

    var_raw_a = emb_a_oof.var(axis=0).sum()
    var_raw_b = emb_b_oof.var(axis=0).sum()
    print(f"[Raw-Emb-A] no-TD (All-OOF): total_var={var_raw_a:.4f}")
    print(f"[Raw-Emb-B] no-DT (All-OOF): total_var={var_raw_b:.4f}")

    pca_a = PCA(n_components=N_PCA, random_state=RANDOM_SEED)
    emb_a_pca = pca_a.fit_transform(emb_a_oof)
    var_a = pca_a.explained_variance_ratio_.sum()
    print(f"\n[PCA-A] no-TD All-OOF: {OUT_DIM}d → {N_PCA}d, cumulative var={var_a:.4f}")
    for i, v in enumerate(pca_a.explained_variance_ratio_):
        print(f"  PC{i:2d}: {v:.4f}")

    pca_b = PCA(n_components=N_PCA, random_state=RANDOM_SEED)
    emb_b_pca = pca_b.fit_transform(emb_b_oof)
    var_b = pca_b.explained_variance_ratio_.sum()
    print(f"\n[PCA-B] no-DT All-OOF: {OUT_DIM}d → {N_PCA}d, cumulative var={var_b:.4f}")
    for i, v in enumerate(pca_b.explained_variance_ratio_):
        print(f"  PC{i:2d}: {v:.4f}")

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

    print(f"\n[Save] → {out_path}  ({merged.shape[0]} x {merged.shape[1]})")

    # ════════════════════════════════════════════════════
    #  Leakage Check
    # ════════════════════════════════════════════════════
    print("\n[Leak Check] Single-feature AUC (v6 All-OOF):")
    from sklearn.metrics import roc_auc_score
    df_chk = pd.read_csv(out_path)
    y_dt = df_chk["is_drug_target"].values
    y_dg = df_chk["is_disease_gene"].values

    max_auc = {"DT": 0.0, "DG": 0.0}
    for check_name, check_cols, check_label in [
        ("gat_nodg_* vs DT", cols_a, y_dt),
        ("gat_nodt_* vs DG", cols_b, y_dg),
    ]:
        print(f"\n  {check_name}:")
        for c in check_cols:
            x = df_chk[c].fillna(0).values
            a = roc_auc_score(check_label, x)
            if a < 0.5:
                a = roc_auc_score(check_label, -x)
            tag = " ⚠️ LEAK" if a > 0.90 else (" ◆ HIGH" if a > 0.80 else "")
            if a > 0.70:
                print(f"    {c}: AUC={a:.4f}{tag}")
            label_key = "DT" if check_label is y_dt else "DG"
            max_auc[label_key] = max(max_auc[label_key], a)

    print(f"\n  Max single-feature AUC: DT={max_auc['DT']:.4f}  DG={max_auc['DG']:.4f}")
    if max_auc["DT"] < 0.90 and max_auc["DG"] < 0.90:
        print("  ✓ No single-dimension leakage (all < 0.90)")
    else:
        print("  ⚠️ Some dimensions still show high AUC — check above")

    print("\n[Done] GAT v6 All-OOF embeddings — uniform scale, no full model.")

if __name__ == "__main__":
    main()