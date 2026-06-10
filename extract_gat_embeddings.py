# -*- coding: utf-8 -*-
"""
Extract GAT 64-dim gene embeddings (z_dict["gene"]) from trained HeteroGAT model.
Output: gat_emb_pca_0 ... gat_emb_pca_14 (15 dims, PCA-reduced)
"""

import os, sys, warnings
import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from torch_geometric.data import HeteroData
from sklearn.decomposition import PCA
from sklearn.model_selection import KFold

warnings.filterwarnings("ignore")

# ── Paths ──
GAT_DIR = r"C:\Users\Jy-Mentor-7\Desktop\GAT"
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RF_DIR = r"C:\Users\Jy-Mentor-7\Desktop\随机森林"
FEATURE_TABLE = os.path.join(RF_DIR, "gene_features_table_with_gat.csv")
OUT_EMBEDDINGS = os.path.join(RF_DIR, "gat_gene_embeddings_64d.npy")
OUT_GENE_LIST = os.path.join(RF_DIR, "gat_gene_list_emb.npy")

os.makedirs(RF_DIR, exist_ok=True)

# ── Model config (must match gat_hetero_link_prediction.py) ──
HIDDEN_DIM = 128
OUT_DIM = 64
GAT_HEADS = 4
DROPOUT = 0.3
LR = 1e-3
WEIGHT_DECAY = 1e-5
EPOCHS = 300
PATIENCE = 40
PPI_THRESHOLD = 700
ENSEMBLE_RUNS = 3
N_CV_FOLDS = 5
RANDOM_SEED = 42

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"[Init] Device: {device}")

import random
def set_seed(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)

set_seed(RANDOM_SEED)

# ── 1. Data Loading ──

def load_drug_targets(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Not found: {path}")
    for enc in ("utf-8", "gbk", "latin1"):
        try:
            with open(path, "r", encoding=enc) as fh:
                lines = [ln.strip() for ln in fh if ln.strip()]
            break
        except UnicodeDecodeError:
            continue
    return list(dict.fromkeys([g.upper() for g in lines]))

def load_disease_genes(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Not found: {path}")
    for enc in ("utf-8", "gbk", "latin1"):
        try:
            with open(path, "r", encoding=enc) as fh:
                lines = [ln.strip() for ln in fh if ln.strip()]
            break
        except UnicodeDecodeError:
            continue
    if "," in lines[0] or "\t" in lines[0]:
        genes = []
        for ln in lines:
            parts = ln.replace("\t", ",").split(",")
            gene = parts[0].strip().upper()
            if gene and all(c.isalnum() or c == "-" for c in gene):
                genes.append(gene)
        return list(dict.fromkeys(genes))
    else:
        return list(dict.fromkeys([g.upper() for g in lines if g]))

def load_ppi(path):
    ppi_edges = []
    for enc in ("utf-8", "gbk", "latin1"):
        try:
            with open(path, "r", encoding=enc) as fh:
                for line in fh:
                    line = line.strip()
                    if not line or "protein" in line.lower():
                        continue
                    parts = line.replace("\t", ",").split(",")
                    if len(parts) >= 3:
                        gene_a = parts[0].strip().upper()
                        gene_b = parts[1].strip().upper()
                        try:
                            score = int(parts[2].strip())
                        except ValueError:
                            continue
                        if score >= PPI_THRESHOLD:
                            ppi_edges.append((gene_a, gene_b))
            break
        except UnicodeDecodeError:
            continue
    print(f"[Load] PPI edges: {len(ppi_edges)} (score >= {PPI_THRESHOLD})")
    return ppi_edges

def load_gene_features(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Not found: {path}")
    df = pd.read_csv(path)
    first_col = df.columns[0]
    if first_col.lower() in ("gene", "gene_symbol", "protein"):
        gene_names = df[first_col].str.upper().tolist()
        feat = df.iloc[:, 1:].values.astype(np.float32)
    else:
        gene_names = [f"GENE_{i}" for i in range(len(df))]
        feat = df.values.astype(np.float32)
    print(f"[Load] Gene features: {len(gene_names)} x {feat.shape[1]}")
    return feat, gene_names

def load_drug_fingerprint(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Not found: {path}")
    df = pd.read_csv(path)
    if df.shape[0] == 1:
        arr = df.values.astype(np.float32)
    elif df.shape[1] == 1:
        arr = df.T.values.astype(np.float32)
    else:
        arr = df.iloc[0:1, :].values.astype(np.float32)
    print(f"[Load] Drug fingerprint: {arr.shape}")
    return arr

def load_all_genes(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Not found: {path}")
    for enc in ("utf-8", "gbk", "latin1"):
        try:
            with open(path, "r", encoding=enc) as fh:
                lines = [ln.strip() for ln in fh if ln.strip()]
            if lines[0].upper().replace('"', '') in ("GENE", "GENE_SYMBOL", "GENES"):
                lines = lines[1:]
            break
        except UnicodeDecodeError:
            continue
    return list(dict.fromkeys([g.upper() for g in lines]))


# ── 2. Build HeteroGraph ──

def build_graph(drug_target_genes, disease_genes_raw, ppi_edges,
                gene_features_arr, gene_feature_names, drug_fp, all_genes):

    gene_feat_set = set(gene_feature_names)
    all_gene_candidates = set(all_genes) if all_genes else set()
    all_gene_candidates |= gene_feat_set | set(drug_target_genes) | set(disease_genes_raw)
    for a, b in ppi_edges:
        all_gene_candidates.add(a)
        all_gene_candidates.add(b)

    all_genes_list = sorted(all_gene_candidates & gene_feat_set)
    if len(all_genes_list) == 0:
        raise ValueError("No genes with features!")

    gene_to_idx = {g: i for i, g in enumerate(all_genes_list)}
    n_genes = len(all_genes_list)
    feat_dim = gene_features_arr.shape[1]
    print(f"[Build] Genes: {n_genes}, feat_dim: {feat_dim}")

    gene_feat_dict = dict(zip(gene_feature_names, gene_features_arr))
    gene_feat = np.zeros((n_genes, feat_dim), dtype=np.float32)
    for i, g in enumerate(all_genes_list):
        if g in gene_feat_dict:
            gene_feat[i] = gene_feat_dict[g]

    drug_feat = drug_fp.reshape(1, -1).astype(np.float32)

    disease_in_graph = [g for g in disease_genes_raw if g in gene_to_idx]
    if disease_in_graph:
        d_idx = [gene_to_idx[g] for g in disease_in_graph]
        disease_feat = np.mean(gene_feat[d_idx], axis=0, keepdims=True)
    else:
        disease_feat = np.mean(gene_feat, axis=0, keepdims=True)

    data = HeteroData()
    data["drug"].x = torch.from_numpy(drug_feat)
    data["gene"].x = torch.from_numpy(gene_feat)
    data["disease"].x = torch.from_numpy(disease_feat)

    dt_src, dt_dst = [], []
    for g in drug_target_genes:
        if g in gene_to_idx:
            dt_src.append(0)
            dt_dst.append(gene_to_idx[g])
    data["drug", "targets", "gene"].edge_index = torch.tensor([dt_src, dt_dst], dtype=torch.long)

    td_src, td_dst = [], []
    for g in disease_genes_raw:
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

    print(f"[Build] DT:{len(dt_src)} TD:{len(td_src)} PPI:{len(ppi_src)}")
    return data, gene_to_idx, all_genes_list


# ── 3. HeteroGAT Model ──

from torch_geometric.nn import GATConv, HeteroConv

class HeteroGAT(torch.nn.Module):
    def __init__(self, metadata, drug_dim, gene_dim, disease_dim,
                 hidden_dim, out_dim, heads, dropout):
        super().__init__()
        self.dropout_rate = dropout
        self.drug_proj = torch.nn.Linear(drug_dim, hidden_dim)
        self.gene_proj = torch.nn.Linear(gene_dim, hidden_dim)
        self.disease_proj = torch.nn.Linear(disease_dim, hidden_dim)

        self.conv1 = HeteroConv({
            ("drug", "targets", "gene"):
                GATConv((-1, -1), hidden_dim // heads, heads=heads, concat=True,
                        dropout=dropout, add_self_loops=False),
            ("gene", "associated_with", "disease"):
                GATConv((-1, -1), hidden_dim // heads, heads=heads, concat=True,
                        dropout=dropout, add_self_loops=False),
            ("gene", "interacts", "gene"):
                GATConv((-1, -1), hidden_dim // heads, heads=heads, concat=True,
                        dropout=dropout, add_self_loops=True),
        }, aggr="mean")

        self.conv2 = HeteroConv({
            ("drug", "targets", "gene"):
                GATConv((-1, -1), out_dim, heads=1, concat=False,
                        dropout=dropout, add_self_loops=False),
            ("gene", "associated_with", "disease"):
                GATConv((-1, -1), out_dim, heads=1, concat=False,
                        dropout=dropout, add_self_loops=False),
            ("gene", "interacts", "gene"):
                GATConv((-1, -1), out_dim, heads=1, concat=False,
                        dropout=dropout, add_self_loops=True),
        }, aggr="mean")

        self.dt_decoder = torch.nn.Bilinear(out_dim, out_dim, 1)
        self.td_decoder = torch.nn.Bilinear(out_dim, out_dim, 1)
        self.drug_out = torch.nn.Linear(hidden_dim, out_dim)

    def forward(self, x_dict, edge_index_dict):
        x_dict = {
            "drug": self.drug_proj(x_dict["drug"]),
            "gene": self.gene_proj(x_dict["gene"]),
            "disease": self.disease_proj(x_dict["disease"]),
        }
        out1 = self.conv1(x_dict, edge_index_dict)
        for k in x_dict:
            if k not in out1:
                out1[k] = x_dict[k]
        out1 = {k: v.relu() for k, v in out1.items()}
        out1 = {k: F.dropout(v, p=self.dropout_rate, training=self.training)
                for k, v in out1.items()}
        out2 = self.conv2(out1, edge_index_dict)
        for k in out1:
            if k not in out2:
                out2[k] = out1[k]
        out2["drug"] = self.drug_out(out2["drug"])
        return out2

    def decode_drug_target(self, z_dict, edge_index):
        return self.dt_decoder(
            z_dict["drug"][edge_index[0]],
            z_dict["gene"][edge_index[1]],
        ).squeeze(-1)

    def decode_target_disease(self, z_dict, edge_index):
        return self.td_decoder(
            z_dict["gene"][edge_index[0]],
            z_dict["disease"][edge_index[1]],
        ).squeeze(-1)


# ── 4. Negative Sampling ──

class NegSampler:
    def __init__(self, candidates, n_genes, seed=42):
        self.candidates = torch.from_numpy(candidates)
        self.n_genes = n_genes
        self.rng = torch.Generator()
        self.rng.manual_seed(seed)

    def sample(self, pos_edge, is_dt):
        n = pos_edge.size(1)
        idx = torch.randint(0, len(self.candidates), (n,), generator=self.rng)
        neg = self.candidates[idx].to(pos_edge.device)
        if is_dt:
            return torch.stack([pos_edge[0], neg])
        else:
            return torch.stack([neg, pos_edge[1]])


# ── 5. PPI Neighbor Set ──

def build_ppi_neighbors(ppi_ei, n_genes):
    neighbors = {i: set() for i in range(n_genes)}
    src = ppi_ei[0].tolist()
    dst = ppi_ei[1].tolist()
    for s, d in zip(src, dst):
        neighbors[s].add(d)
        neighbors[d].add(s)
    return neighbors


# ── 6. Training ──

def train_model(model, data, full_dt, full_td, dt_train, td_train,
                dt_sampler, td_sampler, n_genes, epochs=300, patience=40):
    opt = torch.optim.Adam(model.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)
    best_loss, wait = float("inf"), 0

    for ep in range(epochs):
        model.train()

        pos_dt = full_dt[:, dt_train]
        neg_dt = dt_sampler.sample(pos_dt, is_dt=True)
        dt_edge = torch.cat([pos_dt, neg_dt], dim=1).to(device)
        dt_label = torch.cat([torch.ones(pos_dt.size(1)), torch.zeros(neg_dt.size(1))]).to(device)

        pos_td = full_td[:, td_train]
        neg_td = td_sampler.sample(pos_td, is_dt=False)
        td_edge = torch.cat([pos_td, neg_td], dim=1).to(device)
        td_label = torch.cat([torch.ones(pos_td.size(1)), torch.zeros(neg_td.size(1))]).to(device)

        z = model(data.x_dict, data.edge_index_dict)

        dt_logit = model.decode_drug_target(z, dt_edge)
        td_logit = model.decode_target_disease(z, td_edge)

        loss = F.binary_cross_entropy_with_logits(dt_logit, dt_label) + \
               F.binary_cross_entropy_with_logits(td_logit, td_label)

        opt.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        opt.step()

        if ep % 50 == 0:
            print(f"  Epoch {ep:3d}: loss={loss.item():.4f}")

        if loss.item() < best_loss - 1e-4:
            best_loss = loss.item()
            wait = 0
        else:
            wait += 1
            if wait >= patience:
                print(f"  Early stop at {ep}")
                break


# ── 7. Main Extraction ──

def main():
    print("=" * 60)
    print("  GAT Gene Embedding Extractor (Subgraph OOF-CV)")
    print("=" * 60)

    drug_target_genes = load_drug_targets(os.path.join(GAT_DIR, "drug_targets.txt"))
    disease_genes_raw = load_disease_genes(os.path.join(GAT_DIR, "disease_genes.txt"))
    ppi_edges = load_ppi(os.path.join(GAT_DIR, "ppi_subgraph.csv"))
    gene_feat, gene_names = load_gene_features(os.path.join(SCRIPT_DIR, "subgraph_embeddings.csv"))
    drug_fp = load_drug_fingerprint(os.path.join(GAT_DIR, "drug_fingerprint.csv"))
    all_genes = load_all_genes(os.path.join(GAT_DIR, "subgraph_genes.txt"))

    data, gene_to_idx, gene_list = build_graph(
        drug_target_genes, disease_genes_raw, ppi_edges,
        gene_feat, gene_names, drug_fp, all_genes,
    )
    n_genes = len(gene_list)
    data = data.to(device)

    full_dt = data["drug", "targets", "gene"].edge_index
    full_td = data["gene", "associated_with", "disease"].edge_index
    full_ppi = data["gene", "interacts", "gene"].edge_index
    ppi_src = full_ppi[0].cpu().numpy()
    ppi_dst = full_ppi[1].cpu().numpy()

    ppi_nb = build_ppi_neighbors(full_ppi, n_genes)
    dt_target_idx = [gene_to_idx[g] for g in drug_target_genes if g in gene_to_idx]
    td_disease_idx = [gene_to_idx[g] for g in disease_genes_raw if g in gene_to_idx]

    dt_cand = set()
    for idx in dt_target_idx: dt_cand.update(ppi_nb.get(idx, set()))
    dt_cand -= set(dt_target_idx)
    dt_cand_arr = np.array(sorted(dt_cand), dtype=np.int64)

    td_cand = set()
    for idx in td_disease_idx: td_cand.update(ppi_nb.get(idx, set()))
    td_cand -= set(td_disease_idx)
    td_cand_arr = np.array(sorted(td_cand), dtype=np.int64)

    dt_gene_per_edge = full_dt[1].cpu().numpy()
    td_gene_per_edge = full_td[0].cpu().numpy()

    all_labeled_genes = sorted(set(dt_target_idx + td_disease_idx))
    print(f"[CV] Labeled genes: {len(all_labeled_genes)}, folds: {N_CV_FOLDS}, ensemble: {ENSEMBLE_RUNS}")
    print(f"[CV] Total models to train: {ENSEMBLE_RUNS * N_CV_FOLDS}")

    oof_embeddings = np.zeros((n_genes, OUT_DIM), dtype=np.float32)
    oof_counts = np.zeros(n_genes, dtype=np.int32)

    for ensemble_run in range(ENSEMBLE_RUNS):
        rs = RANDOM_SEED + ensemble_run
        print(f"\n{'='*50}")
        print(f"  Ensemble Run {ensemble_run+1}/{ENSEMBLE_RUNS} (seed={rs})")
        print(f"{'='*50}")

        kf = KFold(n_splits=N_CV_FOLDS, shuffle=True, random_state=rs)

        for fold_i, (train_pos, test_pos) in enumerate(kf.split(all_labeled_genes)):
            test_gene_set = set(all_labeled_genes[i] for i in test_pos)
            print(f"\n  Fold {fold_i+1}/{N_CV_FOLDS}: train_g={len(train_pos)}, test_g={len(test_pos)}")

            dt_train_mask = np.array([g not in test_gene_set for g in dt_gene_per_edge], dtype=bool)
            td_train_mask = np.array([g not in test_gene_set for g in td_gene_per_edge], dtype=bool)
            ppi_train_mask = np.array([s not in test_gene_set and d not in test_gene_set
                                       for s, d in zip(ppi_src, ppi_dst)], dtype=bool)

            dt_train = torch.from_numpy(np.where(dt_train_mask)[0])
            td_train = torch.from_numpy(np.where(td_train_mask)[0])
            ppi_train = torch.from_numpy(np.where(ppi_train_mask)[0])

            n_ppi_removed = len(ppi_train_mask) - ppi_train_mask.sum()
            print(f"    DT: {len(dt_train)}/{len(dt_train_mask)}, "
                  f"TD: {len(td_train)}/{len(td_train_mask)}, "
                  f"PPI: {len(ppi_train)}/{len(ppi_train_mask)} (-{n_ppi_removed})")

            fold_data = HeteroData()
            fold_data["drug"].x = data["drug"].x
            fold_data["gene"].x = data["gene"].x
            fold_data["disease"].x = data["disease"].x
            fold_data["drug", "targets", "gene"].edge_index = full_dt[:, dt_train]
            fold_data["gene", "associated_with", "disease"].edge_index = full_td[:, td_train]
            fold_data["gene", "interacts", "gene"].edge_index = full_ppi[:, ppi_train]
            fold_data = fold_data.to(device)

            dt_cand_filt = np.array([c for c in dt_cand_arr if c not in test_gene_set], dtype=np.int64)
            td_cand_filt = np.array([c for c in td_cand_arr if c not in test_gene_set], dtype=np.int64)

            dt_sampler = NegSampler(dt_cand_filt, n_genes, rs + fold_i)
            td_sampler = NegSampler(td_cand_filt, n_genes, rs + fold_i)

            set_seed(rs + fold_i)
            model = HeteroGAT(
                metadata=fold_data.metadata(),
                drug_dim=data["drug"].x.size(1),
                gene_dim=data["gene"].x.size(1),
                disease_dim=data["disease"].x.size(1),
                hidden_dim=HIDDEN_DIM, out_dim=OUT_DIM,
                heads=GAT_HEADS, dropout=DROPOUT,
            ).to(device)

            f_dt = fold_data["drug", "targets", "gene"].edge_index
            f_td = fold_data["gene", "associated_with", "disease"].edge_index
            dt_all = torch.arange(f_dt.size(1))
            td_all = torch.arange(f_td.size(1))

            train_model(model, fold_data, f_dt, f_td, dt_all, td_all,
                        dt_sampler, td_sampler, n_genes, epochs=EPOCHS, patience=PATIENCE)

            model.eval()
            with torch.no_grad():
                z = model(fold_data.x_dict, fold_data.edge_index_dict)
                emb = z["gene"].cpu().numpy()

            for g_idx in test_gene_set:
                oof_embeddings[g_idx] += emb[g_idx]
                oof_counts[g_idx] += 1

    for i in range(n_genes):
        if oof_counts[i] > 0:
            oof_embeddings[i] /= oof_counts[i]

    n_oof = int((oof_counts > 0).sum())
    print(f"\n{'='*50}")
    print(f"[OOF] Labeled covered: {int((oof_counts[all_labeled_genes] > 0).sum())}/{len(all_labeled_genes)}")

    # Train final model on FULL graph for unlabeled gene embeddings (no label leakage possible)
    unlabeled_mask = oof_counts == 0
    n_unlabeled = int(unlabeled_mask.sum())
    if n_unlabeled > 0:
        print(f"\n  [FullGraph] Training for {n_unlabeled} unlabeled genes ...")
        set_seed(RANDOM_SEED + 999)
        dt_full = torch.arange(full_dt.size(1))
        td_full = torch.arange(full_td.size(1))
        dt_sampler_f = NegSampler(dt_cand_arr, n_genes, RANDOM_SEED + 999)
        td_sampler_f = NegSampler(td_cand_arr, n_genes, RANDOM_SEED + 999)

        full_model = HeteroGAT(
            metadata=data.metadata(),
            drug_dim=data["drug"].x.size(1),
            gene_dim=data["gene"].x.size(1),
            disease_dim=data["disease"].x.size(1),
            hidden_dim=HIDDEN_DIM, out_dim=OUT_DIM,
            heads=GAT_HEADS, dropout=DROPOUT,
        ).to(device)

        train_model(full_model, data, full_dt, full_td, dt_full, td_full,
                    dt_sampler_f, td_sampler_f, n_genes, epochs=EPOCHS, patience=PATIENCE)

        full_model.eval()
        with torch.no_grad():
            z = full_model(data.x_dict, data.edge_index_dict)
            full_emb = z["gene"].cpu().numpy()
        oof_embeddings[unlabeled_mask] = full_emb[unlabeled_mask]
        print(f"  [FullGraph] Filled {n_unlabeled} unlabeled gene embeddings")

    n_final = int((oof_embeddings.sum(axis=1) != 0).sum())
    print(f"[Final] Non-zero embeddings: {n_final}/{n_genes}")

    emb_mean = oof_embeddings.astype(np.float32)

    np.save(OUT_EMBEDDINGS, emb_mean)
    np.save(OUT_GENE_LIST, np.array(gene_list))
    print(f"[Save] Raw → {OUT_EMBEDDINGS}")

    pca = PCA(n_components=15, random_state=RANDOM_SEED)
    emb_pca = pca.fit_transform(emb_mean)
    var_explained = pca.explained_variance_ratio_.sum()
    print(f"[PCA] 64d → 15d, var={var_explained:.3f}")

    pca_cols = [f"gat_emb_pca_{i}" for i in range(15)]
    df_emb = pd.DataFrame(emb_pca, columns=pca_cols, index=gene_list)
    df_emb.index.name = "gene_symbol"
    df_emb = df_emb.reset_index()
    df_emb["gene_symbol"] = df_emb["gene_symbol"].str.upper()

    ft = pd.read_csv(FEATURE_TABLE)
    print(f"[Merge] FT: {ft.shape[0]} × {ft.shape[1]}, GAT: {df_emb.shape[0]} × {df_emb.shape[1]}")

    ft["gene_symbol"] = ft["gene_symbol"].str.upper()
    for c in pca_cols:
        if c in ft.columns:
            ft = ft.drop(columns=[c])
    merged = ft.merge(df_emb, on="gene_symbol", how="left")
    for c in pca_cols:
        merged[c] = merged[c].fillna(0.0)

    n_matched = (merged[pca_cols].sum(axis=1) != 0).sum()
    print(f"[Merge] {n_matched}/{merged.shape[0]} genes with GAT embedding")

    out_path = os.path.join(RF_DIR, "gene_features_table_with_gat_emb.csv")
    merged.to_csv(out_path, index=False, encoding="utf-8-sig")
    print(f"[Save] → {out_path}  ({merged.shape[0]} × {merged.shape[1]})")

    print("\nPCA variance:")
    for i, vr in enumerate(pca.explained_variance_ratio_):
        print(f"  gat_emb_pca_{i:2d}: {vr:.4f}")

    return merged


if __name__ == "__main__":
    main()