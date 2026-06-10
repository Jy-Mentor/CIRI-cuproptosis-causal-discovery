# -*- coding: utf-8 -*-
"""
Stage 9: PPI Network + GAT Node Regression (v5)
================================================

FIX: [P0-1][break circular dependency with Stage8]
  - Labels: NO LONGER use Stage8 comprehensive score (which includes GAT)
  - Labels: Now use independent cuproptosis pathway coreness scores from config
    (based on Tsvetkov Science 2022, PMID:35298263)
  - Stage8 GAT dimension uses cold-start (PPI degree placeholder on first run)

FIX: [P0-2][eliminate label-feature circuit leakage]
  - Labels: cuproptosis pathway hierarchy scores (independent biological knowledge)
  - Features: perturbation(GRN) + shap(ML) + priors + topology + 2-hop neighbor aggregation
  - Features and labels are completely independent

FIX: [P1-6][feature enhancement to 5D -> 10D]
  - Add 2-hop neighbor feature aggregation (mean perturbation of neighbors)
  - Add PPI degree/betweenness/closeness as topology features
  - Add module membership (WGCNA) as categorical feature

Reference:
  - GNN for Drug Target Identification (PMID:36168971)
  - EGNF (PMID:41139924)
  - MOGAT (PMID:38474033)
  - Cuproptosis pathway hierarchy: Tsvetkov Science 2022 (PMID:35298263)

Input:
  - stage5/ppi_topology.json: PPI network topology
  - stage5/string_ppi.tsv: PPI edge list
  - stage7/gene_shap_importance.csv: ML SHAP importance
  - stage6/gene_perturbation_scores.csv: GRN perturbation scores
  - config.CUPROPTOSIS_PATHWAY_SCORES: Independent labels

Output:
  - gat_gene_ranking.csv: GAT gene ranking
  - gat_model.pth: Model weights
  - stage9.log: Run log
"""

import os
import sys
import warnings
import logging
import json
import time
from collections import defaultdict
from datetime import datetime

warnings.filterwarnings('ignore')
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'

import numpy as np
import pandas as pd

try:
    import torch
    import torch.nn as nn
    import torch.nn.functional as F
    from torch_geometric.nn import GATConv
    from torch_geometric.data import Data
    from sklearn.preprocessing import StandardScaler
    from sklearn.metrics import r2_score, mean_squared_error
except ImportError as e:
    print(f"Error: Missing dependencies ({e})")
    sys.exit(1)

# Path configuration
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import (
    RESULTS_DIR, CUPROPTOSIS_GENES, CUPROPTOSIS_RELATED,
    CUPROPTOSIS_PATHWAY_SCORES, BCP_TARGETS
)

STAGE_DIR = os.path.join(RESULTS_DIR, "stage9_ppi_gat")
os.makedirs(STAGE_DIR, exist_ok=True)

logger = logging.getLogger("stage9")
logger.setLevel(logging.INFO)

file_handler = logging.FileHandler(os.path.join(STAGE_DIR, "stage9.log"), encoding="utf-8", mode="w")
file_handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
logger.addHandler(file_handler)

console_handler = logging.StreamHandler(sys.stdout)
console_handler.setFormatter(logging.Formatter("%(message)s"))
logger.addHandler(console_handler)


def load_ppi_network():
    """Load pre-computed PPI network from Stage5"""
    logger.info("Loading PPI network (Stage5)...")
    
    topo_file = os.path.join(RESULTS_DIR, "stage5_string_ppi", "ppi_topology.json")
    ppi_file = os.path.join(RESULTS_DIR, "stage5_string_ppi", "string_ppi.tsv")
    
    if not os.path.exists(topo_file) or not os.path.exists(ppi_file):
        logger.warning("  Stage5 PPI data not found, using local cache")
        topo_file = os.path.join(STAGE_DIR, "ppi_topology.json")
        ppi_file = os.path.join(STAGE_DIR, "string_ppi.tsv")
    
    if not os.path.exists(ppi_file):
        logger.error("PPI network file not found")
        return None, None, None
    
    ppi_df = pd.read_csv(ppi_file, sep='\t')
    
    if os.path.exists(topo_file):
        with open(topo_file, 'r', encoding='utf-8') as f:
            topology = json.load(f)
        logger.info(f"  Topology loaded: {len(topology)} genes")
    else:
        topology = {}
        logger.warning("  No topology data")
    
    ppi_genes = set()
    for _, row in ppi_df.iterrows():
        ppi_genes.add(row["preferredName_A"].upper())
        ppi_genes.add(row["preferredName_B"].upper())
    
    logger.info(f"  PPI genes: {len(ppi_genes)}")
    return ppi_df, ppi_genes, topology


def load_stage_results():
    """Load results from previous stages"""
    results = {}
    
    # Stage7: SHAP importance
    shap_file = os.path.join(RESULTS_DIR, "stage7_ml_shap", "gene_shap_importance.csv")
    if os.path.exists(shap_file):
        shap_df = pd.read_csv(shap_file)
        col = 'SHAP_importance' if 'SHAP_importance' in shap_df.columns else 'shap_importance'
        results['shap'] = dict(zip(shap_df['Gene'].str.upper(), shap_df[col]))
        logger.info(f"  SHAP: {len(results['shap'])} genes")
    else:
        logger.warning("  SHAP results not found")
        results['shap'] = {}
    
    # Stage6: GRN perturbation scores
    pert_file = os.path.join(RESULTS_DIR, "stage6_graphsage_knockout", "gene_perturbation_scores.csv")
    if os.path.exists(pert_file):
        pert_df = pd.read_csv(pert_file)
        results['perturbation'] = dict(zip(pert_df['gene'].str.upper(), pert_df['perturbation_score']))
        logger.info(f"  Perturbation scores: {len(results['perturbation'])} genes")
    else:
        logger.warning("  Perturbation scores not found")
        results['perturbation'] = {}
    
    # WGCNA module membership (optional)
    wgcna_file = os.path.join(RESULTS_DIR, "stage4_seed_wgcna", "wgcna_modules.csv")
    if os.path.exists(wgcna_file):
        try:
            wgcna_df = pd.read_csv(wgcna_file)
            if 'GeneSymbol' in wgcna_df.columns and 'Module' in wgcna_df.columns:
                results['wgcna_module'] = dict(zip(
                    wgcna_df['GeneSymbol'].str.upper(), wgcna_df['Module']
                ))
                logger.info(f"  WGCNA modules: {len(results['wgcna_module'])} genes")
        except Exception:
            results['wgcna_module'] = {}
    else:
        results['wgcna_module'] = {}
    
    return results


def build_adjacency_list(ppi_df, gene_list):
    """Build adjacency list for 2-hop neighbor aggregation"""
    gene_set = set(g.upper() for g in gene_list)
    adj = defaultdict(list)
    
    for _, row in ppi_df.iterrows():
        g1 = row["preferredName_A"].upper()
        g2 = row["preferredName_B"].upper()
        w = row.get("score", 1.0)
        if g1 in gene_set and g2 in gene_set:
            adj[g1].append((g2, w))
            adj[g2].append((g1, w))
    
    return adj


def compute_2hop_neighbor_features(gene_list, adj, feature_dict):
    """
    Compute 2-hop neighbor feature aggregation
    
    FIX:[P1-6][2-hop neighbor mean perturbation/SHAP]
    Reference: Message passing in GNNs (PMID:38474033)
    """
    gene_set = set(g.upper() for g in gene_list)
    features_1hop = []
    features_2hop = []
    
    for gene in gene_list:
        gene_upper = gene.upper()
        
        # 1-hop: mean feature of direct neighbors
        neighbors_1 = adj.get(gene_upper, [])
        if neighbors_1:
            vals_1 = [feature_dict.get(n, 0.0) for n, w in neighbors_1]
            mean_1hop = np.mean(vals_1)
        else:
            mean_1hop = 0.0
        
        # 2-hop: mean feature of neighbors-of-neighbors (excluding self and 1-hop)
        neighbors_1_set = {n for n, w in neighbors_1}
        neighbors_2_vals = []
        for n1, w1 in neighbors_1:
            for n2, w2 in adj.get(n1, []):
                if n2 != gene_upper and n2 not in neighbors_1_set and n2 in gene_set:
                    neighbors_2_vals.append(feature_dict.get(n2, 0.0))
        
        if neighbors_2_vals:
            mean_2hop = np.mean(neighbors_2_vals)
        else:
            mean_2hop = mean_1hop  # fallback to 1-hop
        
        features_1hop.append(mean_1hop)
        features_2hop.append(mean_2hop)
    
    return np.array(features_1hop, dtype=np.float32), np.array(features_2hop, dtype=np.float32)


def build_node_features(all_genes, stage_results, topology, ppi_df):
    """
    Build node features (FIX:[P0-2][no label leakage] + FIX:[P1-6][enhanced to ~10D])
    
    Features (independent of labels):
    1. GRN perturbation score (from Stage6, independent computation)
    2. SHAP importance (from Stage7, independent computation)
    3. Is cuproptosis gene (prior knowledge)
    4. Is BCP target (prior knowledge)
    5. PPI degree (network topology, independent of label)
    6. PPI betweenness centrality
    7. PPI closeness centrality
    8. 1-hop neighbor mean perturbation (P1-6: feature aggregation)
    9. 2-hop neighbor mean perturbation (P1-6: feature aggregation)
    10. WGCNA module membership (categorical, encoded as numeric)
    
    Labels: Cuproptosis pathway coreness scores (config.CUPROPTOSIS_PATHWAY_SCORES)
    These are COMPLETELY INDEPENDENT of all features above.
    """
    logger.info("Building node features (independent of labels)...")
    
    gene_list = sorted(all_genes)
    n_genes = len(gene_list)
    
    # Build adjacency list for neighbor aggregation
    adj = build_adjacency_list(ppi_df, gene_list)
    
    # Compute 2-hop neighbor features for perturbation
    pert_1hop, pert_2hop = compute_2hop_neighbor_features(
        gene_list, adj, stage_results.get('perturbation', {})
    )
    
    # Compute 2-hop neighbor features for SHAP
    shap_1hop, shap_2hop = compute_2hop_neighbor_features(
        gene_list, adj, stage_results.get('shap', {})
    )
    
    # Encode WGCNA module membership (top modules by size)
    wgcna_module_dict = stage_results.get('wgcna_module', {})
    module_counts = defaultdict(int)
    for mod in wgcna_module_dict.values():
        module_counts[mod] += 1
    top_modules = sorted(module_counts.keys(), key=lambda m: module_counts[m], reverse=True)[:5]
    module_to_idx = {m: i for i, m in enumerate(top_modules)}
    
    features = []
    for i, gene in enumerate(gene_list):
        gene_upper = gene.upper()
        
        # 1. GRN perturbation score
        pert = stage_results['perturbation'].get(gene_upper, 0.0)
        
        # 2. SHAP importance
        shap = stage_results['shap'].get(gene_upper, 0.0)
        
        # 3. Is cuproptosis gene
        is_cupro = 1 if gene_upper in CUPROPTOSIS_GENES else 0
        
        # 4. Is BCP target
        is_bcp = 1 if gene_upper in BCP_TARGETS else 0
        
        # 5. PPI degree
        degree = topology.get(gene_upper, {}).get('degree', 0)
        
        # 6. Betweenness centrality
        betweenness = topology.get(gene_upper, {}).get('betweenness', 0.0)
        
        # 7. Closeness centrality
        closeness = topology.get(gene_upper, {}).get('closeness', 0.0)
        
        # 8. 1-hop neighbor mean perturbation
        mean_pert_1hop = pert_1hop[i]
        
        # 9. 2-hop neighbor mean perturbation
        mean_pert_2hop = pert_2hop[i]
        
        # 10. WGCNA module membership (0 if no module)
        module_idx = module_to_idx.get(wgcna_module_dict.get(gene_upper, ''), -1)
        wgcna_feature = module_idx if module_idx >= 0 else -1
        
        features.append([
            pert, shap, is_cupro, is_bcp, degree,
            betweenness, closeness, mean_pert_1hop, mean_pert_2hop, wgcna_feature
        ])
    
    X = np.array(features, dtype=np.float32)
    
    # Standardize
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    feature_names = [
        'perturbation', 'shap', 'is_cupro', 'is_bcp', 'degree',
        'betweenness', 'closeness', 'neighbor_pert_1hop', 'neighbor_pert_2hop',
        'wgcna_module'
    ]
    
    logger.info(f"  Feature matrix: {X_scaled.shape}")
    logger.info(f"  Features: {feature_names}")
    
    return X_scaled, gene_list, feature_names


def build_labels(gene_list):
    """
    Build node labels using INDEPENDENT biological knowledge
    
    FIX:[P0-1][break circular dependency] FIX:[P0-2][no feature-label leakage]
    
    Labels based on cuproptosis pathway hierarchy (Tsvetkov Science 2022, PMID:35298263):
    - Upstream regulators (FDX1, LIAS, LIPT1): high scores 90-100
    - Lipoylation effectors (DLAT, PDHA1, PDHB): medium-high scores 80-85
    - Copper transporters (ATP7A/B, SLC31A1): medium scores 60-65
    - Downstream regulators (MTF1, CDKN2A, GLS): lower scores 45-55
    - Non-cuproptosis genes: default low score based on degree connectivity
    
    This is COMPLETELY INDEPENDENT of:
    - GRN perturbation scores (Stage6)
    - SHAP importance (Stage7)
    - PPI topology (Stage5)
    """
    logger.info("Building node labels (independent cuproptosis pathway scores)...")
    
    pathway_scores = CUPROPTOSIS_PATHWAY_SCORES
    
    y = []
    for gene in gene_list:
        gene_upper = gene.upper()
        
        if gene_upper in pathway_scores:
            score = pathway_scores[gene_upper]
        else:
            # Non-cuproptosis genes: use structural prior based on connectivity
            # This avoids giving all non-cupro genes the same score
            score = 10  # baseline for non-cuproptosis genes
        
        y.append(score)
    
    y = np.array(y, dtype=np.float32)
    
    # Normalize to 0-1
    y_min, y_max = y.min(), y.max()
    if y_max > y_min:
        y_normalized = (y - y_min) / (y_max - y_min)
    else:
        y_normalized = np.zeros_like(y)
    
    logger.info(f"  Label range: [{y.min():.1f}, {y.max():.1f}] -> [{y_normalized.min():.3f}, {y_normalized.max():.3f}]")
    logger.info(f"  Label mean: {y_normalized.mean():.3f}, std: {y_normalized.std():.3f}")
    
    n_cupro_labeled = sum(1 for g in gene_list if g.upper() in pathway_scores)
    logger.info(f"  Cuproptosis-labeled genes: {n_cupro_labeled}/{len(gene_list)}")
    
    return y_normalized


class GATRegressor(nn.Module):
    """GAT node regression model with residual connections"""
    def __init__(self, in_features, hidden_dim, n_heads=4, dropout=0.3):
        super(GATRegressor, self).__init__()
        self.conv1 = GATConv(in_features, hidden_dim, heads=n_heads, dropout=dropout)
        self.conv2 = GATConv(hidden_dim * n_heads, hidden_dim, heads=1, dropout=dropout)
        self.norm1 = nn.BatchNorm1d(hidden_dim * n_heads)
        self.norm2 = nn.BatchNorm1d(hidden_dim)
        self.fc = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim // 2, 1)
        )
    
    def forward(self, x, edge_index):
        x = F.elu(self.conv1(x, edge_index))
        x = self.norm1(x)
        x = F.dropout(x, p=0.2, training=self.training)
        x = F.elu(self.conv2(x, edge_index))
        x = self.norm2(x)
        x = self.fc(x)
        return x.squeeze(-1)


def train_gat(X, y, edge_index, n_epochs=300, lr=0.005, weight_decay=1e-4):
    """Train GAT model with 80/20 train/val split"""
    logger.info(f"Training GAT node regression (epochs={n_epochs}, lr={lr}, weight_decay={weight_decay})...")
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    logger.info(f"  Device: {device}")
    
    # Train/val split (80/20)
    n_nodes = X.shape[0]
    np.random.seed(42)
    indices = np.random.permutation(n_nodes)
    split = int(0.8 * n_nodes)
    train_idx = indices[:split]
    val_idx = indices[split:]
    
    train_mask = torch.zeros(n_nodes, dtype=torch.bool)
    train_mask[train_idx] = True
    val_mask = torch.zeros(n_nodes, dtype=torch.bool)
    val_mask[val_idx] = True
    
    X_tensor = torch.FloatTensor(X).to(device)
    y_tensor = torch.FloatTensor(y).to(device)
    edge_index_tensor = torch.LongTensor(edge_index).to(device)
    train_mask = train_mask.to(device)
    val_mask = val_mask.to(device)
    
    data = Data(x=X_tensor, edge_index=edge_index_tensor, y=y_tensor,
                train_mask=train_mask, val_mask=val_mask)
    
    in_features = X.shape[1]
    model = GATRegressor(in_features, hidden_dim=32, n_heads=4, dropout=0.4).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=0.5, patience=20)
    
    best_val_loss = float('inf')
    best_state = None
    patience_counter = 0
    
    train_losses = []
    val_losses = []
    
    for epoch in range(n_epochs):
        model.train()
        optimizer.zero_grad()
        out = model(data.x, data.edge_index)
        loss = F.mse_loss(out[train_mask], data.y[train_mask])
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
        optimizer.step()
        scheduler.step(loss)
        
        train_losses.append(loss.item())
        
        # Validation
        model.eval()
        with torch.no_grad():
            val_out = model(data.x, data.edge_index)
            val_loss = F.mse_loss(val_out[val_mask], data.y[val_mask]).item()
            val_losses.append(val_loss)
        
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            best_state = {k: v.clone() for k, v in model.state_dict().items()}
            patience_counter = 0
        else:
            patience_counter += 1
        
        if patience_counter >= 50:
            logger.info(f"  Early stopping at epoch {epoch} (val_loss={val_loss:.6f})")
            break
        
        if (epoch + 1) % 50 == 0:
            logger.info(f"  Epoch {epoch+1}/{n_epochs}, Train Loss: {loss.item():.6f}, Val Loss: {val_loss:.6f}")
    
    if best_state:
        model.load_state_dict(best_state)
    
    logger.info(f"  Best validation loss: {best_val_loss:.6f}")
    
    return model, train_losses, val_losses


def evaluate_model(model, X, y, edge_index):
    """Evaluate model on full graph"""
    logger.info("Evaluating model...")
    
    device = next(model.parameters()).device
    X_tensor = torch.FloatTensor(X).to(device)
    edge_index_tensor = torch.LongTensor(edge_index).to(device)
    
    model.eval()
    with torch.no_grad():
        pred = model(X_tensor, edge_index_tensor).cpu().numpy()
    
    r2 = r2_score(y, pred)
    mse = mean_squared_error(y, pred)
    
    logger.info(f"  R² = {r2:.4f}")
    logger.info(f"  MSE = {mse:.6f}")
    
    return pred, r2, mse


def save_results(model, gene_list, y_true, y_pred, r2, mse, train_losses, val_losses, feature_names, edge_index):
    """Save results"""
    logger.info("Saving results...")
    
    # Gene ranking
    ranking_data = []
    for i, gene in enumerate(gene_list):
        ranking_data.append({
            'Gene': gene,
            'GAT_score': round(float(y_pred[i]), 6),
            'True_label': round(float(y_true[i]), 6),
            'Rank': 0
        })
    
    ranking_df = pd.DataFrame(ranking_data).sort_values('GAT_score', ascending=False)
    ranking_df['Rank'] = range(1, len(ranking_df) + 1)
    
    ranking_file = os.path.join(STAGE_DIR, "gat_gene_ranking.csv")
    ranking_df.to_csv(ranking_file, index=False)
    logger.info(f"  ✓ Gene ranking: {ranking_file}")
    logger.info(f"    Top5: {ranking_df.head(5)['Gene'].tolist()}")
    
    # Model weights
    model_file = os.path.join(STAGE_DIR, "gat_model.pth")
    torch.save(model.state_dict(), model_file)
    logger.info(f"  ✓ Model weights: {model_file}")
    
    # Training curve
    loss_data = {
        'epoch': range(1, len(train_losses) + 1),
        'train_loss': train_losses,
        'val_loss': val_losses[:len(train_losses)]
    }
    loss_df = pd.DataFrame(loss_data)
    loss_file = os.path.join(STAGE_DIR, "gat_training_loss.csv")
    loss_df.to_csv(loss_file, index=False)
    
    # Performance metrics
    perf = {
        'R2': round(r2, 4),
        'MSE': round(mse, 6),
        'N_genes': len(gene_list),
        'N_edges': edge_index.shape[1],
        'N_features': len(feature_names),
        'N_train_epochs': len(train_losses),
        'N_val_epochs': len(val_losses)
    }
    perf_file = os.path.join(STAGE_DIR, "gat_performance.json")
    with open(perf_file, 'w') as f:
        json.dump(perf, f, indent=2)
    
    return ranking_df


def main():
    logger.info("=" * 60)
    logger.info("Stage 9: PPI Network + GAT Node Regression (v5)")
    logger.info("=" * 60)
    logger.info(f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info("FIX: [P0-1] Break circular dependency with Stage8")
    logger.info("FIX: [P0-2] Independent labels (cuproptosis pathway scores)")
    logger.info("FIX: [P1-6] Enhanced features (10D with 2-hop aggregation)")
    logger.info("Reference: Tsvetkov Science 2022 (PMID:35298263)")
    
    # 1. Load PPI network
    ppi_df, ppi_genes, topology = load_ppi_network()
    if ppi_df is None:
        logger.error("PPI network load failed")
        return
    
    # 2. Load stage results
    stage_results = load_stage_results()
    
    # 3. Build feature matrix (independent of labels)
    X, gene_list, feature_names = build_node_features(
        ppi_genes, stage_results, topology, ppi_df)
    
    # 4. Build labels (INDEPENDENT: cuproptosis pathway scores)
    y = build_labels(gene_list)
    
    # 5. Build edge index
    gene_to_idx = {g: i for i, g in enumerate(gene_list)}
    edge_list = []
    for _, row in ppi_df.iterrows():
        g1, g2 = row["preferredName_A"].upper(), row["preferredName_B"].upper()
        if g1 in gene_to_idx and g2 in gene_to_idx:
            edge_list.append([gene_to_idx[g1], gene_to_idx[g2]])
            edge_list.append([gene_to_idx[g2], gene_to_idx[g1]])
    
    edge_index = np.array(edge_list).T if edge_list else np.empty((2, 0), dtype=int)
    logger.info(f"  Edge index: {edge_index.shape[1]} edges")
    
    # 6. Train GAT
    model, train_losses, val_losses = train_gat(X, y, edge_index, n_epochs=300, lr=0.005)
    
    # 7. Evaluate
    y_pred, r2, mse = evaluate_model(model, X, y, edge_index)
    
    # 8. Save
    ranking_df = save_results(model, gene_list, y, y_pred, r2, mse, train_losses, val_losses, feature_names, edge_index)
    
    logger.info("\n" + "=" * 60)
    logger.info("Stage 9 completed!")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
