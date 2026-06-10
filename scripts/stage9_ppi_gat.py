# -*- coding: utf-8 -*-
"""
Stage 9: PPI Network + GAT Node Regression (v9)
================================================
FIX: [v9][无偏泛靶点发现管线 + 铜死亡模块整合]
  - Labels: 使用 Stage8 纯数据驱动综合得分 (GRN perturbation + ML SHAP + PPI topology)
  - Stage8 v9 不含 GAT 维度 → 无循环依赖 (cold-start resolved)
  - 删除 is_cupro / is_bcp 特征 (10D → 8D)
  - 整合铜死亡模块分析结果作为验证层

FIX: [P1-6][feature enhancement to 5D -> 10D]
  - Add 2-hop neighbor feature aggregation (mean perturbation of neighbors)
  - Add PPI degree/betweenness/closeness as topology features
  - Add module membership (WGCNA) as categorical feature

Reference:
  - GNN for Drug Target Identification (PMID:36168971)
  - EGNF (PMID:41139924)
  - MOGAT (PMID:38474033)
  - Tsvetkov P, et al. Science 2022 (PMID:35298263) - 铜死亡核心基因集
  - Zhou J, et al. Transl Stroke Res 2026 (PMID:41673363) - 单细胞铜死亡CIRI

Input:
  - stage5/ppi_topology.json: PPI network topology
  - stage5/string_ppi.tsv: PPI edge list
  - stage7/gene_shap_importance.csv: ML SHAP importance
  - stage6/gene_perturbation_scores.csv: GRN perturbation scores
  - stage8/core_targets.csv: Data-driven comprehensive scores (labels)
  - cuproptosis_gsva/cuproptosis_gsva_scores.csv: 铜死亡通路评分 (验证用)
  - cuproptosis_gsea/cuproptosis_gsea_summary.csv: GSEA结果 (验证用)

Output:
  - gat_gene_ranking.csv: GAT gene ranking
  - gat_model.pth: Model weights
  - stage9.log: Run log
  - cuproptosis_validation.csv: 铜死亡基因在GAT中的排名验证
"""

import os
import sys
import warnings
import logging
import time
from collections import defaultdict
from datetime import datetime
from typing import Dict, List, Optional, Tuple, Set

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

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import RESULTS_DIR
from scripts.data_manager import StageDataManager

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


# ============================================================
# 铜死亡核心基因集 (Tsvetkov Science 2022)
# ============================================================
CUPROPTOSIS_CORE_GENES = [
    "FDX1", "LIAS", "LIPT1", "DLAT", "PDHA1", "PDHB",
    "MTF1", "GLS", "CDKN2A", "SLC31A1", "ATP7A", "ATP7B",
    "DLD", "DBT", "DLST", "PDHA2", "GCSH"
]

CUPROPTOSIS_EXTENDED_GENES = CUPROPTOSIS_CORE_GENES + [
    "ATOX1", "COX17", "CCS", "COX11", "SCO1", "SCO2",
    "STEAP1", "STEAP2", "STEAP3", "STEAP4",
    "CP", "COMMD1", "MT1A", "MT2A"
]


def load_ppi_network():
    """使用StageDataManager加载PPI网络数据"""
    logger.info("Loading PPI network (Stage5)...")
    
    dm = StageDataManager()
    
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
            import json
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
    """使用StageDataManager加载跨阶段结果"""
    dm = StageDataManager()
    results = {}
    
    # Stage7: SHAP importance
    results['shap'] = dm.load_csv_as_dict(
        'stage7_ml_shap',
        'gene_shap_importance.csv',
        key_col='Gene',
        value_col='SHAP_importance',
        default={}
    )
    logger.info(f"  SHAP: {len(results['shap'])} genes")
    
    # Stage6: GRN perturbation scores
    results['perturbation'] = dm.load_csv_as_dict(
        'stage6_sctenifold_knockout',
        'gene_perturbation_scores.csv',
        key_col='gene',
        value_col='perturbation_score',
        default={}
    )
    logger.info(f"  Perturbation scores: {len(results['perturbation'])} genes")
    
    # WGCNA module membership (optional)
    wgcna_df = dm.load_csv('stage4_seed_wgcna', 'wgcna_modules.csv')
    if wgcna_df is not None and not wgcna_df.empty:
        try:
            if 'Module' in wgcna_df.columns:
                wgcna_col = 'ProbeID' if 'ProbeID' in wgcna_df.columns else 'GeneSymbol'
                results['wgcna_module'] = dict(zip(
                    wgcna_df[wgcna_col].str.upper(), wgcna_df['Module']
                ))
                logger.info(f"  WGCNA modules: {len(results['wgcna_module'])} genes (col: {wgcna_col})")
        except Exception:
            results['wgcna_module'] = {}
    else:
        results['wgcna_module'] = {}
    
    return results


def build_labels(gene_list, stage_results=None, ppi_df=None):
    """使用StageDataManager构建节点标签"""
    logger.info("Building node labels (v9: Stage8 data-driven comprehensive)...")
    
    dm = StageDataManager()
    core_df = dm.load_csv('stage8_final_targets', 'core_targets.csv')
    
    if core_df is None or core_df.empty:
        logger.error(f"Stage8 core_targets.csv not found or empty")
        logger.error("Please run Stage8 first")
        sys.exit(1)
    
    if 'Comprehensive' not in core_df.columns:
        logger.error("core_targets.csv missing 'Comprehensive' column")
        sys.exit(1)
    
    gene_to_score = dict(zip(core_df['Gene'].str.upper(), core_df['Comprehensive']))
    logger.info(f"  Loaded {len(gene_to_score)} gene scores from Stage8")
    
    y = []
    n_missing = 0
    for gene in gene_list:
        gene_upper = gene.upper()
        score = gene_to_score.get(gene_upper, 0.0)
        if score == 0.0 and gene_upper not in gene_to_score:
            n_missing += 1
        y.append(score)
    
    y = np.array(y, dtype=np.float32)
    
    y_min, y_max = y.min(), y.max()
    if y_max > y_min:
        y_normalized = (y - y_min) / (y_max - y_min)
    else:
        y_normalized = np.zeros_like(y)
    
    logger.info(f"  Label range: [{y.min():.2f}, {y.max():.2f}] -> [{y_normalized.min():.3f}, {y_normalized.max():.3f}]")
    logger.info(f"  Label mean: {y_normalized.mean():.3f}, std: {y_normalized.std():.3f}")
    logger.info(f"  Genes in PPI but not in Stage8: {n_missing}/{len(gene_list)}")
    
    hist, _ = np.histogram(y_normalized, bins=10, range=(0, 1))
    logger.info(f"  Label distribution (10 bins): {hist}")
    
    return y_normalized


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
    Build node features (v9: 无偏8D + 铜死亡验证)
    
    Features (v9 data-driven, no priors):
    1. GRN perturbation score (Stage6, independent computation)
    2. SHAP importance (Stage7, independent computation)
    3. PPI degree (network topology)
    4. PPI betweenness centrality
    5. PPI closeness centrality
    6. 1-hop neighbor mean perturbation (message passing)
    7. 2-hop neighbor mean perturbation (message passing)
    8. WGCNA module membership (categorical, encoded as numeric)
    
    Labels: Stage8 v9 纯数据驱动综合得分
    """
    logger.info("Building node features (v9: 8D data-driven)...")
    
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
        
        pert = stage_results['perturbation'].get(gene_upper, 0.0)
        shap = stage_results['shap'].get(gene_upper, 0.0)
        
        degree = topology.get(gene_upper, {}).get('degree', 0)
        betweenness = topology.get(gene_upper, {}).get('betweenness', 0.0)
        closeness = topology.get(gene_upper, {}).get('closeness', 0.0)
        
        mean_pert_1hop = pert_1hop[i]
        mean_pert_2hop = pert_2hop[i]
        
        module_idx = module_to_idx.get(wgcna_module_dict.get(gene_upper, ''), -1)
        wgcna_feature = module_idx if module_idx >= 0 else -1
        
        features.append([
            pert, shap, degree,
            betweenness, closeness, mean_pert_1hop, mean_pert_2hop, wgcna_feature
        ])
    
    X = np.array(features, dtype=np.float32)
    
    # Standardize
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    feature_names = [
        'perturbation', 'shap', 'degree',
        'betweenness', 'closeness', 'neighbor_pert_1hop', 'neighbor_pert_2hop',
        'wgcna_module'
    ]
    
    logger.info(f"  Feature matrix: {X_scaled.shape}")
    logger.info(f"  Features: {feature_names}")
    
    return X_scaled, gene_list, feature_names



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


def validate_cuproptosis_genes(gene_list: List[str], 
                                y_pred: np.ndarray,
                                y_true: np.ndarray,
                                output_dir: str) -> pd.DataFrame:
    """
    验证铜死亡基因在GAT排名中的位置
    
    目的: 检查数据驱动管线是否能自然发现铜死亡相关基因
    """
    logger.info("Validating cuproptosis genes in GAT ranking...")
    
    # 创建排名
    ranking = pd.DataFrame({
        'Gene': gene_list,
        'GAT_score': y_pred,
        'True_label': y_true
    }).sort_values('GAT_score', ascending=False)
    ranking['Rank'] = range(1, len(ranking) + 1)
    
    # 检查铜死亡基因
    cupro_set = set(g.upper() for g in CUPROPTOSIS_CORE_GENES)
    cupro_in_ppi = [g for g in cupro_set if g in [x.upper() for x in gene_list]]
    
    validation_results = []
    
    for gene in cupro_in_ppi:
        gene_row = ranking[ranking['Gene'].str.upper() == gene]
        if not gene_row.empty:
            rank = gene_row['Rank'].iloc[0]
            score = gene_row['GAT_score'].iloc[0]
            percentile = (1 - rank / len(ranking)) * 100
            
            validation_results.append({
                'Gene': gene,
                'GAT_Rank': rank,
                'Total_Genes': len(ranking),
                'Percentile': percentile,
                'GAT_Score': score,
                'True_Label': gene_row['True_label'].iloc[0],
                'In_Top_10pct': rank <= len(ranking) * 0.1,
                'In_Top_5pct': rank <= len(ranking) * 0.05
            })
    
    validation_df = pd.DataFrame(validation_results).sort_values('GAT_Rank')
    
    logger.info(f"  铜死亡基因在PPI中: {len(cupro_in_ppi)}/{len(cupro_set)}")
    logger.info(f"  平均排名: {validation_df['GAT_Rank'].mean():.0f} / {len(ranking)}")
    logger.info(f"  前10%基因数: {validation_df['In_Top_10pct'].sum()}")
    logger.info(f"  前5%基因数: {validation_df['In_Top_5pct'].sum()}")
    
    # 保存验证结果
    validation_file = os.path.join(output_dir, "cuproptosis_validation.csv")
    validation_df.to_csv(validation_file, index=False)
    logger.info(f"  验证结果已保存: {validation_file}")
    
    return validation_df


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
    
    # 铜死亡基因验证
    validate_cuproptosis_genes(gene_list, y_pred, y_true, STAGE_DIR)
    
    return ranking_df


def main():
    logger.info("=" * 60)
    logger.info("Stage 9: PPI Network + GAT Node Regression (v9)")
    logger.info("=" * 60)
    logger.info(f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info("FIX: [v9] 无偏数据驱动管线 + 铜死亡模块整合")
    logger.info("Features: 8D (perturbation, shap, degree, betweenness, closeness,")
    logger.info("         neighbor_pert_1hop, neighbor_pert_2hop, wgcna_module)")
    logger.info("Labels: Stage8 v9 comprehensive (GRN+ML+PPI, no priors)")
    logger.info("Validation: Cuproptosis gene ranking in GAT output")
    
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
    
    # 4. Build labels (v9: Stage8 data-driven comprehensive, no cupro/BCP priors)
    y = build_labels(gene_list, stage_results, ppi_df)
    
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
