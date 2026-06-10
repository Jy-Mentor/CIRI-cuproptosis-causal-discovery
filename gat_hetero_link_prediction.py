# -*- coding: utf-8 -*-
"""
异构图注意力网络 (HeteroGAT) 用于 β-石竹烯-CIRI 桥梁靶点预测
v3: MRHormer TRLA (Sigmoid Attention) + MLGANN Multi-layer Pooling + TaRGET II 多源特征

架构:
  1. SigmoidAttnConv: 替换 GATv2Conv 的 Softmax 为 Sigmoid 注意力
     → 类型特定投影矩阵 + 独立邻居重要性
  2. MLGANN 多层池化: 对 conv1/conv2 输出做可学习注意力加权融合
  3. 多源特征增强: RNA FPKM + ATAC narrowPeak (TaRGET II, 小鼠→人类映射)

解码器: MLP Decoder
负采样: NegativeSampler 预缓存 + PPI 邻居约束
CV 防泄漏: 每折仅用训练集疾病基因构建疾病节点特征

数据文件:
  GAT/drug_targets.txt         → β-石竹烯靶点
  GAT/disease_genes.txt        → CIRI 疾病相关基因
  GAT/ppi_subgraph.csv         → STRING PPI 子图 (471K 边, Tab分隔)
  GAT/subgraph_genes.txt       → 子图基因列表 (15,648)
  GAT/drug_fingerprint.csv     → ECFP4指纹 (1 × 1024)
  GAT/subgraph_embeddings.csv  → 蛋白质嵌入 (15,648 × 1072)
  toxirna_enhanced_features.csv → TaRGET II FPKM+ATAC (N × 32, PCA降维)
"""

import os
import math
import logging
import warnings
import random
import copy
import numpy as np
import gc
import time
import pandas as pd
import torch
import torch.nn.functional as F
from torch.nn import Linear, Parameter, LayerNorm
from torch_geometric.nn import HeteroConv, MessagePassing
from torch_geometric.data import HeteroData
from collections import defaultdict
from sklearn.model_selection import KFold
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")

# ============================================================================
# Logging 配置
# ============================================================================
logger = logging.getLogger("HeteroGAT")
logger.setLevel(logging.INFO)
if not logger.handlers:
    _formatter = logging.Formatter(
        "[%(asctime)s] %(message)s", datefmt="%H:%M:%S")
    # 控制台 handler (stderr)
    _ch = logging.StreamHandler()
    _ch.setLevel(logging.INFO)
    _ch.setFormatter(_formatter)
    _ch.flush = lambda: _ch.stream.flush()  # 强制刷新
    logger.addHandler(_ch)
    # 文件 handler (实时写入, 绕过 PowerShell 缓冲)
    _log_dir = os.path.dirname(os.path.abspath(__file__))
    _fh = logging.FileHandler(
        os.path.join(_log_dir, "training_progress.log"), mode="a", encoding="utf-8")
    _fh.setLevel(logging.INFO)
    _fh.setFormatter(_formatter)
    logger.addHandler(_fh)
    # 确保每次日志输出都刷新
    logger.handlers[0].flush()


def drop_edge(edge_index, p=0.1):
    """训练时随机丢弃 PPI 边 (DropEdge)，缓解过平滑和过拟合"""
    if p <= 0:
        return edge_index
    mask = torch.rand(edge_index.size(1), device=edge_index.device) > p
    return edge_index[:, mask]


# ============================================================================
# 0. 配置与随机种子
# ============================================================================
GAT_DATA_DIR = r"C:\Users\Jy-Mentor-7\Desktop\GAT"
GAT_EXT_DIR = r"D:\反向网络药理学\GAT拓展维度"
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT = os.path.join(SCRIPT_DIR, "top20_bridge_genes.csv")

HIDDEN_DIM = 128
OUT_DIM = 64
GAT_HEADS = 4
DROPOUT = 0.5
LR = 1e-3
WEIGHT_DECAY = 5e-4
EPOCHS = 200
PATIENCE = 30
N_FOLDS = 5
N_LAYERS = 2               # GNN层数: 2层最佳平衡点（>2导致基因过平滑, Li et al. 2024）
VAL_RATIO = 0.15
TD_LOSS_WEIGHT = 2.0
PPI_SCORE_THRESHOLD = 700
RANDOM_SEED = 42
ENSEMBLE_RUNS = 5
SWA_START = 120
DROPEDGE_P = 0.45          # 训练时随机丢弃45% PPI边 (Rong et al. 2020 ICLR: DropEdge=0.3~0.5最优)
                            # 原0.15过小: avg_degree ~60, 15%裁剪后仍有~51邻居, 过平滑不减
                            # 45%裁剪后 ~33邻居, 有效缓解过平滑
DROPEDGE_P_MIN = 0.05      # 余弦退火DropEdge最低值 (前期多用边学拓扑)
RANK_WEIGHT = 0.5
RANK_MARGIN = 0.3
FOCAL_GAMMA = 1.5
FOCAL_ALPHA = 0.25
LOGIT_CLAMP = 10.0
SIGMOID_TEMP = 5.0
NEG_DEGREE_POWER = 0.5     # 节点度数感知负采样: 平滑系数 (Cappelletti 2024)
EARLY_STOP_HITS_K = 100    # 早停指标包含Hits@K, 避免AUROC过拟合信号

# ---- GPU 资源限制（缓解系统卡顿） ----
# 注释掉显存上限限制: 50xx 系驱动已知 set_per_process_memory_fraction 可能导致 CUDA 内核卡死,
# 且实际显存 8151 MiB 已由 PyTorch AMP 自动管理, 无需手动限制.
# GPU_MEM_FRACTION = 0.75     # 取消, 让 CUDA driver 按需分配
CUDA_BENCHMARK = False      # 禁用 CuDNN autotune（降低 GPU 利用率尖峰）
PIN_MEMORY = False          # 禁用 pin_memory（减少内存占用）

# ---- 多源异构图增强 ----
USE_TF_EDGES = False
USE_METHYLATION_FEATURES = False
USE_PATHWAY_EDGES = False
USE_CTD_EDGES = True           # Step 3: + CTD gene-level features (16维)
# Step 1 结果 (基线+CTD直接边合并): DT AUC=0.7699, DT AP=0.7807, TD AUC=0.6671, TD AP=0.6496, P@20=0.77
CTD_MERGE_TARGETS = False       # Step 2: 不合并靶点，保持chemical节点独立，通过化学物节点增强药物表示
CTD_GENE_FEATURES = True        # Step 3: 启用基因级CTD聚合特征 (16维), 拼接到基因初始特征
CTD_TOP_CHEMICALS = 10          # Step 2: 减少化学物数量避免信息泛滥 (25→10)
CTD_FOCUSED_EDGES = True        # Step 2: 仅保留drug_targets/disease_genes范围内的chem→gene边
PATHWAY_PCA_DIM = 256
USE_COEXP_EDGES = False

# ---- 消融实验模式 (已完成, 关闭) ----
ABLATION_MODE = False

ABLATION_CONFIGS = [
    # ========== 递增实验: 从基线逐项添加 (Addition) ==========
    {"name": "Baseline (PPI+DT+GD only)",      "tf": False, "meth": False, "pw": False, "ctd": False, "coexp": False},
    {"name": "+ CTD (direct edges)",            "tf": False, "meth": False, "pw": False, "ctd": True,  "coexp": False},
    {"name": "+ TF (TRRUST)",                   "tf": True,  "meth": False, "pw": False, "ctd": False, "coexp": False},
    {"name": "+ Co-expression (ARCHS4)",        "tf": False, "meth": False, "pw": False, "ctd": False, "coexp": True},
    {"name": "+ Pathway (KEGG/Reactome)",       "tf": False, "meth": False, "pw": True,  "ctd": False, "coexp": False},
    {"name": "+ Methylation (EWAS Atlas)",      "tf": False, "meth": True,  "pw": False, "ctd": False, "coexp": False},
    {"name": "Full Model (all components)",     "tf": True,  "meth": True,  "pw": True,  "ctd": True,  "coexp": True},

    # ========== 消融实验: 从全模型中逐一移除 (Ablation) ==========
    {"name": "Full - CTD",                      "tf": True,  "meth": True,  "pw": True,  "ctd": False, "coexp": True},
    {"name": "Full - TF",                       "tf": False, "meth": True,  "pw": True,  "ctd": True,  "coexp": True},
    {"name": "Full - Co-expression",            "tf": True,  "meth": True,  "pw": True,  "ctd": True,  "coexp": False},
    {"name": "Full - Pathway",                  "tf": True,  "meth": True,  "pw": False, "ctd": True,  "coexp": True},
    {"name": "Full - Methylation",              "tf": True,  "meth": False, "pw": True,  "ctd": True,  "coexp": True},
]

# ---- 解析度冗余缓存（INT8 令牌化缓存） ----

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
USE_AMP = False  # AMP 在 float16 下产生 NaN (SigmoidAttnConv 消息传递数值溢出), 使用 float32 训练
USE_COMPILE = False  # Windows 无 Triton，AMP 已提供主要加速

# ---- GPU 资源限制（缓解系统卡顿） ----
if device.type == "cuda":
    # 不再设置显存上限: 旧版 PyTorch 的 set_per_process_memory_fraction 
    # 在 50xx 系列驱动 + CUDA 12.8 环境中可能导致 CUDA kernel 挂起.
    # PyTorch AMP + cuDNN 自动调优已足够管理 ~8GB 显存.
    torch.backends.cudnn.benchmark = CUDA_BENCHMARK
    logger.info(f"[GPU] Memory limit: none (PyTorch AMP manages ~8GB), cudnn.benchmark={CUDA_BENCHMARK}")


def set_seed(seed):
    """全局确定性种子: Python + NumPy + PyTorch + CUDA (禁用TF32)"""
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    # TF32 在 Ampere 架构 (RTX 50xx) 上可安全提速 ~1.5x，对精度影响极小
    # torch.backends.cuda.matmul.allow_tf32 = False
    # torch.backends.cudnn.allow_tf32 = False


set_seed(RANDOM_SEED)
logger.info(f"[Config] Device: {device}, Hidden: {HIDDEN_DIM}, Out: {OUT_DIM}, "
      f"Heads: {GAT_HEADS}, Layers: {N_LAYERS}, Folds: {N_FOLDS}, PPI threshold: {PPI_SCORE_THRESHOLD}, "
      f"Ensemble: {ENSEMBLE_RUNS} runs")


# ============================================================================
# 1. 数据加载（针对实际文件格式）
# ============================================================================

def load_drug_targets(path):
    """
    读取 drug_targets.txt — 无 header 的基因列表
    实际格式: 第一行是 LYN（被 pandas 误读为列名），后续为基因名
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"drug_targets.txt not found: {path}")
    for enc in ("utf-8", "gbk", "latin1"):
        try:
            with open(path, "r", encoding=enc) as fh:
                lines = [ln.strip() for ln in fh if ln.strip()]
            break
        except UnicodeDecodeError:
            continue
    genes = [g.upper() for g in lines]
    unique_genes = list(dict.fromkeys(genes))
    logger.info(f"[Load] drug_targets: {len(unique_genes)} unique genes from {path}")
    return unique_genes


def load_disease_genes(path):
    """
    读取疾病基因文件 — 支持 CSV (gene_symbol,source) 和纯文本两种格式
    仅保留全大写字母/数字组成的 human-like 基因名
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"disease genes file not found: {path}")
    for enc in ("utf-8", "gbk", "latin1"):
        try:
            with open(path, "r", encoding=enc) as fh:
                lines = [ln.strip() for ln in fh if ln.strip()]
            break
        except UnicodeDecodeError:
            continue
    # 判断是否为 CSV 格式（包含逗号）
    if "," in lines[0] or "\t" in lines[0]:
        genes = []
        for ln in lines:
            parts = ln.replace("\t", ",").split(",")
            gene = parts[0].strip().upper()
            if gene and gene not in ("GENE_SYMBOL", "GENE"):
                genes.append(gene)
    else:
        genes = [g.upper() for g in lines if g.upper() not in ("GENE_SYMBOL", "GENE")]
    unique_genes = list(dict.fromkeys(genes))
    logger.info(f"[Load] disease_genes raw: {len(unique_genes)} (may include non-human)")
    return unique_genes


def load_ppi(path, score_thresh=PPI_SCORE_THRESHOLD):
    """
    读取 ppi.csv — Tab 分隔 (gene_a\\tgene_b\\tcombined_score)
    实际格式: 三列 Tab 分隔, header 为 gene_a / gene_b / combined_score
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"ppi.csv not found: {path}")
    for enc in ("utf-8", "gbk", "latin1"):
        try:
            df = pd.read_csv(path, sep="\t", encoding=enc)
            break
        except (UnicodeDecodeError, pd.errors.ParserError):
            continue
    if df.shape[1] == 1:
        df = pd.read_csv(path, sep="\t", encoding="utf-8")
    col_a, col_b = df.columns[0], df.columns[1]
    score_col = df.columns[2] if len(df.columns) >= 3 else None
    df = df.dropna(subset=[col_a, col_b])
    df[col_a] = df[col_a].astype(str).str.strip().str.upper()
    df[col_b] = df[col_b].astype(str).str.strip().str.upper()
    if score_col is not None:
        df[score_col] = pd.to_numeric(df[score_col], errors="coerce")
        df = df[df[score_col] >= score_thresh]
    edges = list(zip(df[col_a], df[col_b]))
    logger.info(f"[Load] PPI: {len(edges)} edges (score >= {score_thresh}) from {path}")
    return edges


def load_gene_features(path):
    """
    读取 subgraph_embeddings.csv — 蛋白质语言模型嵌入 (15,467 × 1025)
    第一列 gene_symbol, 后续列 feat_0001 ~ feat_1024
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"subgraph_embeddings.csv not found: {path}")
    for enc in ("utf-8", "gbk", "latin1"):
        try:
            df = pd.read_csv(path, encoding=enc, index_col=0)
            break
        except (UnicodeDecodeError, pd.errors.ParserError):
            continue
    df.index = df.index.astype(str).str.strip().str.upper()
    features = df.apply(pd.to_numeric, errors="coerce").fillna(0.0).values.astype(np.float32)
    gene_list = df.index.tolist()
    logger.info(f"[Load] gene features: {features.shape} (dim={features.shape[1]}) from {path}")
    return features, gene_list


def load_drug_fingerprint(path):
    """
    读取 drug_fingerprint.csv — ECFP4 指纹 (1 × 1024)
    列名: fp_0 ~ fp_1023
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"drug_fingerprint.csv not found: {path}")
    for enc in ("utf-8", "gbk", "latin1"):
        try:
            df = pd.read_csv(path, encoding=enc)
            break
        except (UnicodeDecodeError, pd.errors.ParserError):
            continue
    fp = df.apply(pd.to_numeric, errors="coerce").fillna(0).values.astype(np.float32)
    if fp.ndim == 1:
        fp = fp.reshape(1, -1)
    logger.info(f"[Load] drug fingerprint: {fp.shape} from {path}")
    return fp


def load_all_genes(path):
    """
    读取 all_genes.txt — 每行一个 gene_symbol（有 header 'gene'）
    """
    if not os.path.exists(path):
        logger.info(f"[Warn] all_genes.txt not found: {path}")
        return None
    for enc in ("utf-8", "gbk", "latin1"):
        try:
            with open(path, "r", encoding=enc) as fh:
                lines = [ln.strip() for ln in fh if ln.strip()]
            if lines[0].upper().replace('"', '') in ("GENE", "GENE_SYMBOL", "GENES"):
                lines = lines[1:]
            break
        except UnicodeDecodeError:
            continue
    genes = list(dict.fromkeys([g.upper() for g in lines]))
    logger.info(f"[Load] all_genes: {len(genes)} genes from {path}")
    return genes


def load_toxirna_features(path):
    """
    读取 toxirna_enhanced_features.csv — TaRGET II FPKM+ATAC PCA 特征
    第一列 gene_symbol (index), 后续列 toxirna_000 ~ toxirna_031
    返回: DataFrame (gene_symbol → features)
    """
    if not os.path.exists(path):
        logger.info(f"[Warn] toxirna features not found: {path}")
        return None
    for enc in ("utf-8", "gbk", "latin1"):
        try:
            df = pd.read_csv(path, encoding=enc, index_col=0)
            break
        except (UnicodeDecodeError, pd.errors.ParserError):
            continue
    df.index = df.index.astype(str).str.strip().str.upper()
    logger.info(f"[Load] toxirna features: {df.shape} from {path}")
    return df  # Return DataFrame with gene_symbol index for easy lookup


def load_tf_edges(path):
    """加载 TF-靶基因调控边 (TRRUST) — 无 header: TF\\tGene\\tRegulationType"""
    if not os.path.exists(path):
        logger.info(f"[Warn] TF edges not found: {path}")
        return [], {}
    tf_edges = []
    tf_nodes = {}
    with open(path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                tf_name = parts[0].upper()
                gene_name = parts[1].upper()
                tf_edges.append((tf_name, gene_name))
                if tf_name not in tf_nodes:
                    tf_nodes[tf_name] = {'degree': 0, 'activation': 0, 'repression': 0, 'unknown': 0}
                tf_nodes[tf_name]['degree'] += 1
                if len(parts) >= 3:
                    reg = parts[2].lower()
                    if 'activation' in reg:
                        tf_nodes[tf_name]['activation'] += 1
                    elif 'repression' in reg:
                        tf_nodes[tf_name]['repression'] += 1
                    else:
                        tf_nodes[tf_name]['unknown'] += 1
    logger.info(f"[Load] TF edges: {len(tf_edges)} edges, {len(tf_nodes)} TFs")
    return tf_edges, tf_nodes


def load_methylation_features(path, gene_to_idx):
    """聚合甲基化特征到基因级别 (来自 EWAS Atlas)"""
    if not os.path.exists(path):
        logger.info(f"[Warn] Methylation edges not found: {path}")
        return None
    gene_data = defaultdict(list)
    with open(path) as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                gene = parts[0].upper()
                try:
                    beta = float(parts[2])
                except ValueError:
                    continue
                status = parts[3]
                gene_data[gene].append((beta, status))

    n_genes = len(gene_to_idx)
    # 20 维甲基化特征 (10→20):
    #   mean_beta, std_beta, min_beta, max_beta, n_cpgs_log,
    #   hypo_ratio, hyper_ratio, normal_ratio, median_beta, beta_range,
    #   skewness, kurtosis, q25, q75, iqr,
    #   extreme_hypo_ratio (beta<0.1), extreme_hyper_ratio (beta>0.9),
    #   mean_deviation, coef_variation, entropy_approx
    METH_DIM = 20
    meth_feat = np.zeros((n_genes, METH_DIM), dtype=np.float32)
    matched = 0
    for gene, idx in gene_to_idx.items():
        if gene in gene_data:
            betas = [d[0] for d in gene_data[gene]]
            statuses = [d[1] for d in gene_data[gene]]
            betas_arr = np.array(betas)
            n = len(betas)
            hypo = sum(1 for s in statuses if 'hypo' in s.lower())
            hyper = sum(1 for s in statuses if 'hyper' in s.lower())
            normal = n - hypo - hyper
            extreme_hypo = int(np.sum(betas_arr < 0.1))
            extreme_hyper = int(np.sum(betas_arr > 0.9))

            # 偏度 & 峰度 (使用无偏估计)
            mean_b = np.mean(betas_arr)
            std_b = np.std(betas_arr)
            if std_b > 1e-8:
                skew = np.mean(((betas_arr - mean_b) / std_b) ** 3)
                kurt = np.mean(((betas_arr - mean_b) / std_b) ** 4) - 3  # 超值峰度
            else:
                skew, kurt = 0.0, 0.0

            # 分位数
            q25 = np.percentile(betas_arr, 25)
            q75 = np.percentile(betas_arr, 75)

            # 熵近似 (beta 值离散化)
            hist, _ = np.histogram(betas_arr, bins=10, range=(0, 1))
            hist = hist / max(n, 1)
            entropy = -np.sum(hist[hist > 0] * np.log(hist[hist > 0] + 1e-10))

            meth_feat[idx] = [
                mean_b, std_b, np.min(betas_arr), np.max(betas_arr),
                np.log1p(n),
                hypo / max(n, 1), hyper / max(n, 1), normal / max(n, 1),
                np.median(betas_arr), np.max(betas_arr) - np.min(betas_arr),
                skew, kurt, q25, q75, q75 - q25,
                extreme_hypo / max(n, 1), extreme_hyper / max(n, 1),
                np.mean(np.abs(betas_arr - mean_b)),  # mean absolute deviation
                std_b / max(mean_b, 1e-8),  # coefficient of variation
                entropy,
            ]
            matched += 1
    logger.info(f"[Load] Methylation features: {matched}/{n_genes} genes, {METH_DIM} dims")
    return meth_feat


def load_pathway_edges(path):
    """加载基因-通路边 (KEGG/Reactome) — 无 header: Gene\\tPathway"""
    if not os.path.exists(path):
        logger.info(f"[Warn] Pathway edges not found: {path}")
        return [], {}
    pw_edges = []
    pw_nodes = {}
    with open(path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                gene = parts[0].upper()
                pathway = parts[1]
                pw_edges.append((gene, pathway))
                pw_nodes[pathway] = pw_nodes.get(pathway, 0) + 1
    logger.info(f"[Load] Pathway edges: {len(pw_edges)} edges, {len(pw_nodes)} pathways")
    return pw_edges, pw_nodes


def load_pathway_features(path):
    """加载通路 PCA 特征 + 名称列表 (若存在)"""
    if not os.path.exists(path):
        logger.info(f"[Warn] Pathway features not found: {path}")
        return None, None
    pwy_feat = np.load(path)
    logger.info(f"[Load] Pathway features: {pwy_feat.shape}")

    # 尝试从同伴文件加载通路名称列表 (pathway_pca_names.txt 或 .npy)
    name_path = path.replace(".npy", "_names.txt")
    pw_names = None
    if os.path.exists(name_path):
        with open(name_path, "r") as f:
            pw_names = [line.strip() for line in f if line.strip()]
        logger.info(f"[Load] Pathway names: {len(pw_names)} from {name_path}")
    else:
        # 尝试 .npy 格式  (list/array of strings)
        name_npy = path.replace(".npy", "_names.npy")
        if os.path.exists(name_npy):
            pw_names = list(np.load(name_npy, allow_pickle=True))
            logger.info(f"[Load] Pathway names: {len(pw_names)} from {name_npy}")

    if pw_names is not None and len(pw_names) != pwy_feat.shape[0]:
        logger.warning(f"[Warn] Pathway names count ({len(pw_names)}) != features rows ({pwy_feat.shape[0]}), discarding names")
        pw_names = None

    return pwy_feat, pw_names  # (n_pathways, pca_dim), pathway_names or None


# CTD 交互类型分组 (20 维 chemical 特征，更细粒度)
CTD_ACTION_GROUPS = [
    ("inc_expr",     ["increases^expression"]),
    ("dec_expr",     ["decreases^expression"]),
    ("aff_expr",     ["affects^expression"]),
    ("inc_react",    ["increases^reaction"]),
    ("dec_react",    ["decreases^reaction"]),
    ("aff_react",    ["affects^reaction"]),
    ("aff_bind",     ["affects^binding"]),
    ("inc_act",      ["increases^activity"]),
    ("dec_act",      ["decreases^activity"]),
    ("aff_loc",      ["affects^localization", "affects^transport", "affects^uptake"]),
    ("methyl",       ["increases^methylation", "decreases^methylation", "affects^methylation"]),
    ("phospho",      ["increases^phosphorylation", "decreases^phosphorylation", "affects^phosphorylation"]),
    ("abundance",    ["increases^abundance", "decreases^abundance"]),
    ("deg",          ["increases^degradation", "decreases^degradation", "affects^degradation"]),
    ("stab",         ["increases^stability", "decreases^stability"]),
    ("sec",          ["increases^secretion", "decreases^secretion"]),
    ("cotreatment",  ["affects^cotreatment"]),
    ("metabolic",    ["increases^metabolic processing", "decreases^metabolic processing",
                      "affects^metabolic processing", "increases^chemical synthesis",
                      "decreases^chemical synthesis"]),
    ("hydroxyl",     ["increases^hydroxylation", "decreases^hydroxylation", "affects^hydroxylation"]),
    ("other",        []),  # 兜底
]
CTD_CHEM_FEAT_DIM = len(CTD_ACTION_GROUPS) + 4  # 20 groups + n_genes + n_pmids + n_pmids_log + n_actions_total = 24

# CTD 直接互作过滤: 排除 expression/abundance/response_to_substance 等间接调控
# expression 占 CTD 行为的 55.1%, 是最大的间接调控噪声源
CTD_DIRECT_ONLY = True
CTD_INDIRECT_ACTION_TYPES = {"expression", "abundance", "response to substance"}
CTD_GENE_FEAT_DIM = 16  # 基因级 CTD 聚合特征维度 (10→16)


def load_ctd_edges(path, gene_set, top_k=50):
    """从 CTD 加载化学-基因互作边，提取丰富的化学物特征 + 基因级 CTD 特征

    返回:
      ctd_edges: [(chemical, gene), ...]
      ctd_chem_nodes: {chemical: {'n_genes': N, 'chem_id': ID, 'features': np.array(24,)}}
      ctd_gene_features: np.array (n_genes_in_set, 16) — 基因级 CTD 聚合特征
    """
    if not os.path.exists(path):
        logger.info(f"[Warn] CTD file not found: {path}")
        return [], {}, None

    # ---- 第一遍扫描: 收集化学物-基因边 + 交互类型 + PubMed ----
    chem_genes = defaultdict(set)          # chemical → set of genes
    chem_names = {}                         # chemical → ChemicalID
    chem_actions = defaultdict(lambda: defaultdict(int))  # chemical → {action_type: count}
    chem_pmids = defaultdict(set)           # chemical → set of PubMed IDs
    gene_chem_actions = defaultdict(lambda: defaultdict(int))  # gene → {action_type: count}
    gene_chem_count = defaultdict(int)      # gene → n_chemicals
    gene_total_actions = defaultdict(int)   # gene → total action count
    pair_has_direct = defaultdict(bool)     # (chem, gene) → has at least one direct action

    action_to_group = {}
    for gi, (grp_name, actions) in enumerate(CTD_ACTION_GROUPS):
        for a in actions:
            action_to_group[a] = gi

    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 10:
                continue
            gene = parts[3].upper()  # GeneSymbol
            if gene not in gene_set:
                continue
            chem = parts[0]
            chem_names[chem] = parts[1]  # ChemicalID

            # ---- 直接互作过滤: 检查该边是否有至少一个直接行为 ----
            interaction_actions = parts[9]  # e.g. "affects^reaction|increases^expression"
            pair_is_direct = False
            for act in interaction_actions.split('|'):
                act = act.strip()
                if '^' in act:
                    action_type = act.split('^')[-1].strip()
                    if action_type not in CTD_INDIRECT_ACTION_TYPES:
                        pair_is_direct = True
                        break  # 至少有一个直接行为即可
            if CTD_DIRECT_ONLY and not pair_is_direct:
                continue  # 纯间接边, 跳过（如纯 expression/abundance 变化）
            pair_has_direct[(chem, gene)] = pair_is_direct

            chem_genes[chem].add(gene)

            # 交互类型统计 (仅统计直接行为, 或全量)
            for act in interaction_actions.split('|'):
                act = act.strip()
                if CTD_DIRECT_ONLY and '^' in act:
                    action_type = act.split('^')[-1].strip()
                    if action_type in CTD_INDIRECT_ACTION_TYPES:
                        continue  # 跳过间接行为，不纳入化学物特征统计
                grp_idx = action_to_group.get(act, len(CTD_ACTION_GROUPS) - 1)
                chem_actions[chem][grp_idx] += 1
                gene_chem_actions[gene][grp_idx] += 1
                gene_total_actions[gene] += 1

            # PubMed IDs
            if len(parts) >= 11:
                for pmid in parts[10].split('|'):
                    pmid = pmid.strip()
                    if pmid:
                        chem_pmids[chem].add(pmid)

            gene_chem_count[gene] += 1

    # ---- 取 top_k 化学物 (按关联基因数排序) ----
    top_chems = sorted(chem_genes.items(), key=lambda x: len(x[1]), reverse=True)[:top_k]

    n_groups = len(CTD_ACTION_GROUPS)
    ctd_edges = []
    ctd_chem_nodes = {}
    for chem, genes in top_chems:
        feat = np.zeros(CTD_CHEM_FEAT_DIM, dtype=np.float32)
        # 交互类型组计数 (dim 0 ~ n_groups-1)
        total_actions = 0
        for gi in range(n_groups):
            count = float(chem_actions[chem].get(gi, 0))
            feat[gi] = count
            total_actions += count
        feat[n_groups] = float(len(genes))              # n_genes
        feat[n_groups + 1] = float(len(chem_pmids[chem]))  # n_pmids
        feat[n_groups + 2] = np.log1p(len(chem_pmids[chem]))  # n_pmids_log
        feat[n_groups + 3] = float(total_actions)         # n_actions_total

        ctd_chem_nodes[chem] = {
            'n_genes': len(genes),
            'chem_id': chem_names.get(chem, ''),
            'features': feat,
        }
        for gene in genes:
            ctd_edges.append((chem, gene))

    logger.info(f"[Load] CTD edges: {len(ctd_edges)} edges, "
                f"{len(ctd_chem_nodes)} chemicals (top {top_k}), "
                f"chem feat dim={CTD_CHEM_FEAT_DIM}"
                + (f", DIRECT_ONLY (excluded: expression/abundance/response)" if CTD_DIRECT_ONLY else ""))
    if CTD_DIRECT_ONLY:
        total_pairs = len(pair_has_direct)
        direct_pairs = sum(1 for v in pair_has_direct.values() if v)
        logger.info(f"[Load] CTD direct-only filter: {direct_pairs}/{total_pairs} pairs have direct actions"
                    f" ({100*direct_pairs/max(1,total_pairs):.1f}%)")

    # ---- 构建基因级 CTD 聚合特征 (16 dims) ----
    gene_list = sorted(gene_set)
    ctd_gene_feat = np.zeros((len(gene_list), CTD_GENE_FEAT_DIM), dtype=np.float32)
    ctd_gene_matched = 0
    for i, gene in enumerate(gene_list):
        if gene in gene_chem_count:
            n_chem = gene_chem_count[gene]
            n_actions = gene_total_actions[gene]
            ctd_gene_feat[i, 0] = float(n_chem)               # 总化学物数
            ctd_gene_feat[i, 1] = np.log1p(n_chem)            # log 变换
            ctd_gene_feat[i, 2] = float(n_actions)            # 总交互事件数
            ctd_gene_feat[i, 3] = np.log1p(n_actions)         # log 变换
            ctd_gene_feat[i, 4] = n_actions / max(n_chem, 1)   # 平均每个化学物的交互数
            # 前 11 个交互类型组归一化（前 11 组 + 位置5作为多样性指标）
            for gi in range(min(n_groups, 11)):
                ctd_gene_feat[i, 5 + gi] = float(gene_chem_actions[gene].get(gi, 0)) / max(n_actions, 1)
            ctd_gene_matched += 1

    logger.info(f"[Load] CTD gene-level features: {ctd_gene_matched}/{len(gene_list)} genes, "
                f"{CTD_GENE_FEAT_DIM} dims")
    return ctd_edges, ctd_chem_nodes, ctd_gene_feat


def load_coexp_edges(path):
    """加载 ARCHS4 共表达边 (|correlation| > 0.7)

    返回:
      coexp_edges: [(gene1, gene2, corr), ...]  — 带权重的边
    """
    if not os.path.exists(path):
        logger.info(f"[Warn] Co-expression edges not found: {path}")
        return []
    edges = []
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                try:
                    corr = float(parts[2])
                except ValueError:
                    continue
                edges.append((parts[0].upper(), parts[1].upper(), corr))
    logger.info(f"[Load] Co-expression edges: {len(edges)} (ARCHS4, |corr|>0.7)")
    return edges


# ============================================================================
# 2. 异构图构建
# ============================================================================

def build_hetero_data(
    drug_target_genes,
    disease_genes_raw,
    ppi_edges,
    gene_features_arr,
    gene_feature_names,
    drug_fingerprint_arr,
    all_genes=None,
    toxirna_df=None,
    tf_edges=None, tf_nodes=None,
    meth_path=None, use_methylation=True,
    pw_edges=None, pw_nodes=None,
    pw_features=None, pw_feature_names=None,
    ctd_edges=None, ctd_chem_nodes=None,
    ctd_gene_feat=None,
    coexp_edges=None, use_coexp=True,
):
    """
    构建 PyG HeteroData 异构图

    节点类型: drug (1个), gene (N个), disease (1个)
              + tf (转录因子), pathway (通路), chemical (化学物) [可选]
    边类型:
      - ('drug', 'targets', 'gene')
      - ('gene', 'associated_with', 'disease')
      - ('gene', 'interacts', 'gene')   [PPI]
      - ('gene', 'coexpressed_with', 'gene')  [ARCHS4 共表达, 可选]
      - ('tf', 'regulates', 'gene') / ('gene', 'regulated_by', 'tf')  [可选]
      - ('gene', 'in_pathway', 'pathway') / ('pathway', 'has_gene', 'gene')  [可选]
      - ('chemical', 'interacts_with', 'gene') / ('gene', 'targeted_by_chemical', 'chemical')  [可选]

    toxirna_df: TaRGET II FPKM+ATAC 特征 DataFrame (gene_symbol → features)
                若提供，将拼接到基因特征上
    meth_path: 甲基化数据文件路径 — 在 gene_to_idx 确定后内部加载
    pw_features: 通路 PCA 特征 (n_pathways, PCA_DIM) — 作为通路节点特征,
                 同时均值聚合到基因特征
    tf_edges/tf_nodes: TF-靶基因调控边
    pw_edges/pw_nodes: 基因-通路边
    ctd_edges/ctd_chem_nodes: CTD 化学-基因互作边
    ctd_gene_feat: 基因级 CTD 聚合特征 (n_genes_in_set, CTD_GENE_FEAT_DIM)
    coexp_edges: ARCHS4 共表达边 [(gene1, gene2, corr), ...]
    """
    # ---- 确定基因全集 ----
    gene_feat_set = set(gene_feature_names)
    if all_genes is not None:
        all_gene_candidates = set(all_genes)
    else:
        all_gene_candidates = (
            gene_feat_set
            | set(drug_target_genes)
            | set(disease_genes_raw)
        )
    for a, b in ppi_edges:
        all_gene_candidates.add(a)
        all_gene_candidates.add(b)

    # 仅保留在 subgraph_embeddings 中有特征的基因（有特征才能参与GAT学习）
    all_genes_list = sorted(all_gene_candidates & gene_feat_set)
    if len(all_genes_list) == 0:
        raise ValueError("No genes with features in the graph!")

    gene_to_idx = {g: i for i, g in enumerate(all_genes_list)}
    n_genes = len(all_genes_list)
    feat_dim = gene_features_arr.shape[1]
    logger.info(f"[Build] Graph genes (with features): {n_genes}")

    # ---- 基因特征矩阵 ----
    gene_feat_dict = dict(zip(gene_feature_names, gene_features_arr))
    gene_feat = np.zeros((n_genes, feat_dim), dtype=np.float32)
    matched = 0
    for i, g in enumerate(all_genes_list):
        if g in gene_feat_dict:
            gene_feat[i] = gene_feat_dict[g]
            matched += 1
    logger.info(f"[Build] Gene features assigned: {matched}/{n_genes}")

    # ---- 拼接 TaRGET II 多源特征 ----
    if toxirna_df is not None and len(toxirna_df) > 0:
        toxirna_feat_dict = {g.upper(): toxirna_df.loc[g.upper()].values
                            for g in toxirna_df.index if g.upper() in gene_to_idx}
        toxirna_dim = toxirna_df.shape[1]
        toxirna_extra = np.zeros((n_genes, toxirna_dim), dtype=np.float32)
        toxirna_matched = 0
        for i, g in enumerate(all_genes_list):
            if g in toxirna_feat_dict:
                toxirna_extra[i] = toxirna_feat_dict[g]
                toxirna_matched += 1

        # StandardScaler: 统一 toxirna 特征尺度，避免与 SapBERT 嵌入（~1072维）
        # 直接拼接时量纲差异导致训练不稳定
        scaler = StandardScaler()
        toxirna_extra = scaler.fit_transform(toxirna_extra).astype(np.float32)

        gene_feat = np.concatenate([gene_feat, toxirna_extra], axis=1)
        logger.info(f"[Build] TaRGET II features appended (StandardScaler): "
              f"{toxirna_matched}/{n_genes} genes, "
              f"+{toxirna_dim} dims → total {gene_feat.shape[1]} dims")
    else:
        logger.info(f"[Build] No TaRGET II features (file not found, will skip)")

    # ---- 拼接甲基化特征 (EWAS Atlas, 10 dims) ----
    if meth_path is not None and use_methylation:
        meth_features = load_methylation_features(meth_path, gene_to_idx)
        if meth_features is not None and meth_features.shape[0] == n_genes:
            gene_feat = np.concatenate([gene_feat, meth_features], axis=1)
            logger.info(f"[Build] Methylation features appended: +{meth_features.shape[1]} dims → total {gene_feat.shape[1]} dims")
        else:
            logger.info(f"[Build] No methylation features (will skip)")
    else:
        logger.info(f"[Build] No methylation features (will skip)")

    # ---- 拼接通路 PCA 特征 (均值聚合到基因) ----
    if pw_features is not None and pw_edges and len(pw_edges) > 0:
        # 构建基因→通路映射, 取通路特征的均值
        pw_name_to_idx = {}
        for pw_name, _ in pw_nodes.items():
            if pw_name not in pw_name_to_idx:
                pw_name_to_idx[pw_name] = len(pw_name_to_idx)

        # ---- 名称显式对齐: pw_features 行顺序通过 pw_feature_names 确定 ----
        pw_feat_aligned = np.zeros((len(pw_name_to_idx), pw_features.shape[1]), dtype=np.float32)
        pw_matched = 0
        if pw_feature_names is not None and len(pw_feature_names) == pw_features.shape[0]:
            # 可靠方案: 通过名称列表逐行映射
            feat_name_to_row = {name: i for i, name in enumerate(pw_feature_names)}
            for pw_name, pi in pw_name_to_idx.items():
                row = feat_name_to_row.get(pw_name)
                if row is not None:
                    pw_feat_aligned[pi] = pw_features[row]
                    pw_matched += 1
        else:
            # 回退: 按 pw_nodes 迭代顺序与 pw_features 行顺序对齐 (有风险, 发出警告)
            if pw_feature_names is None:
                logger.warning("[Build] Pathway feature names not provided; "
                               "assuming pw_features row order matches pw_nodes iteration order.")
            for pw_name, pi in pw_name_to_idx.items():
                if pi < pw_features.shape[0]:
                    pw_feat_aligned[pi] = pw_features[pi]
                    pw_matched += 1
        logger.info(f"[Build] Pathway features aligned: {pw_matched}/{len(pw_name_to_idx)}")

        gene_to_pw = defaultdict(list)
        for gene, pw_name in pw_edges:
            if gene in gene_to_idx and pw_name in pw_name_to_idx:
                gene_to_pw[gene].append(pw_name_to_idx[pw_name])

        pw_agg = np.zeros((n_genes, pw_features.shape[1]), dtype=np.float32)
        pw_agg_matched = 0
        for i, g in enumerate(all_genes_list):
            if g in gene_to_pw:
                pw_indices = gene_to_pw[g]
                pw_agg[i] = pw_feat_aligned[pw_indices].mean(axis=0)
                pw_agg_matched += 1

        gene_feat = np.concatenate([gene_feat, pw_agg], axis=1)
        logger.info(f"[Build] Pathway PCA features appended: {pw_agg_matched}/{n_genes} genes, "
              f"+{pw_features.shape[1]} dims → total {gene_feat.shape[1]} dims")
    else:
        logger.info(f"[Build] No pathway features (will skip)")
        pw_feat_aligned = None  # 后续通路节点特征使用

    # ---- 拼接基因级 CTD 特征 (化学-基因互作聚合, 10 dims) ----
    if ctd_gene_feat is not None:
        # ctd_gene_feat 按 sorted(ctd_gene_set) 排列, 直接按基因名对齐到 all_genes_list
        # ctd_gene_set 是 gene_feature_names ∪ drug_targets ∪ disease_genes 的超集
        ctd_gene_set = gene_feat_set | set(drug_target_genes) | set(disease_genes_raw)
        sorted_ctd_genes = sorted(ctd_gene_set)
        if ctd_gene_feat.shape[0] == len(sorted_ctd_genes):
            ctd_gene_dict = {g: ctd_gene_feat[i] for i, g in enumerate(sorted_ctd_genes)}
            ctd_gene_agg = np.zeros((n_genes, CTD_GENE_FEAT_DIM), dtype=np.float32)
            ctd_matched = 0
            for i, g in enumerate(all_genes_list):
                if g in ctd_gene_dict:
                    ctd_gene_agg[i] = ctd_gene_dict[g]
                    ctd_matched += 1
            gene_feat = np.concatenate([gene_feat, ctd_gene_agg], axis=1)
            logger.info(f"[Build] CTD gene features appended: {ctd_matched}/{n_genes} genes, "
                  f"+{CTD_GENE_FEAT_DIM} dims → total {gene_feat.shape[1]} dims")
        else:
            logger.info(f"[Build] CTD gene feature shape mismatch: "
                  f"{ctd_gene_feat.shape[0]} vs {len(sorted_ctd_genes)} (will skip)")
    else:
        logger.info(f"[Build] No CTD gene features (will skip)")

    # ---- 药物特征 ----
    drug_feat = drug_fingerprint_arr.reshape(1, -1)

    # ---- 疾病特征：仅用匹配到的疾病基因均值 ----
    disease_in_graph = [g for g in disease_genes_raw if g in gene_to_idx]
    if len(disease_in_graph) > 0:
        d_indices = [gene_to_idx[g] for g in disease_in_graph]
        disease_feat = np.mean(gene_feat[d_indices], axis=0, keepdims=True)
    else:
        disease_feat = np.mean(gene_feat, axis=0, keepdims=True)
    logger.info(f"[Build] Disease feature from {len(disease_in_graph)} matched genes")

    # ---- 构建 HeteroData ----
    data = HeteroData()
    data["drug"].x = torch.from_numpy(drug_feat)
    data["gene"].x = torch.from_numpy(gene_feat)
    data["disease"].x = torch.from_numpy(disease_feat)

    # ---- 药物-靶点边（仅保留在图中且在 subgraph 有特征的基因） ----
    dt_src, dt_dst = [], []
    for g in drug_target_genes:
        if g in gene_to_idx:
            dt_src.append(0)
            dt_dst.append(gene_to_idx[g])
    data["drug", "targets", "gene"].edge_index = torch.tensor(
        [dt_src, dt_dst], dtype=torch.long
    )
    logger.info(f"[Build] drug->gene edges: {len(dt_src)}")

    # ---- 反向边: gene → targeted_by → drug (药物节点接收靶基因聚合信息, 解决药物孤立) ----
    data["gene", "targeted_by", "drug"].edge_index = torch.tensor(
        [dt_dst, dt_src], dtype=torch.long
    )
    logger.info(f"[Build] gene->drug reverse edges: {len(dt_src)}")

    # ---- 靶点-疾病边 ----
    td_src, td_dst = [], []
    for g in disease_genes_raw:
        if g in gene_to_idx:
            td_src.append(gene_to_idx[g])
            td_dst.append(0)
    data["gene", "associated_with", "disease"].edge_index = torch.tensor(
        [td_src, td_dst], dtype=torch.long
    )
    logger.info(f"[Build] gene->disease edges: {len(td_src)}")

    # ---- PPI 边（仅保留两端均在图中的边） ----
    ppi_src, ppi_dst = [], []
    for a, b in ppi_edges:
        if a in gene_to_idx and b in gene_to_idx:
            ppi_src.append(gene_to_idx[a])
            ppi_dst.append(gene_to_idx[b])
    data["gene", "interacts", "gene"].edge_index = torch.tensor(
        [ppi_src, ppi_dst], dtype=torch.long
    )
    # 双向 PPI
    rev = torch.stack([data["gene", "interacts", "gene"].edge_index[1],
                       data["gene", "interacts", "gene"].edge_index[0]], dim=0)
    data["gene", "interacts", "gene"].edge_index = torch.cat([
        data["gene", "interacts", "gene"].edge_index, rev
    ], dim=1)
    logger.info(f"[Build] PPI edges (bidirectional): {len(ppi_src) * 2}")

    # ---- 共表达边 (ARCHS4, |corr| > 0.7) ----
    if coexp_edges and len(coexp_edges) > 0 and use_coexp:
        coexp_src, coexp_dst, coexp_weights = [], [], []
        for g1, g2, corr in coexp_edges:
            if g1 in gene_to_idx and g2 in gene_to_idx:
                coexp_src.append(gene_to_idx[g1])
                coexp_dst.append(gene_to_idx[g2])
                coexp_weights.append(abs(corr))  # 使用绝对值作为权重
        if coexp_src:
            data["gene", "coexpressed_with", "gene"].edge_index = torch.tensor(
                [coexp_src, coexp_dst], dtype=torch.long)
            # 存储边权重 (用于 SigmoidAttnConv 的 edge_attr)
            data["gene", "coexpressed_with", "gene"].edge_weight = torch.tensor(
                coexp_weights, dtype=torch.float32)
            # 双向
            rev_coexp = torch.stack([
                data["gene", "coexpressed_with", "gene"].edge_index[1],
                data["gene", "coexpressed_with", "gene"].edge_index[0]], dim=0)
            data["gene", "coexpressed_with", "gene"].edge_index = torch.cat([
                data["gene", "coexpressed_with", "gene"].edge_index, rev_coexp
            ], dim=1)
            data["gene", "coexpressed_with", "gene"].edge_weight = torch.cat([
                data["gene", "coexpressed_with", "gene"].edge_weight,
                data["gene", "coexpressed_with", "gene"].edge_weight
            ], dim=0)
        logger.info(f"[Build] Co-expression edges (bidirectional): {len(coexp_src) * 2}")
    else:
        logger.info(f"[Build] No co-expression edges (will skip)")

    # ---- TF 节点 + 边 ----
    if tf_edges and len(tf_edges) > 0 and tf_nodes:
        tf_name_to_idx = {}
        tf_feat_list = []
        for tf_name, info in tf_nodes.items():
            if tf_name not in tf_name_to_idx:
                tf_name_to_idx[tf_name] = len(tf_name_to_idx)
                tf_feat_list.append([
                    float(info['degree']), float(info['activation']),
                    float(info['repression']), float(info['unknown'])
                ])
        tf_feat = np.array(tf_feat_list, dtype=np.float32)
        tf_feat = StandardScaler().fit_transform(tf_feat)

        tf_reg_src, tf_reg_dst = [], []
        for tf_name, gene_name in tf_edges:
            if tf_name in tf_name_to_idx and gene_name in gene_to_idx:
                tf_reg_src.append(tf_name_to_idx[tf_name])
                tf_reg_dst.append(gene_to_idx[gene_name])

        data["tf"].x = torch.from_numpy(tf_feat)
        if tf_reg_src:
            data["tf", "regulates", "gene"].edge_index = torch.tensor(
                [tf_reg_src, tf_reg_dst], dtype=torch.long)
            data["gene", "regulated_by", "tf"].edge_index = torch.tensor(
                [tf_reg_dst, tf_reg_src], dtype=torch.long)
        logger.info(f"[Build] TF nodes: {len(tf_name_to_idx)}, "
              f"TF→gene edges: {len(tf_reg_src)}")
    else:
        logger.info(f"[Build] No TF edges (will skip)")

    # ---- Pathway 节点 + 边 ----
    if pw_edges and len(pw_edges) > 0 and pw_nodes:
        pw_name_to_idx = {}
        for pw_name in pw_nodes:
            if pw_name not in pw_name_to_idx:
                pw_name_to_idx[pw_name] = len(pw_name_to_idx)

        # 通路节点特征 (PCA 维度)
        if pw_feat_aligned is not None:
            pw_node_feat = pw_feat_aligned
        else:
            # 回退: 随机正交初始化 (保证所有节点非零且公平)
            # 避免 np.eye(n, d) 在 n > d 时后 n-d 行全为零的问题
            n_pw = len(pw_name_to_idx)
            pw_node_feat = np.random.randn(n_pw, PATHWAY_PCA_DIM).astype(np.float32) * 0.02
            logger.info(f"[Build] Pathway features fallback: random init ({n_pw}x{PATHWAY_PCA_DIM})")

        pw_src, pw_dst = [], []
        for gene, pw_name in pw_edges:
            if gene in gene_to_idx and pw_name in pw_name_to_idx:
                pw_src.append(gene_to_idx[gene])
                pw_dst.append(pw_name_to_idx[pw_name])

        data["pathway"].x = torch.from_numpy(pw_node_feat)
        if pw_src:
            data["gene", "in_pathway", "pathway"].edge_index = torch.tensor(
                [pw_src, pw_dst], dtype=torch.long)
            data["pathway", "has_gene", "gene"].edge_index = torch.tensor(
                [pw_dst, pw_src], dtype=torch.long)
        logger.info(f"[Build] Pathway nodes: {len(pw_name_to_idx)}, "
              f"gene→pathway edges: {len(pw_src)}")
    else:
        logger.info(f"[Build] No pathway edges (will skip)")

    # ---- Chemical 节点 + 边 (CTD) ----
    if ctd_edges and len(ctd_edges) > 0 and ctd_chem_nodes:
        chem_name_to_idx = {}
        chem_feat_list = []
        for chem_name, info in ctd_chem_nodes.items():
            if chem_name not in chem_name_to_idx:
                chem_name_to_idx[chem_name] = len(chem_name_to_idx)
                # 使用 load_ctd_edges 提取的 16 维丰富特征
                chem_feat_list.append(info['features'])
        chem_feat = np.array(chem_feat_list, dtype=np.float32)
        # StandardScaler 安全防护: 单样本或零方差列时回退
        if chem_feat.shape[0] > 1:
            scaler = StandardScaler()
            chem_feat = scaler.fit_transform(chem_feat)
            if np.any(np.isnan(chem_feat)) or np.any(np.isinf(chem_feat)):
                logger.warning("[Build] StandardScaler produced NaN/Inf on chem features, falling back to raw")
                chem_feat = np.array(chem_feat_list, dtype=np.float32)
        else:
            logger.info("[Build] Single chemical node, skipping StandardScaler")

        chem_src, chem_dst = [], []
        # Step 2 聚焦过滤: 仅保留 drug_targets ∪ disease_genes 范围内的 chem→gene 边
        focused_gene_set = set(drug_target_genes) | set(disease_genes_raw) if CTD_FOCUSED_EDGES else None
        for chem_name, gene_name in ctd_edges:
            if chem_name in chem_name_to_idx and gene_name in gene_to_idx:
                if focused_gene_set is not None and gene_name not in focused_gene_set:
                    continue
                chem_src.append(chem_name_to_idx[chem_name])
                chem_dst.append(gene_to_idx[gene_name])

        data["chemical"].x = torch.from_numpy(chem_feat)
        if chem_src:
            data["chemical", "interacts_with", "gene"].edge_index = torch.tensor(
                [chem_src, chem_dst], dtype=torch.long)
            data["gene", "targeted_by_chemical", "chemical"].edge_index = torch.tensor(
                [chem_dst, chem_src], dtype=torch.long)
        logger.info(f"[Build] Chemical nodes: {len(chem_name_to_idx)}, "
              f"chem→gene edges: {len(chem_src)}"
              + (f" (focused: drug+disease genes only)" if CTD_FOCUSED_EDGES else ""))
    else:
        logger.info(f"[Build] No CTD edges (will skip)")

    logger.info(f"[Build] HeteroData: drug=1, gene={n_genes}, disease=1 nodes")
    logger.info(f"[Build] Edge types: {list(data.edge_types)}")
    return data, gene_to_idx


# ============================================================================
# 3. 高效负采样
# ============================================================================

class NegativeSampler:
    """
    节点度数感知负采样器 (Degree-Aware Negative Sampling)

    依据 Cappelletti et al. (Bioinform Adv 2024, PMID:38577542):
    低度节点更可能是真负样本→赋予更高采样权重；
    高度节点可能是潜在假阴性（未发现的关联）→赋予更低采样权重。

    同时保留 NSCaching (Zhang et al., EMNLP 2019) 的预缓存思想：
    用 numpy 预计算权重池 + shuffle 顺序取用，复杂度 O(K)。

    参数:
        candidate_indices: 候选负样本索引 (PPI 1-hop 邻居排除后)
        n_genes: 总基因数
        degree_power: 度数幂权重 (0=均匀, >0=低度偏好, <0=高度偏好)
    """
    def __init__(self, candidate_indices, n_genes, degree_power=NEG_DEGREE_POWER,
                 ppi_node_degrees=None):
        """
        参数:
            candidate_indices: 候选负样本索引
            n_genes: 总基因数
            degree_power: 度数幂权重，0=均匀采样，>0=低度偏好
            ppi_node_degrees: (n_genes,) 或 None — PPI图真实度数 (Cappelletti 2024)
        """
        if candidate_indices is not None and len(candidate_indices) > 0:
            self.pool = np.asarray(candidate_indices, dtype=np.int64)
        else:
            self.pool = np.arange(n_genes, dtype=np.int64)

        # 节点度数感知权重: 使用真实PPI度数 (Cappelletti 2024, PMID:38577542)
        if degree_power > 0 and len(self.pool) > 1:
            if ppi_node_degrees is not None and len(ppi_node_degrees) == n_genes:
                candidate_degrees = ppi_node_degrees[self.pool]
                inv_degree = 1.0 / (1.0 + candidate_degrees)
            else:
                # 回退: 基于 pool 索引的近似（v1 兼容）
                degree_order = np.argsort(np.argsort(self.pool))
                inv_degree = 1.0 / (1.0 + degree_order / max(1, len(self.pool)))
            weights = inv_degree ** degree_power
            weights /= weights.sum()
            self.weights = weights
            self.use_weighted = True
        else:
            self.weights = None
            self.use_weighted = False

        self._reshuffle()

    def _reshuffle(self):
        """基于权重重新排列池"""
        if self.use_weighted and len(self.pool) > 1:
            # 按权重采样排列: 高权重(低度节点)先出现
            sampled_idx = np.random.choice(
                len(self.pool), size=len(self.pool), replace=False,
                p=self.weights
            )
            self.pool = self.pool[sampled_idx]
        else:
            np.random.shuffle(self.pool)
        self.pointer = 0

    def sample(self, num):
        if num > len(self.pool):
            return np.random.choice(self.pool, num, replace=True)
        if self.pointer + num > len(self.pool):
            self._reshuffle()
        batch = self.pool[self.pointer:self.pointer + num]
        self.pointer += num
        return batch


def build_ppi_neighbor_set(ppi_edge_index, gene_range):
    """
    从 PPI 边索引构建邻接表，用于约束负采样范围。
    返回每个基因的 PPI 1-hop 邻居集合。
    """
    neighbors = {i: set() for i in range(gene_range)}
    src, dst = ppi_edge_index[0].tolist(), ppi_edge_index[1].tolist()
    for s, d in zip(src, dst):
        if s < gene_range and d < gene_range:
            neighbors[s].add(d)
    return neighbors


def focal_bce_with_logits(logits, labels, gamma=1.5, alpha=0.25):
    """
    Focal BCE Loss — 自动聚焦难分类样本

    float16 安全: 内部用 float32 计算避免 sigmoid 饱和导致 -log(0)=inf→NaN
    """
    # 钳制 logits 防止极端值
    logits = torch.clamp(logits, -LOGIT_CLAMP, LOGIT_CLAMP)
    # float16 安全: 转为 float32 计算 loss, 再转回原始 dtype
    logits_f32 = logits.float()
    labels_f32 = labels.float()
    bce = F.binary_cross_entropy_with_logits(logits_f32, labels_f32, reduction="none")
    pt = torch.where(labels_f32 == 1, torch.sigmoid(logits_f32), 1 - torch.sigmoid(logits_f32))
    focal_weight = (1 - pt) ** gamma
    alpha_weight = torch.where(labels_f32 == 1, alpha, 1 - alpha)
    loss = (alpha_weight * focal_weight * bce).mean()
    return loss.to(logits.dtype)


def safe_ranking_loss(pos_scores, neg_scores, margin=0.3):
    """
    轻量 pairwise ranking loss — 每个正样本与其对应的负样本比较
    比 batch-all 更稳定，比 hard-mining 更简单
    """
    n_pos, n_neg = pos_scores.size(0), neg_scores.size(0)
    if n_pos == 0 or n_neg == 0:
        return torch.tensor(0.0, device=pos_scores.device)
    n = min(n_pos, n_neg)
    loss = torch.relu(margin - (pos_scores[:n] - neg_scores[:n]))
    return loss.mean()


def sample_negative_edges_via_set(pos_edge_index, n_genes, sampler, is_dt):
    """
    从 NegativeSampler 中采样与 pos 不重叠的负边。
    完全独立于 HeteroData，无依赖泄漏。
    """
    pos_set = set(zip(pos_edge_index[0].tolist(), pos_edge_index[1].tolist()))
    neg_src, neg_dst = [], []
    max_attempts = pos_edge_index.size(1) * 5

    while len(neg_src) < pos_edge_index.size(1) and max_attempts > 0:
        max_attempts -= 1
        batch = sampler.sample(min(pos_edge_index.size(1) * 2, 1000))
        for t in batch:
            if len(neg_src) >= pos_edge_index.size(1):
                break
            key = (0, int(t)) if is_dt else (int(t), 0)
            if key not in pos_set:
                neg_src.append(key[0])
                neg_dst.append(key[1])
                pos_set.add(key)

    return torch.tensor([neg_src, neg_dst], dtype=torch.long, device=device)


# ============================================================================
# 3b. SWA 公共函数 + DropEdge 调度
# ============================================================================

def update_swa_model(swa_model, model, swa_n):
    """
    更新 SWA (Stochastic Weight Averaging) 滑动平均

    在 SWA_START 轮后开始累积权重，
    训练结束后将累积和除以 swa_n 得到平均权重。
    """
    if swa_model is None:
        swa_model = copy.deepcopy(model)
        swa_n = 1
    else:
        swa_n += 1
        swa_sd = swa_model.state_dict()
        model_sd = model.state_dict()
        for key in swa_sd:
            if swa_sd[key].dtype in (torch.float32, torch.float64, torch.float16):
                swa_sd[key].data.add_(model_sd[key].data)
    return swa_model, swa_n


def finalize_swa(swa_model, swa_n):
    """训练结束后，将 SWA 累积和除以计数，返回平均模型"""
    if swa_model is None or swa_n <= 1:
        return swa_model
    swa_sd = swa_model.state_dict()
    for key in swa_sd:
        if swa_sd[key].dtype in (torch.float32, torch.float64, torch.float16):
            swa_sd[key].data.div_(swa_n)
    return swa_model


def ensemble_rrf_fusion(all_run_scores, n_runs, k_rrf=60):
    """
    RRF (Reciprocal Rank Fusion) 集成融合
    - 对每次运行分别计算 DT/TD 排名 → RRF
    - 返回融合后的 DataFrame
    """
    all_dt_ranks = []
    all_td_ranks = []
    rrf_list = []

    for run_df in all_run_scores:
        dt_rank = run_df["drug_target_score"].rank(ascending=False)
        td_rank = run_df["target_disease_score"].rank(ascending=False)
        rrf = 1.0 / (k_rrf + dt_rank) + 1.0 / (k_rrf + td_rank)
        all_dt_ranks.append(dt_rank.values)
        all_td_ranks.append(td_rank.values)
        rrf_list.append(rrf.values)

    merged = all_run_scores[0][["gene_symbol"]].copy()
    merged["drug_target_score"] = np.mean(
        [r["drug_target_score"].values for r in all_run_scores], axis=0)
    merged["target_disease_score"] = np.mean(
        [r["target_disease_score"].values for r in all_run_scores], axis=0)
    merged["combined_score"] = np.mean(
        [r["combined_score"].values for r in all_run_scores], axis=0)
    merged["dt_mean_rank"] = np.mean(all_dt_ranks, axis=0)
    merged["td_mean_rank"] = np.mean(all_td_ranks, axis=0)
    merged["rrf_score"] = np.mean(rrf_list, axis=0)

    merged = merged.sort_values("rrf_score", ascending=False).reset_index(drop=True)
    return merged


def compute_ensemble_jaccard(all_run_scores, n_runs, k_rrf=60, top_k=20):
    """
    计算集成运行间 Top-K Jaccard 相似度
    衡量排名稳定性 (>0.7 = 稳定)
    """
    top_sets = []
    for run_df in all_run_scores:
        dt_rank = run_df["drug_target_score"].rank(ascending=False)
        td_rank = run_df["target_disease_score"].rank(ascending=False)
        rrf = 1.0 / (k_rrf + dt_rank) + 1.0 / (k_rrf + td_rank)
        top_sets.append(
            set(run_df.iloc[np.argsort(-rrf.values)[:top_k]]["gene_symbol"]))

    jaccards = []
    for i in range(n_runs):
        for j in range(i + 1, n_runs):
            si, sj = top_sets[i], top_sets[j]
            jaccards.append(len(si & sj) / max(len(si | sj), 1))

    return np.mean(jaccards) if jaccards else 1.0


def get_dropedge_p(epoch, epochs=200,
                   min_p=DROPEDGE_P_MIN, max_p=DROPEDGE_P):
    """
    余弦退火 DropEdge 概率调度 (Rong et al. 2020 ICLR)

    训练初期大丢弃率(抗过平滑, 学习粗粒度拓扑),
    训练后期小丢弃率(保留全图结构, 学习细粒度模式).

    范围: [min_p, max_p], cosine 从 max_p 降为 min_p.
    """
    return min_p + 0.5 * (max_p - min_p) * (
        1 + math.cos(math.pi * epoch / epochs)
    )


def get_edge_index_dict_for_fold(data, dt_train_idx, td_train_idx):
    """从 HeteroData 提取仅含训练边的 edge_index_dict（消除验证边消息泄漏）"""
    edge_dict = {}
    for etype in data.edge_types:
        edge_dict[etype] = data[etype].edge_index.clone()
    edge_dict[("drug", "targets", "gene")] = \
        data["drug", "targets", "gene"].edge_index[:, dt_train_idx].clone()
    # 反向边: gene→drug 也必须裁剪验证边（否则药物嵌入泄漏验证集靶点信息）
    rev_full = data["gene", "targeted_by", "drug"].edge_index
    edge_dict[("gene", "targeted_by", "drug")] = torch.stack([
        rev_full[0, dt_train_idx], rev_full[1, dt_train_idx]
    ], dim=0).clone()
    edge_dict[("gene", "associated_with", "disease")] = \
        data["gene", "associated_with", "disease"].edge_index[:, td_train_idx].clone()
    return edge_dict


# ============================================================================
# 4. SigmoidAttnConv (MRHormer TRLA) — 替换 GATv2Conv 的 Softmax
# ============================================================================

class SigmoidAttnConv(MessagePassing):
    """
    MRHormer TRLA 注意力卷积 — Sigmoid 替代 Softmax

    核心改进（相对于 GATv2Conv）：
    1. 类型特定投影: 源节点和目标节点独立 Linear 投影
    2. Sigmoid 注意力: 每个邻居独立计算重要性，不做归一化
       → 避免 Softmax 在异质特征下过度集中于某一类型节点
    3. 保留多头机制和自环（add_self_loops）

    Reference: MRHormer (TRLA), BioRxiv 2024
    """
    def __init__(self, in_channels, out_channels, heads=1, concat=True,
                 dropout=0.0, add_self_loops=True, bias=True, use_edge_attr=False):
        # in_channels can be int or tuple (src_dim, dst_dim)
        if isinstance(in_channels, int):
            in_src = in_dst = in_channels
        else:
            in_src, in_dst = in_channels
            if in_src == -1:
                in_src = in_dst
            if in_dst == -1:
                in_dst = in_src

        super().__init__(aggr='add', node_dim=0)

        self.in_src = in_src
        self.in_dst = in_dst
        self.out_channels = out_channels
        self.heads = heads
        self.concat = concat
        self.dropout = dropout
        self.add_self_loops = add_self_loops
        self.use_edge_attr = use_edge_attr

        # 类型特定投影矩阵
        self.lin_src = Linear(in_src, heads * out_channels, bias=bias)
        self.lin_dst = Linear(in_dst, heads * out_channels, bias=bias)

        # 注意力参数 (每个头独立)
        # LayerScale 初始化 (Apple ICLR 2025): 用 0.1 缩放避免早期大注意力范数
        self.att_src = Parameter(torch.empty(1, heads, out_channels))
        self.att_dst = Parameter(torch.empty(1, heads, out_channels))

        # 边权重编码器: 将 1 维 edge_weight 映射为多头注意力偏置
        if use_edge_attr:
            self.edge_encoder = Linear(1, heads, bias=False)

        if bias and concat:
            self.bias = Parameter(torch.empty(heads * out_channels))
        elif bias and not concat:
            self.bias = Parameter(torch.empty(out_channels))
        else:
            self.register_parameter('bias', None)

        # ---- 注意力存储（可解释性） ----
        self._store_attention = False     # 推理时设为 True 启用保存
        self._saved_edge_index = None     # 当前前向的边索引
        self._saved_alphas = None         # [E, heads] 每条边的多头注意力系数

        self.reset_parameters()

    def reset_parameters(self):
        super().reset_parameters()
        self.lin_src.reset_parameters()
        self.lin_dst.reset_parameters()
        # LayerScale 初始化: 缩放因子 0.1 防止早期注意力范数过大 (Apple, ICLR 2025)
        torch.nn.init.xavier_uniform_(self.att_src)
        torch.nn.init.xavier_uniform_(self.att_dst)
        with torch.no_grad():
            self.att_src.mul_(0.1)
            self.att_dst.mul_(0.1)
        if self.bias is not None:
            torch.nn.init.zeros_(self.bias)

    def enable_attention_store(self, enable=True):
        """启用/禁用注意力系数存储（推理时用于可解释性分析）"""
        self._store_attention = enable

    def forward(self, x, edge_index, edge_weight=None):
        # x 可以是 tuple (x_src, x_dst) 用于二分图
        if isinstance(x, tuple):
            x_src, x_dst = x
        else:
            x_src = x_dst = x

        # 保存当前边索引用于可解释性
        if self._store_attention:
            self._saved_edge_index = edge_index.detach().cpu()
            self._saved_alphas = None

        # 类型特定投影
        x_src = self.lin_src(x_src).view(-1, self.heads, self.out_channels)
        x_dst = self.lin_dst(x_dst).view(-1, self.heads, self.out_channels)

        # 添加自环（兼容 PyG add_remaining_self_loops 语义）
        # 仅对尚未有自环的节点添加，避免重复边
        if self.add_self_loops:
            if isinstance(edge_index, tuple):
                raise NotImplementedError("Self-loops not supported for bipartite graphs")
            num_nodes = x_dst.size(0)
            # 检测已有自环，避免重复
            existing_self_loops = (edge_index[0] == edge_index[1])
            has_self_loop = torch.zeros(num_nodes, dtype=torch.bool, device=edge_index.device)
            has_self_loop[edge_index[0, existing_self_loops]] = True
            missing = torch.arange(num_nodes, device=edge_index.device)[~has_self_loop]
            if missing.numel() > 0:
                loop_index = torch.stack([missing, missing], dim=0)
                edge_index = torch.cat([edge_index, loop_index], dim=1)
                # 自环边用 1.0 权重填充
                if edge_weight is not None:
                    loop_weight = torch.ones(missing.numel(), device=edge_weight.device, dtype=edge_weight.dtype)
                    edge_weight = torch.cat([edge_weight, loop_weight], dim=0)
                # 更新保存的边索引以包含自环
                if self._store_attention:
                    self._saved_edge_index = edge_index.detach().cpu()

        # 将 1-D edge_weight 转为 2-D [E, 1] 供 message 使用
        if edge_weight is not None and edge_weight.dim() == 1:
            edge_weight = edge_weight.unsqueeze(-1)

        # 消息传递
        out = self.propagate(edge_index, x=(x_src, x_dst), size=None, edge_weight=edge_weight)

        if self.concat:
            out = out.view(-1, self.heads * self.out_channels)
        else:
            out = out.mean(dim=1)

        if self.bias is not None:
            out = out + self.bias

        return out

    def message(self, x_i, x_j, index, ptr, size_i, edge_weight=None):
        """
        Sigmoid 注意力计算

        x_i: 目标节点嵌入 [E, heads, out_channels]
        x_j: 源节点嵌入   [E, heads, out_channels]
        edge_weight: 边先验权重 [E, 1] 或 None (用于共表达相关系数等)
        """
        # 注意力系数: LeakyReLU((x_i·att_src + x_j·att_dst))
        alpha_src = (x_i * self.att_src).sum(dim=-1)  # [E, heads]
        alpha_dst = (x_j * self.att_dst).sum(dim=-1)  # [E, heads]
        alpha = alpha_src + alpha_dst

        # 边权重编码: 将先验权重 (如 ARCHS4 相关系数) 融入注意力偏置
        if edge_weight is not None and self.use_edge_attr:
            edge_bias = self.edge_encoder(edge_weight)  # [E, heads]
            alpha = alpha + edge_bias

        # SIGMOID 替代 Softmax — 每个邻居独立计算重要性
        alpha = torch.sigmoid(alpha)
        alpha = F.dropout(alpha, p=self.dropout, training=self.training)

        # [可解释性] 保存注意力系数（仅在推理时）
        if self._store_attention and not self.training:
            self._saved_alphas = alpha.detach().cpu()  # [E, heads]

        # 加权消息
        return x_j * alpha.unsqueeze(-1)


# ============================================================================
# 4b. 异构图 GAT 模型 (MRHormer + MLGANN)
# ============================================================================

class HeteroGAT(torch.nn.Module):
    """
    异构图 GAT 模型 — MRHormer TRLA + MLGANN 多层注意力池化

    改进 (v4):
      - SigmoidAttnConv: 替换 GATv2Conv 的 Softmax 为 Sigmoid 注意力
        → 类型特定投影矩阵 + 独立邻居重要性 (避免异质特征下注意集中)
      - MLGANN 多层池化: 对所有 GNN 层输出做可学习注意力加权融合
        → 保留浅层局部信息 + 深层全局信息
      - 支持 2~3 层 GNN (N_LAYERS): 3层捕获更高级邻域关系 (MSH-DTI 2024)
      - LayerNorm 残差连接: 稳定深层特征分布 (Apple ICLR 2025)
    """
    def __init__(self, metadata, drug_dim, gene_dim, disease_dim,
                 hidden_dim, out_dim, heads, dropout, n_layers=N_LAYERS):
        super().__init__()
        self.drug_proj = torch.nn.Linear(drug_dim, hidden_dim)
        self.gene_proj = torch.nn.Linear(gene_dim, hidden_dim)
        self.disease_proj = torch.nn.Linear(disease_dim, hidden_dim)
        self.tf_proj = torch.nn.Linear(4, hidden_dim)           # TF: degree, activation, repression, unknown
        self.pathway_proj = torch.nn.Linear(PATHWAY_PCA_DIM, hidden_dim)
        self.chemical_proj = torch.nn.Linear(CTD_CHEM_FEAT_DIM, hidden_dim)  # chemical: 24 维
        self.dropout_rate = dropout
        self.drug_out = torch.nn.Linear(hidden_dim, out_dim)
        self.out_dim = out_dim
        self.n_layers = n_layers

        # 残差连接投影层 + LayerNorm
        self.gene_res = torch.nn.Linear(hidden_dim, out_dim)
        self.disease_res = torch.nn.Linear(hidden_dim, out_dim)
        self.drug_res = torch.nn.Linear(hidden_dim, out_dim)
        self.tf_res = torch.nn.Linear(hidden_dim, out_dim)
        self.pathway_res = torch.nn.Linear(hidden_dim, out_dim)
        self.chemical_res = torch.nn.Linear(hidden_dim, out_dim)
        self.gene_norm = LayerNorm(out_dim)
        self.disease_norm = LayerNorm(out_dim)
        self.drug_norm = LayerNorm(out_dim)
        self.tf_norm = LayerNorm(out_dim)
        self.pathway_norm = LayerNorm(out_dim)
        self.chemical_norm = LayerNorm(out_dim)

        hid_out = hidden_dim // heads
        conv_dropout = dropout * 0.6

        conv1_in = hidden_dim
        conv1_out = heads * hid_out
        conv2_in = conv1_out
        conv2_out = out_dim
        conv3_in = out_dim
        conv3_out = out_dim

        # ---- MLGANN: 多层注意力池化 ----
        # 为所有 GNN 层输出保留独立投影 (MLGANN, Lu et al., Sci Rep 2024)
        num_layers_to_pool = n_layers  # 所有层都参与池化
        self.layer1_proj = torch.nn.ModuleDict()
        self.layer_attn = torch.nn.ModuleDict()
        self.mlgann_norm = torch.nn.ModuleDict()
        self.layer_scalar = torch.nn.ParameterDict()

        for node_type in ["drug", "gene", "disease", "tf", "pathway", "chemical"]:
            self.layer1_proj[node_type] = torch.nn.Linear(conv1_out, out_dim)
            self.layer_attn[node_type] = torch.nn.Linear(out_dim, 1)
            self.mlgann_norm[node_type] = LayerNorm(out_dim)
            self.layer_scalar[node_type] = torch.nn.Parameter(torch.zeros(num_layers_to_pool))

        # ---- Conv1: SigmoidAttnConv (MRHormer TRLA) ----
        self.conv1 = HeteroConv({
            ("drug", "targets", "gene"):
                SigmoidAttnConv((conv1_in, conv1_in), hid_out, heads=heads, concat=True,
                                dropout=conv_dropout, add_self_loops=False),
            ("gene", "targeted_by", "drug"):   # 反向边: 药物接收靶基因聚合
                SigmoidAttnConv((conv1_in, conv1_in), hid_out, heads=heads, concat=True,
                                dropout=conv_dropout, add_self_loops=False),
            ("gene", "associated_with", "disease"):
                SigmoidAttnConv((conv1_in, conv1_in), hid_out, heads=heads, concat=True,
                                dropout=conv_dropout, add_self_loops=False),
            ("gene", "interacts", "gene"):
                SigmoidAttnConv(conv1_in, hid_out, heads=heads, concat=True,
                                dropout=conv_dropout, add_self_loops=True),
            # TF 调控
            ("tf", "regulates", "gene"):
                SigmoidAttnConv((conv1_in, conv1_in), hid_out, heads=heads, concat=True,
                                dropout=conv_dropout, add_self_loops=False),
            ("gene", "regulated_by", "tf"):
                SigmoidAttnConv((conv1_in, conv1_in), hid_out, heads=heads, concat=True,
                                dropout=conv_dropout, add_self_loops=False),
            # 共表达
            ("gene", "coexpressed_with", "gene"):
                SigmoidAttnConv(conv1_in, hid_out, heads=heads, concat=True,
                                dropout=conv_dropout, add_self_loops=True, use_edge_attr=True),
            # 通路
            ("gene", "in_pathway", "pathway"):
                SigmoidAttnConv((conv1_in, conv1_in), hid_out, heads=heads, concat=True,
                                dropout=conv_dropout, add_self_loops=False),
            ("pathway", "has_gene", "gene"):
                SigmoidAttnConv((conv1_in, conv1_in), hid_out, heads=heads, concat=True,
                                dropout=conv_dropout, add_self_loops=False),
            # 化学物
            ("chemical", "interacts_with", "gene"):
                SigmoidAttnConv((conv1_in, conv1_in), hid_out, heads=heads, concat=True,
                                dropout=conv_dropout, add_self_loops=False),
            ("gene", "targeted_by_chemical", "chemical"):
                SigmoidAttnConv((conv1_in, conv1_in), hid_out, heads=heads, concat=True,
                                dropout=conv_dropout, add_self_loops=False),
        }, aggr="mean")

        # ---- Conv2: SigmoidAttnConv (1 head, no concat) ----
        self.conv2 = HeteroConv({
            ("drug", "targets", "gene"):
                SigmoidAttnConv((conv2_in, conv2_in), conv2_out, heads=1, concat=False,
                                dropout=conv_dropout, add_self_loops=False),
            ("gene", "targeted_by", "drug"):
                SigmoidAttnConv((conv2_in, conv2_in), conv2_out, heads=1, concat=False,
                                dropout=conv_dropout, add_self_loops=False),
            ("gene", "associated_with", "disease"):
                SigmoidAttnConv((conv2_in, conv2_in), conv2_out, heads=1, concat=False,
                                dropout=conv_dropout, add_self_loops=False),
            ("gene", "interacts", "gene"):
                SigmoidAttnConv(conv2_in, conv2_out, heads=1, concat=False,
                                dropout=conv_dropout, add_self_loops=True),
            # TF 调控
            ("tf", "regulates", "gene"):
                SigmoidAttnConv((conv2_in, conv2_in), conv2_out, heads=1, concat=False,
                                dropout=conv_dropout, add_self_loops=False),
            ("gene", "regulated_by", "tf"):
                SigmoidAttnConv((conv2_in, conv2_in), conv2_out, heads=1, concat=False,
                                dropout=conv_dropout, add_self_loops=False),
            # 共表达
            ("gene", "coexpressed_with", "gene"):
                SigmoidAttnConv(conv2_in, conv2_out, heads=1, concat=False,
                                dropout=conv_dropout, add_self_loops=True, use_edge_attr=True),
            # 通路
            ("gene", "in_pathway", "pathway"):
                SigmoidAttnConv((conv2_in, conv2_in), conv2_out, heads=1, concat=False,
                                dropout=conv_dropout, add_self_loops=False),
            ("pathway", "has_gene", "gene"):
                SigmoidAttnConv((conv2_in, conv2_in), conv2_out, heads=1, concat=False,
                                dropout=conv_dropout, add_self_loops=False),
            # 化学物
            ("chemical", "interacts_with", "gene"):
                SigmoidAttnConv((conv2_in, conv2_in), conv2_out, heads=1, concat=False,
                                dropout=conv_dropout, add_self_loops=False),
            ("gene", "targeted_by_chemical", "chemical"):
                SigmoidAttnConv((conv2_in, conv2_in), conv2_out, heads=1, concat=False,
                                dropout=conv_dropout, add_self_loops=False),
        }, aggr="mean")

        # ---- Conv3 (可选, N_LAYERS=3 时启用): 深层特征提取 (MSH-DTI 2024) ----
        if n_layers >= 3:
            self.conv3 = HeteroConv({
                ("drug", "targets", "gene"):
                    SigmoidAttnConv((conv3_in, conv3_in), conv3_out, heads=1, concat=False,
                                    dropout=conv_dropout, add_self_loops=False),
                ("gene", "targeted_by", "drug"):   # 反向边: 药物接收靶基因聚合 (v3.1 补全)
                    SigmoidAttnConv((conv3_in, conv3_in), conv3_out, heads=1, concat=False,
                                    dropout=conv_dropout, add_self_loops=False),
                ("gene", "associated_with", "disease"):
                    SigmoidAttnConv((conv3_in, conv3_in), conv3_out, heads=1, concat=False,
                                    dropout=conv_dropout, add_self_loops=False),
                ("gene", "interacts", "gene"):
                    SigmoidAttnConv(conv3_in, conv3_out, heads=1, concat=False,
                                    dropout=conv_dropout, add_self_loops=True),
                # TF 调控
                ("tf", "regulates", "gene"):
                    SigmoidAttnConv((conv3_in, conv3_in), conv3_out, heads=1, concat=False,
                                    dropout=conv_dropout, add_self_loops=False),
                ("gene", "regulated_by", "tf"):
                    SigmoidAttnConv((conv3_in, conv3_in), conv3_out, heads=1, concat=False,
                                    dropout=conv_dropout, add_self_loops=False),
                # 共表达
                ("gene", "coexpressed_with", "gene"):
                    SigmoidAttnConv(conv3_in, conv3_out, heads=1, concat=False,
                                    dropout=conv_dropout, add_self_loops=True),
                # 通路
                ("gene", "in_pathway", "pathway"):
                    SigmoidAttnConv((conv3_in, conv3_in), conv3_out, heads=1, concat=False,
                                    dropout=conv_dropout, add_self_loops=False),
                ("pathway", "has_gene", "gene"):
                    SigmoidAttnConv((conv3_in, conv3_in), conv3_out, heads=1, concat=False,
                                    dropout=conv_dropout, add_self_loops=False),
                # 化学物
                ("chemical", "interacts_with", "gene"):
                    SigmoidAttnConv((conv3_in, conv3_in), conv3_out, heads=1, concat=False,
                                    dropout=conv_dropout, add_self_loops=False),
                ("gene", "targeted_by_chemical", "chemical"):
                    SigmoidAttnConv((conv3_in, conv3_in), conv3_out, heads=1, concat=False,
                                    dropout=conv_dropout, add_self_loops=False),
            }, aggr="mean")
            # Conv3 残差连接
            self.conv3_norm = LayerNorm(out_dim)
        else:
            self.conv3 = None

        self.dt_decoder = torch.nn.Sequential(
            torch.nn.Linear(out_dim * 2, out_dim),
            torch.nn.ReLU(),
            torch.nn.Dropout(dropout * 0.5),
            torch.nn.Linear(out_dim, 1),
        )
        self.td_decoder = torch.nn.Sequential(
            torch.nn.Linear(out_dim * 2, out_dim),
            torch.nn.ReLU(),
            torch.nn.Dropout(dropout * 0.5),
            torch.nn.Linear(out_dim, 1),
        )

    def forward(self, x_dict, edge_index_dict, edge_weight_dict=None):
        """
        MRHormer + MLGANN 前向传播 (支持 N_LAYERS)

        流程:
        1. 投影到 hidden_dim
        2. Conv1 (SigmoidAttnConv, 多头, concat) → ReLU → Dropout
        3. Conv2 (SigmoidAttnConv, 1头, no concat) → LayerNorm 残差
        4. [可选] Conv3 (SigmoidAttnConv, 更深层特征) → LayerNorm 残差
        5. MLGANN 多层注意力池化: 所有层 ← 逐节点注意力 + 逐类型标量

        edge_weight_dict: 可选, 边先验权重 (如 ARCHS4 共表达相关系数)
        """
        x_proj = {
            "drug": self.drug_proj(x_dict["drug"]),
            "gene": self.gene_proj(x_dict["gene"]),
            "disease": self.disease_proj(x_dict["disease"]),
        }
        # 动态添加新节点类型投影
        if "tf" in x_dict:
            x_proj["tf"] = self.tf_proj(x_dict["tf"])
        if "pathway" in x_dict:
            x_proj["pathway"] = self.pathway_proj(x_dict["pathway"])
        if "chemical" in x_dict:
            x_proj["chemical"] = self.chemical_proj(x_dict["chemical"])

        # ---- Layer 1: SigmoidAttnConv ----
        ew_dict = edge_weight_dict if edge_weight_dict is not None else {}
        out1 = self.conv1(x_proj, edge_index_dict, edge_weight_dict=ew_dict)
        for k in x_proj:
            if k not in out1:
                out1[k] = x_proj[k]
        out1 = {k: v.relu() for k, v in out1.items()}
        out1 = {k: F.dropout(v, p=self.dropout_rate, training=self.training)
                for k, v in out1.items()}

        # ---- Layer 2: SigmoidAttnConv ----
        out2 = self.conv2(out1, edge_index_dict, edge_weight_dict=ew_dict)
        for k in out1:
            if k not in out2:
                out2[k] = out1[k]

        # LayerNorm 残差连接: 稳定融合特征分布
        if "gene" in out2 and "gene" in x_proj:
            out2["gene"] = self.gene_norm(
                out2["gene"] + self.gene_res(x_proj["gene"]))
        if "disease" in out2 and "disease" in x_proj:
            out2["disease"] = self.disease_norm(
                out2["disease"] + self.disease_res(x_proj["disease"]))
        # 药物残差: 药物节点通过 gene→drug 反向边接收靶基因聚合后, 加入残差连接
        if "drug" in out2 and "drug" in x_proj:
            out2["drug"] = self.drug_norm(
                out2["drug"] + self.drug_res(x_proj["drug"]))
        # TF 残差
        if "tf" in out2 and "tf" in x_proj:
            out2["tf"] = self.tf_norm(
                out2["tf"] + self.tf_res(x_proj["tf"]))
        # Pathway 残差
        if "pathway" in out2 and "pathway" in x_proj:
            out2["pathway"] = self.pathway_norm(
                out2["pathway"] + self.pathway_res(x_proj["pathway"]))
        # Chemical 残差
        if "chemical" in out2 and "chemical" in x_proj:
            out2["chemical"] = self.chemical_norm(
                out2["chemical"] + self.chemical_res(x_proj["chemical"]))

        # ---- Layer 3 (可选): 更深层特征提取 (MSH-DTI 2024) ----
        if self.conv3 is not None:
            out3 = self.conv3(out2, edge_index_dict, edge_weight_dict=ew_dict)
            for k in out2:
                if k not in out3:
                    out3[k] = out2[k]
            # Conv3 残差连接 (从 out2 跳接)
            if "gene" in out3:
                out3["gene"] = self.conv3_norm(out3["gene"] + out2["gene"])
            if "disease" in out3:
                out3["disease"] = self.conv3_norm(out3["disease"] + out2["disease"])
            if "drug" in out3:
                out3["drug"] = self.conv3_norm(out3["drug"] + out2["drug"])
            if "tf" in out3:
                out3["tf"] = self.conv3_norm(out3["tf"] + out2["tf"])
            if "pathway" in out3:
                out3["pathway"] = self.conv3_norm(out3["pathway"] + out2["pathway"])
            if "chemical" in out3:
                out3["chemical"] = self.conv3_norm(out3["chemical"] + out2["chemical"])
        else:
            out3 = None

        # ---- MLGANN 多层注意力池化 ----
        # 将所有 GNN 层输出投影到 out_dim 并收集
        layer_outputs = {"out1": out1, "out2": out2}
        if out3 is not None:
            layer_outputs["out3"] = out3

        # 投影层1输出到 out_dim
        layer1_adapted = {}
        for k in out1:
            if k in out2 and k in self.layer1_proj:
                layer1_adapted[k] = self.layer1_proj[k](out1[k])

        # 融合所有层
        z_dict = {}
        for k in out2:
            # 收集该节点类型的所有层输出
            layers_k = {}
            if k in layer1_adapted:
                layers_k["l1"] = layer1_adapted[k]  # [N, out_dim]
            else:
                layers_k["l1"] = None
            layers_k["l2"] = out2[k]                 # [N, out_dim]
            if out3 is not None and k in out3:
                layers_k["l3"] = out3[k]
            else:
                layers_k["l3"] = None

            # LayerNorm + 注意力 + 标量融合
            processed = []
            scalars = F.softmax(self.layer_scalar[k], dim=0)  # [n_layers]
            idx = 0
            for lkey in ["l1", "l2", "l3"]:
                feat = layers_k[lkey]
                if feat is None:
                    continue
                feat_norm = self.mlgann_norm[k](feat)  # LayerNorm 稳定分布
                attn = torch.sigmoid(self.layer_attn[k](feat_norm))
                processed.append((attn, scalars[idx], feat_norm))
                idx += 1

            if len(processed) == 1:
                z_dict[k] = processed[0][2]  # 只有一层, 直接返回
            else:
                # 加权融合: 逐节点注意力 × 逐类型标量
                total_weight = sum(a * s for a, s, _ in processed) + 1e-8
                z_dict[k] = sum(
                    (a * s / total_weight) * f for a, s, f in processed
                )

        return z_dict

    @torch.no_grad()
    def diagnose_oversmoothing(self, x_dict, edge_index_dict):
        """
        过平滑诊断: 计算每层基因嵌入的标准差和余弦相似度。

        原理 (Li et al. 2024, Over-Smoothing in GNNs):
          - 层间标准差逐层下降 → 过平滑
          - 正/负样本余弦相似度趋同 → 区分度降低

        返回:
          layer_stds: {l1_std, l2_std} — 各层基因嵌入标准差
          pos_neg_cos: 正/负样本对余弦相似度对比
        """
        self.eval()
        x_proj = {
            "drug": self.drug_proj(x_dict["drug"]),
            "gene": self.gene_proj(x_dict["gene"]),
            "disease": self.disease_proj(x_dict["disease"]),
        }
        if "tf" in x_dict:
            x_proj["tf"] = self.tf_proj(x_dict["tf"])
        if "pathway" in x_dict:
            x_proj["pathway"] = self.pathway_proj(x_dict["pathway"])
        if "chemical" in x_dict:
            x_proj["chemical"] = self.chemical_proj(x_dict["chemical"])
        out1 = self.conv1(x_proj, edge_index_dict)
        for k in x_proj:
            if k not in out1:
                out1[k] = x_proj[k]

        out2 = self.conv2(out1, edge_index_dict)
        for k in out1:
            if k not in out2:
                out2[k] = out1[k]

        # 层间标准差诊断
        gene_std_l1 = float(out1["gene"].std().cpu())
        gene_std_l2 = float(out2["gene"].std().cpu())
        drug_std_l1 = float(out1.get("drug", x_proj["drug"]).std().cpu())
        drug_std_l2 = float(out2["drug"].std().cpu()) if "drug" in out2 else float(x_proj["drug"].std().cpu())

        logger.info(f"[Oversmoothing] gene_std: L1={gene_std_l1:.4f}, L2={gene_std_l2:.4f}, "
              f"ratio L2/L1={gene_std_l2/gene_std_l1:.3f}")
        logger.info(f"[Oversmoothing] drug_std: L1={drug_std_l1:.4f}, L2={drug_std_l2:.4f}")

        if gene_std_l2 / gene_std_l1 < 0.5:
            logger.warning(f"[Oversmoothing] ⚠ gene std dropped >50% across layers — severe oversmoothing")
        elif gene_std_l2 / gene_std_l1 < 0.7:
            logger.warning(f"[Oversmoothing] ⚠ gene std dropped ~30% — moderate oversmoothing")
        else:
            logger.info(f"[Oversmoothing] ✓ gene std stable — no significant oversmoothing")

        return gene_std_l1, gene_std_l2

    def decode_drug_target(self, z_dict, edge_index):
        feat = torch.cat([
            z_dict["drug"][edge_index[0]],
            z_dict["gene"][edge_index[1]],
        ], dim=-1)
        return self.dt_decoder(feat).squeeze(-1)

    def decode_target_disease(self, z_dict, edge_index):
        feat = torch.cat([
            z_dict["gene"][edge_index[0]],
            z_dict["disease"][edge_index[1]],
        ], dim=-1)
        return self.td_decoder(feat).squeeze(-1)

    def enable_attention_store(self, enable=True):
        """递归启用所有 SigmoidAttnConv 的注意力存储（推理时可解释性）"""
        for conv_name in ['conv1', 'conv2', 'conv3']:
            conv = getattr(self, conv_name, None)
            if conv is None:
                continue
            # HeteroConv.convs 是 ModuleDict[edge_type -> SigmoidAttnConv]
            for etype, subconv in conv.convs.items():
                if hasattr(subconv, 'enable_attention_store'):
                    subconv.enable_attention_store(enable)

    def get_attention_weights(self):
        """
        收集所有存储的注意力权重
        返回: {edge_type: {layer_name: (edge_index, alphas)}}
          - edge_index: [2, E] CPU 张量
          - alphas: [E, heads] CPU 张量, 每个头的注意力系数
        """
        attention = {}
        for layer_name in ['conv1', 'conv2', 'conv3']:
            conv = getattr(self, layer_name, None)
            if conv is None:
                continue
            for etype, subconv in conv.convs.items():
                if (hasattr(subconv, '_saved_alphas')
                        and subconv._saved_alphas is not None):
                    if etype not in attention:
                        attention[etype] = {}
                    attention[etype][layer_name] = (
                        subconv._saved_edge_index,
                        subconv._saved_alphas
                    )
        return attention


# ============================================================================
# 5. 训练 & 评估
# ============================================================================

def _get_coexp_edge_weight(hetero_data):
    """从 HeteroData 提取共表达边权重字典 (用于 SigmoidAttnConv 的 edge_weight)"""
    et = ("gene", "coexpressed_with", "gene")
    try:
        if et in hetero_data.edge_types and hasattr(hetero_data[et], "edge_weight"):
            w = hetero_data[et].edge_weight
            if w is not None and w.numel() > 0:
                return {et: w}
    except Exception:
        pass
    return None


def train_epoch(model, train_data, full_data_dt, full_data_td,
                dt_train_idx, td_train_idx, n_genes, optimizer,
                dt_sampler, td_sampler, scaler=None, dropedge_p=None):
    """
    单轮训练：Focal BCE + Safe Ranking Loss + DropEdge
    无验证边泄漏，无数值溢出风险

    dropedge_p: 若为 None 则使用全局 DROPEDGE_P，否则使用调度值
    """
    model.train()

    # ---- 药物-靶点 ----
    pos_dt = full_data_dt[:, dt_train_idx]
    neg_dt = sample_negative_edges_via_set(pos_dt, n_genes, dt_sampler, is_dt=True)
    dt_edge = torch.cat([pos_dt, neg_dt], dim=1)
    dt_label = torch.cat([
        torch.ones(pos_dt.size(1)), torch.zeros(neg_dt.size(1))
    ], dim=0).to(device)

    # ---- 靶点-疾病 ----
    pos_td = full_data_td[:, td_train_idx]
    neg_td = sample_negative_edges_via_set(pos_td, n_genes, td_sampler, is_dt=False)
    td_edge = torch.cat([pos_td, neg_td], dim=1)
    td_label = torch.cat([
        torch.ones(pos_td.size(1)), torch.zeros(neg_td.size(1))
    ], dim=0).to(device)

    # ---- DropEdge + Forward ----
    de_p = dropedge_p if dropedge_p is not None else DROPEDGE_P
    train_edge_dict = dict(train_data.edge_index_dict)
    train_edge_dict["gene", "interacts", "gene"] = drop_edge(
        train_data["gene", "interacts", "gene"].edge_index, de_p)

    if scaler is not None:
        optimizer.zero_grad()
        with torch.cuda.amp.autocast():
            z_dict = model(train_data.x_dict, train_edge_dict,
                           edge_weight_dict=_get_coexp_edge_weight(train_data))
            dt_logit = model.decode_drug_target(z_dict, dt_edge)
            td_logit = model.decode_target_disease(z_dict, td_edge)
            loss = (
                focal_bce_with_logits(dt_logit, dt_label, FOCAL_GAMMA, FOCAL_ALPHA) +
                TD_LOSS_WEIGHT * focal_bce_with_logits(td_logit, td_label, FOCAL_GAMMA, FOCAL_ALPHA) +
                RANK_WEIGHT * safe_ranking_loss(td_logit[:pos_td.size(1)], td_logit[pos_td.size(1):], RANK_MARGIN)
            )
        if torch.isnan(loss) or torch.isinf(loss):
            logger.warning(f"  [NaN/Inf] loss={loss.item()}, skipping update")
            return float('nan')
        scaler.scale(loss).backward()
        scaler.unscale_(optimizer)
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        scaler.step(optimizer)
        scaler.update()
    else:
        z_dict = model(train_data.x_dict, train_edge_dict,
                       edge_weight_dict=_get_coexp_edge_weight(train_data))
        dt_logit = model.decode_drug_target(z_dict, dt_edge)
        td_logit = model.decode_target_disease(z_dict, td_edge)
        loss = (
            focal_bce_with_logits(dt_logit, dt_label, FOCAL_GAMMA, FOCAL_ALPHA) +
            TD_LOSS_WEIGHT * focal_bce_with_logits(td_logit, td_label, FOCAL_GAMMA, FOCAL_ALPHA) +
            RANK_WEIGHT * safe_ranking_loss(td_logit[:pos_td.size(1)], td_logit[pos_td.size(1):], RANK_MARGIN)
        )
        if torch.isnan(loss) or torch.isinf(loss):
            logger.warning(f"  [NaN/Inf] loss={loss.item()}, skipping update")
            return float('nan')
        optimizer.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()

    return loss.item()


@torch.no_grad()
def evaluate(model, train_data, full_data_dt, full_data_td,
             dt_val_idx, td_val_idx, n_genes,
             dt_sampler, td_sampler,
             cached_neg_dt=None, cached_neg_td=None,
             dt_train_idx=None):
    """
    评估：使用训练子图做消息传递，验证边仅作为打分目标（不参与消息传递）

    支持缓存负边 (cached_neg_dt, cached_neg_td) 避免重复采样。
    若提供 dt_train_idx，则同时计算训练集 DT AUROC（复用同一次 forward pass）。

    train_data:  仅含训练边的子图，用于 forward 获取无泄漏嵌入
    full_data_dt: 原始完整 drug-target 边索引（验证正边来源）
    full_data_td: 原始完整 gene-disease 边索引（验证正边来源）
    """
    model.eval()
    t0 = time.time()
    z_dict = model(train_data.x_dict, train_data.edge_index_dict,
                   edge_weight_dict=_get_coexp_edge_weight(train_data))
    if device.type == 'cuda':
        torch.cuda.synchronize()
    logger.info(f"  [Eval] forward pass took {time.time()-t0:.1f}s")
    import sys; sys.stderr.flush()  # 强制刷新

    # ---- 药物-靶点 ----
    pos_dt = full_data_dt[:, dt_val_idx]
    if cached_neg_dt is not None:
        neg_dt = cached_neg_dt
    else:
        neg_dt = sample_negative_edges_via_set(
            pos_dt, n_genes, dt_sampler, is_dt=True,
        )
    dt_edge = torch.cat([pos_dt, neg_dt], dim=1)
    dt_label = torch.cat([
        torch.ones(pos_dt.size(1)), torch.zeros(neg_dt.size(1))
    ], dim=0).cpu().numpy()
    dt_score = torch.nan_to_num(
        torch.sigmoid(model.decode_drug_target(z_dict, dt_edge) / SIGMOID_TEMP), nan=0.5).cpu().numpy()
    dt_auroc = roc_auc_score(dt_label, dt_score) if len(np.unique(dt_label)) > 1 else 0.5
    dt_auprc = average_precision_score(dt_label, dt_score)
    logger.info(f"  [Eval] DT scoring done: AUC={dt_auroc:.4f}, AP={dt_auprc:.4f}"); import sys; sys.stderr.flush()

    # ---- 训练集 DT AUROC (复用 z_dict, 无需额外 forward pass) ----
    train_dt_auroc = None
    if dt_train_idx is not None:
        pos_tr = full_data_dt[:, dt_train_idx]
        n_tr = pos_tr.size(1)
        # 向量化负采样: 从基因池中排除正样本, 批量选择 (避免 while 循环重试)
        pos_tr_cpu = pos_tr.cpu()
        pos_dst_set = set(pos_tr_cpu[1].tolist())  # 只有 drug=0, 所以只收集 dst
        neg_candidates = np.setdiff1d(np.arange(n_genes), np.array(sorted(pos_dst_set)))
        neg_tr_dst = neg_candidates[np.random.choice(len(neg_candidates), size=n_tr, replace=True)]
        neg_tr = torch.stack([
            torch.zeros(n_tr, dtype=torch.long),
            torch.tensor(neg_tr_dst, dtype=torch.long)
        ], dim=0).to(device)
        tr_pos_score = torch.nan_to_num(
            torch.sigmoid(model.decode_drug_target(z_dict, pos_tr) / SIGMOID_TEMP), nan=0.5).cpu().numpy()
        tr_neg_score = torch.nan_to_num(
            torch.sigmoid(model.decode_drug_target(z_dict, neg_tr) / SIGMOID_TEMP), nan=0.5).cpu().numpy()
        tr_scores = np.concatenate([tr_pos_score, tr_neg_score])
        tr_labels = np.array([1]*n_tr + [0]*n_tr)
        train_dt_auroc = roc_auc_score(tr_labels, tr_scores) if len(np.unique(tr_labels)) > 1 else 0.5

    # ---- 靶点-疾病 ----
    pos_td = full_data_td[:, td_val_idx]
    if cached_neg_td is not None:
        neg_td = cached_neg_td
    else:
        neg_td = sample_negative_edges_via_set(
            pos_td, n_genes, td_sampler, is_dt=False,
        )
    td_edge = torch.cat([pos_td, neg_td], dim=1)
    td_label = torch.cat([
        torch.ones(pos_td.size(1)), torch.zeros(neg_td.size(1))
    ], dim=0).cpu().numpy()
    td_score = torch.nan_to_num(
        torch.sigmoid(model.decode_target_disease(z_dict, td_edge) / SIGMOID_TEMP), nan=0.5).cpu().numpy()
    td_auroc = roc_auc_score(td_label, td_score) if len(np.unique(td_label)) > 1 else 0.5
    td_auprc = average_precision_score(td_label, td_score)
    logger.info(f"  [Eval] TD scoring done: AUC={td_auroc:.4f}, AP={td_auprc:.4f}"); import sys; sys.stderr.flush()

    # ---- Precision@K / MRR / Hits@K (诊断排序质量) ----
    n_pos_td = pos_td.size(1)
    order = np.argsort(td_score)[::-1]
    td_p20 = td_label[order[:20]].sum() / 20.0

    # MRR
    pos_ranks = []
    for i in range(n_pos_td):
        rank = np.where(order == i)[0]
        if len(rank) > 0:
            pos_ranks.append(1.0 / (rank[0] + 1))
    td_mrr = np.mean(pos_ranks) if pos_ranks else 0.0

    # Hits@K
    hits_k = {}
    for k in [5, 10, 50, 100]:
        hits_k[k] = td_label[order[:k]].sum() / min(k, n_pos_td)

    return dt_auroc, dt_auprc, td_auroc, td_auprc, td_p20, td_mrr, hits_k, train_dt_auroc


# ============================================================================
# 5b. 训练集 AUROC 快速计算（过拟合监控）
# ============================================================================

@torch.no_grad()
def fast_train_auroc(model, train_data, full_dt_ei, dt_train_idx, n_genes, device):
    """
    快速计算 DT 训练集 AUROC：正边 = 训练 DT 边，负边 = 随机采样等量。

    用于监控 val-train gap，识别过拟合。
    不缓存负边，每次重新采样（轻量级，仅数百条）。
    负采样与 evaluate 一致：使用 np.setdiff1d 批量排除正样本。
    """
    model.eval()
    pos_ei = full_dt_ei[:, dt_train_idx]                     # [2, n_pos]
    n_pos = pos_ei.size(1)

    # 批量负采样: 从基因池中排除正样本, 向量化选择 (与 evaluate 一致)
    pos_ei_cpu = pos_ei.cpu()
    pos_dst_set = set(pos_ei_cpu[1].tolist())
    neg_candidates = np.setdiff1d(np.arange(n_genes), np.array(sorted(pos_dst_set)))
    neg_dst_np = neg_candidates[np.random.choice(len(neg_candidates), size=n_pos, replace=True)]

    neg_src = torch.zeros(n_pos, dtype=torch.long, device=device)
    neg_dst = torch.from_numpy(neg_dst_np).long().to(device)
    neg_ei = torch.stack([neg_src, neg_dst], dim=0)

    z_dict = model(train_data.x_dict, train_data.edge_index_dict,
                   edge_weight_dict=_get_coexp_edge_weight(train_data))
    print(f"  [fast_auroc] forward pass done, computing scores...", flush=True)
    pos_score = torch.sigmoid(model.decode_drug_target(z_dict, pos_ei) / SIGMOID_TEMP)
    neg_score = torch.sigmoid(model.decode_drug_target(z_dict, neg_ei) / SIGMOID_TEMP)
    scores = torch.cat([pos_score, neg_score]).cpu().numpy()
    labels = np.array([1]*n_pos + [0]*n_pos)
    auroc = roc_auc_score(labels, scores)
    return auroc, scores


# ============================================================================
# 6. 交叉验证
# ============================================================================

def cross_validate(data, n_genes, drug_target_genes, disease_genes_raw,
                  gene_to_idx, n_folds=N_FOLDS, epochs=EPOCHS):
    """
    5折交叉验证，关键改进：
      1. 每折仅用训练集疾病基因构建疾病特征（无数据泄漏）
      2. 使用 NegativeSampler 高效缓存
      3. 报告 AUROC / AUPRC / Precision@20 / MRR
    """
    dt_idx = torch.arange(data["drug", "targets", "gene"].edge_index.size(1))
    td_idx = torch.arange(data["gene", "associated_with", "disease"].edge_index.size(1))

    logger.info(f"[CV] DT edges: {len(dt_idx)}, TD edges: {len(td_idx)}")

    # ---- 构建 NegSampler 候选池 (PPI 1-hop 邻居约束) ----
    ppi_ei = data["gene", "interacts", "gene"].edge_index
    ppi_neighbors = build_ppi_neighbor_set(ppi_ei, n_genes)
    # 计算PPI图真实节点度数（用于度数感知负采样, Cappelletti 2024）
    ppi_node_degrees = torch.bincount(ppi_ei[0], minlength=n_genes).cpu().numpy().astype(np.float32)
    dt_target_indices = [gene_to_idx[g] for g in drug_target_genes if g in gene_to_idx]
    td_disease_indices = [gene_to_idx[g] for g in disease_genes_raw if g in gene_to_idx]

    dt_candidates = set()
    for idx in dt_target_indices: dt_candidates.update(ppi_neighbors.get(idx, set()))
    dt_candidates -= set(dt_target_indices)
    dt_candidates_arr = np.array(sorted(dt_candidates), dtype=np.int64)

    td_candidates = set()
    for idx in td_disease_indices: td_candidates.update(ppi_neighbors.get(idx, set()))
    td_candidates -= set(td_disease_indices)
    td_candidates_arr = np.array(sorted(td_candidates), dtype=np.int64)

    logger.info(f"[CV] DT hard-neg pool: {len(dt_candidates_arr)}, TD: {len(td_candidates_arr)}")

    if len(dt_idx) < n_folds or len(td_idx) < n_folds:
        k = min(len(dt_idx), len(td_idx), n_folds)
        if k < 2:
            logger.warning("[CV] Too few edges, skipping CV")
            return []
        n_folds = k

    kf_dt = KFold(n_splits=n_folds, shuffle=True, random_state=RANDOM_SEED)
    kf_td = KFold(n_splits=n_folds, shuffle=True, random_state=RANDOM_SEED)
    dt_splits = list(kf_dt.split(dt_idx))
    td_splits = list(kf_td.split(td_idx))

    # 基因索引 → 嵌入数组映射
    gene_feat_full = data["gene"].x.cpu().numpy()
    # 保存完整边索引（用于验证边位置引用）
    full_dt_ei = data["drug", "targets", "gene"].edge_index
    full_td_ei = data["gene", "associated_with", "disease"].edge_index
    full_ppi_ei = data["gene", "interacts", "gene"].edge_index

    cv_metrics = []
    for fold in range(n_folds):
        logger.info(f"\n{'='*50}")
        logger.info(f"[CV] Fold {fold + 1}/{n_folds}")
        logger.info(f"{'='*50}")

        dt_train, dt_val = dt_splits[fold]
        td_train, td_val = td_splits[fold]

        # ---- 每折动态疾病特征：仅用训练集疾病基因 ----
        td_train_edges = full_td_ei[:, td_train]
        train_disease_nodes = set(td_train_edges[0].tolist())
        if len(train_disease_nodes) > 0:
            disease_feat_fold = np.mean(
                gene_feat_full[list(train_disease_nodes)], axis=0, keepdims=True
            )
        else:
            disease_feat_fold = np.mean(gene_feat_full, axis=0, keepdims=True)

        # ---- 构建仅含训练边的子图（消除验证边消息泄漏） ----
        train_data = HeteroData()
        train_data["drug"].x = data["drug"].x.clone()
        train_data["gene"].x = data["gene"].x.clone()
        train_data["disease"].x = torch.from_numpy(disease_feat_fold).to(device)

        train_data["drug", "targets", "gene"].edge_index = full_dt_ei[:, dt_train].clone()
        # 反向边: gene→drug (与正向边同步, 药物接收靶基因聚合)
        rev_gene_drug = torch.stack([
            full_dt_ei[1, dt_train],   # gene index
            full_dt_ei[0, dt_train],   # drug index (0)
        ], dim=0)
        train_data["gene", "targeted_by", "drug"].edge_index = rev_gene_drug.clone()
        train_data["gene", "associated_with", "disease"].edge_index = full_td_ei[:, td_train].clone()
        train_data["gene", "interacts", "gene"].edge_index = full_ppi_ei.clone()

        # ---- 复制新节点类型和边 (TF/Pathway/Chemical, 不涉及标签, 无需裁剪) ----
        for nt in ["tf", "pathway", "chemical"]:
            if nt in data.node_types:
                train_data[nt].x = data[nt].x.clone()
        for et in data.edge_types:
            if et not in train_data.edge_types:
                train_data[et].edge_index = data[et].edge_index.clone()

        logger.info(f"[CV] Train subgraph: DT={dt_train.shape[0]}, TD={td_train.shape[0]}, "
              f"PPI={full_ppi_ei.size(1)} edges, pool DT={len(dt_candidates_arr)} TD={len(td_candidates_arr)}")

        dt_sampler = NegativeSampler(dt_candidates_arr, n_genes,
                                       ppi_node_degrees=ppi_node_degrees)
        td_sampler = NegativeSampler(td_candidates_arr, n_genes,
                                       ppi_node_degrees=ppi_node_degrees)

        # 预缓存验证负边（避免每次 evaluate 重复采样）
        pos_dt_val = full_dt_ei[:, dt_val]
        cached_neg_dt = sample_negative_edges_via_set(
            pos_dt_val, n_genes, dt_sampler, is_dt=True).to(device)
        pos_td_val = full_td_ei[:, td_val]
        cached_neg_td = sample_negative_edges_via_set(
            pos_td_val, n_genes, td_sampler, is_dt=False).to(device)
        logger.info(f"[CV] Cached neg dt={cached_neg_dt.shape[1]}, "
                    f"td={cached_neg_td.shape[1]} for fast evaluation")

        model = HeteroGAT(
            metadata=data.metadata(),
            drug_dim=data["drug"].x.size(1),
            gene_dim=data["gene"].x.size(1),
            disease_dim=data["disease"].x.size(1),
            hidden_dim=HIDDEN_DIM, out_dim=OUT_DIM,
            heads=GAT_HEADS, dropout=DROPOUT, n_layers=N_LAYERS,
        ).to(device)

        n_params = sum(p.numel() for p in model.parameters())
        logger.info(f"[CV] Model created: {n_params:,} params, device={device}")

        if USE_COMPILE:
            try:
                model = torch.compile(model, dynamic=True)
            except Exception:
                pass

        optimizer = torch.optim.AdamW(
            model.parameters(), lr=LR, weight_decay=WEIGHT_DECAY
        )
        scaler = torch.cuda.amp.GradScaler() if USE_AMP else None

        best_val = float("inf")
        patience_cnt = 0
        fold_best = {}
        swa_model = None
        swa_n = 0
        nan_streak = 0
        nan_resets_left = 3  # 最多重初始化 3 次

        # 快速冒烟测试: 首次 train + eval
        logger.info(f"[CV] Starting training loop ({epochs} epochs, eval every 50)...")

        for epoch in range(epochs):
            # 余弦退火 DropEdge 调度
            de_p = get_dropedge_p(epoch, epochs)

            t0 = time.time()
            loss = train_epoch(model, train_data, full_dt_ei, full_td_ei,
                               dt_train, td_train, n_genes,
                               optimizer, dt_sampler, td_sampler, scaler,
                               dropedge_p=de_p)
            if epoch == 0:
                logger.info(f"  [Debug] train_epoch took {time.time()-t0:.1f}s")

            # ---- NaN 恢复: 连续 NaN 则重初始化 ----
            if math.isnan(loss):
                nan_streak += 1
                if nan_streak >= 5 and nan_resets_left > 0:
                    nan_resets_left -= 1
                    logger.warning(f"[CV] NaN streak={nan_streak}, reinitializing model ({nan_resets_left} resets left)...")
                    del model
                    torch.cuda.empty_cache()
                    model = HeteroGAT(
                        metadata=data.metadata(),
                        drug_dim=data["drug"].x.size(1),
                        gene_dim=data["gene"].x.size(1),
                        disease_dim=data["disease"].x.size(1),
                        hidden_dim=HIDDEN_DIM, out_dim=OUT_DIM,
                        heads=GAT_HEADS, dropout=DROPOUT, n_layers=N_LAYERS,
                    ).to(device)
                    optimizer = torch.optim.AdamW(
                        model.parameters(), lr=LR, weight_decay=WEIGHT_DECAY
                    )
                    scaler = torch.cuda.amp.GradScaler() if USE_AMP else None
                    nan_streak = 0
                    best_val = float("inf")
                    patience_cnt = 0
                    swa_model = None
                    swa_n = 0
                    continue
                elif nan_streak >= 10:
                    logger.error(f"[CV] NaN persists after 3 resets, skipping fold {fold+1}")
                    fold_best = dict(fold=fold+1, epoch=0,
                                     dt_auroc=0.5, dt_auprc=0.5,
                                     td_auroc=0.5, td_auprc=0.5,
                                     td_p20=0.0, td_mrr=0.0,
                                     hits_k={5:0.0, 10:0.0, 50:0.0, 100:0.0})
                    break
                continue
            else:
                nan_streak = 0  # 正常 loss, 重置 NaN 计数器

            # SWA: 从 SWA_START 开始维护权重滑动平均
            if epoch >= SWA_START:
                swa_model, swa_n = update_swa_model(swa_model, model, swa_n)

            if (epoch + 1) % 50 == 0 or epoch == 0:
                dt_auroc, dt_auprc, td_auroc, td_auprc, td_p20, td_mrr, hits_k, train_dt_auroc = evaluate(
                    model, train_data, full_dt_ei, full_td_ei,
                    dt_val, td_val, n_genes,
                    dt_sampler, td_sampler,
                    cached_neg_dt=cached_neg_dt, cached_neg_td=cached_neg_td,
                    dt_train_idx=dt_train,
                )
                # 简化早停: 纯 AUROC 均值 (用户建议, 避免 Hits@K 在正样本多时饱和)
                val_metric = -(dt_auroc + td_auroc) / 2.0
                # --- 监控过拟合: 训练集 AUROC (复用 evaluate 的 z_dict, 无需额外 forward) ---
                overfit_gap = dt_auroc - train_dt_auroc if train_dt_auroc is not None else 0.0
                logger.info(f"  Epoch {epoch+1:3d}/{epochs} | Loss: {loss:.4f} | "
                      f"DT AUC: val={dt_auroc:.4f}/tr={train_dt_auroc:.4f} | "
                      f"TD AUC: {td_auroc:.4f}/{td_auprc:.4f} "
                      f"P@20: {td_p20:.3f} MRR: {td_mrr:.3f} "
                      f"gap={overfit_gap:.3f}")
                if overfit_gap < -0.1:
                    logger.warning(f"  ⚠ 过拟合: val-train DT gap={overfit_gap:.3f}")
                if val_metric < best_val:
                    best_val = val_metric
                    patience_cnt = 0
                    fold_best = dict(fold=fold+1, epoch=epoch+1,
                                     dt_auroc=dt_auroc, dt_auprc=dt_auprc,
                                     td_auroc=td_auroc, td_auprc=td_auprc,
                                     td_p20=td_p20, td_mrr=td_mrr,
                                     hits_k=hits_k)
                else:
                    patience_cnt += 1
                    if patience_cnt >= PATIENCE:
                        logger.info(f"  [Early stop] epoch {epoch+1}")
                        break

        cv_metrics.append(fold_best)
        logger.info(f"[CV] Fold {fold+1} best: DT AUC={fold_best['dt_auroc']:.4f}, "
              f"TD AUC={fold_best['td_auroc']:.4f}, "
              f"MRR={fold_best['td_mrr']:.3f}, "
              f"H@5={fold_best['hits_k'][5]:.2f} H@10={fold_best['hits_k'][10]:.2f} "
              f"H@50={fold_best['hits_k'][50]:.2f} H@100={fold_best['hits_k'][100]:.2f}")

        # ---- 过平滑诊断 ----
        model.eval()
        try:
            model.diagnose_oversmoothing(data.x_dict, data.edge_index_dict)
        except Exception as e:
            logger.warning(f"[Diagnose] Oversmoothing check skipped: {e}")

        if device.type == "cuda":
            torch.cuda.empty_cache()

    return cv_metrics


# ============================================================================
# 7. 全数据训练 + 桥梁靶点预测
# ============================================================================

def _single_final_train(model, train_data, full_dt_ei, full_td_ei,
                       dt_train, td_train, dt_val, td_val,
                       dt_sampler, td_sampler, n_genes, epochs):
    """单次最终训练（SWA + 早停），返回训练好的模型"""
    if USE_COMPILE:
        try:
            model = torch.compile(model, dynamic=True)
        except Exception:
            pass
    optimizer = torch.optim.AdamW(model.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)
    scaler = torch.cuda.amp.GradScaler() if USE_AMP else None

    # 预缓存验证负边
    pos_dt_val = full_dt_ei[:, dt_val]
    cached_neg_dt = sample_negative_edges_via_set(
        pos_dt_val, n_genes, dt_sampler, is_dt=True).to(device)
    pos_td_val = full_td_ei[:, td_val]
    cached_neg_td = sample_negative_edges_via_set(
        pos_td_val, n_genes, td_sampler, is_dt=False).to(device)

    best_val = float("inf")
    patience_cnt = 0
    swa_model = None
    swa_n = 0

    for epoch in range(epochs):
        # 余弦退火 DropEdge 调度
        de_p = get_dropedge_p(epoch, epochs)

        loss = train_epoch(model, train_data, full_dt_ei, full_td_ei,
                           dt_train, td_train, n_genes,
                           optimizer, dt_sampler, td_sampler, scaler,
                           dropedge_p=de_p)

        if epoch >= SWA_START:
            swa_model, swa_n = update_swa_model(swa_model, model, swa_n)

        if (epoch + 1) % 50 == 0 or epoch == 0:
            dt_auroc, dt_auprc, td_auroc, td_auprc, td_p20, td_mrr, hits_k, _ = evaluate(
                model, train_data, full_dt_ei, full_td_ei,
                dt_val, td_val, n_genes, dt_sampler, td_sampler,
                cached_neg_dt=cached_neg_dt, cached_neg_td=cached_neg_td)
            # 纯 AUROC 早停，与 cross_validate 一致
            val_metric = -(dt_auroc + td_auroc) / 2.0
            if val_metric < best_val:
                best_val = val_metric
                patience_cnt = 0
            else:
                patience_cnt += 1
                if patience_cnt >= PATIENCE:
                    break

    # 返回 SWA 平均后的模型
    return finalize_swa(swa_model, swa_n) or model


def train_final_and_predict(data, n_genes, gene_list, gene_to_idx,
                           drug_target_genes, disease_genes_raw):
    """
    3次集成训练 + 打分平均，确保排名稳定可复现

    每次训练使用 seed = RANDOM_SEED + run_idx，
    所有随机源（Python / NumPy / PyTorch / CUDA）确定性锁定。
    最终 combined_score = mean(combined_score_1, combined_score_2, combined_score_3)
    """
    logger.info(f"\n{'='*50}")
    logger.info(f"[Final] Ensemble training ({ENSEMBLE_RUNS} runs, deterministic seeds)...")
    logger.info(f"{'='*50}")

    # ---- 候选池 & 已知基因（与训练无关，预计算一次） ----
    ppi_ei = data["gene", "interacts", "gene"].edge_index
    ppi_neighbors = build_ppi_neighbor_set(ppi_ei, n_genes)
    # 计算PPI图真实节点度数（用于度数感知负采样, Cappelletti 2024）
    ppi_node_degrees = torch.bincount(ppi_ei[0], minlength=n_genes).cpu().numpy().astype(np.float32)
    dt_target_indices = [gene_to_idx[g] for g in drug_target_genes if g in gene_to_idx]
    td_disease_indices = [gene_to_idx[g] for g in disease_genes_raw if g in gene_to_idx]

    dt_candidates = set()
    for idx in dt_target_indices: dt_candidates.update(ppi_neighbors.get(idx, set()))
    dt_candidates -= set(dt_target_indices)
    dt_candidates_arr = np.array(sorted(dt_candidates), dtype=np.int64)

    td_candidates = set()
    for idx in td_disease_indices: td_candidates.update(ppi_neighbors.get(idx, set()))
    td_candidates -= set(td_disease_indices)
    td_candidates_arr = np.array(sorted(td_candidates), dtype=np.int64)

    known_dt_genes = set(drug_target_genes)
    known_td_genes = set(disease_genes_raw)

    # ---- 验证集划分 ----
    dt_all = torch.arange(data["drug", "targets", "gene"].edge_index.size(1))
    td_all = torch.arange(data["gene", "associated_with", "disease"].edge_index.size(1))
    full_dt_ei = data["drug", "targets", "gene"].edge_index
    full_td_ei = data["gene", "associated_with", "disease"].edge_index
    full_ppi_ei = data["gene", "interacts", "gene"].edge_index

    n_dt_val = max(1, int(len(dt_all) * VAL_RATIO))
    n_td_val = max(1, int(len(td_all) * VAL_RATIO))
    dt_perm = torch.randperm(len(dt_all), generator=torch.Generator().manual_seed(RANDOM_SEED))
    td_perm = torch.randperm(len(td_all), generator=torch.Generator().manual_seed(RANDOM_SEED))
    dt_train_idx = dt_all[dt_perm[n_dt_val:]]
    dt_val_idx = dt_all[dt_perm[:n_dt_val]]
    td_train_idx = td_all[td_perm[n_td_val:]]
    td_val_idx = td_all[td_perm[:n_td_val]]
    logger.info(f"[Final] DT: {len(dt_train_idx)} train / {len(dt_val_idx)} val, "
          f"TD: {len(td_train_idx)} train / {len(td_val_idx)} val")

    # ---- 构建仅含训练边的子图（最终验证也需隔离验证边） ----
    td_train_edges = full_td_ei[:, td_train_idx]
    train_disease_nodes = set(td_train_edges[0].tolist())
    gene_feat_full = data["gene"].x.cpu().numpy()
    if len(train_disease_nodes) > 0:
        disease_feat = np.mean(gene_feat_full[list(train_disease_nodes)], axis=0, keepdims=True)
    else:
        disease_feat = np.mean(gene_feat_full, axis=0, keepdims=True)

    final_train_data = HeteroData()
    final_train_data["drug"].x = data["drug"].x.clone()
    final_train_data["gene"].x = data["gene"].x.clone()
    final_train_data["disease"].x = torch.from_numpy(disease_feat).to(device)
    final_train_data["drug", "targets", "gene"].edge_index = full_dt_ei[:, dt_train_idx].clone()
    # 反向边: gene→drug (与正向边同步)
    rev_gene_drug = torch.stack([
        full_dt_ei[1, dt_train_idx],
        full_dt_ei[0, dt_train_idx],
    ], dim=0)
    final_train_data["gene", "targeted_by", "drug"].edge_index = rev_gene_drug.clone()
    final_train_data["gene", "associated_with", "disease"].edge_index = full_td_ei[:, td_train_idx].clone()
    final_train_data["gene", "interacts", "gene"].edge_index = full_ppi_ei.clone()

    # ---- 复制新节点类型和边 (TF/Pathway/Chemical, 不涉及标签, 无需裁剪) ----
    for nt in ["tf", "pathway", "chemical"]:
        if nt in data.node_types:
            final_train_data[nt].x = data[nt].x.clone()
    for et in data.edge_types:
        if et not in final_train_data.edge_types:
            final_train_data[et].edge_index = data[et].edge_index.clone()
            # 复制共表达边权重
            if et == ("gene", "coexpressed_with", "gene") and hasattr(data[et], "edge_weight"):
                final_train_data[et].edge_weight = data[et].edge_weight.clone()

    # ---- 集成训练 & 打分 ----
    all_run_scores = []
    final_trained_model = None  # 保存最终模型用于注意力分析

    for run_idx in range(ENSEMBLE_RUNS):
        # 每轮启动前清理显存 + 堆内存，避免跨轮累积
        if run_idx > 0:
            torch.cuda.empty_cache()
            gc.collect()

        run_seed = RANDOM_SEED + run_idx + 100
        set_seed(run_seed)
        logger.info(f"\n  --- Ensemble Run {run_idx+1}/{ENSEMBLE_RUNS} (seed={run_seed}) ---")

        dt_sampler = NegativeSampler(dt_candidates_arr, n_genes,
                                   ppi_node_degrees=ppi_node_degrees)
        td_sampler = NegativeSampler(td_candidates_arr, n_genes,
                                   ppi_node_degrees=ppi_node_degrees)

        model = HeteroGAT(
            metadata=data.metadata(),
            drug_dim=data["drug"].x.size(1),
            gene_dim=data["gene"].x.size(1),
            disease_dim=data["disease"].x.size(1),
            hidden_dim=HIDDEN_DIM, out_dim=OUT_DIM,
            heads=GAT_HEADS, dropout=DROPOUT, n_layers=N_LAYERS,
        ).to(device)

        # 捕获 SWA 平均后的模型（否则 SWA 权重被丢弃）
        model = _single_final_train(model, final_train_data, full_dt_ei, full_td_ei,
                                    dt_train_idx, td_train_idx, dt_val_idx, td_val_idx,
                                    dt_sampler, td_sampler, n_genes, EPOCHS)

        # 对新颖基因评分：使用全图 forward（不存在泄漏——评分目标不在已知边中）
        model.eval()
        # 保存最后一次集成的模型用于注意力分析
        if run_idx == ENSEMBLE_RUNS - 1:
            final_trained_model = copy.deepcopy(model)
        # NaN 检查：验证模型权重是否完好
        has_nan = False
        for name, p in model.named_parameters():
            if torch.isnan(p).any() or torch.isinf(p).any():
                logger.info(f"  [WARN] NaN/Inf in {name}, skipping this run")
                has_nan = True
                break
        if has_nan:
            logger.error("[ERROR] NaN detected in model weights during ensemble run. Aborting.")
            if device.type == "cuda":
                torch.cuda.empty_cache()
            return None

        with torch.no_grad():
            z_dict = model(data.x_dict, data.edge_index_dict,
                           edge_weight_dict=_get_coexp_edge_weight(data))

        run_results = []
        for gene_symbol in gene_list:
            if gene_symbol in known_dt_genes or gene_symbol in known_td_genes:
                continue
            gi = gene_to_idx[gene_symbol]
            dt_logit = model.decode_drug_target(z_dict,
                torch.tensor([[0], [gi]], dtype=torch.long, device=device))
            td_logit = model.decode_target_disease(z_dict,
                torch.tensor([[gi], [0]], dtype=torch.long, device=device))
            drug_score = torch.nan_to_num(torch.sigmoid(dt_logit / SIGMOID_TEMP), nan=0.5).item()
            disease_score = torch.nan_to_num(torch.sigmoid(td_logit / SIGMOID_TEMP), nan=0.5).item()
            run_results.append(dict(
                gene_symbol=gene_symbol,
                drug_target_score=drug_score,
                target_disease_score=disease_score,
                combined_score=drug_score * disease_score,
            ))
        all_run_scores.append(pd.DataFrame(run_results))
        if device.type == "cuda":
            torch.cuda.empty_cache()

    # ---- 打分平均 & RRF 排序 ----
    # RRF (Reciprocal Rank Fusion): 对排名而非原始分数取平均，对极端值鲁棒
    # refactored into a public function to eliminate code duplication
    merged = ensemble_rrf_fusion(all_run_scores, ENSEMBLE_RUNS, k_rrf=60)

    # 计算排名稳定性 (Top 20 跨运行的 Jaccard 相似度)
    mean_jaccard = compute_ensemble_jaccard(all_run_scores, ENSEMBLE_RUNS, k_rrf=60)
    logger.info(f"\n[Ensemble] Top-20 Jaccard (RRF): {mean_jaccard:.3f} "
          f"(1.0 = identical, >0.7 = stable)")

    # 保存
    merged["drug_target_score"] = merged["drug_target_score"].round(6)
    merged["target_disease_score"] = merged["target_disease_score"].round(6)
    merged["combined_score"] = merged["combined_score"].round(6)
    merged["rrf_score"] = merged["rrf_score"].round(6)
    merged["dt_mean_rank"] = merged["dt_mean_rank"].round(1)
    merged["td_mean_rank"] = merged["td_mean_rank"].round(1)

    all_output = os.path.join(SCRIPT_DIR, "all_bridge_genes.csv")
    # 列顺序: gene_symbol, drug_target_score, target_disease_score, combined_score, rrf_score, dt_mean_rank, td_mean_rank
    merged.to_csv(all_output, index=False)

    top20 = merged.head(20)
    try:
        top20.to_csv(OUTPUT, index=False)
    except PermissionError:
        logger.info(f"[Warn] Cannot write {OUTPUT} (file locked by IDE)")

    logger.info(f"[Output] Top 20 → {OUTPUT}")
    logger.info(f"[Output] All {len(merged)} genes → {all_output}")
    logger.info(f"[Info]  Skipped {len(known_dt_genes | known_td_genes)} known genes, "
          f"scored {len(merged)} novel candidates")
    logger.info(f"\n{'='*50}")
    logger.info(f"[Result] Top 20 Novel Bridge Targets (RRF, {ENSEMBLE_RUNS}-run ensemble)")
    logger.info(f"{'='*50}")
    for _, row in top20.iterrows():
        logger.info(f"  {row['gene_symbol']:15s} | DT={row['drug_target_score']:.4f} "
              f"| TD={row['target_disease_score']:.4f} "
              f"| RRF={row['rrf_score']:.4f}")
    logger.info(f"{'='*50}")
    return top20, final_trained_model


# ============================================================================
# 8a. 注意力权重可解释性分析
# ============================================================================

@torch.no_grad()
def analyze_attention_weights(model, data, gene_to_idx, gene_list_sorted,
                               top_n_bridge=20, top_k_neighbors=10):
    """
    提取并分析异构图注意力权重，解释桥梁基因预测的驱动因素。

    流程:
    1. 启用注意力存储 → 全图推理
    2. 收集每层/每种边类型的注意力系数 alpha
    3. 对 top-N 桥梁基因, 提取其 top-K 高注意力邻居
    4. 保存注意力分析结果到 CSV

    输出文件:
      attention_analysis/attention_per_gene.csv  — 每个基因的 top 邻居
      attention_analysis/attention_summary.csv    — 按边类型/层汇总
    """
    out_dir = os.path.join(SCRIPT_DIR, "attention_analysis")
    os.makedirs(out_dir, exist_ok=True)

    # 反转索引: idx → gene_symbol
    idx_to_gene = {v: k for k, v in gene_to_idx.items()}
    n_genes = len(gene_to_idx)

    # 1. 启用存储 + 前向传播
    model.eval()
    model.enable_attention_store(True)
    _ = model(data.x_dict, data.edge_index_dict,
             edge_weight_dict=_get_coexp_edge_weight(data))  # 前向, 触发注意力保存
    attn_weights = model.get_attention_weights()
    model.enable_attention_store(False)

    if not attn_weights:
        logger.warning("[Attention] No attention weights collected")
        return

    # 2. 汇总统计
    logger.info(f"\n{'='*60}")
    logger.info(f"  注意力权重可解释性分析")
    logger.info(f"{'='*60}")
    logger.info(f"  收集到 {len(attn_weights)} 种边类型的注意力:")
    for etype, layers in attn_weights.items():
        for lname, (ei, alp) in layers.items():
            e_str = f"({' → '.join(etype)})"
            n_edges = ei.shape[1]
            mean_attn = alp.mean().item()
            logger.info(f"    {lname:6s} {e_str:40s} {n_edges:7d} edges, "
                  f"mean α={mean_attn:.4f}")

    # 3. 读取桥梁基因排名
    bridge_path = os.path.join(SCRIPT_DIR, "all_bridge_genes.csv")
    if not os.path.exists(bridge_path):
        logger.warning(f"[Attention] 未找到 {bridge_path}, 跳过基因级分析")
        return
    bridge_df = pd.read_csv(bridge_path)
    top_genes = bridge_df.head(top_n_bridge)["gene_symbol"].tolist()
    logger.info(f"\n  Top-{top_n_bridge} 桥梁基因注意力分析:")

    # 4. 逐基因提取 top-K 注意力边
    per_gene_rows = []
    for gene_symbol in top_genes:
        if gene_symbol not in gene_to_idx:
            continue
        gi = gene_to_idx[gene_symbol]
        logger.info(f"  {gene_symbol} (idx={gi}):")

        for etype, layers in attn_weights.items():
            etype_str = " → ".join(etype)
            for lname, (ei, alp) in layers.items():
                # 找到与该基因相连的边
                if etype == ("gene", "interacts", "gene"):
                    mask = (ei[0] == gi) | (ei[1] == gi)
                elif etype == ("drug", "targets", "gene"):
                    mask = (ei[1] == gi)  # 基因是目标
                elif etype == ("gene", "associated_with", "disease"):
                    mask = (ei[0] == gi)  # 基因是源
                elif etype == ("gene", "targeted_by", "drug"):
                    mask = (ei[0] == gi)  # 基因→药物 (基因是源, 药物接收聚合)
                else:
                    continue

                n_conn = mask.sum().item()
                if n_conn == 0:
                    continue

                conn_alp = alp[mask]  # [n_conn, heads]
                conn_ei = ei[:, mask]  # [2, n_conn]
                attn_mean = conn_alp.mean(dim=1).numpy()  # [n_conn]
                top_k = min(top_k_neighbors, n_conn)
                top_idx = np.argsort(attn_mean)[::-1][:top_k]

                for rank, idx_pos in enumerate(top_idx, 1):
                    src, dst = conn_ei[0, idx_pos].item(), conn_ei[1, idx_pos].item()
                    if etype == ("gene", "interacts", "gene"):
                        neighbor = idx_to_gene.get(dst if src == gi else src, f"NODE_{dst}")
                    elif etype == ("drug", "targets", "gene"):
                        neighbor = "BCP(drug)" if src == 0 else idx_to_gene.get(src, f"NODE_{src}")
                    elif etype == ("gene", "associated_with", "disease"):
                        neighbor = "CIRI(disease)" if dst == 0 else idx_to_gene.get(dst, f"NODE_{dst}")
                    elif etype == ("gene", "targeted_by", "drug"):
                        neighbor = "BCP(drug→gene)"
                    else:
                        neighbor = f"NODE_{dst}"

                    per_gene_rows.append(dict(
                        gene_symbol=gene_symbol,
                        edge_type=etype_str,
                        layer=lname,
                        rank=rank,
                        neighbor=neighbor,
                        attention_score=round(float(attn_mean[idx_pos]), 6),
                        src_idx=src,
                        dst_idx=dst,
                    ))

                top_score = attn_mean[top_idx[0]] if len(top_idx) > 0 else 0
                logger.info(f"    {lname:6s} {etype_str:40s} "
                      f"{n_conn} edges, top α={top_score:.4f}")

    # 5. 保存
    attn_df = pd.DataFrame(per_gene_rows)
    attn_csv = os.path.join(out_dir, "attention_per_gene.csv")
    attn_df.to_csv(attn_csv, index=False)
    logger.info(f"\n  [Output] 基因级注意力 → {attn_csv}")

    # 6. 边类型汇总
    summary_rows = []
    for etype, layers in attn_weights.items():
        etype_str = " → ".join(etype)
        for lname, (ei, alp) in layers.items():
            summary_rows.append(dict(
                edge_type=etype_str,
                layer=lname,
                n_edges=ei.shape[1],
                n_heads=alp.shape[1],
                mean_alpha=round(float(alp.mean().item()), 6),
                std_alpha=round(float(alp.std().item()), 6),
                max_alpha=round(float(alp.max().item()), 6),
                min_alpha=round(float(alp.min().item()), 6),
            ))
    summary_df = pd.DataFrame(summary_rows)
    summary_csv = os.path.join(out_dir, "attention_summary.csv")
    summary_df.to_csv(summary_csv, index=False)
    logger.info(f"  [Output] 注意力汇总 → {summary_csv}")
    logger.info(f"{'='*60}\n")


# ============================================================================
# 8. 主流程
# ============================================================================

def _load_base_data():
    """加载所有配置共享的基础数据（缓存一次，避免重复IO）"""
    drug_target_genes = load_drug_targets(
        os.path.join(GAT_DATA_DIR, "drug_targets.txt"))
    disease_genes_raw = load_disease_genes(
        os.path.join(GAT_DATA_DIR, "disease_genes.txt"))
    ppi_edges = load_ppi(
        os.path.join(GAT_DATA_DIR, "ppi_subgraph.csv"))
    gene_features_arr, gene_feature_names = load_gene_features(
        os.path.join(GAT_DATA_DIR, "subgraph_embeddings.csv"))
    drug_fingerprint_arr = load_drug_fingerprint(
        os.path.join(GAT_DATA_DIR, "drug_fingerprint.csv"))
    all_genes = load_all_genes(
        os.path.join(GAT_DATA_DIR, "subgraph_genes.txt"))
    toxirna_df = load_toxirna_features(
        os.path.join(SCRIPT_DIR, "toxirna_enhanced_features.csv"))
    return dict(
        drug_target_genes=drug_target_genes,
        disease_genes_raw=disease_genes_raw,
        ppi_edges=ppi_edges,
        gene_features_arr=gene_features_arr,
        gene_feature_names=gene_feature_names,
        drug_fingerprint_arr=drug_fingerprint_arr,
        all_genes=all_genes,
        toxirna_df=toxirna_df,
    )


def run_single_experiment(exp_name, cfg, base_data, skip_attention=False):
    """
    运行单次消融实验。

    参数:
        exp_name: 实验名称 (用于日志 & 结果记录)
        cfg: dict with keys tf, meth, pw, ctd, coexp
        base_data: _load_base_data() 返回的共享数据
        skip_attention: 消融模式跳过注意力分析 (加速)

    返回:
        dict: {name, dt_auroc, dt_auprc, td_auroc, td_auprc, td_p20, td_mrr, top1_bridge}
        or None on failure
    """
    use_tf   = cfg["tf"]
    use_meth = cfg["meth"]
    use_pw   = cfg["pw"]
    use_ctd  = cfg["ctd"]
    use_coexp = cfg["coexp"]

    logger.info(f"\n{'#'*70}")
    logger.info(f"#  Experiment: {exp_name}")
    logger.info(f"#  Config: TF={use_tf} Meth={use_meth} PW={use_pw} CTD={use_ctd} Coexp={use_coexp}")
    logger.info(f"{'#'*70}")

    # 重置随机种子 — 确保每个实验可从相同基线开始
    set_seed(RANDOM_SEED)

    # ---- 提取共享数据 ----
    drug_target_genes = base_data["drug_target_genes"]
    disease_genes_raw  = base_data["disease_genes_raw"]
    ppi_edges          = base_data["ppi_edges"]
    gene_features_arr  = base_data["gene_features_arr"]
    gene_feature_names = base_data["gene_feature_names"]
    drug_fingerprint_arr = base_data["drug_fingerprint_arr"]
    all_genes          = base_data["all_genes"]
    toxirna_df         = base_data["toxirna_df"]

    # ---- 加载多源异构图数据 ----
    tf_edges, tf_nodes = None, None
    if use_tf:
        tf_edges, tf_nodes = load_tf_edges(
            os.path.join(GAT_EXT_DIR, "tf_target_edges.txt"))

    meth_path = None
    if use_meth:
        meth_path = os.path.join(GAT_EXT_DIR, "gene_methylation_edges.txt")

    pw_edges, pw_nodes = None, None
    pw_features, pw_feature_names = None, None
    if use_pw:
        pw_edges, pw_nodes = load_pathway_edges(
            os.path.join(GAT_EXT_DIR, "gene_pathway_edges.txt"))
        pw_features, pw_feature_names = load_pathway_features(
            os.path.join(GAT_EXT_DIR, "cache", "pathway_features_pca256.npy"))

    ctd_edges, ctd_chem_nodes = None, None
    ctd_gene_feat = None
    if use_ctd:
        ctd_gene_set = set(gene_feature_names) | set(drug_target_genes) | set(disease_genes_raw)
        ctd_edges, ctd_chem_nodes, ctd_gene_feat = load_ctd_edges(
            os.path.join(GAT_EXT_DIR, "CTD_chem_gene_ixns.tsv"),
            ctd_gene_set, top_k=CTD_TOP_CHEMICALS)

        # Step 1 模式: 将CTD直接互作基因合并到drug_targets (不引入chemical节点, 不加CTD特征)
        if CTD_MERGE_TARGETS and ctd_edges:
            ctd_target_genes = set(g for _, g in ctd_edges)
            n_original = len(drug_target_genes)
            drug_target_genes = list(dict.fromkeys(drug_target_genes + list(ctd_target_genes)))
            n_new = len(drug_target_genes) - n_original
            logger.info(f"[CTD Merge] Added {n_new} CTD direct target genes to drug_targets "
                        f"({n_original} → {len(drug_target_genes)})")
            # 清空 chemical 节点/边/特征, 避免 build_hetero_data 创建 chemical 节点或拼接特征
            ctd_edges = None
            ctd_chem_nodes = None
            ctd_gene_feat = None

    # Step 2 模式: 不合并靶点，保持chemical节点独立 (不引入gene-level CTD特征)
    if use_ctd and not CTD_GENE_FEATURES:
        ctd_gene_feat = None

    coexp_edges = None
    if use_coexp:
        coexp_edges = load_coexp_edges(
            os.path.join(GAT_EXT_DIR, "gene_coexp_edges.txt"))

    # ---- [数据层验证] ----
    ctd_has_features = use_ctd and CTD_GENE_FEATURES
    ctd_has_chem = use_ctd and not CTD_MERGE_TARGETS
    data_dim_summary = {
        "gene_base": gene_features_arr.shape[1],
        "drug_ecfp4": drug_fingerprint_arr.shape[1],
        "toxirna": toxirna_df.shape[1] if toxirna_df is not None else 0,
        "methylation": 20 if use_meth else 0,
        "pathway_pca": PATHWAY_PCA_DIM if use_pw else 0,
        "ctd_gene": CTD_GENE_FEAT_DIM if ctd_has_features else 0,
        "gene_final": (gene_features_arr.shape[1]
                       + (toxirna_df.shape[1] if toxirna_df is not None else 0)
                       + (20 if use_meth else 0)
                       + (PATHWAY_PCA_DIM if use_pw else 0)
                       + (CTD_GENE_FEAT_DIM if ctd_has_features else 0)),
        "tf_dim": 4,
        "chem_dim": CTD_CHEM_FEAT_DIM if ctd_has_chem else 0,
    }
    logger.info(f"\n{'='*60}")
    logger.info(f"  [{exp_name}] 数据层维度验证")
    logger.info(f"{'='*60}")
    for k, v in data_dim_summary.items():
        logger.info(f"  {k:>12s}: {v:5d} dims")
    logger.info(f"  {'genes':>12s}: {len(all_genes):5d}")
    logger.info(f"  {'drug_targets':>12s}: {len(drug_target_genes):5d}")
    logger.info(f"  {'disease_genes':>12s}: {len(disease_genes_raw):5d}")
    if tf_edges:   logger.info(f"  {'tf_edges':>12s}: {len(tf_edges):5d}")
    if pw_edges:   logger.info(f"  {'pw_edges':>12s}: {len(pw_edges):5d}")
    if ctd_edges:  logger.info(f"  {'ctd_edges':>12s}: {len(ctd_edges):5d}")
    if coexp_edges: logger.info(f"  {'coexp_edges':>12s}: {len(coexp_edges):5d}")
    logger.info(f"{'='*60}\n")

    # ---- 构建异构图 ----
    data, gene_to_idx = build_hetero_data(
        drug_target_genes, disease_genes_raw, ppi_edges,
        gene_features_arr, gene_feature_names, drug_fingerprint_arr,
        all_genes, toxirna_df,
        tf_edges=tf_edges, tf_nodes=tf_nodes,
        meth_path=meth_path, use_methylation=use_meth,
        pw_edges=pw_edges, pw_nodes=pw_nodes,
        pw_features=pw_features, pw_feature_names=pw_feature_names,
        ctd_edges=ctd_edges, ctd_chem_nodes=ctd_chem_nodes,
        ctd_gene_feat=ctd_gene_feat,
        coexp_edges=coexp_edges, use_coexp=use_coexp,
    )

    n_genes = data["gene"].x.size(0)
    gene_list_sorted = list(gene_to_idx.keys())
    data = data.to(device)

    # ---- 5折交叉验证 ----
    cv_metrics = cross_validate(data, n_genes, drug_target_genes,
                                disease_genes_raw, gene_to_idx, n_folds=N_FOLDS)

    # ---- 全数据训练 + 预测 ----
    train_result = train_final_and_predict(data, n_genes, gene_list_sorted, gene_to_idx,
                                          drug_target_genes, disease_genes_raw)
    if train_result is None or train_result[0] is None:
        logger.critical(f"[{exp_name}] Ensemble training failed. Skipping.")
        return None

    top20, final_model = train_result

    # ---- 注意力分析 (仅非消融模式) ----
    if not skip_attention and final_model is not None:
        logger.info(f"\n{'='*60}")
        logger.info(f"  [{exp_name}] 注意力可解释性分析...")
        logger.info(f"{'='*60}")
        data = data.to(device)
        try:
            analyze_attention_weights(
                final_model, data, gene_to_idx, gene_list_sorted,
                top_n_bridge=20, top_k_neighbors=10,
            )
        except Exception as e:
            logger.warning(f"[{exp_name}] 注意力分析异常: {e}")

    # ---- 汇总 ----
    result = {"name": exp_name}
    if cv_metrics:
        result["dt_auroc"] = round(np.mean([m["dt_auroc"] for m in cv_metrics]), 4)
        result["dt_auroc_std"] = round(np.std([m["dt_auroc"] for m in cv_metrics]), 4)
        result["dt_auprc"] = round(np.mean([m["dt_auprc"] for m in cv_metrics]), 4)
        result["dt_auprc_std"] = round(np.std([m["dt_auprc"] for m in cv_metrics]), 4)
        result["td_auroc"] = round(np.mean([m["td_auroc"] for m in cv_metrics]), 4)
        result["td_auroc_std"] = round(np.std([m["td_auroc"] for m in cv_metrics]), 4)
        result["td_auprc"] = round(np.mean([m["td_auprc"] for m in cv_metrics]), 4)
        result["td_auprc_std"] = round(np.std([m["td_auprc"] for m in cv_metrics]), 4)
        result["td_p20"] = round(np.mean([m.get("td_p20", 0) for m in cv_metrics]), 4)
        result["td_mrr"] = round(np.mean([m.get("td_mrr", 0) for m in cv_metrics]), 4)
        result["top1_bridge"] = top20.iloc[0]["gene_symbol"]
        result["top1_score"] = round(top20.iloc[0]["combined_score"], 4)

        logger.info(f"\n{'='*60}")
        logger.info(f"  [{exp_name}] Complete")
        logger.info(f"  DT AUROC: {result['dt_auroc']:.4f} ± {result['dt_auroc_std']:.4f}")
        logger.info(f"  DT AUPRC: {result['dt_auprc']:.4f} ± {result['dt_auprc_std']:.4f}")
        logger.info(f"  TD AUROC: {result['td_auroc']:.4f} ± {result['td_auroc_std']:.4f}")
        logger.info(f"  TD AUPRC: {result['td_auprc']:.4f} ± {result['td_auprc_std']:.4f}")
        logger.info(f"  TD P@20:  {result['td_p20']:.4f}")
        logger.info(f"  TD MRR:   {result['td_mrr']:.4f}")
        logger.info(f"  Top-1 Bridge: {result['top1_bridge']} (score={result['top1_score']:.4f})")
        logger.info(f"{'='*60}")
    else:
        result["dt_auroc"] = result["td_auroc"] = 0
        result["dt_auprc"] = result["td_auprc"] = 0
        result["top1_bridge"] = "N/A"
        logger.info(f"\n  [{exp_name}] Complete (CV skipped)")

    # 释放 GPU 显存
    del data, final_model
    torch.cuda.empty_cache()
    return result


def main():
    global USE_TF_EDGES, USE_METHYLATION_FEATURES, USE_PATHWAY_EDGES, USE_CTD_EDGES, USE_COEXP_EDGES

    logger.info("=" * 60)
    logger.info("  HeteroGAT v3: BCP x CIRI Bridge Target Prediction")
    logger.info("  MRHormer TRLA + MLGANN Pooling + TaRGET II Features")
    logger.info("=" * 60)

    # ---- 加载共享基础数据 ----
    logger.info("[Init] Loading base data (shared across all experiments)...")
    base_data = _load_base_data()

    if ABLATION_MODE:
        # ====================================================================
        # 消融实验模式: 遍历所有配置, 逐步增加 + 严格消融
        # ====================================================================
        logger.info(f"\n{'#'*70}")
        logger.info(f"#  ABLATION STUDY MODE")
        logger.info(f"#  {len(ABLATION_CONFIGS)} experiments: "
                     f"{sum(1 for c in ABLATION_CONFIGS if not c['name'].startswith('Full -'))} "
                     f"additions + "
                     f"{sum(1 for c in ABLATION_CONFIGS if c['name'].startswith('Full -'))} "
                     f"ablations")
        logger.info(f"#  Each: {N_FOLDS}-fold CV + {ENSEMBLE_RUNS} ensemble runs")
        logger.info(f"{'#'*70}\n")

        all_results = []
        for i, cfg in enumerate(ABLATION_CONFIGS):
            logger.info(f"\n{'▶'*30}")
            logger.info(f"▶  [{i+1}/{len(ABLATION_CONFIGS)}] Running: {cfg['name']}")
            logger.info(f"{'▶'*30}")

            # 同步模块级标志 (build_hetero_data 内部不再依赖它们, 但日志/其他函数可能引用)
            USE_TF_EDGES             = cfg["tf"]
            USE_METHYLATION_FEATURES = cfg["meth"]
            USE_PATHWAY_EDGES        = cfg["pw"]
            USE_CTD_EDGES            = cfg["ctd"]
            USE_COEXP_EDGES          = cfg["coexp"]

            result = run_single_experiment(cfg["name"], cfg, base_data, skip_attention=True)
            if result is not None:
                all_results.append(result)

            # 每次实验后清理
            gc.collect()
            torch.cuda.empty_cache()

        # ---- 生成消融实验结果对比表 ----
        if all_results:
            df = pd.DataFrame(all_results)
            # 添加 Δ 列: 相对基线的变化
            baseline = next((r for r in all_results if "Baseline" in r["name"]), None)
            if baseline:
                for key in ["dt_auroc", "dt_auprc", "td_auroc", "td_auprc"]:
                    df[f"Δ_{key}"] = df[key] - baseline[key]

            logger.info(f"\n{'#'*70}")
            logger.info(f"#  ABLATION STUDY RESULTS")
            logger.info(f"{'#'*70}")
            logger.info(f"\n{df.to_string(index=False)}")

            # 保存结果
            result_csv = os.path.join(SCRIPT_DIR, "ablation_results.csv")
            df.to_csv(result_csv, index=False)
            logger.info(f"\n[Output] Ablation results saved to: {result_csv}")

            # Markdown 对比表
            md_lines = []
            md_lines.append("# Ablation Study Results\n")
            md_lines.append("| Experiment | DT AUROC | DT AUPRC | TD AUROC | TD AUPRC | TD P@20 | TD MRR | Top-1 Bridge | Δ DT AUC |")
            md_lines.append("|:-----------|:--------:|:--------:|:--------:|:--------:|:-------:|:------:|:------------:|:--------:|")
            for _, row in df.iterrows():
                delta = f"{row.get('Δ_dt_auroc', 0):+.4f}" if 'Δ_dt_auroc' in row else "-"
                md_lines.append(
                    f"| {row['name']} | {row['dt_auroc']:.4f} | {row['dt_auprc']:.4f} | "
                    f"{row['td_auroc']:.4f} | {row['td_auprc']:.4f} | "
                    f"{row['td_p20']:.4f} | {row['td_mrr']:.4f} | {row['top1_bridge']} | {delta} |")
            md_lines.append("")
            md_lines.append(f"*Baseline: {baseline['name']}*")
            md_path = os.path.join(SCRIPT_DIR, "ablation_results.md")
            with open(md_path, "w", encoding="utf-8") as f:
                f.write("\n".join(md_lines))
            logger.info(f"[Output] Markdown table saved to: {md_path}")
        else:
            logger.critical("[ABLATION] No experiments completed successfully.")

    else:
        # ---- 单次运行模式 (向后兼容) ----
        cfg = {"tf": USE_TF_EDGES, "meth": USE_METHYLATION_FEATURES,
               "pw": USE_PATHWAY_EDGES, "ctd": USE_CTD_EDGES,
               "coexp": USE_COEXP_EDGES}
        run_single_experiment("Single Run", cfg, base_data, skip_attention=False)


if __name__ == "__main__":
    main()
