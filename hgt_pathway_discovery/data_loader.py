# -*- coding: utf-8 -*-
"""数据加载模块：统一加载所有输入数据（基因特征、边列表、通路、药物等）。"""

from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set

import numpy as np
import pandas as pd
import torch


# ============================================================================
# 文本/边列表加载
# ============================================================================

def load_txt(path: Path, skip_header_if_colname: bool = True) -> List[str]:
    """加载文本行列表，自动去重和跳过表头。"""
    if not path or not path.exists():
        return []
    for enc in ("utf-8", "gbk", "latin1"):
        try:
            with open(path, "r", encoding=enc) as f:
                lines = [ln.strip() for ln in f if ln.strip()]
            break
        except (UnicodeDecodeError, OSError):
            continue
    else:
        return []
    if skip_header_if_colname and lines and (
        lines[0].lower().startswith("gene") or lines[0].lower().startswith("symbol")
    ):
        lines = lines[1:]
    return list(dict.fromkeys(g.upper() for g in lines))


def load_edge_list(path: Optional[Path], sep: Optional[str] = None,
                   col_pair: Tuple[int, int] = (0, 1)) -> List[Tuple[str, str]]:
    """加载边列表（支持 .txt 和 .csv）。"""
    if not path or not path.exists():
        return []
    if sep is None:
        sep = "\t" if path.suffix == ".txt" else ","
    raw: List[Tuple[str, str]] = []
    for enc in ("utf-8", "gbk", "latin1"):
        try:
            with open(path, "r", encoding=enc) as f:
                for ln in f:
                    ln = ln.strip()
                    if not ln or ln.startswith("#"):
                        continue
                    parts = ln.split(sep)
                    if len(parts) > max(col_pair):
                        raw.append((parts[col_pair[0]].upper().strip(),
                                    parts[col_pair[1]].strip()))
            break
        except (UnicodeDecodeError, OSError):
            continue
    return raw


# ============================================================================
# 特征加载
# ============================================================================

def load_gene_features(path: Path) -> Tuple[np.ndarray, List[str]]:
    """加载增强基因特征矩阵 (CSV: 行=基因, 列=特征)。"""
    df = pd.read_csv(path, index_col=0)
    df.index = df.index.astype(str).str.strip().str.upper()
    arr = df.apply(pd.to_numeric, errors="coerce").fillna(0.0).values.astype(np.float32)
    gene_names = df.index.tolist()
    return arr, gene_names


def load_drug_fingerprint(path: Path) -> np.ndarray:
    """加载药物指纹特征。"""
    df = pd.read_csv(path)
    arr = df.apply(pd.to_numeric, errors="coerce").fillna(0).values.astype(np.float32)
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    return arr


def load_disease_features(path: Path) -> Optional[np.ndarray]:
    """加载疾病特征 (npy)。"""
    if not path.exists():
        return None
    arr = np.load(str(path))
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    elif arr.ndim == 2 and arr.shape[0] > 1:
        arr = np.mean(arr, axis=0, keepdims=True)
    return arr


def load_pathway_features(path: Path) -> np.ndarray:
    """加载通路特征 (npy)。"""
    arr = np.load(str(path))
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    return arr


# ============================================================================
# PPI 加载 (含置信度分数)
# ============================================================================

def load_ppi(path: Path, score_thresh: int = 700,
             bridge_set: Optional[Set[str]] = None,
             max_edges: int = 80000,
             subsample: bool = True) -> List[Tuple[str, str, float]]:
    """加载 PPI 边，返回 (gene_a, gene_b, combined_score)。"""
    df = pd.read_csv(path, sep="\t", encoding="utf-8")
    col_a, col_b = df.columns[0], df.columns[1]
    score_col = df.columns[2] if len(df.columns) >= 3 else None

    df = df.dropna(subset=[col_a, col_b])
    df[col_a] = df[col_a].astype(str).str.strip().str.upper()
    df[col_b] = df[col_b].astype(str).str.strip().str.upper()

    if score_col is not None:
        df[score_col] = pd.to_numeric(df[score_col], errors="coerce")
        df = df[df[score_col] >= score_thresh]
        df = df.sort_values(score_col, ascending=False)

    if subsample and bridge_set and len(bridge_set) > 0:
        mask = df[col_a].isin(bridge_set) | df[col_b].isin(bridge_set)
        df_bridge = df[mask]
        df_other = df[~mask]
        if len(df_bridge) > max_edges:
            df_bridge = df_bridge.head(max_edges)
        n_left = max(0, max_edges - len(df_bridge))
        df_other = df_other.head(n_left) if n_left > 0 else df_other.iloc[:0]
        df = pd.concat([df_bridge, df_other], ignore_index=True)
    elif subsample and len(df) > max_edges:
        df = df.head(max_edges)

    if score_col is not None:
        edges = list(zip(df[col_a], df[col_b], df[score_col]))
    else:
        edges = [(a, b, 0.0) for a, b in zip(df[col_a], df[col_b])]
    return edges


# ============================================================================
# 通路相关
# ============================================================================

def load_pathway_list(path: Path) -> List[str]:
    """加载通路名称列表。"""
    for enc in ("utf-8", "gbk", "latin1"):
        try:
            df = pd.read_csv(path, encoding=enc)
            break
        except UnicodeDecodeError:
            continue
    name_col = "pathway_name" if "pathway_name" in df.columns else df.columns[0]
    return df[name_col].tolist()


def load_gene_pathway_edges(path: Path) -> List[Tuple[str, str]]:
    """加载基因-通路关联边（保留通路名原始大小写）。"""
    edges: List[Tuple[str, str]] = []
    for enc in ("utf-8", "gbk", "latin1"):
        try:
            with open(path, "r", encoding=enc) as f:
                for ln in f:
                    ln = ln.strip()
                    if not ln:
                        continue
                    parts = ln.split("\t")
                    if len(parts) >= 2:
                        edges.append((parts[0].upper(), parts[1].strip()))
            break
        except (UnicodeDecodeError, OSError):
            continue
    return edges


def parse_reactome_hierarchy(gmt_path: Path) -> List[Tuple[str, str]]:
    """从 Reactome GMT 文件提取通路层级关系（父-子通路）。"""
    if not gmt_path.exists():
        return []

    import re
    hierarchy: List[Tuple[str, str]] = []
    pathway_data: Dict[str, Tuple[str, str, Set[str]]] = {}

    with open(gmt_path, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            name = parts[0].strip()
            reactome_id = parts[1].strip() if len(parts) > 1 else ""
            genes = set(g.strip().upper() for g in parts[2:] if g.strip())
            pathway_data[name] = (reactome_id, name, genes)

    pathway_names = sorted(pathway_data.keys())

    # 基于基因重叠推断层级：若通路 A 的基因大部分是通路 B 的子集，A 可能是 B 的子通路
    for i, name_a in enumerate(pathway_names):
        genes_a = pathway_data[name_a][2]
        if len(genes_a) < 3:
            continue
        for name_b in pathway_names[i + 1:]:
            genes_b = pathway_data[name_b][2]
            if len(genes_b) < 3:
                continue
            overlap = len(genes_a & genes_b)
            if overlap == 0:
                continue
            ratio_a = overlap / len(genes_a)
            ratio_b = overlap / len(genes_b)
            # A 几乎完全被 B 包含 → A 是 B 的子集
            if ratio_a > 0.9 and ratio_b < 0.9:
                hierarchy.append((name_b, name_a))  # parent → child
            elif ratio_b > 0.9 and ratio_a < 0.9:
                hierarchy.append((name_a, name_b))  # parent → child

    # 去重
    hierarchy = list(set(hierarchy))
    return hierarchy


# ============================================================================
# Bridge 基因
# ============================================================================

def load_bridge_genes(path: Path, alt_path: Optional[Path] = None) -> List[str]:
    """加载桥梁基因列表。"""
    if path.exists():
        df = pd.read_csv(path)
        col = "gene_symbol" if "gene_symbol" in df.columns else df.columns[0]
        genes = df[col].astype(str).str.strip().str.upper().tolist()
        return genes
    if alt_path and alt_path.exists():
        df = pd.read_csv(alt_path)
        genes = df["gene_symbol"].astype(str).str.strip().str.upper().tolist()
        return genes
    return []


# ============================================================================
# miRNA 边
# ============================================================================

def load_mirna_edges(path: Optional[Path]) -> List[Tuple[str, str]]:
    """加载 miRNA-靶基因边（若文件存在）。"""
    if not path or not path.exists():
        return []
    return load_edge_list(path, sep="\t", col_pair=(0, 1))


# ============================================================================
# 疾病-通路边
# ============================================================================

def load_disease_pathway_edges(path: Optional[Path]) -> List[Tuple[str, str]]:
    """加载疾病-通路关联边（DisGeNET/KEGG格式，若文件存在）。"""
    if not path or not path.exists():
        return []
    return load_edge_list(path, sep="\t", col_pair=(0, 1))


# ============================================================================
# 聚合加载
# ============================================================================

def load_all_data(data_dir: Path, gat_data_dir: Path, cache_dir: Path,
                  bridge_genes_path: Path, config: object) -> Dict:
    """一次性加载所有输入数据，返回字典。"""
    result = {}

    # 基因特征
    result["gene_feat_arr"], result["gene_feat_names"] = load_gene_features(
        cache_dir / "enhanced_gene_features.csv")

    # 通路特征
    pw_path = data_dir / "pathway_features.npy"
    result["pathway_feat_arr"] = load_pathway_features(pw_path)

    # 药物指纹
    result["drug_fp_arr"] = load_drug_fingerprint(gat_data_dir / "drug_fingerprint.csv")

    # 疾病特征
    result["disease_feat_arr"] = load_disease_features(data_dir / "disease_features.npy")

    # 通路列表
    result["pathway_names"] = load_pathway_list(cache_dir / "pathway_nodes.csv")

    # Bridge 基因
    result["bridge_genes"] = load_bridge_genes(
        bridge_genes_path, alt_path=gat_data_dir / "all_bridge_genes.csv")

    # 边数据
    bridge_set = set(result["bridge_genes"]) if result["bridge_genes"] else None
    result["ppi_edges"] = load_ppi(
        gat_data_dir / "ppi_subgraph.csv",
        score_thresh=config.graph.ppi_score_threshold,
        bridge_set=bridge_set,
        max_edges=config.graph.ppi_max_edges,
        subsample=config.graph.subsample_ppi,
    )

    result["coexp_edges"] = load_edge_list(data_dir / "gene_coexp_edges.txt", sep="\t")
    result["tf_edges"] = load_edge_list(data_dir / "tf_target_edges.txt", sep="\t")
    result["gene_pathway_edges"] = load_gene_pathway_edges(data_dir / "gene_pathway_edges.txt")
    result["all_genes_list"] = load_txt(gat_data_dir / "subgraph_genes.txt")

    meth_path = data_dir / "gene_methylation_edges.txt"
    result["methyl_edges"] = load_edge_list(meth_path, sep=",", col_pair=(0, 1)) if meth_path.exists() else None

    # miRNA (可选)
    mirna_path = config.paths.mirna_target
    if mirna_path:
        result["mirna_edges"] = load_mirna_edges(Path(mirna_path))
    else:
        result["mirna_edges"] = None

    # Reactome 层级
    gmt_path = Path(config.paths.reactome_gmt) if config.paths.reactome_gmt else None
    if gmt_path and gmt_path.exists():
        result["pathway_hierarchy"] = parse_reactome_hierarchy(gmt_path)
    else:
        result["pathway_hierarchy"] = []

    # 疾病-通路关联 (DisGeNET/KEGG, 可选)
    dp_path = data_dir / "disease_pathway_edges.txt"
    result["disease_pathway_edges"] = load_disease_pathway_edges(dp_path) if dp_path.exists() else []

    return result