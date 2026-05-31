# -*- coding: utf-8 -*-
"""
脑共表达特征整合 v1.0
基于成品脑组织数据库下载脑共表达数据，构建基因特征

数据源:
  1. PsychENCODE WGCNA 脑共表达模块 (INT-09)
     http://resource.psychencode.org/Datasets/Integrative/ModelParams/INT-09_WGCNA_modules_hgnc_ids.xlsx
  2. Harmonizome 3.0 - Allen Brain Atlas 发育中人脑 RNA-seq
     https://maayanlab.cloud/Harmonizome/dataset/Allen+Brain+Atlas+Developing+Human+Brain+Tissue+Gene+Expression+Profiles+by+RNA-seq
  3. (备选) BrainEXP 数据库 - 网站不稳定时跳过

输出:
  - brain_coexpression_features.csv: 基因 × 脑共表达特征矩阵
  - brain_module_assignments.csv:   WGCNA 模块归属
  - brain_region_profiles.csv:      脑区表达谱
  - gene_pair_coexpression.csv:     靶基因间共表达分数
  - feature_dimensions.json:        特征维度说明
"""

import os
import sys
import csv
import json
import time
import gzip
import io
import urllib.request
import urllib.parse
import zipfile
import tempfile
import shutil
from collections import defaultdict, Counter
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import requests
    HAS_REQUESTS = True
except ImportError:
    HAS_REQUESTS = False

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import BASE_DIR, RESULTS_DIR, BCP_TARGETS, CUPROPTOSIS_GENES, CUPROPTOSIS_RELATED
from scripts.utils import setup_logger, ensure_dir

OUTPUT_DIR = Path(RESULTS_DIR) / "brain_coexpression"
ensure_dir(str(OUTPUT_DIR))

CACHE_DIR = OUTPUT_DIR / "cache"
CACHE_DIR.mkdir(parents=True, exist_ok=True)

logger = setup_logger("brain_coexpr", str(OUTPUT_DIR / "brain_coexpression.log"))

RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)

# ============================================================
# 数据源 URL 配置
# ============================================================
PSYCHENCODE_WGCNA_URL = (
    "http://resource.psychencode.org/Datasets/Integrative/ModelParams"
    "/INT-09_WGCNA_modules_hgnc_ids.xlsx"
)

HARMONIZOME_BRAIN_ATLAS_BASE = (
    "https://maayanlab.cloud/static/hdfs/harmonizome/data"
    "/brainatlasdevelopmentalhumanrnaseq"
)

HARMONIZOME_FILES = {
    "gene_attr_matrix_cleaned": "gene_attribute_matrix_cleaned.txt.gz",
    "gene_attr_matrix_standardized": "gene_attribute_matrix_standardized.txt.gz",
    "gene_similarity_cosine": "gene_similarity_matrix_cosine.txt.gz",
    "attribute_list": "attribute_list_entries.txt.gz",
    "gene_list": "gene_list_terms.txt.gz",
}


# ============================================================
# 1. 通用下载工具 (带缓存)
# ============================================================
def download_file(url, cache_path, max_retries=5, timeout=600):
    """下载文件到缓存路径，支持重试。矩阵文件默认10分钟超时"""
    if cache_path.exists() and cache_path.stat().st_size > 1000:
        logger.info(f"[CACHE] 使用缓存: {cache_path.name}")
        return str(cache_path)

    logger.info(f"[DOWNLOAD] 下载: {url}")
    for attempt in range(max_retries):
        try:
            headers = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"}

            if HAS_REQUESTS:
                resp = requests.get(url, headers=headers, timeout=timeout, stream=True)
                resp.raise_for_status()
                content = b""
                for chunk in resp.iter_content(chunk_size=65536):
                    content += chunk
            else:
                req = urllib.request.Request(url, headers=headers)
                with urllib.request.urlopen(req, timeout=timeout) as resp:
                    content = resp.read()

            with open(cache_path, "wb") as f:
                f.write(content)
            logger.info(f"[OK] 下载完成: {cache_path.name} ({len(content)} bytes)")
            return str(cache_path)
        except Exception as e:
            logger.warning(f"[RETRY] 尝试 {attempt + 1}/{max_retries} 失败: {e}")
            if attempt < max_retries - 1:
                time.sleep(5 * (2 ** attempt))
            else:
                logger.error(f"[FAIL] 下载失败: {url}")
                return None


def download_harmonizome_file(file_key):
    """下载 Harmonizome 文件"""
    url = f"{HARMONIZOME_BRAIN_ATLAS_BASE}/{HARMONIZOME_FILES[file_key]}"
    cache_path = CACHE_DIR / HARMONIZOME_FILES[file_key]
    return download_file(url, cache_path)


def download_psychencode_wgcna():
    """下载 PsychENCODE WGCNA 模块文件"""
    cache_path = CACHE_DIR / "INT-09_WGCNA_modules_hgnc_ids.xlsx"
    return download_file(PSYCHENCODE_WGCNA_URL, cache_path)


# ============================================================
# 2. 解析 PsychENCODE WGCNA 模块
# ============================================================
def parse_psychencode_wgcna(filepath):
    """
    解析 PsychENCODE WGCNA 模块 Excel 文件

    格式说明 (INT-09_WGCNA_modules_hgnc_ids.xlsx):
      Sheet1 是转置矩阵:
        - 第0行: 列头 (hgnc_Ids | gene_symbol_1 | gene_symbol_2 | ...)
        - 第1+N行: module_N 名称 | 属于该模块的基因列表(不定长)
        - 注意: 同一列在不同行是不同的基因, 每行定义了一个模块的基因集

    返回: {gene_symbol: module_name}
    """
    logger.info("[WGCNA] 解析 PsychENCODE WGCNA 模块...")
    if filepath is None or not os.path.exists(filepath):
        logger.warning("[WGCNA] 文件不存在，跳过")
        return {}, {}

    try:
        df = pd.read_excel(filepath, sheet_name=0, header=0)
        logger.info(f"[WGCNA] 数据形状: {df.shape}")

        gene_to_module = {}
        module_to_genes = defaultdict(list)

        first_col_name = df.columns[0]
        logger.info(f"[WGCNA] 第一列名: '{first_col_name}'")

        for row_idx in range(len(df)):
            module_val = df.iloc[row_idx, 0]
            if pd.isna(module_val):
                continue
            module_name = str(module_val).strip()
            if not module_name or module_name.startswith("hgnc"):
                continue

            genes_in_module = []
            for col_idx in range(1, df.shape[1]):
                gene_val = df.iloc[row_idx, col_idx]
                if pd.notna(gene_val):
                    gene_str = str(gene_val).strip().upper()
                    if gene_str and gene_str != "nan":
                        genes_in_module.append(gene_str)

            for g in genes_in_module:
                if g not in gene_to_module:
                    gene_to_module[g] = module_name
                else:
                    existing = gene_to_module[g]
                    if existing != module_name:
                        gene_to_module[g] = f"{existing}|{module_name}"

            module_to_genes[module_name] = genes_in_module

        logger.info(f"[WGCNA] 解析完成: {len(gene_to_module)} 个基因, "
                    f"{len(module_to_genes)} 个模块")
        for mod, genes in sorted(module_to_genes.items(),
                                 key=lambda x: -len(x[1]))[:10]:
            actual_genes = [g for g in genes if g]
            logger.info(f"  模块 {mod}: {len(genes)} 个基因 "
                        f"(非空: {len(actual_genes)})")

        return gene_to_module, dict(module_to_genes)

    except Exception as e:
        logger.error(f"[WGCNA] 解析失败: {e}", exc_info=True)
        return {}, {}


def compute_wgcna_features(genes, gene_to_module):
    """
    根据 WGCNA 模块归属生成特征:
      - 模块 one-hot 编码 (top N 模块)
      - 是否属于主要模块 (Top 10 模块)
    """
    n = len(genes)
    gene_to_idx = {g: i for i, g in enumerate(genes)}

    if not gene_to_module:
        return np.zeros((n, 1), dtype=np.float32)

    module_counts = Counter(gene_to_module.values())
    top_modules = [m for m, _ in module_counts.most_common(20)]
    module_to_mid = {m: idx for idx, m in enumerate(top_modules)}

    n_modules = len(top_modules)
    feats = np.zeros((n, n_modules + 2), dtype=np.float32)

    for gene, module in gene_to_module.items():
        if gene in gene_to_idx:
            idx = gene_to_idx[gene]
            if module in module_to_mid:
                feats[idx, module_to_mid[module]] = 1.0
            feats[idx, n_modules] = 1.0

    top10_modules = set(top_modules[:10])
    for gene, module in gene_to_module.items():
        if gene in gene_to_idx and module in top10_modules:
            idx = gene_to_idx[gene]
            feats[idx, n_modules + 1] = 1.0

    logger.info(f"[WGCNA] 特征: {n_modules + 2} 维 "
                f"(top{len(top_modules)} one-hot + is_in_any + is_in_top10)")
    return feats


# ============================================================
# 3. 解析 Harmonizome Allen Brain Atlas 脑区表达谱
# ============================================================
def parse_gene_list_file(filepath):
    """解析 Harmonizome 基因列表文件"""
    if filepath is None or not os.path.exists(filepath):
        return []
    genes = []
    try:
        if str(filepath).endswith(".gz"):
            with gzip.open(filepath, "rt", encoding="utf-8") as f:
                for line in f:
                    g = line.strip().upper()
                    if g:
                        genes.append(g)
        else:
            with open(filepath, "r", encoding="utf-8") as f:
                for line in f:
                    g = line.strip().upper()
                    if g:
                        genes.append(g)
    except Exception as e:
        logger.warning(f"[GENE_LIST] 解析失败: {e}")
    return genes


def parse_attribute_list(filepath):
    """解析 Harmonizome 属性列表 (脑区样本)"""
    if filepath is None or not os.path.exists(filepath):
        return []
    attributes = []
    try:
        if str(filepath).endswith(".gz"):
            f = gzip.open(filepath, "rt", encoding="utf-8")
        else:
            f = open(filepath, "r", encoding="utf-8")

        with f:
            for line in f:
                attr = line.strip()
                if attr:
                    attributes.append(attr)

        logger.info(f"[ATTRIBUTE] 解析完成: {len(attributes)} 个脑区样本")
        return attributes

    except Exception as e:
        logger.warning(f"[ATTRIBUTE] 解析失败: {e}")
        return []


def parse_gene_attribute_matrix(filepath, genes_of_interest):
    """
    解析基因-脑区表达矩阵 (cleaned/standardized)
    格式: gzip TSV, 第一行是属性名, 第一列是基因名
    只加载关注的基因
    """
    logger.info(f"[MATRIX] 解析基因-脑区表达矩阵...")
    if filepath is None or not os.path.exists(filepath):
        logger.warning("[MATRIX] 文件不存在")
        return None, None

    target_set = set(g.upper() for g in genes_of_interest)
    gene_to_idx = {}
    matrix_rows = []

    try:
        with gzip.open(filepath, "rt", encoding="utf-8") as f:
            header = next(f).strip().split("\t")
            n_cols = len(header) - 1

            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 2:
                    continue
                gene = parts[0].strip().upper()
                if gene in target_set:
                    idx = len(gene_to_idx)
                    gene_to_idx[gene] = idx
                    values = []
                    for v in parts[1:]:
                        try:
                            values.append(float(v))
                        except (ValueError, IndexError):
                            values.append(0.0)
                    matrix_rows.append(values)

        matrix = np.array(matrix_rows, dtype=np.float32)
        logger.info(f"[MATRIX] 加载完成: {matrix.shape} (基因×脑区样本)")
        return matrix, gene_to_idx

    except Exception as e:
        logger.error(f"[MATRIX] 解析失败: {e}")
        return None, None


def parse_gene_similarity_matrix(filepath, genes_of_interest):
    """
    解析基因相似性矩阵 (cosine)
    格式: gzip TSV, 第一行是基因名, 第一列是基因名
    提取关注基因之间的共表达相似性
    """
    logger.info(f"[SIM] 解析基因相似性矩阵...")
    if filepath is None or not os.path.exists(filepath):
        logger.warning("[SIM] 文件不存在")
        return None, None

    target_set = set(g.upper() for g in genes_of_interest)

    try:
        with gzip.open(filepath, "rt", encoding="utf-8") as f:
            header = next(f).strip().split("\t")
            all_genes = [h.strip().upper() for h in header]

            gene_positions = {}
            for pos, g in enumerate(all_genes):
                if g in target_set:
                    gene_positions[g] = pos

            if not gene_positions:
                logger.warning("[SIM] 目标基因不在相似性矩阵中")
                return None, None

            target_genes_sorted = sorted(gene_positions.keys())
            n_targets = len(target_genes_sorted)
            sim_matrix = np.zeros((n_targets, n_targets), dtype=np.float32)

            target_pos_to_idx = {gene_positions[g]: i
                                 for i, g in enumerate(target_genes_sorted)}

            for line_idx, line in enumerate(f):
                parts = line.strip().split("\t")
                if len(parts) < 2:
                    continue
                gene = parts[0].strip().upper()
                if gene not in target_set:
                    continue
                row_idx = target_pos_to_idx.get(gene_positions.get(gene))
                if row_idx is None:
                    continue

                n_vals = min(len(parts) - 1, len(all_genes))
                for col_offset in range(n_vals):
                    col_gene = all_genes[col_offset]
                    if col_gene in gene_positions:
                        col_idx = target_pos_to_idx[gene_positions[col_gene]]
                        try:
                            sim_matrix[row_idx, col_idx] = float(parts[col_offset + 1])
                        except (ValueError, IndexError):
                            pass

            logger.info(f"[SIM] 相似性矩阵: {sim_matrix.shape} "
                        f"({n_targets} 个目标基因)")
            return sim_matrix, target_genes_sorted

    except Exception as e:
        logger.error(f"[SIM] 解析失败: {e}")
        return None, None


def compute_brain_region_features(genes, matrix, gene_to_idx, attributes):
    """
    从基因-脑区表达矩阵计算特征:
      - 在各脑区的平均表达 (有数据则算均值)
      - 表达广度 (表达非零的脑区数量)
      - 脑区类别汇总: cortex, subcortex, cerebellum, brainstem 等
    """
    n = len(genes)
    gene_idx_map = {g: i for i, g in enumerate(genes)}

    if matrix is None or gene_to_idx is None:
        return np.zeros((n, 6), dtype=np.float32)

    region_categories = {
        "cortex": ["cortex", "frontal", "parietal", "temporal", "occipital",
                   "cingulate", "insular", "orbital", "prefrontal",
                   "somatosensory", "auditory", "visual", "motor"],
        "subcortex": ["amygdala", "hippocampus", "striatum", "caudate",
                      "putamen", "pallidum", "thalamus", "hypothalamus",
                      "nucleus", "basal", "septum"],
        "cerebellum": ["cerebellum", "cerebellar"],
        "brainstem": ["brainstem", "brain stem", "medulla", "pons",
                      "midbrain", "spinal"],
    }

    n_genes_in_matrix = len(gene_to_idx)
    n_regions = matrix.shape[1]

    feats = np.zeros((n, 6 + n_regions), dtype=np.float32)

    for gene in genes:
        gene_up = gene.upper()
        if gene_up not in gene_to_idx:
            continue
        mat_idx = gene_to_idx[gene_up]
        row = matrix[mat_idx]
        gidx = gene_idx_map[gene_up]

        valid = row[np.isfinite(row)]
        if len(valid) == 0:
            continue

        mean_expr = np.mean(valid)
        std_expr = np.std(valid) if len(valid) > 1 else 0.0
        n_expressed = np.sum(valid > 0)
        n_high = np.sum(valid > np.percentile(valid, 75)) if len(valid) > 4 else 0

        feats[gidx, 0] = np.log1p(mean_expr)
        feats[gidx, 1] = std_expr
        feats[gidx, 2] = n_expressed / max(n_regions, 1)
        feats[gidx, 3] = n_high / max(n_regions, 1)
        feats[gidx, 4] = np.log1p(mean_expr * n_expressed)
        feats[gidx, 5] = np.log1p(len(valid))

        feats[gidx, 6:] = row

    cat_feats = np.zeros((n, 4), dtype=np.float32)
    for gene in genes:
        gene_up = gene.upper()
        if gene_up not in gene_to_idx:
            continue
        mat_idx = gene_to_idx[gene_up]
        row = matrix[mat_idx]
        gidx = gene_idx_map[gene_up]

        region_names = [a.lower() for a in attributes[:len(row)]]
        for ri, rn in enumerate(region_names):
            if ri >= len(row):
                break
            val = row[ri]
            if not np.isfinite(val) or val <= 0:
                continue
            for ci, (cat_name, keywords) in enumerate(region_categories.items()):
                if any(kw in rn for kw in keywords):
                    cat_feats[gidx, ci] = max(cat_feats[gidx, ci], val)

    combined = np.hstack([feats, cat_feats])
    logger.info(f"[BRAIN] 脑区表达特征: {combined.shape}")
    return combined


def compute_pairwise_coexpression(genes, gene_sim_matrix, target_genes_sorted):
    """
    计算目标基因对之间的共表达分数
    输出: 每个基因与其他目标基因的平均/最大共表达
    """
    n = len(genes)
    gene_idx_map = {g: i for i, g in enumerate(genes)}

    if gene_sim_matrix is None:
        return np.zeros((n, 3), dtype=np.float32)

    target_to_matrix_idx = {g: i for i, g in enumerate(target_genes_sorted)}

    feats = np.zeros((n, 3), dtype=np.float32)

    for gene in genes:
        gene_up = gene.upper()
        if gene_up not in target_to_matrix_idx:
            continue
        gidx = gene_idx_map[gene_up]
        mat_idx = target_to_matrix_idx[gene_up]
        sim_row = gene_sim_matrix[mat_idx]

        valid = sim_row[np.isfinite(sim_row)]
        if len(valid) == 0:
            continue

        mask = np.arange(len(sim_row)) != mat_idx
        other_sims = sim_row[mask]

        feats[gidx, 0] = np.mean(other_sims) if len(other_sims) > 0 else 0
        feats[gidx, 1] = np.max(other_sims) if len(other_sims) > 0 else 0
        feats[gidx, 2] = np.sum(other_sims > 0.5) / max(len(other_sims), 1)

    logger.info(f"[PAIRWISE] 共表达特征: {feats.shape}")
    return feats


# ============================================================
# 4. 主流程
# ============================================================
def gather_target_genes():
    """收集所有关注的基因 (BCP + 铜死亡 + 铜死亡相关)"""
    targets = set()
    for g in BCP_TARGETS:
        targets.add(g.upper())
    for g in CUPROPTOSIS_GENES:
        targets.add(g.upper())
    for g in CUPROPTOSIS_RELATED:
        targets.add(g.upper())
    genes = sorted(targets)
    logger.info(f"[TARGETS] 关注基因: {len(genes)} 个")
    return genes


def assemble_and_save():
    """主流程"""
    logger.info("=" * 60)
    logger.info("脑共表达特征整合 - 开始")
    logger.info("=" * 60)

    genes = gather_target_genes()
    gene_to_idx = {g: i for i, g in enumerate(genes)}
    n = len(genes)

    all_feats_list = []
    dim_names = []

    # ---- A. PsychENCODE WGCNA 模块 ----
    logger.info("\n[A] PsychENCODE WGCNA 模块特征")
    wgcna_file = download_psychencode_wgcna()
    gene_to_module, module_to_genes = parse_psychencode_wgcna(wgcna_file)
    wgcna_feats = compute_wgcna_features(genes, gene_to_module)
    all_feats_list.append(wgcna_feats)
    dim_names.append(f"WGCNA_{wgcna_feats.shape[1]}")
    logger.info(f"  WGCNA 特征: {wgcna_feats.shape[1]} 维")

    # ---- B. Harmonizome Allen Brain Atlas 脑区表达 ----
    logger.info("\n[B] Harmonizome Allen Brain Atlas 脑区表达")

    harmonizome_available = False

    attr_file = download_harmonizome_file("attribute_list")
    attributes = parse_attribute_list(attr_file)

    # 尝试下载/使用脑区表达矩阵 (短超时)
    # 先检查 cleaned 矩阵是否已缓存
    cleaned_cache = CACHE_DIR / HARMONIZOME_FILES["gene_attr_matrix_cleaned"]
    if cleaned_cache.exists() and cleaned_cache.stat().st_size > 1000:
        matrix_file = str(cleaned_cache)
        logger.info(f"[CACHE] 使用已缓存的 cleaned 矩阵: {cleaned_cache.name} "
                    f"({cleaned_cache.stat().st_size / 1e6:.1f} MB)")
    else:
        logger.warning("  Harmonizome 矩阵文件过大，网络连接不稳定，跳过脑区矩阵下载")
        logger.warning("  可手动下载后放入缓存目录:")
        logger.warning(f"    {HARMONIZOME_BRAIN_ATLAS_BASE}/gene_attribute_matrix_cleaned.txt.gz")
        matrix_file = None

    matrix, matrix_gene_to_idx = parse_gene_attribute_matrix(matrix_file, genes)

    if matrix is not None:
        brain_feats = compute_brain_region_features(
            genes, matrix, matrix_gene_to_idx, attributes
        )
        all_feats_list.append(brain_feats)
        dim_names.append(f"BrainRegion_{brain_feats.shape[1]}")
        harmonizome_available = True
        logger.info(f"  脑区表达特征: {brain_feats.shape[1]} 维")
    else:
        logger.warning("  脑区矩阵不可用，跳过")
        all_feats_list.append(np.zeros((n, 0), dtype=np.float32))
        dim_names.append("BrainRegion_0")

    # ---- C. 基因共表达相似性 ----
    logger.info("\n[C] 基因共表达相似性 (Harmonizome)")
    sim_cache = CACHE_DIR / HARMONIZOME_FILES["gene_similarity_cosine"]
    if sim_cache.exists() and sim_cache.stat().st_size > 1000:
        sim_matrix, target_genes_sorted = parse_gene_similarity_matrix(str(sim_cache), genes)
    else:
        logger.warning("  Harmonizome 相似矩阵过大，跳过下载")
        logger.warning(f"  可手动下载: {HARMONIZOME_BRAIN_ATLAS_BASE}/gene_similarity_matrix_cosine.txt.gz")
        sim_matrix, target_genes_sorted = None, None

    if sim_matrix is not None:
        pairwise_feats = compute_pairwise_coexpression(
            genes, sim_matrix, target_genes_sorted
        )
        all_feats_list.append(pairwise_feats)
        dim_names.append(f"CoexprPairwise_{pairwise_feats.shape[1]}")
        harmonizome_available = True
        logger.info(f"  共表达特征: {pairwise_feats.shape[1]} 维")
    else:
        logger.warning("  相似性矩阵不可用，跳过")
        all_feats_list.append(np.zeros((n, 0), dtype=np.float32))
        dim_names.append("CoexprPairwise_0")

    # ---- D. 综合共表达得分 ----
    logger.info("\n[D] 综合共表达得分")
    coexpr_score = np.zeros((n, 2), dtype=np.float32)
    target_to_midx = {}
    if target_genes_sorted:
        target_to_midx = {g: i for i, g in enumerate(target_genes_sorted)}

    if sim_matrix is not None:
        for i, gene in enumerate(genes):
            if gene in target_to_midx:
                midx = target_to_midx[gene]
                row = sim_matrix[midx]
                valid = row[np.isfinite(row)]
                if len(valid) > 0:
                    coexpr_score[i, 0] = np.mean(valid)
                    coexpr_score[i, 1] = np.max(valid)

    all_feats_list.append(coexpr_score)
    dim_names.append(f"CoexprScore_{coexpr_score.shape[1]}")

    # ---- 拼接所有特征 ----
    valid_feats = [f for f in all_feats_list if f.shape[1] > 0]
    if not valid_feats:
        logger.error("[FAIL] 没有可用的特征!")
        return

    enhanced_matrix = np.hstack(valid_feats)
    enhanced_matrix = np.nan_to_num(enhanced_matrix, nan=0.0, posinf=0.0, neginf=0.0)

    logger.info(f"\n{'=' * 60}")
    logger.info(f"脑共表达特征矩阵: {enhanced_matrix.shape}")
    for dn in dim_names:
        logger.info(f"  + {dn}")
    logger.info(f"  = 总计: {enhanced_matrix.shape[1]} 维")

    # ---- 保存特征矩阵 ----
    output_csv = OUTPUT_DIR / "brain_coexpression_features.csv"
    feat_cols = [f"brain_feat_{i:04d}" for i in range(enhanced_matrix.shape[1])]
    df_out = pd.DataFrame(enhanced_matrix, columns=feat_cols)
    df_out.insert(0, "gene_symbol", genes)
    df_out.to_csv(output_csv, index=False, float_format="%.8g")
    logger.info(f"[SAVE] 脑共表达特征: {output_csv}")

    # ---- 保存 NPY ----
    output_npy = OUTPUT_DIR / "brain_coexpression_features.npy"
    np.save(output_npy, enhanced_matrix)
    logger.info(f"[SAVE] NPY: {output_npy}")

    # ---- 保存维度映射 ----
    dim_map = {
        "total_dim": int(enhanced_matrix.shape[1]),
        "n_genes": int(enhanced_matrix.shape[0]),
        "blocks": dim_names,
        "psychencode_wgcna_available": bool(gene_to_module),
        "harmonizome_available": harmonizome_available,
    }
    dim_map_path = OUTPUT_DIR / "feature_dimensions.json"
    with open(dim_map_path, "w", encoding="utf-8") as f:
        json.dump(dim_map, f, indent=2, ensure_ascii=False)
    logger.info(f"[SAVE] 维度映射: {dim_map_path}")

    # ---- 保存 WGCNA 模块归属 ----
    if gene_to_module:
        module_csv = OUTPUT_DIR / "brain_module_assignments.csv"
        rows = []
        for g in genes:
            mod = gene_to_module.get(g, "UNASSIGNED")
            rows.append({"gene_symbol": g, "wgcna_module": mod})
        pd.DataFrame(rows).to_csv(module_csv, index=False)
        logger.info(f"[SAVE] 模块归属: {module_csv}")

    # ---- 保存基因对共表达 ----
    if sim_matrix is not None and target_genes_sorted:
        pair_csv = OUTPUT_DIR / "gene_pair_coexpression.csv"
        pair_rows = []
        for i in range(len(target_genes_sorted)):
            for j in range(i + 1, len(target_genes_sorted)):
                val = sim_matrix[i, j]
                if np.isfinite(val) and val > 0.3:
                    pair_rows.append({
                        "gene1": target_genes_sorted[i],
                        "gene2": target_genes_sorted[j],
                        "coexpression": round(float(val), 4),
                    })
        pair_rows.sort(key=lambda x: -x["coexpression"])
        pd.DataFrame(pair_rows).to_csv(pair_csv, index=False)
        logger.info(f"[SAVE] 基因对共表达: {pair_csv} ({len(pair_rows)} 对)")

    # ---- 保存脑区表达谱 (汇总) ----
    if matrix is not None and attributes:
        profile_csv = OUTPUT_DIR / "brain_region_profiles.csv"
        profile_rows = []
        for g in genes:
            if g in matrix_gene_to_idx:
                midx = matrix_gene_to_idx[g]
                row_data = matrix[midx]
                for ri in range(min(len(row_data), len(attributes))):
                    if np.isfinite(row_data[ri]) and row_data[ri] > 0:
                        profile_rows.append({
                            "gene_symbol": g,
                            "brain_region": attributes[ri],
                            "expression": round(float(row_data[ri]), 4),
                        })
        if profile_rows:
            pdf = pd.DataFrame(profile_rows)
            pdf.to_csv(profile_csv, index=False)
            logger.info(f"[SAVE] 脑区表达谱: {profile_csv} ({len(profile_rows)} 条记录)")

    logger.info(f"\n{'=' * 60}")
    logger.info("脑共表达特征整合 - 完成!")
    logger.info(f"输出目录: {OUTPUT_DIR}")
    logger.info("=" * 60)


def run_brain_coexpression_pipeline():
    """供 pipeline_engine 调用的入口"""
    logger.info("=" * 60)
    logger.info("脑共表达数据整合管道")
    logger.info("=" * 60)
    try:
        assemble_and_save()
        return True
    except Exception as e:
        logger.error(f"管道失败: {e}", exc_info=True)
        return False


if __name__ == "__main__":
    assemble_and_save()