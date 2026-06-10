# -*- coding: utf-8 -*-
"""
从 STRING v12.0 预训练蛋白质嵌入中提取指定基因的特征向量

STRONG> 输入:
  1. subgraph_genes.txt           — 基因 Symbol 列表（一行一个，首行可能有标题）
  2. 9606.protein.aliases.v12.0.txt  — STRING 别名文件（TSV: protein_id, alias, source）
  3. 9606.protein.sequence.embeddings.v12.0.h5  — STRING 嵌入 HDF5 文件

STRONG> 输出:
  subgraph_embeddings.csv — 第一列 gene_symbol，后续列 feat_0001, feat_0002, ...
  failed_genes.txt        — 未成功提取特征的基因列表

STRONG> HDF5 结构 (v12.0):
  - /proteins:    (N_proteins,)   bytes → STRING 蛋白 ID (如 9606.ENSP00000000233)
  - /embeddings:  (N_proteins, 1024)  float16 → 逐行对应 proteins

STRONG> 依赖:
  pip install h5py numpy pandas tqdm
"""

import os
import sys
import logging
import warnings
from typing import Optional

import numpy as np
import pandas as pd
import h5py
from tqdm import tqdm

warnings.filterwarnings("ignore")
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

# ==================== 配置 ====================

GENE_LIST_FILE = r"C:\Users\Jy-Mentor-7\Desktop\subgraph_genes.txt"
ALIASES_FILE = r"C:\Users\Jy-Mentor-7\Desktop\9606蛋白质\人靶点\9606.protein.aliases.v12.0.txt"
HDF5_FILE = r"C:\Users\Jy-Mentor-7\Desktop\9606蛋白质\9606.protein.sequence.embeddings.v12.0.h5"
OUTPUT_CSV = "subgraph_embeddings.csv"
FAILED_TXT = "failed_genes.txt"

ALIAS_SOURCE_FOR_SYMBOL = "Ensembl_HGNC_symbol"
SKIP_HEADER = True

# ==================== 基因列表读取 ====================

def read_gene_list(filepath: str, skip_header: bool = True) -> list[str]:
    """读取基因Symbol列表，每行一个，转为大写"""
    if not os.path.exists(filepath):
        logger.error("文件不存在: %s", filepath)
        sys.exit(1)
    genes = []
    with open(filepath, "r", encoding="utf-8") as f:
        lines = f.readlines()
    start = 1 if skip_header else 0
    for line in lines[start:]:
        s = line.strip()
        if s and not s.startswith("#"):
            genes.append(s.upper())
    logger.info("读取到 %d 个基因", len(genes))
    return genes


# ==================== 别名映射构建 ====================

def build_alias_mapping(
    aliases_file: str,
    source_filter: str = "Ensembl_HGNC_symbol",
) -> dict[str, str]:
    """
    从 STRING 别名文件中构建 Gene Symbol → STRING 蛋白 ID 的映射

    STRING 别名文件格式（TSV，无压缩）:
      #string_protein_id\talias\tsource
      9606.ENSP00000000233\tARF5\tEnsembl_HGNC_symbol
      ...

    返回: {GENE_SYMBOL_UPPER: protein_id}
    一个 Symbol 可能对应多个蛋白 ID，取第一个遇到的。
    """
    logger.info("正在构建别名映射 (source=%s) ...", source_filter)
    symbol_to_protein: dict[str, str] = {}

    if not os.path.exists(aliases_file):
        logger.error("别名文件不存在: %s", aliases_file)
        sys.exit(1)

    with open(aliases_file, "r", encoding="utf-8") as f:
        total_lines = sum(1 for _ in f)

    with open(aliases_file, "r", encoding="utf-8") as f:
        for line in tqdm(f, total=total_lines, desc="解析别名文件"):
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            protein_id, alias, source = parts[0], parts[1], parts[2]
            if source == source_filter:
                gene_upper = alias.upper()
                if gene_upper not in symbol_to_protein:
                    symbol_to_protein[gene_upper] = protein_id

    logger.info("别名映射完成: %d 个唯一 Gene Symbol", len(symbol_to_protein))
    return symbol_to_protein


# ==================== HDF5 嵌入提取 ====================

def extract_embeddings(
    genes: list[str],
    symbol_to_protein: dict[str, str],
    hdf5_file: str,
) -> tuple[dict[str, np.ndarray], list[str]]:
    """
    从 HDF5 中提取基因对应的嵌入向量

    参数:
      genes: 基因 Symbol 列表（大写）
      symbol_to_protein: Symbol → STRING protein ID 映射
      hdf5_file: HDF5 文件路径

    返回:
      (features_dict, failed_genes)
    """
    logger.info("正在加载 HDF5 嵌入文件 ...")
    with h5py.File(hdf5_file, "r") as h5:

        # 解析 HDF5 结构: 可能有两种格式
        # 格式A (v12.0 标准): /proteins (N,) + /embeddings (N, D)
        # 格式B (旧版): 键为蛋白 ID，值为嵌入

        if "proteins" in h5 and "embeddings" in h5:
            # 格式A: 矩阵格式
            proteins_ds = h5["proteins"]
            embeddings_ds = h5["embeddings"]
            n_proteins = proteins_ds.shape[0]
            embed_dim = embeddings_ds.shape[1]
            logger.info("HDF5 格式: 矩阵型 | proteins=%d | embeddings=(%d, %d) | dtype=%s",
                        n_proteins, n_proteins, embed_dim, embeddings_ds.dtype)

            # 构建 protein_id → row_index 映射
            protein_to_idx: dict[str, int] = {}
            logger.info("正在构建蛋白ID索引 ...")
            for i in tqdm(range(n_proteins), desc="索引蛋白ID"):
                pid = proteins_ds[i]
                if isinstance(pid, bytes):
                    pid = pid.decode("utf-8")
                protein_to_idx[pid] = i

            # 提取嵌入
            features: dict[str, np.ndarray] = {}
            failed: list[str] = []

            logger.info("正在提取嵌入向量 ...")
            for gene in tqdm(genes, desc="提取嵌入"):
                emb = None

                # 策略1: 通过 Symbol → protein_id 映射
                protein_id = symbol_to_protein.get(gene)
                if protein_id and protein_id in protein_to_idx:
                    idx = protein_to_idx[protein_id]
                    emb = embeddings_ds[idx, :].astype(np.float32)
                elif protein_id and protein_id not in protein_to_idx:
                    # 尝试变体格式
                    for variant in _protein_id_variants(protein_id):
                        if variant in protein_to_idx:
                            idx = protein_to_idx[variant]
                            emb = embeddings_ds[idx, :].astype(np.float32)
                            break

                # 策略2: 直接查找基因名作为蛋白ID（备选）
                if emb is None:
                    if gene in protein_to_idx:
                        idx = protein_to_idx[gene]
                        emb = embeddings_ds[idx, :].astype(np.float32)
                    else:
                        for variant in _protein_id_variants(gene):
                            if variant in protein_to_idx:
                                idx = protein_to_idx[variant]
                                emb = embeddings_ds[idx, :].astype(np.float32)
                                break

                if emb is not None:
                    features[gene] = emb
                else:
                    failed.append(gene)

        else:
            # 格式B: 键-值格式（旧版兼容）
            logger.info("HDF5 格式: 键-值型（旧版兼容模式）")
            keys = list(h5.keys())
            embed_dim = h5[keys[0]].shape[0] if len(h5[keys[0]].shape) == 1 else h5[keys[0]].shape[-1]
            logger.info("蛋白数: %d | 嵌入维度: %d", len(keys), embed_dim)

            features = {}
            failed = []

            for gene in tqdm(genes, desc="提取嵌入"):
                emb = None
                protein_id = symbol_to_protein.get(gene)

                if protein_id and protein_id in h5:
                    emb = h5[protein_id][:].astype(np.float32).flatten()
                elif protein_id:
                    for variant in _protein_id_variants(protein_id):
                        if variant in h5:
                            emb = h5[variant][:].astype(np.float32).flatten()
                            break

                if emb is None and gene in h5:
                    emb = h5[gene][:].astype(np.float32).flatten()

                if emb is not None:
                    features[gene] = emb
                else:
                    failed.append(gene)

    logger.info("提取完成: 成功 %d | 失败 %d", len(features), len(failed))
    return features, failed


def _protein_id_variants(protein_id: str) -> list[str]:
    """生成 STRING 蛋白 ID 的可能变体"""
    variants = [protein_id]
    # 9606.ENSP00000000233 ↔ ENSP00000000233 (去掉物种前缀)
    if protein_id.startswith("9606."):
        variants.append(protein_id[5:])
    if not protein_id.startswith("9606."):
        variants.append("9606." + protein_id)
    return variants


# ==================== 主流程 ====================

def main():
    logger.info("=" * 60)
    logger.info("STRING v12 蛋白质嵌入提取")
    logger.info("=" * 60)

    # 1. 读取基因列表
    genes = read_gene_list(GENE_LIST_FILE, skip_header=SKIP_HEADER)
    n_input = len(genes)

    # 2. 构建别名映射
    symbol_to_protein = build_alias_mapping(ALIASES_FILE, ALIAS_SOURCE_FOR_SYMBOL)

    # 3. 提取嵌入
    features, failed = extract_embeddings(genes, symbol_to_protein, HDF5_FILE)

    # 4. 保存失败基因
    if failed:
        with open(FAILED_TXT, "w", encoding="utf-8") as f:
            for g in failed:
                f.write(g + "\n")
        logger.warning("失败基因已保存: %s (%d 个)", FAILED_TXT, len(failed))

    if not features:
        logger.error("未成功提取任何基因特征，请检查输入文件")
        sys.exit(1)

    # 5. 组装特征矩阵并按输入顺序排列
    ordered_genes = [g for g in genes if g in features]
    feature_matrix = np.array([features[g] for g in ordered_genes], dtype=np.float32)
    embed_dim = feature_matrix.shape[1]

    # 6. 保存 CSV
    columns = ["gene_symbol"] + [f"feat_{i+1:04d}" for i in range(embed_dim)]
    df = pd.DataFrame(feature_matrix, columns=columns[1:])
    df.insert(0, "gene_symbol", ordered_genes)
    df.to_csv(OUTPUT_CSV, index=False, encoding="utf-8")

    file_size_mb = os.path.getsize(OUTPUT_CSV) / (1024 * 1024)
    n_success = len(ordered_genes)
    n_failed = len(failed)

    logger.info("=" * 60)
    logger.info("结果已保存: %s", OUTPUT_CSV)
    logger.info("  输入基因数:  %d", n_input)
    logger.info("  成功提取:    %d (%.1f%%)", n_success, 100 * n_success / n_input)
    logger.info("  失败/未找到: %d (%.1f%%)", n_failed, 100 * n_failed / n_input if n_input else 0)
    logger.info("  特征维度:    %d", embed_dim)
    logger.info("  文件大小:    %.1f MB", file_size_mb)
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
