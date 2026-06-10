# -*- coding: utf-8 -*-
"""
大规模批量蛋白质特征嵌入脚本
输入: C:/Users/Jy-Mentor-7/Desktop/subgraph_genes.txt (5723个基因)
输出: subgraph_embeddings.csv

核心优化:
  1. ESM-2分批推理 (batch_size=500) → GPU显存友好
  2. UniProt映射本地缓存 → 避免重复API调用
  3. 每批自动保存checkpoint → 断点续跑
  4. 进度条 + 耗时预估

依赖:
  pip install torch transformers pandas numpy scikit-learn requests tqdm
"""

import os
import sys
import time
import json
import logging
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Optional

os.environ.setdefault("HF_ENDPOINT", "https://hf-mirror.com")

import numpy as np
import pandas as pd
import requests
from tqdm import tqdm

warnings.filterwarnings("ignore")
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

# ==================== 配置区 ====================

GENE_LIST_FILE = r"C:\Users\Jy-Mentor-7\Desktop\subgraph_genes.txt"
OUTPUT_CSV = "subgraph_embeddings.csv"
CHECKPOINT_FILE = "subgraph_checkpoint.json"
CACHE_FILE = "uniprot_cache.json"

ESM2_MODEL_NAME = "facebook/esm2_t6_8M_UR50D"
BATCH_SIZE = 500
MAX_SEQUENCE_LENGTH = 1022
REQUEST_TIMEOUT = 30
REQUEST_DELAY = 0.15
SKIP_HEADER = True


# ==================== 缓存管理 ====================

class UniProtCache:
    def __init__(self, cache_file: str):
        self.cache_file = cache_file
        self.data: dict[str, dict] = {}
        self._load()

    def _load(self):
        if os.path.exists(self.cache_file):
            with open(self.cache_file, "r", encoding="utf-8") as f:
                self.data = json.load(f)
            logger.info("加载UniProt缓存: %d 条记录", len(self.data))

    def save(self):
        with open(self.cache_file, "w", encoding="utf-8") as f:
            json.dump(self.data, f, ensure_ascii=False, indent=1)
        logger.debug("缓存已保存: %s", self.cache_file)

    def get(self, gene: str) -> Optional[dict]:
        return self.data.get(gene.upper())

    def set(self, gene: str, uniprot_id: Optional[str], sequence: Optional[str]):
        self.data[gene.upper()] = {
            "uniprot_id": uniprot_id,
            "sequence": sequence,
        }
        if len(self.data) % 500 == 0:
            self.save()

    def save_force(self):
        self.save()


# ==================== 基因读取 ====================

def read_gene_list(filepath: str, skip_header: bool = True) -> list[str]:
    if not os.path.exists(filepath):
        logger.error("文件不存在: %s", filepath)
        sys.exit(1)
    genes = []
    with open(filepath, "r", encoding="utf-8") as f:
        lines = f.readlines()
    start = 1 if skip_header else 0
    for line in lines[start:]:
        line = line.strip()
        if line and not line.startswith("#"):
            genes.append(line.upper())
    logger.info("读取到 %d 个基因 (文件: %s)", len(genes), filepath)
    return genes


# ==================== UniProt映射 & 序列获取 ====================

def map_symbol_to_uniprot(gene_symbol: str) -> Optional[str]:
    url = "https://mygene.info/v3/query"
    params = {
        "q": gene_symbol,
        "scopes": "symbol",
        "fields": "uniprot",
        "species": "human",
        "size": 3,
    }
    try:
        resp = requests.get(url, params=params, timeout=REQUEST_TIMEOUT)
        resp.raise_for_status()
        data = resp.json()
        hits = data.get("hits", [])
        for hit in hits:
            uniprot_info = hit.get("uniprot", {})
            if isinstance(uniprot_info, dict):
                swissprot = uniprot_info.get("Swiss-Prot")
                if swissprot:
                    return swissprot[0] if isinstance(swissprot, list) else swissprot
            elif isinstance(uniprot_info, list):
                for entry in uniprot_info:
                    if isinstance(entry, dict) and "Swiss-Prot" in entry:
                        val = entry["Swiss-Prot"]
                        return val[0] if isinstance(val, list) else val
        for hit in hits:
            uniprot_info = hit.get("uniprot", {})
            if isinstance(uniprot_info, dict):
                for key in ("TrEMBL",):
                    val = uniprot_info.get(key)
                    if val:
                        return val[0] if isinstance(val, list) else val
        return None
    except requests.exceptions.RequestException:
        return None


def fetch_protein_sequence(uniprot_id: str) -> Optional[str]:
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    try:
        resp = requests.get(url, timeout=REQUEST_TIMEOUT)
        resp.raise_for_status()
        lines = resp.text.strip().splitlines()
        seq_lines = [l.strip() for l in lines if not l.startswith(">") and l.strip()]
        return "".join(seq_lines) or None
    except requests.exceptions.RequestException:
        return None


def resolve_single_gene(gene: str) -> tuple[str, Optional[str]]:
    uniprot_id = map_symbol_to_uniprot(gene)
    seq = None
    if uniprot_id:
        seq = fetch_protein_sequence(uniprot_id)
    time.sleep(REQUEST_DELAY)
    return gene, seq


def resolve_genes(
    genes: list[str],
    cache: UniProtCache,
) -> tuple[dict[str, str], list[str]]:
    gene_to_seq: dict[str, str] = {}
    failed: list[str] = []

    unresolved = [g for g in genes if cache.get(g) is None]
    resolved = [g for g in genes if cache.get(g) is not None]

    logger.info("缓存命中: %d / %d", len(resolved), len(genes))

    for g in resolved:
        entry = cache.get(g)
        if entry and entry.get("sequence"):
            gene_to_seq[g] = entry["sequence"]
        elif entry and entry.get("uniprot_id"):
            seq = fetch_protein_sequence(entry["uniprot_id"])
            if seq:
                gene_to_seq[g] = seq
                cache.set(g, entry["uniprot_id"], seq)
            else:
                failed.append(g)
        else:
            failed.append(g)

    logger.info("并发API查询: %d 个基因 (workers=20)", len(unresolved))

    with ThreadPoolExecutor(max_workers=20) as executor:
        futures = {executor.submit(resolve_single_gene, g): g for g in unresolved}
        for future in tqdm(as_completed(futures), total=len(unresolved), desc="UniProt映射+序列"):
            gene = futures[future]
            try:
                _, seq = future.result()
                if seq:
                    gene_to_seq[gene] = seq
                    entry = cache.get(gene)
                    uniprot_id = entry["uniprot_id"] if entry else None
                    cache.set(gene, uniprot_id, seq)
                else:
                    failed.append(gene)
            except Exception as e:
                logger.warning("  %s 查询异常: %s", gene, e)
                failed.append(gene)

    cache.save_force()
    logger.info("序列获取完成: 成功 %d | 失败 %d", len(gene_to_seq), len(failed))
    return gene_to_seq, failed


# ==================== ESM-2 分批推理 ====================

def load_esm2_model(model_name: str):
    logger.info("加载ESM-2模型: %s ...", model_name)
    try:
        from transformers import AutoModel, AutoTokenizer
    except ImportError:
        logger.error("请安装: pip install transformers torch")
        sys.exit(1)
    try:
        tokenizer = AutoTokenizer.from_pretrained(model_name)
        model = AutoModel.from_pretrained(model_name)
    except Exception as e:
        raise RuntimeError(f"模型加载失败: {e}")

    hidden_size = model.config.hidden_size
    logger.info("模型: %s | 维度: %d | 层数: %d",
                model_name, hidden_size, model.config.num_hidden_layers)
    model.eval()
    return model, tokenizer, hidden_size


def extract_embeddings_batch(
    sequences: list[str],
    gene_names: list[str],
    model,
    tokenizer,
    max_length: int = 1022,
) -> dict[str, np.ndarray]:
    import torch

    # 截断 & tokenize
    truncated = [s[:max_length] for s in sequences]
    inputs = tokenizer(
        truncated,
        return_tensors="pt",
        padding=True,
        truncation=True,
        max_length=max_length + 2,
        add_special_tokens=True,
    )

    with torch.no_grad():
        outputs = model(**inputs)
        last_hidden = outputs.last_hidden_state

    # 对每条序列做平均池化
    results = {}
    for i, gene in enumerate(gene_names):
        seq_len = (inputs["attention_mask"][i] == 1).sum()
        emb = last_hidden[i, :seq_len].mean(dim=0).cpu().numpy().astype(np.float32)
        results[gene] = emb

    return results


# ==================== Checkpoint管理 ====================

def load_checkpoint(checkpoint_file: str) -> dict:
    if os.path.exists(checkpoint_file):
        with open(checkpoint_file, "r", encoding="utf-8") as f:
            data = json.load(f)
        logger.info("加载checkpoint: %d 个基因已处理", len(data.get("done", [])))
        return data
    return {"done": [], "features": {}}


def save_checkpoint(checkpoint_file: str, done_genes: list, features: dict):
    # 只保存基因名列表 + 特征路径（特征本身在CSV中）
    data = {
        "done": done_genes,
        "count": len(done_genes),
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    with open(checkpoint_file, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)


# ==================== 主流程 ====================

def main():
    logger.info("=" * 60)
    logger.info("大规模蛋白质特征嵌入 - ESM-2批处理")
    logger.info("=" * 60)

    # 1. 读取基因列表
    genes = read_gene_list(GENE_LIST_FILE, skip_header=SKIP_HEADER)
    total = len(genes)

    # 2. 加载缓存 & checkpoint
    cache = UniProtCache(CACHE_FILE)
    checkpoint = load_checkpoint(CHECKPOINT_FILE)
    done_genes = checkpoint.get("done", [])

    # 3. 过滤已处理的基因
    remaining = [g for g in genes if g not in done_genes]
    logger.info("待处理: %d | 已完成: %d", len(remaining), len(done_genes))

    if not remaining and len(done_genes) > 0:
        logger.info("所有基因已完成！跳过UniProt/序列查询阶段。")
    else:
        # 4. UniProt映射 + 序列获取
        gene_to_seq, failed_genes = resolve_genes(remaining, cache)
        save_checkpoint(CHECKPOINT_FILE, done_genes, {})

    # 5. 加载序列数据（包含之前已完成的）
    gene_to_seq = {}
    all_failed = []
    for gene in genes:
        entry = cache.get(gene)
        if entry and entry.get("sequence"):
            gene_to_seq[gene] = entry["sequence"]
        else:
            all_failed.append(gene)

    logger.info("可用序列: %d / %d", len(gene_to_seq), total)

    if not gene_to_seq:
        logger.error("无可用序列，终止")
        sys.exit(1)

    # 6. 加载ESM-2模型
    model, tokenizer, hidden_size = load_esm2_model(ESM2_MODEL_NAME)

    # 7. 分批提取特征
    seq_items = [(g, s) for g, s in gene_to_seq.items() if g not in done_genes]
    batches = [
        seq_items[i:i + BATCH_SIZE]
        for i in range(0, len(seq_items), BATCH_SIZE)
    ]
    logger.info("分批推理: %d 批 (batch_size=%d)", len(batches), BATCH_SIZE)

    all_features: dict[str, np.ndarray] = {}

    for batch_idx, batch in enumerate(batches):
        batch_genes = [g for g, s in batch]
        batch_seqs = [s for g, s in batch]
        logger.info("  批 %d/%d (%d 个基因) ...",
                     batch_idx + 1, len(batches), len(batch_genes))

        try:
            batch_features = extract_embeddings_batch(
                batch_seqs, batch_genes, model, tokenizer, MAX_SEQUENCE_LENGTH
            )
            all_features.update(batch_features)

            done_genes.extend(batch_genes)
            save_checkpoint(CHECKPOINT_FILE, done_genes, all_features)
            logger.info("    ✓ 已完成: %d / %d", len(done_genes), total)

        except Exception as e:
            logger.error("    批 %d 失败: %s — 逐个回退...", batch_idx + 1, e)
            for g, s in batch:
                try:
                    from .extract_protein_features_esm2 import extract_esm2_embedding
                except ImportError:
                    pass
                emb = _extract_single_fallback(g, s, model, tokenizer)
                if emb is not None:
                    all_features[g] = emb
                    done_genes.append(g)
                else:
                    all_failed.append(g)
            save_checkpoint(CHECKPOINT_FILE, done_genes, all_features)

    # 8. 组装最终结果
    ordered_genes = [g for g in genes if g in all_features or g in checkpoint.get("done", [])]
    if not ordered_genes:
        logger.error("无任何特征提取成功")
        sys.exit(1)

    # 需要加载之前checkpoint的特征（简化处理：重新计算所有已完成基因的特征）
    # 实际上应该逐步累积，但为简洁，全部重走一遍已完成的
    logger.info("组装最终特征矩阵...")
    final_genes = []
    final_feats = []
    for g in tqdm(genes, desc="组装"):
        if g in all_features:
            final_genes.append(g)
            final_feats.append(all_features[g])
        elif g in gene_to_seq and g not in all_features:
            try:
                emb = _extract_single_fallback(g, gene_to_seq[g], model, tokenizer)
                if emb is not None:
                    final_genes.append(g)
                    final_feats.append(emb)
                    all_features[g] = emb
            except Exception:
                pass

    if not final_genes:
        logger.error("特征矩阵为空")
        sys.exit(1)

    feature_array = np.array(final_feats, dtype=np.float32)
    logger.info("特征矩阵: %s", feature_array.shape)

    # 9. 保存CSV
    dim = feature_array.shape[1]
    columns = ["gene_symbol"] + [f"feat_{i+1:04d}" for i in range(dim)]
    df = pd.DataFrame(feature_array, columns=columns[1:])
    df.insert(0, "gene_symbol", final_genes)
    df.to_csv(OUTPUT_CSV, index=False, encoding="utf-8")
    logger.info("结果保存: %s", OUTPUT_CSV)
    logger.info("  基因: %d | 维度: %d | 大小: %.1f MB",
                len(final_genes), dim, os.path.getsize(OUTPUT_CSV) / 1024 / 1024)

    if all_failed:
        with open("subgraph_failed_genes.txt", "w") as f:
            for g in all_failed:
                f.write(g + "\n")
        logger.warning("失败基因: %d → subgraph_failed_genes.txt", len(all_failed))

    logger.info("全部完成！")


def _extract_single_fallback(
    gene: str, seq: str, model, tokenizer
) -> Optional[np.ndarray]:
    """单个基因的ESM-2推理（回退用）"""
    import torch
    try:
        seq = seq[:MAX_SEQUENCE_LENGTH]
        inputs = tokenizer(seq, return_tensors="pt", truncation=True,
                           max_length=MAX_SEQUENCE_LENGTH + 2)
        with torch.no_grad():
            outputs = model(**inputs)
            emb = outputs.last_hidden_state.squeeze(0).mean(dim=0)
        return emb.cpu().numpy().astype(np.float32)
    except Exception as e:
        logger.error("    %s 推理失败: %s", gene, e)
        return None


if __name__ == "__main__":
    main()