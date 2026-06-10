#!/usr/bin/env python3
"""
批量获取基因注释特征 (GO / KEGG / InterPro)
输出到 local_data/go_terms.tsv, kegg_pathways.tsv, interpro_domains.tsv
数据来源: MyGene.info API (免费, 无需注册)
"""

import os
import sys
import json
import time
import logging
from pathlib import Path
from typing import List, Set, Dict, Tuple
from collections import defaultdict

import requests
import pandas as pd

# ---------------------------------------------------------------------------
# 配置
# ---------------------------------------------------------------------------
BASE_DIR = Path(__file__).parent.resolve()
LOCAL_DATA_DIR = BASE_DIR / "local_data"
PROCESSED_DIR = BASE_DIR / "processed"
LOGS_DIR = BASE_DIR / "logs"

for d in [LOCAL_DATA_DIR, LOGS_DIR]:
    d.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)

MG_URL = "https://mygene.info/v3/query"
BATCH_SIZE = 1000          # MyGene.info 建议每批≤1000
SLEEP_SEC = 0.5            # 礼貌延迟，避免429
FIELDS = "go,pathway.kegg,interpro"


# ---------------------------------------------------------------------------
# 读取基因池
# ---------------------------------------------------------------------------
def load_gene_pool() -> List[str]:
    """从 processed/node_features.csv 读取基因池"""
    nf_path = PROCESSED_DIR / "node_features.csv"
    if not nf_path.exists():
        logger.error(f"找不到基因池文件: {nf_path}，请先运行 --mode data")
        sys.exit(1)
    df = pd.read_csv(nf_path)
    genes = df["GeneSymbol"].dropna().str.strip().str.upper().unique().tolist()
    logger.info(f"基因池大小: {len(genes)}")
    return genes


# ---------------------------------------------------------------------------
# API 查询
# ---------------------------------------------------------------------------
def fetch_batch(genes: List[str]) -> List[dict]:
    """向 MyGene.info 发送一批基因查询"""
    payload = {
        "q": genes,
        "scopes": "symbol",
        "species": "human",
        "fields": FIELDS,
    }
    try:
        r = requests.post(
            MG_URL,
            json=payload,
            headers={"Content-Type": "application/json"},
            timeout=60,
        )
        r.raise_for_status()
        return r.json()
    except requests.exceptions.RequestException as e:
        logger.error(f"API请求失败: {e}")
        return []


def parse_go_terms(hit: dict) -> Set[str]:
    """从MyGene返回解析GO ID集合"""
    terms = set()
    go_data = hit.get("go", {})
    for aspect in ["BP", "CC", "MF"]:
        if aspect in go_data and isinstance(go_data[aspect], list):
            for entry in go_data[aspect]:
                if isinstance(entry, dict) and "id" in entry:
                    terms.add(entry["id"])
    return terms


def parse_kegg_pathways(hit: dict) -> Set[str]:
    """解析KEGG pathway ID集合"""
    terms = set()
    pw = hit.get("pathway", {})
    if isinstance(pw, dict) and "kegg" in pw:
        kegg_list = pw["kegg"]
        if isinstance(kegg_list, list):
            for entry in kegg_list:
                if isinstance(entry, dict):
                    pid = entry.get("id") or entry.get("name")
                    if pid:
                        terms.add(str(pid))
    return terms


def parse_interpro(hit: dict) -> Set[str]:
    """解析InterPro结构域ID集合"""
    terms = set()
    ipr = hit.get("interpro", {})
    if isinstance(ipr, list):
        for entry in ipr:
            if isinstance(entry, dict):
                iid = entry.get("id") or entry.get("accession")
                if iid:
                    terms.add(str(iid))
    return terms


# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------
def main():
    gene_pool = load_gene_pool()
    n_genes = len(gene_pool)
    n_batches = (n_genes + BATCH_SIZE - 1) // BATCH_SIZE

    go_records: List[Tuple[str, str]] = []
    kegg_records: List[Tuple[str, str]] = []
    interpro_records: List[Tuple[str, str]] = []
    missing_genes: List[str] = []

    for i in range(n_batches):
        batch = gene_pool[i * BATCH_SIZE : (i + 1) * BATCH_SIZE]
        logger.info(f"查询批次 {i+1}/{n_batches} ({len(batch)}个基因)...")
        hits = fetch_batch(batch)

        if not hits:
            missing_genes.extend(batch)
            continue

        for gene, hit in zip(batch, hits):
            if not isinstance(hit, dict) or ("notfound" in hit and hit["notfound"]):
                missing_genes.append(gene)
                continue

            # GO
            go_set = parse_go_terms(hit)
            if go_set:
                for t in go_set:
                    go_records.append((gene, t))
            else:
                go_records.append((gene, "NA"))

            # KEGG
            kegg_set = parse_kegg_pathways(hit)
            if kegg_set:
                for t in kegg_set:
                    kegg_records.append((gene, t))
            else:
                kegg_records.append((gene, "NA"))

            # InterPro
            ipr_set = parse_interpro(hit)
            if ipr_set:
                for t in ipr_set:
                    interpro_records.append((gene, t))
            else:
                interpro_records.append((gene, "NA"))

        if i < n_batches - 1:
            time.sleep(SLEEP_SEC)

    # -----------------------------------------------------------------------
    # 保存为 TSV
    # -----------------------------------------------------------------------
    go_df = pd.DataFrame(go_records, columns=["GeneSymbol", "GO_term"])
    go_df = go_df[go_df["GO_term"] != "NA"]
    go_path = LOCAL_DATA_DIR / "go_terms.tsv"
    go_df.to_csv(go_path, sep="\t", index=False)
    logger.info(f"GO terms: {go_df['GeneSymbol'].nunique()}/{n_genes} 个基因有注释, 共 {len(go_df)} 条记录 -> {go_path}")

    kegg_df = pd.DataFrame(kegg_records, columns=["GeneSymbol", "Pathway"])
    kegg_df = kegg_df[kegg_df["Pathway"] != "NA"]
    kegg_path = LOCAL_DATA_DIR / "kegg_pathways.tsv"
    kegg_df.to_csv(kegg_path, sep="\t", index=False)
    logger.info(f"KEGG pathways: {kegg_df['GeneSymbol'].nunique()}/{n_genes} 个基因有注释, 共 {len(kegg_df)} 条记录 -> {kegg_path}")

    ipr_df = pd.DataFrame(interpro_records, columns=["GeneSymbol", "Domain"])
    ipr_df = ipr_df[ipr_df["Domain"] != "NA"]
    ipr_path = LOCAL_DATA_DIR / "interpro_domains.tsv"
    ipr_df.to_csv(ipr_path, sep="\t", index=False)
    logger.info(f"InterPro domains: {ipr_df['GeneSymbol'].nunique()}/{n_genes} 个基因有注释, 共 {len(ipr_df)} 条记录 -> {ipr_path}")

    if missing_genes:
        missing_path = LOGS_DIR / "annotation_missing_genes.log"
        with open(missing_path, "w", encoding="utf-8") as f:
            for g in missing_genes:
                f.write(g + "\n")
        logger.warning(f"无API返回的基因: {len(missing_genes)} 个，已记录到 {missing_path}")

    logger.info("=" * 60)
    logger.info("三个特征文件已就绪，可直接运行主脚本的 --mode data")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
