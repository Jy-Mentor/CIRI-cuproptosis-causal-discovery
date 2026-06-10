# -*- coding: utf-8 -*-
"""
ARCHS4 基因共表达边提取脚本 v2.0
=================================
通过 ARCHS4 API 查询基因-基因共表达相关性，生成异构图所需的共表达边。

API: POST https://maayanlab.cloud/matrixapi/coltop
     Body: {"id": "GENE", "count": 200}
Fallback: POST https://maayanlab.cloud/sigpy/data/correlation

输入:  C:/Users/Jy-Mentor-7/Desktop/GAT/subgraph_genes.txt (~15,648 genes)
输出:  D:/反向网络药理学/GAT拓展维度/gene_coexp_edges.txt

特性:
  - 自动会话管理 (连接复用 + 重试)
  - 断点续传 (Resume)
  - 双API回退 (matrixapi → sigpy)
  - 进度条 + ETA 预估
  - 增量写入 (逐基因追加)
"""

import os
import sys
import time
import json
import random
import argparse
from pathlib import Path

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


ARCHS4_MATRIX_API = "https://maayanlab.cloud/matrixapi/coltop"
ARCHS4_SIGPY_API = "https://maayanlab.cloud/sigpy/data/correlation"

CORRELATION_THRESHOLD = 0.7
TOP_N = 200
GENE_LIST_FILE = r"C:\Users\Jy-Mentor-7\Desktop\GAT\subgraph_genes.txt"
OUTPUT_FILE = r"D:\反向网络药理学\GAT拓展维度\gene_coexp_edges.txt"

REQUEST_TIMEOUT = 60
MAX_RETRIES = 3
RETRY_BACKOFF = 2.0
DELAY_MIN = 0.3
DELAY_MAX = 0.5


def create_session():
    session = requests.Session()
    retry_strategy = Retry(
        total=MAX_RETRIES,
        backoff_factor=RETRY_BACKOFF,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["POST"],
    )
    adapter = HTTPAdapter(max_retries=retry_strategy, pool_connections=10, pool_maxsize=20)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session


def query_matrix_api(session, gene, count):
    payload = {"id": gene, "count": count}
    resp = session.post(ARCHS4_MATRIX_API, json=payload, timeout=REQUEST_TIMEOUT)
    resp.raise_for_status()
    data = resp.json()
    rowids = data.get("rowids", [])
    values = data.get("values", [])
    if not rowids:
        return []
    results = []
    for i in range(len(rowids)):
        if rowids[i].upper() == gene.upper():
            continue
        results.append((rowids[i].upper(), values[i]))
    return results


def query_sigpy_api(session, gene):
    payload = {"gene": gene, "species": "human", "meta": ""}
    resp = session.post(ARCHS4_SIGPY_API, json=payload, timeout=REQUEST_TIMEOUT)
    resp.raise_for_status()
    data = resp.json()
    results = []
    for entry in data.get("positive_correlated_genes", []):
        results.append((entry["gene"].upper(), entry["correlation"]))
    for entry in data.get("negative_correlated_genes", []):
        results.append((entry["gene"].upper(), entry["correlation"]))
    return results


def query_gene(session, gene, top_n):
    last_exc = None
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            results = query_matrix_api(session, gene, top_n)
            return results, None
        except requests.exceptions.HTTPError as e:
            last_exc = e
            if e.response is not None and e.response.status_code == 404:
                return [], None
            if attempt < MAX_RETRIES:
                time.sleep(RETRY_BACKOFF ** attempt)
        except Exception as e:
            last_exc = e
            if attempt < MAX_RETRIES:
                time.sleep(RETRY_BACKOFF ** attempt)

    try:
        results = query_sigpy_api(session, gene)
        return results, None
    except Exception as e:
        return [], str(e)


def load_genes(filepath):
    genes = set()
    with open(filepath, "r") as f:
        first = True
        for line in f:
            if first:
                first = False
                stripped = line.strip().upper()
                if stripped in ("GENE", "GENES", "SYMBOL", "GENE_SYMBOL"):
                    continue
            gene = line.strip().upper()
            if gene:
                genes.add(gene)
    return sorted(genes)


def load_existing_pairs(output_path):
    pairs = set()
    last_gene = None
    if not output_path.exists():
        return pairs, last_gene
    with open(output_path, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                g1, g2 = parts[0].upper(), parts[1].upper()
                pairs.add((g1, g2))
                last_gene = g1
    return pairs, last_gene


def find_resume_idx(gene_list, last_gene):
    if not last_gene:
        return 0
    for i, g in enumerate(gene_list):
        if g >= last_gene:
            return i
    return 0


def main():
    parser = argparse.ArgumentParser(description="ARCHS4 共表达边查询")
    parser.add_argument("-i", "--input", default=GENE_LIST_FILE)
    parser.add_argument("-o", "--output", default=OUTPUT_FILE)
    parser.add_argument("-t", "--threshold", type=float, default=CORRELATION_THRESHOLD)
    parser.add_argument("-n", "--top-n", type=int, default=TOP_N)
    parser.add_argument("--delay-min", type=float, default=DELAY_MIN)
    parser.add_argument("--delay-max", type=float, default=DELAY_MAX)
    parser.add_argument("--resume", action="store_true", help="从已有输出断点续传")
    args = parser.parse_args()

    top_n = args.top_n

    input_path = Path(args.input)
    output_path = Path(args.output)
    threshold = args.threshold

    if not input_path.exists():
        print(f"[ERROR] 输入文件不存在: {input_path}", file=sys.stderr)
        sys.exit(1)

    gene_list = load_genes(input_path)
    gene_set = set(gene_list)
    total = len(gene_list)
    print(f"[INFO] 加载 {total} 个基因")

    output_path.parent.mkdir(parents=True, exist_ok=True)

    written_pairs = set()
    total_edges = 0
    errors = 0
    not_found = 0
    start_idx = 0

    if args.resume:
        written_pairs, last_gene = load_existing_pairs(output_path)
        total_edges = len(written_pairs)
        start_idx = find_resume_idx(gene_list, last_gene)
        print(f"[RESUME] 已有 {total_edges} 条边, 从基因 #{start_idx} ({gene_list[start_idx]}) 续传")
    else:
        with open(output_path, "w") as f:
            f.write(f"# gene1\tgene2\tcorrelation\n")
            f.write(f"# Source: ARCHS4 (https://maayanlab.cloud/archs4/)\n")
            f.write(f"# Threshold: |correlation| > {threshold}\n")
            f.write(f"# Input genes: {total}\n")
            f.write(f"# Top-N per gene: {TOP_N}\n")
            f.write(f"# Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")

    session = create_session()
    start_time = time.time()
    last_report = 0

    for idx in range(start_idx, total):
        gene = gene_list[idx]
        now = time.time()
        if now - last_report >= 10 or idx == start_idx:
            elapsed = now - start_time
            done = idx - start_idx
            rate = done / elapsed if elapsed > 0 else 0
            eta = (total - idx) / rate if rate > 0 else 0
            print(
                f"[{idx + 1}/{total}] {100 * (idx + 1) / total:.1f}% | "
                f"边:{total_edges} 错:{errors} 缺:{not_found} | "
                f"ETA:{eta / 60:.1f}min"
            )
            sys.stdout.flush()
            last_report = now

        results, error = query_gene(session, gene, top_n)

        if error:
            errors += 1
            if errors <= 5:
                print(f"  [WARN] {gene}: {error}")
        elif not results:
            not_found += 1
        else:
            new_edges = []
            for target, corr in results:
                if target not in gene_set:
                    continue
                if abs(corr) <= threshold:
                    continue
                g1, g2 = (gene, target) if gene < target else (target, gene)
                pair = (g1, g2)
                if pair in written_pairs:
                    continue
                written_pairs.add(pair)
                new_edges.append((g1, g2, corr))

            if new_edges:
                with open(output_path, "a") as f:
                    for g1, g2, corr in new_edges:
                        f.write(f"{g1}\t{g2}\t{corr:.6f}\n")
                total_edges += len(new_edges)

        if idx < total - 1:
            time.sleep(random.uniform(args.delay_min, args.delay_max))

    elapsed = time.time() - start_time
    print(f"\n[DONE] 耗时 {elapsed / 60:.1f} 分钟")
    print(f"[RESULT] 共表达边 (|corr|>{threshold}): {total_edges}")
    print(f"[RESULT] 查询基因: {total}")
    print(f"[RESULT] 错误: {errors}, 未找到: {not_found}")
    print(f"[RESULT] 输出: {output_path}")


if __name__ == "__main__":
    main()