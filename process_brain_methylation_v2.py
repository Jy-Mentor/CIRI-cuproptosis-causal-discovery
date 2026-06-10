# -*- coding: utf-8 -*-
"""
脑组织甲基化数据 流式处理 (v2)
================================
来源: brain_methylation_v1.txt (4.89 GB, 1997 samples x ~482K CpG probes)
输出: gene_methylation_edges.txt (基因-甲基化边)

处理策略:
  1. 下载 Illumina 450K/EPIC 探针注释 (CpG→基因)
  2. 逐行流式读取 TXT, 不加载全部到内存
  3. 计算每个 CpG 的平均 beta, 筛选高/低甲基化
  4. 映射到基因并保存边
"""

import os
import sys
import csv
import json
import time
from pathlib import Path
from collections import Counter, defaultdict

OUTPUT_DIR = Path(r"D:\反向网络药理学\GAT拓展维度")
DATA_DIR = Path(r"C:\Users\Jy-Mentor-7\Desktop\GAT")
METH_TXT = DATA_DIR / "brain_methylation_temp" / "brain_methylation" / "brain_methylation_v1.txt"
CPG_MAP_CSV = DATA_DIR / "cpg_gene_map.csv"
OUTPUT_EDGES = OUTPUT_DIR / "gene_methylation_edges.txt"
OUTPUT_STATS = OUTPUT_DIR / "gene_methylation_stats.json"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def load_cpg_gene_map(map_path):
    """从 CSV 加载 CpG→基因映射 (cpg_id,gene)"""
    cpg_gene_map = {}
    print(f"[INFO] 加载 CpG→基因映射: {map_path}")
    if not map_path.exists():
        print(f"[ERROR] 映射文件不存在: {map_path}")
        return cpg_gene_map

    with open(map_path, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        header = next(reader, None)
        if not header or header[0] != 'cpg_id':
            print(f"[WARN] 表头异常: {header}")

        for row in reader:
            if len(row) >= 2:
                cpg = row[0].strip()
                gene = row[1].strip()
                if cpg.startswith('cg') and gene:
                    cpg_gene_map[cpg] = gene

    print(f"[INFO] CpG→基因映射: {len(cpg_gene_map)} 条")
    return cpg_gene_map


def process_methylation_stream(cpg_gene_map):
    """流式处理甲基化 TXT 文件"""
    if not METH_TXT.exists():
        print(f"[ERROR] 未找到文件: {METH_TXT}")
        return [], Counter()

    print(f"\n[PROCESS] 流式读取: {METH_TXT}")
    file_size_gb = METH_TXT.stat().st_size / 1024 / 1024 / 1024
    print(f"  文件大小: {file_size_gb:.2f} GB")

    edges = []
    gene_degree = Counter()
    total_cpg = 0
    mapped_cpg = 0
    unmapped_cpg = 0
    start_time = time.time()
    last_report = time.time()

    with open(METH_TXT, 'r', encoding='utf-8', errors='replace') as f:
        header_line = f.readline().strip()
        tissue_line = f.readline().strip()
        sep = '\t'
        n_samples = len(header_line.split(sep)) - 1
        print(f"  样本数: {n_samples}")

        for line_no, line in enumerate(f, start=3):
            line = line.strip()
            if not line:
                continue

            parts = line.split(sep)
            if len(parts) < 2:
                continue

            cpg = parts[0].strip()
            if not cpg.startswith('cg'):
                continue

            total_cpg += 1

            gene = cpg_gene_map.get(cpg, '')
            if not gene:
                unmapped_cpg += 1
                if total_cpg % 50000 == 0:
                    now = time.time()
                    elapsed = now - start_time
                    speed = total_cpg / elapsed if elapsed > 0 else 0
                    eta = (482421 - total_cpg) / speed if speed > 0 else 0
                    print(f"  已处理 {total_cpg:,} CpG | 已映射 {mapped_cpg:,} | 速度 {speed:.0f} CpG/s | ETA {eta:.0f}s", end='\r')
                continue

            betas = []
            for val in parts[1:]:
                try:
                    b = float(val.strip().strip('"'))
                    if 0 <= b <= 1:
                        betas.append(b)
                except (ValueError, IndexError):
                    continue

            if not betas:
                continue

            avg_beta = sum(betas) / len(betas)

            if avg_beta > 0.7 or avg_beta < 0.3:
                edges.append((gene, cpg, f"{avg_beta:.4f}", "hypermethylated" if avg_beta > 0.7 else "hypomethylated"))
                gene_degree[gene] += 1
                mapped_cpg += 1

            if total_cpg % 50000 == 0:
                now = time.time()
                elapsed = now - start_time
                speed = total_cpg / elapsed if elapsed > 0 else 0
                eta_remaining = (482421 - total_cpg) / speed if speed > 0 else 0
                print(f"  已处理 {total_cpg:,} CpG | 已映射 {mapped_cpg:,} | 速度 {speed:.0f} CpG/s | ETA {eta_remaining:.0f}s", end='\r')

    elapsed = time.time() - start_time
    print(f"\n[STATS] 处理完成! 耗时 {elapsed:.1f}s")
    print(f"  总 CpG: {total_cpg:,}")
    print(f"  映射到基因: {mapped_cpg:,}")
    print(f"  未映射: {unmapped_cpg:,}")
    print(f"  总边数: {len(edges):,}")
    print(f"  唯一基因: {len(gene_degree):,}")

    if gene_degree:
        print(f"  Top-10 甲基化基因: {gene_degree.most_common(10)}")

    return edges, gene_degree, total_cpg, mapped_cpg, unmapped_cpg


def save_edges(edges, output_path):
    """保存基因-甲基化关联边 (去重)"""
    seen = set()
    unique_edges = []
    for e in edges:
        key = (e[0], e[1])
        if key not in seen:
            seen.add(key)
            unique_edges.append(e)

    with open(output_path, 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        writer.writerow(['gene', 'cpg_id', 'avg_beta', 'status'])
        for gene, cpg, beta, status in sorted(unique_edges):
            writer.writerow([gene, cpg, beta, status])

    print(f"[SAVE] 边文件: {output_path} ({len(unique_edges):,} 条)")
    return len(unique_edges)


def save_stats(gene_degree, total_cpg, mapped_cpg, unmapped_cpg, n_edges, output_path):
    """保存统计信息"""
    stats = {
        'total_cpg': total_cpg,
        'genes_mapped': mapped_cpg,
        'cpg_unmapped': unmapped_cpg,
        'total_edges': n_edges,
        'unique_genes': len(gene_degree),
        'unique_cpgs': 0,
        'n_samples': 1997,
        'brain_regions': '25 brain regions',
        'top_methylated_genes': [(k, v) for k, v in gene_degree.most_common(20)],
        'source': 'EWAS Data Hub - Brain Methylation v1',
        'threshold': 'beta > 0.7 (hypermethylated) or beta < 0.3 (hypomethylated)',
    }
    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(stats, f, indent=2, ensure_ascii=False)
    print(f"[SAVE] 统计文件: {output_path}")


if __name__ == "__main__":
    print("=" * 60)
    print("脑组织甲基化数据 流式处理 v2")
    print("=" * 60)

    # Step 1: 加载 CpG→基因映射
    print("\n--- Step 1: 加载 CpG→基因映射 ---")
    cpg_gene_map = load_cpg_gene_map(CPG_MAP_CSV)

    if not cpg_gene_map:
        print("[ERROR] 无法获取 CpG→基因映射")
        sys.exit(1)

    # Step 2: 流式处理
    print("\n--- Step 2: 流式处理甲基化数据 ---")
    edges, gene_degree, total_cpg, mapped_cpg, unmapped_cpg = process_methylation_stream(cpg_gene_map)

    if not edges:
        print("[WARN] 未提取到甲基化边, 输出空文件")
    else:
        print(f"\n[OK] 提取到 {len(edges):,} 条甲基化边")

    # Step 3: 保存结果
    print("\n--- Step 3: 保存结果 ---")
    n_edges = save_edges(edges, OUTPUT_EDGES)
    save_stats(gene_degree, total_cpg, mapped_cpg, unmapped_cpg, n_edges, OUTPUT_STATS)

    print(f"\n{'=' * 60}")
    print("脑组织甲基化数据处理完成!")
    print(f"输出边文件: {OUTPUT_EDGES}")
    print(f"输出统计:  {OUTPUT_STATS}")
    print("=" * 60)