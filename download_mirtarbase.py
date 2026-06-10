# -*- coding: utf-8 -*-
"""
miRTarBase miRNA-靶基因数据 自动化下载 + 处理
===============================================
来源: miRTarBase 2025 (v10.0)
  - https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2025/
  - 380万+ 实验验证的 miRNA-靶基因互作

处理流程:
  1. 尝试自动下载 miRTarBase hsa_MTI.xlsx
  2. 若自动下载失败, 提示手动下载路径
  3. 解析并提取 miRNA-基因关联边

输出: D:\反向网络药理学\GAT拓展维度\gene_mirna_edges.txt
"""

import os
import sys
import csv
import json
import urllib.request
from pathlib import Path
from collections import Counter, defaultdict

# ============================================================
# 0. 配置
# ============================================================
OUTPUT_DIR = Path(r"D:\反向网络药理学\GAT拓展维度")
DATA_DIR = Path(r"C:\Users\Jy-Mentor-7\Desktop\GAT")

# miRTarBase 2025 下载URL (多个备选)
MIRTARBASE_URLS = [
    "https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2025/downloads/hsa_MTI.xlsx",
    "https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2025/MTI.xls",
    "https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2025/downloads/MTI.xls",
]

MIRTARBASE_LOCAL = DATA_DIR / "hsa_MTI.xlsx"

OUTPUT_EDGES = OUTPUT_DIR / "gene_mirna_edges.txt"
OUTPUT_STATS = OUTPUT_DIR / "gene_mirna_stats.json"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# 实验证据类型优先级
EVIDENCE_PRIORITY = {
    'reporter assay': 5, 'western blot': 5, 'qPCR': 4,
    'luciferase reporter assay': 5, 'dual-luciferase reporter assay': 5,
    'microarray': 3, 'pSILAC': 3,
    'clip-seq': 3, 'hits-clip': 3, 'par-clip': 3, 'iclip': 3,
    'clash': 3, 'degradome': 3, 'ngs': 2, 'rna-seq': 2,
    'other': 1,
}


def download_mirtarbase():
    """自动下载 miRTarBase 数据"""
    if MIRTARBASE_LOCAL.exists():
        print(f"[INFO] 已存在: {MIRTARBASE_LOCAL.name} ({MIRTARBASE_LOCAL.stat().st_size/1024/1024:.1f} MB)")
        return str(MIRTARBASE_LOCAL)

    # 尝试备选文件名
    alt_files = [
        DATA_DIR / "miRTarBase_MTI.xls",
        DATA_DIR / "miRTarBase_MTI.csv",
    ]
    for af in alt_files:
        if af.exists():
            print(f"[INFO] 已存在备选文件: {af.name}")
            return str(af)

    # 尝试自动下载
    for url in MIRTARBASE_URLS:
        print(f"[TRY] 尝试下载: {url}")
        try:
            req = urllib.request.Request(url, headers={
                'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
            })
            with urllib.request.urlopen(req, timeout=30) as response:
                total_size = int(response.headers.get('Content-Length', 0))
                print(f"  文件大小: {total_size/1024/1024:.1f} MB")

                with open(MIRTARBASE_LOCAL, 'wb') as f:
                    downloaded = 0
                    while True:
                        chunk = response.read(8192)
                        if not chunk:
                            break
                        f.write(chunk)
                        downloaded += len(chunk)
                        if total_size > 0:
                            pct = downloaded / total_size * 100
                            print(f"\r  进度: {downloaded/1024/1024:.1f}/{total_size/1024/1024:.1f} MB ({pct:.1f}%)", end='')

                print(f"\n[OK] 下载完成: {MIRTARBASE_LOCAL.name}")
                return str(MIRTARBASE_LOCAL)

        except Exception as e:
            print(f"  [FAIL] {e}")
            continue

    print("\n[INFO] 自动下载失败, 请手动下载:")
    print(f"  1. 打开: https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2025/")
    print(f"  2. 点击 Download → Homo sapiens MTI")
    print(f"  3. 保存为: {MIRTARBASE_LOCAL}")
    print(f"  4. 重新运行本脚本")
    return None


def parse_mirtarbase(filepath):
    """解析 miRTarBase 文件, 提取 miRNA-靶基因关联"""
    import pandas as pd

    ext = Path(filepath).suffix.lower()
    print(f"[INFO] 文件格式: {ext}")

    if ext in ('.xls', '.xlsx'):
        df = pd.read_excel(filepath)
    elif ext == '.csv':
        df = pd.read_csv(filepath, sep=None, engine='python', encoding='utf-8', errors='replace')
    else:
        df = pd.read_csv(filepath, sep=None, engine='python', encoding='utf-8', errors='replace')

    print(f"[INFO] 数据维度: {df.shape[0]} 行 × {df.shape[1]} 列")
    print(f"[INFO] 列名: {df.columns.tolist()[:15]}...")

    cols_lower = {c.lower().strip(): c for c in df.columns}

    mirna_col = None
    gene_col = None
    evidence_col = None
    species_col = None

    for key, col in cols_lower.items():
        kl = key.lower()
        if 'mirna' in kl and mirna_col is None:
            mirna_col = col
        elif ('target' in kl and 'gene' in kl) and gene_col is None:
            gene_col = col
        elif 'gene' in kl and ('symbol' in kl or 'name' in kl) and gene_col is None:
            gene_col = col
        elif ('experiment' in kl or 'method' in kl or 'evidence' in kl) and evidence_col is None:
            evidence_col = col
        elif 'species' in kl and species_col is None:
            species_col = col

    if mirna_col is None:
        mirna_col = df.columns[0]
    if gene_col is None:
        gene_col = df.columns[2] if len(df.columns) > 2 else df.columns[1]

    print(f"[INFO] miRNA列: '{mirna_col}'")
    print(f"[INFO] 基因列:   '{gene_col}'")
    print(f"[INFO] 证据列:   '{evidence_col}'")
    print(f"[INFO] 物种列:   '{species_col}'")

    edges = []
    stats = {
        'total_rows': len(df),
        'valid_edges': 0,
        'skipped_non_human': 0,
        'skipped_empty': 0,
        'unique_mirnas': set(),
        'unique_genes': set(),
        'mirna_degree': Counter(),
        'gene_degree': Counter(),
        'evidence_types': Counter(),
    }

    for idx, row in df.iterrows():
        try:
            mirna = str(row[mirna_col]).strip() if mirna_col in df.columns else ''
            gene = str(row[gene_col]).strip().upper() if gene_col in df.columns else ''
        except (ValueError, KeyError):
            stats['skipped_empty'] += 1
            continue

        if not mirna or not gene or gene in ('NA', 'NAN', 'NONE', '', 'N/A'):
            stats['skipped_empty'] += 1
            continue

        # 物种过滤
        if species_col and species_col in df.columns:
            species = str(row[species_col]).lower()
            if 'human' not in species and 'homo' not in species:
                stats['skipped_non_human'] += 1
                continue

        evidence = ''
        if evidence_col and evidence_col in df.columns:
            evidence = str(row[evidence_col]).strip()
            stats['evidence_types'][evidence] += 1

        gene = gene.split('///')[0].split(';')[0].split('/')[0].strip()

        edges.append((gene, mirna, evidence))
        stats['valid_edges'] += 1
        stats['unique_mirnas'].add(mirna)
        stats['unique_genes'].add(gene)
        stats['mirna_degree'][mirna] += 1
        stats['gene_degree'][gene] += 1

        if idx % 500000 == 0 and idx > 0:
            print(f"  已处理 {idx} 行, {stats['valid_edges']} 有效边")

    print(f"\n[STATS] miRTarBase 解析完成:")
    print(f"  总行数:         {stats['total_rows']:,}")
    print(f"  有效边数:       {stats['valid_edges']:,}")
    print(f"  跳过(非人类):   {stats['skipped_non_human']}")
    print(f"  跳过(空值):     {stats['skipped_empty']}")
    print(f"  唯一 miRNA:     {len(stats['unique_mirnas']):,}")
    print(f"  唯一靶基因:     {len(stats['unique_genes']):,}")
    print(f"  Top-5 miRNA:    {stats['mirna_degree'].most_common(5)}")
    print(f"  Top-5 基因:     {stats['gene_degree'].most_common(5)}")

    return edges, stats


def save_edges(edges, output_path):
    """保存 miRNA-靶基因关联边, 去重取最高优先级证据"""
    best_evidence = {}
    for gene, mirna, evidence in edges:
        key = (gene, mirna)
        current_best = best_evidence.get(key, ('Other', 0))
        new_priority = EVIDENCE_PRIORITY.get(evidence.lower(), 1)
        if new_priority > current_best[1]:
            best_evidence[key] = (evidence, new_priority)

    with open(output_path, 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        writer.writerow(['gene_symbol', 'mirna_id', 'evidence_type'])
        for (gene, mirna), (evidence, _) in sorted(best_evidence.items()):
            writer.writerow([gene, mirna, evidence])

    print(f"[SAVE] 边文件: {output_path} ({len(best_evidence):,} 条去重边)")
    return len(best_evidence)


def save_stats(stats, dedup_count, output_path):
    """保存统计信息"""
    stats_out = {
        'total_edges_raw': stats['valid_edges'],
        'total_edges_dedup': dedup_count,
        'unique_mirnas': len(stats['unique_mirnas']),
        'unique_genes': len(stats['unique_genes']),
        'top_mirnas': stats['mirna_degree'].most_common(20),
        'top_targets': stats['gene_degree'].most_common(20),
        'source': 'miRTarBase 2025 (v10.0)',
    }
    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(stats_out, f, indent=2, ensure_ascii=False)
    print(f"[SAVE] 统计文件: {output_path}")


# ============================================================
# MAIN
# ============================================================
if __name__ == "__main__":
    print("=" * 60)
    print("miRTarBase miRNA-靶基因数据 下载 + 处理")
    print("=" * 60)

    # Step 1: 下载
    print("\n--- Step 1: 下载数据 ---")
    filepath = download_mirtarbase()
    if filepath is None:
        print("\n[INFO] 下载失败, 退出")
        sys.exit(1)

    # Step 2: 解析
    print("\n--- Step 2: 解析数据 ---")
    edges, stats = parse_mirtarbase(filepath)

    if len(edges) == 0:
        print("[ERROR] 未提取到任何有效边")
        sys.exit(1)

    # Step 3: 保存
    print("\n--- Step 3: 保存结果 ---")
    dedup_count = save_edges(edges, OUTPUT_EDGES)
    save_stats(stats, dedup_count, OUTPUT_STATS)

    print(f"\n{'='*60}")
    print("miRNA-靶基因关联边提取完成!")
    print(f"输出: {OUTPUT_EDGES}")
    print("=" * 60)