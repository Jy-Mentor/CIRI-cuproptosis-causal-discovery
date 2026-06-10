# -*- coding: utf-8 -*-
"""
从 EWAS Atlas 提取人类基因-甲基化关联 (可选边类型)
=====================================================
来源: EWAS Atlas (https://www.ewasatlas.org/)
  - EWAS Atlas 收录了表观基因组关联研究中甲基化位点-表型/基因关联
  - 下载人类数据集, 筛选 CpG 位点与其调控基因

输入: EWAS_Atlas_human.txt (下载自 https://www.ewasatlas.org/)
输出: gene_methylation_edges.txt (gene_symbol  cpg_id)

备选数据源:
  1. GREAT (http://great.stanford.edu/) - 甲基化峰→基因注释
  2. REDIportal (http://srv00.recas.ba.infn.it/atlas/) - RNA编辑
  3. 自定义甲基化芯片数据 + GREAT注释 (如需处理.dir文件)

作者: 优化版 v2.0
日期: 2026-05-31
"""

import os
import sys
import csv
import json
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

# ============================================================
# 0. 配置
# ============================================================
OUTPUT_DIR = Path(r"D:\反向网络药理学\GAT拓展维度")
DATA_DIR = Path(r"C:\Users\Jy-Mentor-7\Desktop\GAT")

EWAS_ATLAS_URL = "https://www.ewasatlas.org/download"
EWAS_ATLAS_LOCAL = DATA_DIR / "EWAS_Atlas_human.txt"

OUTPUT_EDGES = OUTPUT_DIR / "gene_methylation_edges.txt"
OUTPUT_STATS = OUTPUT_DIR / "gene_methylation_stats.json"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def download_ewas_atlas():
    """下载 EWAS Atlas 数据"""
    if EWAS_ATLAS_LOCAL.exists():
        print(f"[INFO] 已存在本地文件: {EWAS_ATLAS_LOCAL}")
        return str(EWAS_ATLAS_LOCAL)

    print(f"[INFO] EWAS Atlas 需要手动下载")
    print(f"  URL: {EWAS_ATLAS_URL}")
    print(f"  保存到: {EWAS_ATLAS_LOCAL}")
    print(f"\n  步骤:")
    print(f"  1. 访问 https://www.ewasatlas.org/")
    print(f"  2. 点击 'Download' → 选择 'Human' 数据集")
    print(f"  3. 将下载的文件保存到上述路径")
    print(f"  4. 重新运行本脚本")
    return None


def detect_format(filepath):
    """自动检测 EWAS Atlas 文件格式并返回表头映射"""
    with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
        first_line = f.readline().strip()
        sep = '\t' if '\t' in first_line else ','

    print(f"[INFO] 分隔符: {repr(sep)}")

    if first_line.lower().startswith('cpg') or 'cpg' in first_line.lower():
        with open(filepath, 'r', encoding='utf-8') as f:
            header_line = f.readline().strip()
            headers = header_line.split(sep)
            print(f"[INFO] 检测到表头: {headers[:10]}...")

            col_map = {}
            for i, h in enumerate(headers):
                h_lower = h.lower().strip().strip('"').strip("'")
                if 'cpg' in h_lower or 'probe' in h_lower:
                    col_map['cpg'] = i
                elif 'gene' in h_lower or 'symbol' in h_lower or 'target' in h_lower:
                    col_map['gene'] = i
                elif 'p' in h_lower and ('value' in h_lower or 'val' in h_lower):
                    col_map['pvalue'] = i
                elif 'trait' in h_lower or 'phenotype' in h_lower or 'disease' in h_lower:
                    col_map['trait'] = i

            return sep, col_map.to_header if hasattr(col_map, 'to_header') else col_map
    else:
        return sep, {}


def parse_ewas_atlas(filepath):
    """
    解析 EWAS Atlas 文件, 提取基因-甲基化关联

    支持多种格式:
      格式A: CpG  Gene  P-value  Trait  ...
      格式B: probe  gene_symbol  beta  p  ...
      格式C: cpg_id  gene  direction  ...

    输出: gene_symbol  cpg_id 两列
    """
    sep, col_map = detect_format(filepath)

    edges = []
    stats = {
        'total_lines': 0,
        'valid_edges': 0,
        'skipped_no_gene': 0,
        'skipped_no_cpg': 0,
        'unique_genes': set(),
        'unique_cpgs': set(),
        'gene_degree': Counter(),
        'cpg_traits': {},
    }

    with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
        reader = csv.reader(f, delimiter=sep)
        header = next(reader, None)

        for row_num, row in enumerate(reader, start=2):
            stats['total_lines'] += 1

            if col_map:
                gene_idx = col_map.get('gene')
                cpg_idx = col_map.get('cpg')
                trait_idx = col_map.get('trait', -1)
            else:
                gene_idx = 1
                cpg_idx = 0
                trait_idx = -1

            try:
                gene = row[gene_idx].strip().upper() if gene_idx is not None and gene_idx < len(row) else ''
                cpg = row[cpg_idx].strip() if cpg_idx is not None and cpg_idx < len(row) else ''
                trait = row[trait_idx].strip() if trait_idx is not None and trait_idx < len(row) and trait_idx >= 0 else ''
            except (IndexError, ValueError):
                stats['skipped_no_gene'] += 1
                continue

            if not gene or gene in ('NA', 'N/A', '-', ''):
                stats['skipped_no_gene'] += 1
                continue
            if not cpg or cpg in ('NA', 'N/A', '-', ''):
                stats['skipped_no_cpg'] += 1
                continue

            gene = gene.split(';')[0].split(',')[0].split('/')[0].strip()

            edges.append((gene, cpg))
            stats['valid_edges'] += 1
            stats['unique_genes'].add(gene)
            stats['unique_cpgs'].add(cpg)
            stats['gene_degree'][gene] += 1
            if trait:
                stats['cpg_traits'][cpg] = trait

    print(f"\n[STATS] EWAS Atlas 解析完成:")
    print(f"  总行数:         {stats['total_lines']}")
    print(f"  有效边数:       {stats['valid_edges']}")
    print(f"  跳过(无基因):   {stats['skipped_no_gene']}")
    print(f"  跳过(无CpG):    {stats['skipped_no_cpg']}")
    print(f"  唯一基因:       {len(stats['unique_genes'])}")
    print(f"  唯一CpG位点:    {len(stats['unique_cpgs'])}")
    print(f"  Top-10 甲基化多基因: {stats['gene_degree'].most_common(10)}")

    return edges, stats


def save_edges(edges, output_path):
    """保存基因-甲基化关联边"""
    with open(output_path, 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        for gene, cpg in sorted(set(edges)):
            writer.writerow([gene, cpg])
    print(f"[SAVE] 边文件: {output_path} ({len(set(edges))} 条去重边)")


def save_stats(stats, output_path):
    """保存统计信息"""
    stats_out = {
        'total_edges': stats['valid_edges'],
        'unique_genes': len(stats['unique_genes']),
        'unique_cpgs': len(stats['unique_cpgs']),
        'top_methylated_genes': stats['gene_degree'].most_common(20),
        'source': 'EWAS Atlas',
        'url': EWAS_ATLAS_URL,
    }
    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(stats_out, f, indent=2, ensure_ascii=False)
    print(f"[SAVE] 统计文件: {output_path}")


# ============================================================
# MAIN
# ============================================================
if __name__ == "__main__":
    print("=" * 60)
    print("EWAS Atlas 基因-甲基化关联提取")
    print("=" * 60)

    filepath = download_ewas_atlas()
    if filepath is None:
        print("\n[INFO] 跳过甲基化边提取 (文件未就绪)")
        print("[TODO] 请手动下载 EWAS Atlas 数据后重新运行")
        sys.exit(0)

    edges, stats = parse_ewas_atlas(filepath)

    if len(edges) == 0:
        print("[ERROR] 未提取到任何有效边")
        sys.exit(1)

    save_edges(edges, OUTPUT_EDGES)
    save_stats(stats, OUTPUT_STATS)

    print(f"\n{'='*60}")
    print("基因-甲基化关联边提取完成!")
    print(f"输出: {OUTPUT_EDGES}")
    print("=" * 60)