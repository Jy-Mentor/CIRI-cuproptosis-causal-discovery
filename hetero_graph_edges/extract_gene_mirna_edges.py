# -*- coding: utf-8 -*-
"""
从 miRTarBase 提取人类 miRNA-靶基因关联 (可选边类型)
=====================================================
来源: miRTarBase (https://mirtarbase.cuhk.edu.cn/)
  - miRTarBase 收录了经实验验证的 miRNA-靶基因关系
  - 支持多种实验证据: Reporter assay, Western blot, qPCR, NGS 等

输入: miRTarBase 下载的 Excel/CSV/TXT 文件
输出: gene_mirna_edges.txt (gene_symbol  mirna_id  evidence_type)

备选数据源:
  1. TargetScan (http://www.targetscan.org/vert_80/) - 预测保守靶基因
  2. miRDB (https://mirdb.org/) - 机器学习预测
  3. miRWalk (http://mirwalk.umm.uni-heidelberg.de/) - 整合数据库
  4. DIANA-TarBase (https://diana.e-ce.uth.gr/) - 实验验证

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

MIRTARBASE_URL = "https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2025/MTI.xls"
MIRTARBASE_LOCAL_CSV = DATA_DIR / "miRTarBase_MTI.csv"
MIRTARBASE_LOCAL_XLS = DATA_DIR / "miRTarBase_MTI.xls"
MIRTARBASE_LOCAL_HUMAN = DATA_DIR / "hsa_MTI.xlsx"

OUTPUT_EDGES = OUTPUT_DIR / "gene_mirna_edges.txt"
OUTPUT_STATS = OUTPUT_DIR / "gene_mirna_stats.json"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# 实验证据类型优先级 (强 → 弱)
EVIDENCE_PRIORITY = {
    'Reporter assay': 5,
    'Western blot': 5,
    'qPCR': 4,
    'Luciferase reporter assay': 5,
    'Dual-Luciferase reporter assay': 5,
    'Microarray': 3,
    'pSILAC': 3,
    'CLIP-seq': 3,
    'HITS-CLIP': 3,
    'PAR-CLIP': 3,
    'iCLIP': 3,
    'CLASH': 3,
    'Degradome': 3,
    'NGS': 2,
    'RNA-seq': 2,
    'Other': 1,
}


def download_mirtarbase():
    """引导用户下载 miRTarBase 数据"""
    existing = [p for p in [MIRTARBASE_LOCAL_CSV, MIRTARBASE_LOCAL_XLS, MIRTARBASE_LOCAL_HUMAN] if p.exists()]
    if existing:
        print(f"[INFO] 已存在本地文件: {existing[0]}")
        return str(existing[0])

    print("[INFO] miRTarBase 需要手动下载")
    print(f"  URL: {MIRTARBASE_URL}")
    print(f"  保存到: {MIRTARBASE_LOCAL_CSV}")
    print(f"")
    print(f"  步骤:")
    print(f"  1. 访问 https://mirtarbase.cuhk.edu.cn/")
    print(f"  2. 点击 'Download' → 选择最新版本")
    print(f"  3. 下载 'Homo sapiens MTI' 文件 (hsa_MTI.xlsx)")
    print(f"  4. 或将文件另存为 CSV 格式")
    print(f"  5. 保存到上述路径并重新运行本脚本")
    print(f"")
    print(f"  备选: TargetScan 预测数据")
    print(f"    http://www.targetscan.org/vert_80/ → 下载 'Predicted_Targets_Context_Scores.default_predictions.txt'")
    return None


def detect_format(filepath):
    """自动检测 miRTarBase 文件格式"""
    ext = Path(filepath).suffix.lower()

    if ext in ('.xls', '.xlsx'):
        try:
            import openpyxl
            df = pd.read_excel(filepath, nrows=3)
        except ImportError:
            try:
                df = pd.read_excel(filepath, nrows=3, engine='xlrd')
            except ImportError:
                print("[WARN] 缺少 openpyxl/xlrd, 尝试读取为 CSV")
                return None, 'csv'

        cols = df.columns.tolist()
        print(f"[INFO] Excel 列名: {cols[:10]}...")
        return cols, 'excel'
    else:
        with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
            first_line = f.readline().strip()
            sep = '\t' if '\t' in first_line else ','
            print(f"[INFO] 分隔符: {repr(sep)}")
            headers = first_line.split(sep)
            print(f"[INFO] 表头: {headers[:10]}...")
            return headers, sep


def parse_mirtarbase(filepath):
    """
    解析 miRTarBase 文件, 提取 miRNA-靶基因关联

    miRTarBase 常见格式:
      miRTarBase ID  miRNA  Target Gene  Target Gene (Entrez)  Experiments  ...
    """
    import pandas as pd

    ext = Path(filepath).suffix.lower()

    if ext in ('.xls', '.xlsx'):
        df = pd.read_excel(filepath)
    else:
        df = pd.read_csv(filepath, sep=None, engine='python', encoding='utf-8', errors='replace')

    cols_lower = {c.lower(): c for c in df.columns}

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
        elif ('gene' in kl and 'symbol' in kl) and gene_col is None:
            gene_col = col
        elif 'experiment' in kl and evidence_col is None:
            evidence_col = col
        elif 'species' in kl and species_col is None:
            species_col = col

    if mirna_col is None:
        mirna_col = df.columns[0]
    if gene_col is None:
        gene_col = df.columns[2] if len(df.columns) > 2 else df.columns[1]

    print(f"[INFO] miRNA列: '{mirna_col}', 基因列: '{gene_col}', "
          f"证据列: '{evidence_col}', 物种列: '{species_col}'")

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

    for _, row in df.iterrows():
        try:
            mirna = str(row[mirna_col]).strip() if mirna_col in df.columns else ''
            gene = str(row[gene_col]).strip().upper() if gene_col in df.columns else ''
        except (ValueError, KeyError):
            stats['skipped_empty'] += 1
            continue

        if not mirna or not gene or gene in ('NA', 'NAN', 'NONE', ''):
            stats['skipped_empty'] += 1
            continue

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

    print(f"\n[STATS] miRTarBase 解析完成:")
    print(f"  总行数:         {stats['total_rows']}")
    print(f"  有效边数:       {stats['valid_edges']}")
    print(f"  跳过(非人类):   {stats['skipped_non_human']}")
    print(f"  跳过(空值):     {stats['skipped_empty']}")
    print(f"  唯一 miRNA:     {len(stats['unique_mirnas'])}")
    print(f"  唯一靶基因:     {len(stats['unique_genes'])}")
    print(f"  Top-10 miRNA (靶向多基因): {stats['mirna_degree'].most_common(10)}")
    print(f"  Top-10 基因 (被多miRNA靶向): {stats['gene_degree'].most_common(10)}")

    return edges, stats


def save_edges(edges, output_path):
    """保存 miRNA-靶基因关联边, 去重取最高优先级证据"""
    best_evidence = {}
    for gene, mirna, evidence in edges:
        key = (gene, mirna)
        current_best = best_evidence.get(key, ('Other', 0))
        new_priority = EVIDENCE_PRIORITY.get(evidence, 1)
        if new_priority > current_best[1]:
            best_evidence[key] = (evidence, new_priority)

    with open(output_path, 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        for (gene, mirna), (evidence, _) in sorted(best_evidence.items()):
            writer.writerow([gene, mirna, evidence])

    print(f"[SAVE] 边文件: {output_path} ({len(best_evidence)} 条去重边)")
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
        'source': 'miRTarBase',
        'url': MIRTARBASE_URL,
    }
    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(stats_out, f, indent=2, ensure_ascii=False)
    print(f"[SAVE] 统计文件: {output_path}")


# ============================================================
# MAIN
# ============================================================
if __name__ == "__main__":
    print("=" * 60)
    print("miRTarBase miRNA-靶基因关联提取")
    print("=" * 60)

    filepath = download_mirtarbase()
    if filepath is None:
        print("\n[INFO] 跳过 miRNA 边提取 (文件未就绪)")
        print("[TODO] 请手动下载 miRTarBase 数据后重新运行")
        sys.exit(0)

    edges, stats = parse_mirtarbase(filepath)

    if len(edges) == 0:
        print("[ERROR] 未提取到任何有效边")
        sys.exit(1)

    dedup_count = save_edges(edges, OUTPUT_EDGES)
    save_stats(stats, dedup_count, OUTPUT_STATS)

    print(f"\n{'='*60}")
    print("miRNA-靶基因关联边提取完成!")
    print(f"输出: {OUTPUT_EDGES}")
    print("=" * 60)