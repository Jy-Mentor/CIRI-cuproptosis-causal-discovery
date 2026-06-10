# -*- coding: utf-8 -*-
"""
从 TRRUST 数据库提取 TF-靶基因调控边
=====================================================
来源: TRRUST v2 (https://www.grnpedia.org/trrust/)
输入: TRRUST_Network_human.txt (Tab分隔, 约8000条人类调控关系)
输出: tf_target_edges.txt (tf  target  regulation_type)
      tf_target_nodes.csv (tf/target  node_type  regulation_count)

格式说明:
  - 输入列: TF  Target  Mode  PMID
  - Mode: Activation / Repression / Unknown
  - 输出边: TAB分隔, 无表头
  - 输出节点: CSV, 含node_type (TF/Gene) 和 degree

备选数据源:
  1. RegNetwork (https://regnetworkweb.org/) - 整合数据库
  2. ENCODE ChIP-seq - 需额外处理, 不如TRRUST直接
  3. hTFtarget (http://bioinfo.life.hust.edu.cn/hTFtarget/) - ChIP-seq证据

作者: 优化版 v2.0
日期: 2026-05-31
"""

import os
import sys
import csv
import urllib.request
import gzip
import shutil
from collections import Counter, defaultdict
from pathlib import Path

# ============================================================
# 0. 配置
# ============================================================
BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = Path(r"C:\Users\Jy-Mentor-7\Desktop\GAT")
OUTPUT_DIR = Path(r"D:\反向网络药理学\GAT拓展维度")

TRRUST_URL = "https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv"
TRRUST_LOCAL = DATA_DIR / "TRRUST_Network_human.txt"

OUTPUT_EDGES = OUTPUT_DIR / "tf_target_edges.txt"
OUTPUT_NODES = OUTPUT_DIR / "tf_target_nodes.csv"
OUTPUT_STATS = OUTPUT_DIR / "tf_target_stats.json"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def download_trrust():
    """下载 TRRUST 人类调控网络数据"""
    if TRRUST_LOCAL.exists():
        print(f"[INFO] 已存在本地文件: {TRRUST_LOCAL}")
        return str(TRRUST_LOCAL)

    print(f"[DOWNLOAD] 正在从 TRRUST 下载: {TRRUST_URL}")
    try:
        DATA_DIR.mkdir(parents=True, exist_ok=True)
        req = urllib.request.Request(TRRUST_URL, headers={
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
        })
        with urllib.request.urlopen(req, timeout=60) as response:
            content = response.read()
            with open(TRRUST_LOCAL, 'wb') as f:
                f.write(content)
        print(f"[OK] 下载完成: {TRRUST_LOCAL} ({len(content)} bytes)")
        return str(TRRUST_LOCAL)
    except Exception as e:
        print(f"[ERROR] 下载失败: {e}")
        print("[FALLBACK] 请手动下载到以下路径:")
        print(f"          {TRRUST_LOCAL}")
        print(f"          URL: {TRRUST_URL}")
        return None


def detect_separator(filepath):
    """自动检测文件分隔符"""
    with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
        first_line = f.readline().strip()
        if '\t' in first_line:
            return '\t'
        elif ',' in first_line:
            return ','
        elif ' ' in first_line:
            return ' '
        else:
            return '\t'


def parse_trrust(filepath):
    """
    解析 TRRUST 文件, 提取 TF-靶基因调控边

    TRRUST 格式 (Tab分隔):
        TF  Target  Mode  PMID
    """
    sep = detect_separator(filepath)
    print(f"[INFO] 检测到分隔符: {repr(sep)}")

    edges = []
    stats = {
        'total_lines': 0,
        'valid_edges': 0,
        'skipped_empty': 0,
        'skipped_self_loop': 0,
        'regulation_types': Counter(),
        'unique_tfs': set(),
        'unique_targets': set(),
        'tf_degree': Counter(),
        'target_degree': Counter(),
    }

    with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
        reader = csv.reader(f, delimiter=sep)
        for row_idx, row in enumerate(reader, start=1):
            stats['total_lines'] += 1

            if len(row) < 3:
                stats['skipped_empty'] += 1
                continue

            tf = row[0].strip().upper()
            target = row[1].strip().upper()
            mode = row[2].strip() if len(row) > 2 else 'Unknown'

            if not tf or not target:
                stats['skipped_empty'] += 1
                continue

            if tf == target:
                stats['skipped_self_loop'] += 1
                continue

            mode = mode.capitalize()
            if mode not in ('Activation', 'Repression', 'Unknown'):
                mode = 'Unknown'

            edges.append((tf, target, mode))
            stats['valid_edges'] += 1
            stats['regulation_types'][mode] += 1
            stats['unique_tfs'].add(tf)
            stats['unique_targets'].add(target)
            stats['tf_degree'][tf] += 1
            stats['target_degree'][target] += 1

    print(f"\n[STATS] 解析完成:")
    print(f"  总行数:       {stats['total_lines']}")
    print(f"  有效边数:     {stats['valid_edges']}")
    print(f"  跳过空行:     {stats['skipped_empty']}")
    print(f"  跳过自环:     {stats['skipped_self_loop']}")
    print(f"  唯一TF:       {len(stats['unique_tfs'])}")
    print(f"  唯一靶基因:   {len(stats['unique_targets'])}")
    print(f"  调控类型分布: {dict(stats['regulation_types'])}")
    print(f"  Top-10 TF (出度): {stats['tf_degree'].most_common(10)}")

    return edges, stats


def save_edges(edges, output_path):
    """保存边文件: tf  target  regulation_type"""
    with open(output_path, 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        for tf, target, mode in edges:
            writer.writerow([tf, target, mode])
    print(f"[SAVE] 边文件: {output_path} ({len(edges)} 条边)")


def save_nodes(stats, output_path):
    """保存节点文件: node_name  node_type  degree"""
    tfs = set()
    targets = set()
    for tf, target, _ in edges:
        tfs.add(tf)
        targets.add(target)

    all_nodes = []
    for tf in sorted(tfs):
        all_nodes.append({
            'node_name': tf,
            'node_type': 'TF',
            'degree': stats['tf_degree'][tf],
            'role': 'regulator'
        })
    for target in sorted(targets):
        all_nodes.append({
            'node_name': target,
            'node_type': 'Gene',
            'degree': stats['target_degree'][target],
            'role': 'target'
        })

    with open(output_path, 'w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['node_name', 'node_type', 'degree', 'role'])
        writer.writeheader()
        writer.writerows(all_nodes)
    print(f"[SAVE] 节点文件: {output_path} ({len(all_nodes)} 个节点)")


def save_stats(stats, output_path):
    """保存统计信息"""
    import json
    stats_out = {
        'total_edges': stats['valid_edges'],
        'unique_tfs': len(stats['unique_tfs']),
        'unique_targets': len(stats['unique_targets']),
        'regulation_types': dict(stats['regulation_types']),
        'top_tfs': stats['tf_degree'].most_common(20),
        'top_targets': stats['target_degree'].most_common(20),
        'source': 'TRRUST v2',
        'url': TRRUST_URL,
    }
    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(stats_out, f, indent=2, ensure_ascii=False)
    print(f"[SAVE] 统计文件: {output_path}")


def validate_output():
    """验证输出文件完整性"""
    edge_count = 0
    with open(OUTPUT_EDGES, 'r', encoding='utf-8') as f:
        for line in f:
            if line.strip():
                edge_count += 1

    node_count = 0
    with open(OUTPUT_NODES, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        node_count = sum(1 for _ in reader)

    print(f"\n[VALIDATE] 输出验证:")
    print(f"  tf_target_edges.txt:  {edge_count} 条边")
    print(f"  tf_target_nodes.csv:  {node_count} 个节点")

    assert edge_count > 0, "边文件为空!"
    assert node_count > 0, "节点文件为空!"
    print("[PASS] 所有验证通过 ✓")


# ============================================================
# MAIN
# ============================================================
if __name__ == "__main__":
    print("=" * 60)
    print("TRRUST TF-靶基因调控边提取")
    print("=" * 60)

    filepath = download_trrust()
    if filepath is None:
        sys.exit(1)

    edges, stats = parse_trrust(filepath)

    if len(edges) == 0:
        print("[ERROR] 未提取到任何有效边, 请检查输入文件格式")
        sys.exit(1)

    save_edges(edges, OUTPUT_EDGES)
    save_nodes(stats, OUTPUT_NODES)
    save_stats(stats, OUTPUT_STATS)
    validate_output()

    print("\n" + "=" * 60)
    print("TF-靶基因调控边提取完成!")
    print(f"输出: {OUTPUT_EDGES}")
    print(f"输出: {OUTPUT_NODES}")
    print("=" * 60)