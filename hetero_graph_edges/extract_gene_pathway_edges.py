# -*- coding: utf-8 -*-
"""
从 MSigDB / Reactome / KEGG 提取基因-通路关联边
=====================================================
来源 (按优先级):
  1. MSigDB c2.cp (Canonical Pathways) - KEGG + Reactome + BioCarta + PID
  2. Reactome GMT (ReactomePathways.gmt)
  3. KEGG REST API (实时获取)

输入: .gmt 文件 (每行: 通路名  描述  gene1  gene2  ...)
输出: gene_pathway_edges.txt (gene_symbol  pathway_name)
      pathway_nodes.csv (pathway_name  n_genes  source)

下载链接:
  - MSigDB: https://www.gsea-msigdb.org/gsea/msigdb/human/genesets/c2.cp.v2023.2.Hs.symbols.gmt
  - Reactome: https://reactome.org/download/current/ReactomePathways.gmt.zip

作者: 优化版 v2.0
日期: 2026-05-31
"""

import os
import sys
import csv
import urllib.request
import gzip
import zipfile
import io
from collections import Counter, defaultdict
from pathlib import Path

# ============================================================
# 0. 配置
# ============================================================
BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = Path(r"C:\Users\Jy-Mentor-7\Desktop\GAT")
OUTPUT_DIR = Path(r"D:\反向网络药理学\GAT拓展维度")

MSIGDB_URL = (
    "https://www.gsea-msigdb.org/gsea/msigdb/human/genesets/"
    "c2.cp.v2023.2.Hs.symbols.gmt"
)
REACTOME_URL = (
    "https://reactome.org/download/current/ReactomePathways.gmt.zip"
)
KEGG_LIST_URL = "https://rest.kegg.jp/list/pathway/hsa"
KEGG_GET_URL = "https://rest.kegg.jp/get/{}"

MSIGDB_LOCAL = DATA_DIR / "c2.cp.v2023.2.Hs.symbols.gmt"
REACTOME_ZIP = DATA_DIR / "ReactomePathways.gmt.zip"
REACTOME_LOCAL = DATA_DIR / "ReactomePathways.gmt"

OUTPUT_EDGES = OUTPUT_DIR / "gene_pathway_edges.txt"
OUTPUT_NODES = OUTPUT_DIR / "pathway_nodes.csv"
OUTPUT_STATS = OUTPUT_DIR / "gene_pathway_stats.json"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def download_msigdb():
    """下载 MSigDB C2 Canonical Pathways GMT"""
    if MSIGDB_LOCAL.exists():
        size_mb = MSIGDB_LOCAL.stat().st_size / (1024 * 1024)
        print(f"[INFO] 已存在本地文件: {MSIGDB_LOCAL} ({size_mb:.1f} MB)")
        return str(MSIGDB_LOCAL)

    print(f"[DOWNLOAD] 正在从 MSigDB 下载: {MSIGDB_URL}")
    try:
        DATA_DIR.mkdir(parents=True, exist_ok=True)
        req = urllib.request.Request(MSIGDB_URL, headers={
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'
        })
        with urllib.request.urlopen(req, timeout=180) as response:
            content = response.read()
            with open(MSIGDB_LOCAL, 'wb') as f:
                f.write(content)
        size_mb = len(content) / (1024 * 1024)
        print(f"[OK] 下载完成: {MSIGDB_LOCAL} ({size_mb:.1f} MB)")
        return str(MSIGDB_LOCAL)
    except Exception as e:
        print(f"[ERROR] MSigDB 下载失败: {e}")
        return None


def download_reactome():
    """下载 Reactome 通路 GMT"""
    if REACTOME_LOCAL.exists():
        print(f"[INFO] 已存在本地文件: {REACTOME_LOCAL}")
        return str(REACTOME_LOCAL)

    print(f"[DOWNLOAD] 正在从 Reactome 下载: {REACTOME_URL}")
    try:
        DATA_DIR.mkdir(parents=True, exist_ok=True)
        req = urllib.request.Request(REACTOME_URL, headers={
            'User-Agent': 'Mozilla/5.0'
        })
        with urllib.request.urlopen(req, timeout=180) as response:
            content = response.read()
            with zipfile.ZipFile(io.BytesIO(content)) as zf:
                zf.extractall(DATA_DIR)
        print(f"[OK] Reactome 解压完成: {REACTOME_LOCAL}")
        return str(REACTOME_LOCAL)
    except Exception as e:
        print(f"[ERROR] Reactome 下载失败: {e}")
        return None


def parse_gmt(filepath, source_label):
    """
    解析 GMT 文件, 提取基因-通路关联

    GMT 格式:
        pathway_name\tdescription\tgene1\tgene2\tgene3\t...

    返回: (gene, pathway) 列表
    """
    edges = []
    stats = {
        'source': source_label,
        'total_pathways': 0,
        'total_edges': 0,
        'genes_per_pathway': [],
        'pathway_sizes': {},
        'unique_genes': set(),
        'skipped_empty': 0,
    }

    with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
        for line_no, line in enumerate(f, start=1):
            line = line.strip()
            if not line:
                stats['skipped_empty'] += 1
                continue

            parts = line.split('\t')
            if len(parts) < 3:
                stats['skipped_empty'] += 1
                continue

            pathway_name = parts[0].strip()
            genes = [g.strip().upper() for g in parts[2:] if g.strip()]

            if not genes:
                stats['skipped_empty'] += 1
                continue

            stats['total_pathways'] += 1
            stats['pathway_sizes'][pathway_name] = len(genes)
            stats['genes_per_pathway'].append(len(genes))

            for gene in genes:
                edges.append((gene, pathway_name))
                stats['unique_genes'].add(gene)

            stats['total_edges'] += len(genes)

    return edges, stats


def fetch_kegg_pathways():
    """
    从 KEGG REST API 获取人类通路-基因关联 (实时)

    KEGG 的使用限制: 每IP最多10次/秒, 建议分批获取
    返回: (gene, pathway_id|pathway_name) 列表
    """
    print("[KEGG] 正在获取人类通路列表...")
    try:
        req = urllib.request.Request(KEGG_LIST_URL)
        with urllib.request.urlopen(req, timeout=30) as response:
            content = response.read().decode('utf-8')
    except Exception as e:
        print(f"[ERROR] KEGG 通路列表获取失败: {e}")
        return [], {}

    pathway_ids = []
    for line in content.strip().split('\n'):
        if not line.strip():
            continue
        parts = line.split('\t')
        if len(parts) >= 2:
            pathway_id = parts[0].strip().replace('path:', '')
            pathway_name = parts[1].strip()
            pathway_ids.append((pathway_id, pathway_name))

    print(f"[KEGG] 找到 {len(pathway_ids)} 个人类通路")

    edges = []
    stats = {
        'source': 'KEGG',
        'total_pathways': 0,
        'total_edges': 0,
        'pathway_sizes': {},
        'unique_genes': set(),
        'failed_pathways': [],
    }

    import time
    for i, (pid, pname) in enumerate(pathway_ids):
        if i % 5 == 0 and i > 0:
            time.sleep(0.3)

        try:
            url = KEGG_GET_URL.format(pid)
            req = urllib.request.Request(url)
            with urllib.request.urlopen(req, timeout=30) as response:
                text = response.read().decode('utf-8')

            genes = []
            in_gene_section = False
            for line in text.split('\n'):
                if line.startswith('GENE'):
                    in_gene_section = True
                    parts = line.split(';')
                    if len(parts) >= 2:
                        gene_desc = parts[1].strip()
                        gene_symbol = gene_desc.split(' ')[0].strip()
                        if gene_symbol:
                            genes.append(gene_symbol.upper())
                elif line.lstrip() and not line.startswith(' ') and in_gene_section:
                    in_gene_section = False
                elif in_gene_section and line.startswith(' '):
                    parts = line.strip().split(';')
                    if len(parts) >= 2:
                        gene_desc = parts[1].strip()
                        gene_symbol = gene_desc.split(' ')[0].strip()
                        if gene_symbol:
                            genes.append(gene_symbol.upper())

            if genes:
                stats['total_pathways'] += 1
                stats['pathway_sizes'][f"{pid}|{pname}"] = len(genes)
                for gene in set(genes):
                    edges.append((gene, f"{pid}|{pname}"))
                    stats['unique_genes'].add(gene)
                stats['total_edges'] += len(genes)

        except Exception as e:
            stats['failed_pathways'].append(pid)

        if (i + 1) % 50 == 0:
            print(f"[KEGG] 进度: {i+1}/{len(pathway_ids)}")

    print(f"[KEGG] 完成: {stats['total_pathways']} 通路, {stats['total_edges']} 边, "
          f"{len(stats['failed_pathways'])} 个通路失败")
    return edges, stats


def merge_edges(all_edges_list, all_stats_list):
    """合并多个来源的边, 去重, 统计"""
    merged = {}
    merged_stats = {
        'total_edges_dedup': 0,
        'total_edges_raw': 0,
        'unique_genes': set(),
        'unique_pathways': set(),
        'gene_degree': Counter(),
        'pathway_degree': Counter(),
        'sources': {},
    }

    for edges, stats in all_edges_list:
        source = stats['source']
        merged_stats['sources'][source] = {
            'pathways': stats['total_pathways'],
            'edges': stats['total_edges'],
            'genes': len(stats['unique_genes']),
        }
        merged_stats['total_edges_raw'] += stats['total_edges']

        for gene, pathway in edges:
            key = (gene, pathway)
            if key not in merged:
                merged[key] = source
                merged_stats['unique_genes'].add(gene)
                merged_stats['unique_pathways'].add(pathway)
                merged_stats['gene_degree'][gene] += 1
                merged_stats['pathway_degree'][pathway] += 1

    merged_stats['total_edges_dedup'] = len(merged)
    return list(merged.keys()), merged_stats


def save_edges(edges, output_path):
    """保存基因-通路关联边"""
    with open(output_path, 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        for gene, pathway in sorted(edges):
            writer.writerow([gene, pathway])
    print(f"[SAVE] 边文件: {output_path} ({len(edges)} 条边)")


def save_pathway_nodes(edges, output_path):
    """生成通路节点CSV (含统计信息, 特征向量待后续填充)"""
    pathway_info = {}
    for gene, pathway in edges:
        if pathway not in pathway_info:
            pathway_info[pathway] = {'n_genes': 0, 'genes': []}
        pathway_info[pathway]['n_genes'] += 1
        pathway_info[pathway]['genes'].append(gene)

    with open(output_path, 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['pathway_name', 'n_genes', 'source', 'example_genes'])
        for pname in sorted(pathway_info.keys()):
            info = pathway_info[pname]
            source = 'MSigDB' if 'REACTOME' in pname else 'KEGG' if 'KEGG' in pname else 'MSigDB'
            example = ','.join(info['genes'][:5])
            writer.writerow([pname, info['n_genes'], source, example])

    print(f"[SAVE] 通路节点: {output_path} ({len(pathway_info)} 个通路)")


def save_stats(stats, output_path):
    """保存统计信息"""
    import json
    stats_out = {
        'total_edges': stats['total_edges_dedup'],
        'total_edges_raw': stats['total_edges_raw'],
        'unique_genes': len(stats['unique_genes']),
        'unique_pathways': len(stats['unique_pathways']),
        'sources': stats['sources'],
        'top_pathways_by_genes': stats['pathway_degree'].most_common(20),
        'top_genes_by_pathways': stats['gene_degree'].most_common(20),
        'pathway_size_distribution': {
            'min': min(stats['pathway_degree'].values()),
            'max': max(stats['pathway_degree'].values()),
            'mean': sum(stats['pathway_degree'].values()) / len(stats['pathway_degree']),
        },
    }
    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(stats_out, f, indent=2, ensure_ascii=False)
    print(f"[SAVE] 统计文件: {output_path}")


def print_summary(edges, stats):
    """打印摘要"""
    print(f"\n{'='*60}")
    print(f"基因-通路关联提取摘要")
    print(f"{'='*60}")
    print(f"  去重后边数:     {stats['total_edges_dedup']:,}")
    print(f"  原始边数:       {stats['total_edges_raw']:,}")
    print(f"  唯一基因:       {len(stats['unique_genes']):,}")
    print(f"  唯一通路:       {len(stats['unique_pathways']):,}")
    print(f"  各来源统计:")
    for source, s in stats['sources'].items():
        print(f"    {source}: {s['pathways']} 通路, {s['edges']:,} 边, {s['genes']:,} 基因")
    print(f"  Top-10 通路(基因数): {stats['pathway_degree'].most_common(10)}")
    print(f"  Top-10 多通路基因:   {stats['gene_degree'].most_common(10)}")


# ============================================================
# MAIN
# ============================================================
if __name__ == "__main__":
    print("=" * 60)
    print("基因-通路关联边提取")
    print("=" * 60)

    all_edges_list = []

    # 1. MSigDB C2 Canonical Pathways
    msigdb_path = download_msigdb()
    if msigdb_path:
        print(f"\n[MSigDB] 解析 GMT 文件...")
        edges, stats = parse_gmt(msigdb_path, 'MSigDB')
        all_edges_list.append((edges, stats))
        print(f"[MSigDB] {stats['total_pathways']} 通路, {stats['total_edges']:,} 边, "
              f"{len(stats['unique_genes']):,} 基因")

    # 2. Reactome (如果不与MSigDB重复, 可补充)
    reactome_path = download_reactome()
    if reactome_path:
        print(f"\n[Reactome] 解析 GMT 文件...")
        edges, stats = parse_gmt(reactome_path, 'Reactome')
        all_edges_list.append((edges, stats))
        print(f"[Reactome] {stats['total_pathways']} 通路, {stats['total_edges']:,} 边, "
              f"{len(stats['unique_genes']):,} 基因")

    # 3. KEGG REST API (可选, 较慢)
    # 仅当 MSigDB 不可用时启用
    if not msigdb_path and not reactome_path:
        print("\n[KEGG] 回退到 KEGG REST API...")
        edges, stats = fetch_kegg_pathways()
        if edges:
            all_edges_list.append((edges, stats))

    if not all_edges_list:
        print("[ERROR] 未能获取任何数据源, 请检查网络连接")
        print("[MANUAL] 手动下载 MSigDB GMT:")
        print(f"         {MSIGDB_URL}")
        sys.exit(1)

    merged_edges, merged_stats = merge_edges(all_edges_list, [])

    save_edges(merged_edges, OUTPUT_EDGES)
    save_pathway_nodes(merged_edges, OUTPUT_NODES)
    save_stats(merged_stats, OUTPUT_STATS)
    print_summary(merged_edges, merged_stats)

    print(f"\n{'='*60}")
    print("基因-通路关联边提取完成!")
    print(f"输出: {OUTPUT_EDGES}")
    print(f"输出: {OUTPUT_NODES}")
    print("=" * 60)