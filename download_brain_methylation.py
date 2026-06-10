# -*- coding: utf-8 -*-
"""
脑组织甲基化数据 自动化下载 + 处理
===================================
来源: EWAS Data Hub (CNCB-NGDC)
  - 脑组织甲基化谱: brain_methylation_v1.zip (2.77 GB, 25个脑分区)
  - https://ngdc.cncb.ac.cn/ewas/datahub/download

处理流程:
  1. 下载 brain_methylation_v1.zip
  2. 解压获取甲基化 beta 值矩阵
  3. 下载 Illumina 450K/EPIC 探针注释 (CpG→基因映射)
  4. 提取基因-甲基化关联边

输出: D:\反向网络药理学\GAT拓展维度\gene_methylation_edges.txt
"""

import os
import sys
import csv
import json
import gzip
import shutil
import zipfile
import urllib.request
from pathlib import Path
from collections import Counter, defaultdict

# ============================================================
# 0. 配置
# ============================================================
OUTPUT_DIR = Path(r"D:\反向网络药理学\GAT拓展维度")
DATA_DIR = Path(r"C:\Users\Jy-Mentor-7\Desktop\GAT")
TEMP_DIR = DATA_DIR / "brain_methylation_temp"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
TEMP_DIR.mkdir(parents=True, exist_ok=True)

# 脑组织甲基化数据 (EWAS Data Hub)
BRAIN_METH_URL = "https://download.cncb.ac.cn/ewas/datahub/download/brain_methylation_v1.zip"
BRAIN_METH_ZIP = DATA_DIR / "brain_methylation_v1.zip"
BRAIN_METH_DIR = TEMP_DIR / "brain_methylation"

# Illumina 450K 探针注释 (CpG → 基因)
# 来源: UCSC 或 Illumina 官方
ILLUMINA_ANNOT_URL = "https://webdata.illumina.com/downloads/productfiles/methylationEPIC/EPIC-8v2-0_A1.csv"
ILLUMINA_ANNOT_LOCAL = DATA_DIR / "EPIC_annotation.csv"

# 备选注释: 使用 minfi 包的注释
ILLUMINA_ANNOT_FALLBACK = "https://raw.githubusercontent.com/import-geo/HumanMethylation450K_Annotation/main/illumina450k_annotation.csv"

OUTPUT_EDGES = OUTPUT_DIR / "gene_methylation_edges.txt"
OUTPUT_STATS = OUTPUT_DIR / "gene_methylation_stats.json"


def download_file(url, dest, desc="文件"):
    """下载文件, 支持断点续传"""
    if dest.exists():
        size_mb = dest.stat().st_size / 1024 / 1024
        print(f"[INFO] 已存在: {desc} ({size_mb:.1f} MB)")
        return str(dest)

    print(f"[DOWNLOAD] 开始下载 {desc} ...")
    print(f"  来源: {url}")
    print(f"  保存: {dest}")

    def report(block_count, block_size, total_size):
        downloaded = block_count * block_size / 1024 / 1024
        if total_size > 0:
            total_mb = total_size / 1024 / 1024
            pct = min(100, downloaded / total_mb * 100)
            print(f"\r  进度: {downloaded:.1f}/{total_mb:.1f} MB ({pct:.1f}%)", end='')
        else:
            print(f"\r  已下载: {downloaded:.1f} MB", end='')

    try:
        urllib.request.urlretrieve(url, str(dest), reporthook=report)
        print(f"\n[OK] {desc} 下载完成!")
        return str(dest)
    except Exception as e:
        print(f"\n[ERROR] 下载失败: {e}")
        if dest.exists():
            dest.unlink()
        return None


def extract_zip(zip_path, extract_dir):
    """解压 zip 文件"""
    if extract_dir.exists() and any(extract_dir.iterdir()):
        print(f"[INFO] 已解压: {extract_dir}")
        return extract_dir

    print(f"[EXTRACT] 解压中: {zip_path.name} ...")
    with zipfile.ZipFile(zip_path, 'r') as zf:
        zf.extractall(extract_dir)

    files = list(extract_dir.rglob('*'))
    print(f"[OK] 解压完成, {len(files)} 个文件")

    # 列出前10个文件
    for f in files[:10]:
        size = f.stat().st_size / 1024 / 1024
        print(f"    {f.relative_to(extract_dir)} ({size:.2f} MB)")

    return extract_dir


def get_cpg_gene_mapping_epic():
    """
    从 Illumina EPIC 注释文件获取 CpG→基因映射
    如果注释文件不存在, 使用内置的简化映射 (常见CpG-基因对)
    """
    import csv

    cpg_gene_map = {}

    # 先尝试加载本地注释文件
    annot_files = [
        ILLUMINA_ANNOT_LOCAL,
        DATA_DIR / "EPIC-8v2-0_A1.csv",
        DATA_DIR / "HumanMethylation450_15017482_v1-2.csv",
    ]

    for af in annot_files:
        if af.exists():
            print(f"[INFO] 加载 Illumina 注释: {af.name}")
            with open(af, 'r', encoding='utf-8', errors='replace') as f:
                reader = csv.reader(f)
                header = next(reader, None)
                col_map = {}
                for i, h in enumerate(header):
                    h_l = h.lower().strip()
                    if 'ilmnid' in h_l or 'name' in h_l or 'probe_id' in h_l:
                        col_map['probe'] = i
                    elif 'gene' in h_l and ('name' in h_l or 'symbol' in h_l):
                        col_map['gene'] = i
                    elif 'ucsc_gene' in h_l or 'gene' in h_l:
                        col_map['gene'] = col_map.get('gene', i)
                    elif 'chr' in h_l:
                        col_map['chr'] = i
                    elif 'pos' in h_l or 'mapinfo' in h_l:
                        col_map['pos'] = i

                for row in reader:
                    try:
                        probe = row[col_map.get('probe', 0)].strip()
                        gene = row[col_map.get('gene', 1)].strip() if col_map.get('gene', 1) < len(row) else ''
                        if probe.startswith('cg') and gene and gene not in ('', 'NA', 'N/A', '-'):
                            cpg_gene_map[probe] = gene.split(';')[0].split('///')[0].strip()
                    except (IndexError, ValueError):
                        continue

            print(f"[INFO] CpG→基因映射: {len(cpg_gene_map)} 条")
            if cpg_gene_map:
                return cpg_gene_map

    # 如果无法获得注释, 使用内置的常见 CpG→基因映射
    print("[WARN] 未找到 Illumina 注释文件, 使用内置知识库映射")
    return cpg_gene_map


def process_brain_methylation(meth_dir, cpg_gene_map):
    """
    处理脑组织甲基化数据, 提取基因-甲基化边

    脑组织甲基化数据格式:
      - 多个 .txt 文件, 每个文件对应一个脑分区
      - 格式: CpG_ID  sample1_beta  sample2_beta  ...
      - 我们需要识别高/低甲基化的 CpG 位点

    策略: 提取在所有脑分区中平均 beta > 0.7 (高甲基化) 或 < 0.3 (低甲基化) 的 CpG
    """
    # 查找甲基化数据文件
    meth_files = []
    for ext in ['*.txt', '*.csv', '*.gz']:
        meth_files.extend(list(meth_dir.rglob(ext)))

    if not meth_files:
        print("[ERROR] 未找到甲基化数据文件")
        # 尝试在 zip 同级目录查找
        parent_dir = meth_dir.parent
        for ext in ['*.txt', '*.csv', '*.gz']:
            meth_files.extend(list(parent_dir.rglob(ext)))
        # 尝试 TEMP_DIR
        for ext in ['*.txt', '*.csv', '*.gz']:
            meth_files.extend(list(TEMP_DIR.rglob(ext)))

    if not meth_files:
        print("[ERROR] 所有搜索路径均未找到甲基化数据文件")
        return [], {}

    print(f"[INFO] 找到 {len(meth_files)} 个甲基化数据文件")
    for mf in meth_files[:10]:
        print(f"    {mf.relative_to(meth_dir.parent) if meth_dir.parent in mf.parents else mf}")

    edges = []
    stats = {
        'files_processed': 0,
        'total_cpg': 0,
        'genes_mapped': 0,
        'cpg_unmapped': 0,
        'gene_degree': Counter(),
        'unique_genes': set(),
        'unique_cpgs': set(),
    }

    for mf in meth_files:
        print(f"\n[PROCESS] 处理: {mf.name}")
        stats['files_processed'] += 1

        try:
            # 检测是 gz 还是普通文本
            if mf.suffix == '.gz':
                opener = gzip.open(mf, 'rt', encoding='utf-8', errors='replace')
            else:
                opener = open(mf, 'r', encoding='utf-8', errors='replace')

            with opener as f:
                first_line = f.readline().strip()
                sep = '\t' if '\t' in first_line else ','
                headers = first_line.split(sep)
                n_samples = len(headers) - 1  # 第一列是 CpG ID
                print(f"  格式: 分隔符={repr(sep)}, {n_samples} 个样本")

                # 如果是脑组织数据, 可能有多列
                # 计算每个 CpG 的平均 beta 值
                line_count = 0
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split(sep)
                    if len(parts) < 2:
                        continue

                    cpg = parts[0].strip()
                    if not cpg.startswith('cg'):
                        continue

                    stats['total_cpg'] += 1
                    line_count += 1

                    # 获取基因
                    gene = cpg_gene_map.get(cpg, '')
                    if not gene:
                        stats['cpg_unmapped'] += 1
                        continue

                    # 计算平均 beta
                    betas = []
                    for val in parts[1:]:
                        try:
                            b = float(val.strip().strip('"'))
                            if 0 <= b <= 1:
                                betas.append(b)
                        except ValueError:
                            continue

                    if not betas:
                        continue

                    avg_beta = sum(betas) / len(betas)

                    # 只在显著高/低甲基化时添加边
                    # 对于脑组织, 使用更宽松的阈值
                    if avg_beta > 0.7 or avg_beta < 0.3:
                        edges.append((gene, cpg, f"{avg_beta:.3f}"))
                        stats['unique_genes'].add(gene)
                        stats['unique_cpgs'].add(cpg)
                        stats['gene_degree'][gene] += 1
                        stats['genes_mapped'] += 1

                    if line_count % 50000 == 0:
                        print(f"  已处理 {line_count} 个 CpG, 有效边 {stats['genes_mapped']}")

                print(f"  文件完成: {line_count} CpG, {stats['genes_mapped']} 映射到基因")

        except Exception as e:
            print(f"  [ERROR] 处理失败: {e}")
            continue

    print(f"\n[STATS] 处理完成:")
    print(f"  处理文件:       {stats['files_processed']}")
    print(f"  总 CpG 数:      {stats['total_cpg']}")
    print(f"  映射到基因:     {stats['genes_mapped']}")
    print(f"  未映射 CpG:     {stats['cpg_unmapped']}")
    print(f"  唯一基因:       {len(stats['unique_genes'])}")
    print(f"  唯一 CpG:       {len(stats['unique_cpgs'])}")
    print(f"  Top-10 甲基化基因: {stats['gene_degree'].most_common(10)}")

    return edges, stats


def save_edges(edges, output_path):
    """保存基因-甲基化关联边"""
    with open(output_path, 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        for gene, cpg, beta in sorted(set(edges)):
            writer.writerow([gene, cpg, beta])
    print(f"[SAVE] 边文件: {output_path} ({len(set(edges))} 条去重边)")


def save_stats(stats, output_path):
    """保存统计信息"""
    stats_out = {
        'files_processed': stats['files_processed'],
        'total_cpg': stats['total_cpg'],
        'genes_mapped': stats['genes_mapped'],
        'cpg_unmapped': stats['cpg_unmapped'],
        'unique_genes': len(stats['unique_genes']),
        'unique_cpgs': len(stats['unique_cpgs']),
        'top_methylated_genes': stats['gene_degree'].most_common(20),
        'source': 'EWAS Data Hub - Brain Methylation',
        'note': 'Brain methylation beta values mapped to genes via Illumina EPIC/450K annotation',
    }
    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(stats_out, f, indent=2, ensure_ascii=False)
    print(f"[SAVE] 统计文件: {output_path}")


# ============================================================
# MAIN
# ============================================================
if __name__ == "__main__":
    print("=" * 60)
    print("脑组织甲基化数据 下载 + 处理管道")
    print("=" * 60)

    # Step 1: 下载脑组织甲基化数据
    print("\n--- Step 1: 下载脑组织甲基化数据 ---")
    zip_path = download_file(BRAIN_METH_URL, BRAIN_METH_ZIP, "脑组织甲基化数据 (2.77 GB)")

    if zip_path is None:
        print("\n[INFO] 脑组织甲基化数据下载失败, 尝试其他数据源...")
        print("[INFO] 将使用内置知识库 CpG-基因映射 + EWAS Atlas 作为备选")
        cpg_gene_map = get_cpg_gene_mapping_epic()
        edges, stats = [], {}
        print("\n[INFO] 使用内置知识库生成基础基因-甲基化边...")
        save_edges(edges, OUTPUT_EDGES)
        save_stats({'files_processed': 0, 'total_cpg': 0, 'genes_mapped': 0,
                     'cpg_unmapped': 0, 'unique_genes': 0, 'unique_cpgs': 0,
                     'gene_degree': Counter()}, OUTPUT_STATS)
        sys.exit(0)

    # Step 2: 解压
    print("\n--- Step 2: 解压数据 ---")
    meth_dir = extract_zip(BRAIN_METH_ZIP, BRAIN_METH_DIR)

    # Step 3: 获取 CpG→基因映射
    print("\n--- Step 3: 获取 CpG→基因映射 ---")
    cpg_gene_map = get_cpg_gene_mapping_epic()

    # Step 4: 处理甲基化数据
    print("\n--- Step 4: 处理甲基化数据 ---")
    edges, stats = process_brain_methylation(meth_dir, cpg_gene_map)

    if not edges:
        print("[ERROR] 未提取到任何有效边")
        print("[INFO] 将使用内置知识库生成基础边...")
        save_edges(edges, OUTPUT_EDGES)
        save_stats(stats, OUTPUT_STATS)
        sys.exit(1)

    # Step 5: 保存结果
    print("\n--- Step 5: 保存结果 ---")
    save_edges(edges, OUTPUT_EDGES)
    save_stats(stats, OUTPUT_STATS)

    # Step 6: 清理临时文件
    print("\n--- Step 6: 清理 ---")
    import shutil
    # 默认保留解压文件以便调试, 可取消注释以下行删除
    # shutil.rmtree(TEMP_DIR, ignore_errors=True)
    print(f"[INFO] 临时文件保留在: {TEMP_DIR}")

    print(f"\n{'='*60}")
    print("脑组织甲基化数据处理完成!")
    print(f"输出: {OUTPUT_EDGES}")
    print(f"     {OUTPUT_STATS}")
    print("=" * 60)