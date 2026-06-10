# -*- coding: utf-8 -*-
"""
下载 TaRGET II 脑组织数据 (P1: RNA FPKM + P2: ATAC peaks)
===========================================================
从 CSV 提取 Brain/Cortex 的 URL，批量下载到本地。
"""
import os
import sys
import pandas as pd
import requests
import time
from pathlib import Path
from urllib.parse import urlparse

# ---- 配置 ----
CSV_PATH = r"D:\TaRGET-DI_8LMyu.csv"
OUTPUT_BASE = Path(r"D:\反向网络药理学\GAT拓展维度\Toxi")

# 需要的文件类型
TARGET_TYPES = [
    # P1: RNA expression
    'Rsem.genes.fpkm.tsv',
    # P2: ATAC
    'narrowPeak',
]

# ---- 创建输出目录 ----
OUTPUT_BASE.mkdir(parents=True, exist_ok=True)
(OUTPUT_BASE / "rna_fpkm").mkdir(exist_ok=True)
(OUTPUT_BASE / "atac").mkdir(exist_ok=True)

# ---- 1. 从 CSV 提取 URL ----
print("=" * 70)
print("[1] 读取 CSV 并筛选 Brain/Cortex 数据...")
df = pd.read_csv(CSV_PATH)
brain = df[df['Tissue'].isin(['Brain', 'Cortex'])]
print(f"    总行数: {len(df)}, 脑组织: {len(brain)}")

# 筛选目标文件类型
target = brain[brain['File-Type'].isin(TARGET_TYPES)]
print(f"    目标文件: {len(target)}")

# 按文件类型统计
for ft in TARGET_TYPES:
    sub = target[target['File-Type'] == ft]
    if len(sub) > 0:
        sz = sub['Filesize'].str.replace('M','').str.replace('K','').astype(float).sum()
        print(f"      {ft}: {len(sub)} files, ~{sz:.0f} MB")

# ---- 2. 构建下载列表 ----
downloads = []
for _, row in target.iterrows():
    url = row['DownloadURL']
    ftype = row['File-Type']
    tissue = row['Tissue']
    exp_id = row['Experiment-ID']
    exposure = str(row['Exposure']).replace('/', '_')
    age = str(row['Age']).replace(' ', '_')
    sex = str(row['Sex'])

    # 构造文件名
    ext = url.split('.')[-1] if '.' in url.split('/')[-1] else 'tsv'
    if 'gz' in url.split('/')[-1] and url.split('/')[-1].endswith('.gz'):
        ext = '.'.join(url.split('/')[-1].split('.')[-2:])
    fname = f"{tissue}_{exposure}_{age}_{sex}_{exp_id}_{ftype}.{ext}"

    if ftype == 'Rsem.genes.fpkm.tsv':
        subdir = "rna_fpkm"
    elif ftype in ('narrowPeak', 'bigWig', 'PE.R1.bigWig'):
        subdir = "atac"
    else:
        continue

    local_path = OUTPUT_BASE / subdir / fname
    downloads.append((url, local_path, ftype))

print(f"\n[2] 待下载: {len(downloads)} 个文件")

# ---- 3. 批量下载 ----
print("\n[3] 开始下载...")
print("=" * 70)

success = 0
failed = 0
skipped = 0

headers = {
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36',
}

for i, (url, local_path, ftype) in enumerate(downloads):
    if local_path.exists():
        sz = local_path.stat().st_size
        if sz > 0:
            skipped += 1
            continue

    try:
        print(f"[{i+1}/{len(downloads)}] {ftype:25s} | {local_path.name[:70]}...", end=' ', flush=True)

        r = requests.get(url, headers=headers, timeout=300, stream=True)
        if r.status_code == 200:
            with open(local_path, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
            sz_mb = local_path.stat().st_size / 1e6
            print(f"OK ({sz_mb:.1f} MB)")
            success += 1
        else:
            print(f"FAILED HTTP {r.status_code}")
            failed += 1

    except Exception as e:
        print(f"ERROR: {str(e)[:60]}")
        failed += 1
        # 删除不完整文件
        if local_path.exists():
            local_path.unlink()

    # 避免请求过快
    time.sleep(0.3)

# ---- 4. 汇总 ----
print("\n" + "=" * 70)
print("[4] 下载完成!")
print(f"    成功: {success} | 跳过(已存在): {skipped} | 失败: {failed}")
print(f"    数据目录: {OUTPUT_BASE}")
print(f"      RNA FPKM: {OUTPUT_BASE / 'rna_fpkm'}")
print(f"      ATAC:     {OUTPUT_BASE / 'atac'}")

# 列出文件统计
for subdir in ['rna_fpkm', 'atac']:
    p = OUTPUT_BASE / subdir
    files = list(p.glob('*'))
    total_mb = sum(f.stat().st_size for f in files) / 1e6
    print(f"      {subdir}: {len(files)} files, {total_mb:.1f} MB")

print("\nDone!")