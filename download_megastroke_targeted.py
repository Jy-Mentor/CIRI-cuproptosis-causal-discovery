#!/usr/bin/env python3
"""
批量下载目标基因 eQTL 位置对应的 MEGASTROKE SNP
使用批量查询方式提高效率
"""

import pandas as pd
import requests
import time
from pathlib import Path
from collections import defaultdict

# JWT Token
JWT_TOKEN = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiIxNzU3ODgyODc4QHFxLmNvbSIsImlhdCI6MTc3ODMwNDA4MywiZXhwIjoxNzc5NTEzNjgzfQ.ZtcIUEx_xYtrVD_EE-UboKyLlC-lZBq2pjn-iYhJzxocqHdA-02K9n_Qbw-5ngQ07GHjjIYqVtmkZfJ3OJl1yI-tOMBBFzVKe0nkwDcB6-yBgjgBaxVm8vq_pbNrMwy_ZezY5ys9jx7I8T4bZYg9KeUbSwj04OfNP82kGcKXIOErXXVy-Ie3dbUogDRSjnCT-_32yNQxuWpiyYnPWSrWQbQ2HlUQiiDTdFGzWeJJKfSRvjQzdp5g3nccxht0m5A0UsPCdvkyHFEpvPVZ-NpjCjkgy8GbZBv4cmDMSc5JJL6HLO0eV508SRKxMdp-gL6qVdhGiJ2i9XmZQE27-aq_7g"

# 配置
EXPOSURE_DIR = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\exposure_matched\matched_data"
OUTPUT_FILE = r"D:\EQTL\mr_results_megastroke\megastroke_target_genes.csv"

print("="*70)
print("批量下载目标基因的 MEGASTROKE SNP")
print("="*70)

# 步骤 1: 收集所有 SNP 位置
print("\n步骤 1: 从暴露数据收集 SNP 位置...")

positions_by_chr = defaultdict(set)

exposure_files = list(Path(EXPOSURE_DIR).glob("*_exposure.csv"))
print(f"  找到 {len(exposure_files)} 个暴露数据文件")

for exposure_file in exposure_files:
    try:
        df = pd.read_csv(exposure_file)
        
        # 提取染色体和位置
        for _, row in df.iterrows():
            chr_ = str(row['CHR'])
            pos = int(row['BP'])
            positions_by_chr[chr_].add(pos)
        
    except Exception as e:
        print(f"  跳过 {exposure_file.name}: {e}")

# 转换为列表
chr_positions = {}
for chr_, positions in positions_by_chr.items():
    chr_positions[chr_] = sorted(list(positions))
    print(f"  染色体 {chr_}: {len(positions):,} 个位置")

total_positions = sum(len(positions) for positions in chr_positions.values())
print(f"\n  ✓ 共 {total_positions:,} 个唯一位置")

# 步骤 2: 按染色体区域批量下载
print("\n步骤 2: 从 MEGASTROKE API 下载 SNP...")

API_BASE = "https://api.openGWAS.io"
DATASET_ID = "ebi-a-GCST006908"

headers = {
    "Authorization": f"Bearer {JWT_TOKEN}",
    "Content-Type": "application/json"
}

all_snps = []
request_count = 0

# 对每条染色体，合并相邻位置形成连续区域
def merge_regions(positions, window_size=1000000):
    """将相邻的位置合并成连续区域（1MB 窗口）"""
    if not positions:
        return []
    
    regions = []
    positions = sorted(positions)
    region_start = positions[0]
    region_end = positions[0]
    
    for pos in positions[1:]:
        if pos - region_end < window_size:
            # 合并到当前区域
            region_end = pos
        else:
            # 保存当前区域，开始新区域
            regions.append((region_start, region_end))
            region_start = pos
            region_end = pos
    
    # 保存最后一个区域
    regions.append((region_start, region_end))
    
    return regions

for chr_ in sorted(chr_positions.keys(), key=lambda x: int(x) if x.isdigit() else 999):
    positions = chr_positions[chr_]
    
    # 合并区域
    regions = merge_regions(positions, window_size=500000)  # 500kb 窗口
    
    print(f"\n  染色体 {chr_}: {len(regions)} 个区域")
    
    for region_start, region_end in regions:
        try:
            # 查询该区域的 SNP
            response = requests.get(
                f"{API_BASE}/api/v1/associations",
                headers=headers,
                params={
                    "dataset_id": DATASET_ID,
                    "chromosome": f"chr{chr_}",
                    "start": region_start,
                    "end": region_end,
                    "page_size": 10000
                },
                timeout=60
            )
            
            request_count += 1
            
            if response.status_code == 200:
                snps = response.json()
                if len(snps) > 0:
                    all_snps.extend(snps)
                    print(f"    区域 {region_start:,}-{region_end:,}: {len(snps):,} 个 SNP")
            
            # 避免请求过快
            if request_count % 10 == 0:
                time.sleep(1)
            else:
                time.sleep(0.3)
                
        except Exception as e:
            print(f"    区域 {region_start:,}-{region_end:,}: 错误 - {e}")

print(f"\n  ✓ 共下载 {len(all_snps):,} 个 SNP")
print(f"  ✓ 发送 {request_count} 个请求")

# 步骤 3: 保存数据
if len(all_snps) > 0:
    print(f"\n步骤 3: 保存数据...")
    
    # 转换为 DataFrame
    df = pd.DataFrame(all_snps)
    
    # 选择需要的列
    columns_to_keep = ['snp', 'chromosome', 'position', 'effect_allele', 'other_allele', 
                      'eaf', 'beta', 'se', 'pval', 'n']
    
    df = df[columns_to_keep]
    df.columns = ['SNP', 'CHR', 'BP', 'EFFECT_ALLELE', 'OTHER_ALLELE', 
                 'EAF', 'BETA', 'SE', 'PVAL', 'N']
    
    # 添加结局信息
    df['outcome'] = 'Ischemic Stroke'
    df['id.outcome'] = DATASET_ID
    
    # 去重
    original_count = len(df)
    df = df.drop_duplicates(subset=['SNP'])
    print(f"  去重：{original_count:,} -> {len(df):,} 个 SNP")
    
    # 保存
    df.to_csv(OUTPUT_FILE, index=False, encoding='utf-8')
    
    print(f"  ✓ 保存到：{OUTPUT_FILE}")
    print(f"  ✓ 文件大小：{os.path.getsize(OUTPUT_FILE) / 1e6:.2f} MB")
    
    # 显示覆盖的基因
    print(f"\n  覆盖的染色体:")
    for chr_ in sorted(df['CHR'].unique(), key=lambda x: int(str(x)) if str(x).isdigit() else 999):
        count = len(df[df['CHR'] == chr_])
        print(f"    染色体 {chr_}: {count:,} 个 SNP")
else:
    print("\n✗ 未下载任何 SNP")
    print("  请检查:")
    print("  1. JWT Token 是否有效")
    print("  2. 网络连接")
    print("  3. API 限制")

print("\n" + "="*70)
print("完成！")
print("="*70)
