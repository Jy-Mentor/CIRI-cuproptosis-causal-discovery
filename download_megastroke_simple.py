#!/usr/bin/env python3
"""
下载 MEGASTROKE 数据并提取与目标基因相关的 SNP
"""

import requests
import pandas as pd
import os
from datetime import datetime

# JWT Token
JWT_TOKEN = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiIxNzU3ODgyODc4QHFxLmNvbSIsImlhdCI6MTc3ODMwNDA4MywiZXhwIjoxNzc5NTEzNjgzfQ.ZtcIUEx_xYtrVD_EE-UboKyLlC-lZBq2pjn-iYhJzxocqHdA-02K9n_Qbw-5ngQ07GHjjIYqVtmkZfJ3OJl1yI-tOMBBFzVKe0nkwDcB6-yBgjgBaxVm8vq_pbNrMwy_ZezY5ys9jx7I8T4bZYg9KeUbSwj04OfNP82kGcKXIOErXXVy-Ie3dbUogDRSjnCT-_32yNQxuWpiyYnPWSrWQbQ2HlUQiiDTdFGzWeJJKfSRvjQzdp5g3nccxht0m5A0UsPCdvkyHFEpvPVZ-NpjCjkgy8GbZBv4cmDMSc5JJL6HLO0eV508SRKxMdp-gL6qVdhGiJ2i9XmZQE27-aq_7g"

# MEGASTROKE 数据集 ID
DATASET_ID = "ebi-a-GCST006908"
OUTPUT_FILE = r"D:\EQTL\mr_results_megastroke\megastroke_complete.csv"

print("="*70)
print("从 IEU OpenGWAS 下载 MEGASTROKE 数据")
print("="*70)

# API 端点
API_BASE = "https://api.openGWAS.io"

# 设置请求头
headers = {
    "Authorization": f"Bearer {JWT_TOKEN}",
    "Content-Type": "application/json"
}

print(f"\n数据集：{DATASET_ID}")
print(f"输出文件：{OUTPUT_FILE}\n")

# 创建输出目录
os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)

try:
    # 使用分页下载
    all_variants = []
    page = 0
    page_size = 100000  # 每页 10 万条
    
    print("开始下载...")
    print("这可能需要几分钟...\n")
    
    while True:
        try:
            response = requests.get(
                f"{API_BASE}/api/v1/associations",
                headers=headers,
                params={
                    "dataset_id": DATASET_ID,
                    "page": page,
                    "page_size": page_size
                },
                timeout=300
            )
            
            if response.status_code != 200:
                print(f"  ✗ 第 {page} 页下载失败：{response.status_code}")
                break
            
            variants = response.json()
            
            if len(variants) == 0:
                print(f"\n✓ 完成！共下载 {len(all_variants):,} 个 SNP")
                break
            
            all_variants.extend(variants)
            page += 1
            
            if page % 10 == 0:
                print(f"  已下载 {len(all_variants):,} 个 SNP...")
            
            # 避免请求过快
            import time
            time.sleep(0.5)
            
        except Exception as e:
            print(f"  ✗ 错误：{e}")
            break
    
    # 保存数据
    if len(all_variants) > 0:
        print(f"\n保存数据到 {OUTPUT_FILE}...")
        
        # 转换为 DataFrame
        df = pd.DataFrame(all_variants)
        
        # 选择需要的列
        columns_to_keep = ['snp', 'chromosome', 'position', 'effect_allele', 'other_allele', 
                          'eaf', 'beta', 'se', 'pval', 'n']
        
        df = df[columns_to_keep]
        df.columns = ['SNP', 'CHR', 'BP', 'EFFECT_ALLELE', 'OTHER_ALLELE', 
                     'EAF', 'BETA', 'SE', 'PVAL', 'N']
        
        # 添加结局信息
        df['outcome'] = 'Ischemic Stroke'
        df['id.outcome'] = DATASET_ID
        
        # 保存为 CSV
        df.to_csv(OUTPUT_FILE, index=False, encoding='utf-8')
        
        print(f"  ✓ 已保存：{len(df):,} 个 SNP")
        print(f"  ✓ 文件大小：{os.path.getsize(OUTPUT_FILE) / 1e6:.2f} MB")
        
    else:
        print("\n✗ 未下载任何数据")
        
except Exception as e:
    print(f"\n✗ 错误：{e}")
    import traceback
    traceback.print_exc()

print("\n" + "="*70)
print("完成！")
print("="*70)
