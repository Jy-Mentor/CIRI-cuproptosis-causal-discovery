#!/usr/bin/env python3
# ================================================================================
# 从 IEU OpenGWAS 下载 MEGASTROKE 完整数据
# 使用提供的 JWT token 进行认证
# ================================================================================

import os
import sys
import requests
import time

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

# 步骤 1: 获取数据集信息
print("步骤 1: 获取数据集信息...")
try:
    response = requests.get(
        f"{API_BASE}/api/v1/datasets/{DATASET_ID}",
        headers=headers,
        timeout=30
    )
    
    if response.status_code == 200:
        dataset_info = response.json()
        print(f"  ✓ 数据集名称：{dataset_info.get('trait', 'N/A')}")
        print(f"  ✓ 样本量：{dataset_info.get('n', 'N/A'):,}")
        print(f"  ✓ SNP 数量：{dataset_info.get('nsnp', 'N/A'):,}")
    else:
        print(f"  ✗ 获取失败：{response.status_code}")
        print(f"  响应：{response.text[:200]}")
except Exception as e:
    print(f"  ✗ 错误：{e}")

# 步骤 2: 下载完整数据
print("\n步骤 2: 下载完整关联数据...")
print("  这可能需要几分钟...")

try:
    # 使用分页下载
    all_variants = []
    page = 0
    page_size = 100000  # 每页 10 万条
    
    while True:
        response = requests.get(
            f"{API_BASE}/api/v1/associations",
            headers=headers,
            params={
                "dataset_id": DATASET_ID,
                "page": page,
                "page_size": page_size
            },
            timeout=120
        )
        
        if response.status_code != 200:
            print(f"  ✗ 第 {page} 页下载失败：{response.status_code}")
            break
        
        variants = response.json()
        
        if len(variants) == 0:
            print(f"  ✓ 完成！共下载 {len(all_variants):,} 个 SNP")
            break
        
        all_variants.extend(variants)
        page += 1
        
        if page % 10 == 0:
            print(f"    已下载 {len(all_variants):,} 个 SNP...")
        
        # 避免请求过快
        time.sleep(0.5)
    
    # 保存数据
    if len(all_variants) > 0:
        print(f"\n  保存数据到 {OUTPUT_FILE}...")
        
        # 创建输出目录
        os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
        
        # 写入 CSV
        with open(OUTPUT_FILE, 'w', encoding='utf-8') as f:
            # 写入表头
            f.write("SNP,chr,pos.outcome,effect_allele.outcome,other_allele.outcome,eaf.outcome,beta.outcome,se.outcome,pval.outcome,samplesize.outcome\n")
            
            # 写入数据
            for variant in all_variants:
                snp = variant.get('snp', '')
                chr_ = variant.get('chromosome', '')
                pos = variant.get('position', '')
                ea = variant.get('effect_allele', '')
                oa = variant.get('other_allele', '')
                eaf = variant.get('eaf', '')
                beta = variant.get('beta', '')
                se = variant.get('se', '')
                pval = variant.get('pval', '')
                n = variant.get('n', '')
                
                f.write(f"{snp},{chr_},{pos},{ea},{oa},{eaf},{beta},{se},{pval},{n}\n")
        
        print(f"  ✓ 已保存：{len(all_variants):,} 个 SNP")
        print(f"  ✓ 文件大小：{os.path.getsize(OUTPUT_FILE) / 1e6:.2f} MB")
        
except Exception as e:
    print(f"  ✗ 错误：{e}")
    import traceback
    traceback.print_exc()

print("\n" + "="*70)
print("下载完成！")
print("="*70)

print(f"\n下一步:")
print(f"  1. 检查数据文件：{OUTPUT_FILE}")
print(f"  2. 运行 MR 分析：Rscript run_mr_analysis_matched.R")
