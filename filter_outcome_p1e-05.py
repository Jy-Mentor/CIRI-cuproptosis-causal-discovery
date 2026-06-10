#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
为p=1e-05数据过滤结局文件
"""

import pandas as pd
import time
import os

# ============================================
# 参数设置
# ============================================
outcome_file = "D:/EQTL/eqtlgen_ieu_outcome.csv"
output_file = "D:/EQTL/mr_results_p1e-05/outcome_filtered.csv"
chr_pos_file = "D:/EQTL/mr_results_p1e-05/exposure_chr_pos_list.txt"

print("=" * 60)
print("p=1e-05 结局数据过滤")
print("=" * 60)

# ============================================
# 步骤1: 读取chr:pos列表
# ============================================
print("\n步骤1: 读取chr:pos列表...")

with open(chr_pos_file, 'r') as f:
    chr_pos_list = set(line.strip() for line in f if line.strip())

print(f"  需要匹配的chr:pos数: {len(chr_pos_list)}")

# ============================================
# 步骤2: 分块过滤结局数据
# ============================================
print("\n步骤2: 过滤结局数据...")
print(f"  结局文件: {outcome_file}")
print(f"  文件大小: {os.path.getsize(outcome_file) / 1024**3:.2f} GB")

start = time.time()
filtered_chunks = []
chunk_count = 0
matched_rows = 0

# 分块读取结局文件
chunks = pd.read_csv(
    outcome_file,
    chunksize=1000000,  # 每次100万行
    low_memory=False,
    dtype=str  # 避免类型推断
)

for i, chunk in enumerate(chunks):
    chunk_count += 1
    
    # 创建chr:pos列用于匹配
    chunk['chr_pos'] = chunk['chr.outcome'] + ':' + chunk['pos.outcome']
    
    # 过滤匹配的行
    matched = chunk[chunk['chr_pos'].isin(chr_pos_list)]
    
    if len(matched) > 0:
        # 删除临时列
        matched = matched.drop(columns=['chr_pos'])
        filtered_chunks.append(matched)
        matched_rows += len(matched)
        print(f"  块{i}: 匹配到 {len(matched)} 行, 累计 {matched_rows} 行")
    
    if i % 10 == 0 and i > 0:
        elapsed = time.time() - start
        print(f"  已处理 {i}00万行, 耗时 {elapsed:.0f}秒, 累计匹配 {matched_rows} 行")

print(f"\n  总块数: {chunk_count}")

# 合并所有匹配的块
if filtered_chunks:
    result = pd.concat(filtered_chunks, ignore_index=True)
    print(f"  过滤后总行数: {len(result)}")
    
    # 确保输出目录存在
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # 保存结果
    result.to_csv(output_file, index=False)
    print(f"  结果已保存: {output_file}")
else:
    print("  警告: 未匹配到任何行!")
    result = None

elapsed_total = time.time() - start
print(f"\n总耗时: {elapsed_total/60:.1f} 分钟")
print("=" * 60)
