#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
为GTEx v11数据过滤结局文件
提取GTEx暴露SNP的chr:pos，然后过滤结局数据
"""

import pandas as pd
import pyarrow.parquet as pq
import time
import re

# ============================================
# 参数设置
# ============================================
signif_pairs_file = "C:/Users/Jy-Mentor-7/Desktop/生物信息学/ETQL/Whole_Blood.v11.eQTLs.signif_pairs.parquet"
outcome_file = "D:/EQTL/eqtlgen_ieu_outcome.csv"
output_file = "D:/EQTL/mr_results_gtex_v11/outcome_filtered_gtex.csv"

# 目标基因
genes = ["NFKB1", "FDX1", "STAT3", "HIF1A", "HMOX1", "GPX4", "HSPA5", "AGER", "DLAT"]

# ENSG ID到基因名映射
ensg_to_gene = {
    "ENSG00000109320": "NFKB1",
    "ENSG00000204305": "AGER",
    "ENSG00000044574": "HSPA5",
    "ENSG00000137714": "FDX1",
    "ENSG00000150768": "HIF1A",
    "ENSG00000100644": "STAT3",
    "ENSG00000168610": "GPX4",
    "ENSG00000167468": "HMOX1",
    "ENSG00000100292": "DLAT"
}

print("=" * 60)
print("GTEx v11 结局数据过滤")
print("=" * 60)

# ============================================
# 步骤1: 读取signif_pairs.parquet
# ============================================
print("\n步骤1: 读取signif_pairs.parquet...")
start = time.time()

table = pq.read_table(signif_pairs_file)
signif_pairs = table.to_pandas()

print(f"  总行数: {len(signif_pairs)}")
print(f"  列名: {list(signif_pairs.columns)}")
print(f"  耗时: {time.time()-start:.1f}秒")

# ============================================
# 步骤2: 从phenotype_id提取基因名并筛选目标基因
# ============================================
print("\n步骤2: 筛选目标基因SNP...")

# 提取基因基础ID（去除版本号）
signif_pairs['gene_base'] = signif_pairs['phenotype_id'].str.replace(r'\.\d+$', '', regex=True)

# 映射基因名
signif_pairs['gene_name'] = signif_pairs['gene_base'].map(ensg_to_gene)

# 筛选目标基因
target_snps = signif_pairs[signif_pairs['gene_name'].isin(genes)].copy()
print(f"  目标基因SNP数: {len(target_snps)}")
print(f"  找到基因: {target_snps['gene_name'].unique().tolist()}")

# ============================================
# 步骤3: 从variant_id解析chr和pos
# ============================================
print("\n步骤3: 解析variant_id提取chr:pos...")

# variant_id格式: chr1_665098_G_A_b38
def parse_variant(variant_id):
    """解析variant_id获取chr和pos"""
    parts = variant_id.split('_')
    if len(parts) >= 4:
        chr_num = parts[0].replace('chr', '')
        pos = parts[1]
        return f"{chr_num}:{pos}"
    return None

target_snps['chr_pos'] = target_snps['variant_id'].apply(parse_variant)

# 提取唯一的chr:pos用于匹配
unique_chr_pos = set(target_snps['chr_pos'].dropna())
print(f"  唯一chr:pos数: {len(unique_chr_pos)}")

# 保存chr:pos列表用于调试
with open("D:/EQTL/mr_results_gtex_v11/gtex_chr_pos_list.txt", 'w') as f:
    for cp in sorted(unique_chr_pos):
        f.write(f"{cp}\n")
print(f"  chr:pos列表已保存")

# ============================================
# 步骤4: 分块过滤结局数据
# ============================================
print("\n步骤4: 过滤结局数据...")
print(f"  结局文件: {outcome_file}")

start = time.time()
filtered_chunks = []
chunk_count = 0

# 分块读取结局文件
chunks = pd.read_csv(
    outcome_file,
    chunksize=2000000,  # 每次200万行
    low_memory=False,
    dtype=str  # 避免类型推断
)

for i, chunk in enumerate(chunks):
    # 构建chr:pos列用于匹配
    chunk['chr_pos'] = chunk['chr.outcome'] + ':' + chunk['pos.outcome']
    
    # 过滤匹配的行
    matched = chunk[chunk['chr_pos'].isin(unique_chr_pos)]
    
    if len(matched) > 0:
        # 删除临时列
        matched = matched.drop(columns=['chr_pos'])
        filtered_chunks.append(matched)
        print(f"  块{i}: 匹配到 {len(matched)} 行")
    
    chunk_count += 1
    if i % 5 == 0:
        elapsed = time.time() - start
        print(f"  已处理 {i*2}00万行, 耗时 {elapsed:.0f}秒")

print(f"\n  总块数: {chunk_count}")

# 合并所有匹配的块
if filtered_chunks:
    result = pd.concat(filtered_chunks, ignore_index=True)
    print(f"  过滤后总行数: {len(result)}")
    
    # 保存结果
    import os
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    result.to_csv(output_file, index=False)
    print(f"  结果已保存: {output_file}")
else:
    print("  警告: 未匹配到任何行!")

elapsed_total = time.time() - start
print(f"\n总耗时: {elapsed_total/60:.1f} 分钟")
print("=" * 60)
