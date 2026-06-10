#!/usr/bin/env python3
# ================================================================================
# 提取双源 eQTL 数据 - 简化版
# 不依赖 numpy，使用纯 Python + pyarrow
# 参考：GTEx v11 | eQTL Catalogue | TwoSampleMR
# ================================================================================

import os
import sys
from pathlib import Path
import time

try:
    import pyarrow as pa
    import pyarrow.parquet as pq
    print("✓ 使用 pyarrow 处理\n")
except ImportError as e:
    print(f"✗ 错误：{e}")
    sys.exit(1)

# 配置
GTEx_BRAIN_FILE = r"C:\Users\Jy-Mentor-7\Desktop\生物信息学\ETQL\Brain_Cortex.v11.eQTLs.signif_pairs.parquet"
GTEx_BLOOD_FILE = r"C:\Users\Jy-Mentor-7\Desktop\生物信息学\ETQL\Whole_Blood.v11.eQTLs.signif_pairs.parquet"
GTEx_BRAIN_EGENES = r"C:\Users\Jy-Mentor-7\Desktop\生物信息学\ETQL\Brain_Cortex.v11.eGenes.txt"
GTEx_BLOOD_EGENES = r"C:\Users\Jy-Mentor-7\Desktop\生物信息学\ETQL\Whole_Blood.v11.eGenes.txt"
OUTPUT_DIR = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\exposure_dual_source"

print("="*70)
print("双源 eQTL 数据提取 - pyarrow 版")
print("参考：GTEx v11 | eQTL Catalogue | TwoSampleMR")
print("="*70)

os.makedirs(OUTPUT_DIR, exist_ok=True)
print(f"\n输出目录：{OUTPUT_DIR}\n")

start_time = time.time()

# ================================================================================
# 1. 加载基因列表
# ================================================================================
print("步骤 1: 加载基因列表")
print("-"*70)

def load_genes(filename):
    """加载基因列表"""
    genes = set()
    with open(filename, 'r', encoding='utf-8') as f:
        # 跳过标题行
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if parts:
                genes.add(parts[0])  # gene_id 是第一列
    return genes

print("  加载脑皮层 eGenes...")
brain_genes = load_genes(GTEx_BRAIN_EGENES)
print(f"  ✓ 脑皮层：{len(brain_genes):,} 基因")

print("  加载全血 eGenes...")
blood_genes = load_genes(GTEx_BLOOD_EGENES)
print(f"  ✓ 全血：{len(blood_genes):,} 基因")

target_genes = brain_genes | blood_genes
print(f"\n目标总数：{len(target_genes):,}\n")

# ================================================================================
# 2. 流式处理 eQTL 数据
# ================================================================================
print("步骤 2: 流式处理 eQTL 数据")
print("-"*70)

def process_eqtl_streaming(input_file, target_genes, tissue_type, output_dir):
    """流式处理 eQTL 数据"""
    print(f"  处理 {tissue_type}...")
    
    parquet_file = pq.ParquetFile(input_file)
    
    gene_data = {}
    stats = {'processed': 0, 'saved': 0}
    batch_size = 100000
    
    for batch_idx, batch in enumerate(parquet_file.iter_batches(batch_size=batch_size)):
        # 转换为 pandas DataFrame
        batch_df = batch.to_pandas()
        
        # 从 phenotype_id 提取 gene_id
        batch_df['gene_id'] = batch_df['phenotype_id'].str.split('_').str[0]
        
        # 筛选目标基因
        filtered = batch_df[batch_df['gene_id'].isin(target_genes)]
        
        if len(filtered) == 0:
            stats['processed'] += len(batch_df)
            continue
        
        # 解析 variant_id
        variant_parts = filtered['variant_id'].str.split('_', expand=True, n=4)
        filtered['CHR'] = variant_parts[0].str.replace('chr', '', regex=False)
        filtered['BP'] = variant_parts[1].astype(int)
        filtered['REF'] = variant_parts[2]
        filtered['ALT'] = variant_parts[3]
        
        # 按基因分组并保存
        for gene_name, gene_group in filtered.groupby('gene_id'):
            top_eqtl = gene_group.nsmallest(50, 'pval_nominal')
            
            mr_format = pd.DataFrame({
                'SNP': top_eqtl['variant_id'],
                'CHR': top_eqtl['CHR'],
                'BP': top_eqtl['BP'],
                'EFFECT_ALLELE': top_eqtl['ALT'],
                'OTHER_ALLELE': top_eqtl['REF'],
                'BETA': top_eqtl['slope'],
                'SE': top_eqtl['slope_se'],
                'PVAL': top_eqtl['pval_nominal'],
                'EAF': top_eqtl['af'],
                'GENE': gene_name,
                'TISSUE': tissue_type
            })
            
            output_file = os.path.join(output_dir, f"{gene_name}_exposure.csv")
            
            if os.path.exists(output_file):
                existing = pd.read_csv(output_file)
                combined = pd.concat([existing, mr_format], ignore_index=True)
                combined.to_csv(output_file, index=False)
            else:
                mr_format.to_csv(output_file, index=False)
                stats['saved'] += 1
        
        stats['processed'] += len(batch_df)
        
        if (batch_idx + 1) % 10 == 0:
            print(f"    已处理 {stats['processed']:,} 行，保存 {stats['saved']:,} 个基因")
    
    return stats

# 尝试导入 pandas（只在这里使用）
try:
    import pandas as pd
    print("✓ pandas 可用\n")
except:
    print("✗ pandas 不可用，无法继续\n")
    sys.exit(1)

# 处理脑组织数据
brain_stats = process_eqtl_streaming(GTEx_BRAIN_FILE, brain_genes, "Brain_Cortex", OUTPUT_DIR)
print(f"  ✓ 脑皮层：处理 {brain_stats['processed']:,} 行，保存 {brain_stats['saved']:,} 基因\n")

# 处理全血数据
blood_stats = process_eqtl_streaming(GTEx_BLOOD_FILE, blood_genes, "Whole_Blood", OUTPUT_DIR)
print(f"  ✓ 全血：处理 {blood_stats['processed']:,} 行，保存 {blood_stats['saved']:,} 基因\n")

# ================================================================================
# 3. 统计结果
# ================================================================================
elapsed_time = time.time() - start_time

print("="*70)
print("处理完成")
print("="*70)

output_files = [f for f in os.listdir(OUTPUT_DIR) if f.endswith('_exposure.csv')]

print(f"\n总基因数：{len(output_files):,}")
print(f"处理时间：{elapsed_time:.2f} 秒 ({elapsed_time/60:.2f} 分钟)")
print(f"处理速度：{(brain_stats['processed'] + blood_stats['processed']) / elapsed_time / 1000:.1f}K 行/秒\n")

stats = {'both': 0, 'brain_only': 0, 'blood_only': 0}

for filename in output_files:
    filepath = os.path.join(OUTPUT_DIR, filename)
    data = pd.read_csv(filepath, usecols=['TISSUE'])
    tissues = data['TISSUE'].unique()
    
    has_brain = 'Brain_Cortex' in tissues
    has_blood = 'Whole_Blood' in tissues
    
    if has_brain and has_blood:
        stats['both'] += 1
    elif has_brain:
        stats['brain_only'] += 1
    elif has_blood:
        stats['blood_only'] += 1

print("组织分布:")
print(f"  - 双组织都有：{stats['both']:,} 基因")
print(f"  - 仅脑组织：  {stats['brain_only']:,} 基因")
print(f"  - 仅全血：    {stats['blood_only']:,} 基因\n")

with open(os.path.join(OUTPUT_DIR, "stats.txt"), 'w') as f:
    f.write(f"双源 eQTL 整合统计\n")
    f.write(f"参考：GTEx v11 | eQTL Catalogue\n\n")
    f.write(f"总基因数：{len(output_files):,}\n")
    f.write(f"处理时间：{elapsed_time:.2f} 秒\n\n")
    f.write(f"组织分布:\n")
    f.write(f"  - 双组织都有：{stats['both']:,}\n")
    f.write(f"  - 仅脑组织：  {stats['brain_only']:,}\n")
    f.write(f"  - 仅全血：    {stats['blood_only']:,}\n")

print("✓ 统计已保存\n")

if len(output_files) > 0:
    sample = pd.read_csv(os.path.join(OUTPUT_DIR, output_files[0]))
    print(f"样本预览 ({output_files[0]}):")
    print(f"  维度：{sample.shape[0]:,} × {sample.shape[1]}")
    print(f"  列：{', '.join(sample.columns)}\n")

print("="*70)
print("完成！下一步：Rscript run_mr_analysis.R")
print("="*70)
