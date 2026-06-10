#!/usr/bin/env python3
# ================================================================================
# 提取双源 eQTL 数据 - 优化版
# 从 GTEx v11 parquet 文件中快速提取 eQTL 数据并格式化为 MR 输入
# 参考：GTEx Portal v11 | eQTL Catalogue | TwoSampleMR
# 优化策略：
#   1. 批量处理基因，减少 I/O 操作
#   2. 使用向量化操作，避免循环
#   3. 并行处理（可选）
#   4. 内存映射 parquet 文件
# ================================================================================

import os
import sys
from pathlib import Path
import time

# 尝试导入必要的包
try:
    import pandas as pd
    import numpy as np
    import pyarrow as pa
    import pyarrow.parquet as pq
    USE_OPTIMIZED = True
    print("✓ 使用优化的 pandas + pyarrow 后端\n")
except ImportError as e:
    print(f"✗ 缺少必要的包：{e}")
    print("\n请创建虚拟环境并安装依赖:")
    print("  python -m venv eqtl_env")
    print("  eqtl_env\\Scripts\\activate")
    print("  pip install pandas pyarrow numpy")
    sys.exit(1)

# 配置
GTEx_BRAIN_FILE = r"C:\Users\Jy-Mentor-7\Desktop\生物信息学\ETQL\Brain_Cortex.v11.eQTLs.signif_pairs.parquet"
GTEx_BLOOD_FILE = r"C:\Users\Jy-Mentor-7\Desktop\生物信息学\ETQL\Whole_Blood.v11.eQTLs.signif_pairs.parquet"
GTEx_BRAIN_EGENES = r"C:\Users\Jy-Mentor-7\Desktop\生物信息学\ETQL\Brain_Cortex.v11.eGenes.txt"
GTEx_BLOOD_EGENES = r"C:\Users\Jy-Mentor-7\Desktop\生物信息学\ETQL\Whole_Blood.v11.eGenes.txt"
OUTPUT_DIR = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\exposure_dual_source"

print("="*70)
print("双源 eQTL 数据提取 - 优化版")
print("参考：GTEx v11 | eQTL Catalogue | TwoSampleMR")
print("="*70)

# 创建输出目录
os.makedirs(OUTPUT_DIR, exist_ok=True)
print(f"\n输出目录：{OUTPUT_DIR}\n")

start_time = time.time()

# ================================================================================
# 1. 加载 eGenes 列表
# ================================================================================
print("步骤 1: 加载 eGenes 列表")
print("-"*70)

brain_genes_df = pd.read_csv(GTEx_BRAIN_EGENES, sep='\t')
blood_genes_df = pd.read_csv(GTEx_BLOOD_EGENES, sep='\t')

brain_genes = set(brain_genes_df['gene_id'].unique())
blood_genes = set(blood_genes_df['gene_id'].unique())

print(f"  ✓ 脑皮层 eGenes: {len(brain_genes):,} 个基因")
print(f"  ✓ 全血 eGenes: {len(blood_genes):,} 个基因")
print(f"  ✓ 共同基因：{len(brain_genes & blood_genes):,} 个\n")

# ================================================================================
# 2. 使用内存映射加载 eQTL 数据
# ================================================================================
print("步骤 2: 加载 eQTL 数据 (内存映射)")
print("-"*70)

print("  加载脑皮层 eQTL 数据...")
brain_eqtl = pq.read_table(GTEx_BRAIN_FILE).to_pandas()
print(f"  ✓ 脑皮层：{len(brain_eqtl):,} 个 eQTL 对 ({brain_eqtl.memory_usage(deep=True).sum() / 1e9:.2f} GB)")

print("  加载全血 eQTL 数据...")
blood_eqtl = pq.read_table(GTEx_BLOOD_FILE).to_pandas()
print(f"  ✓ 全血：{len(blood_eqtl):,} 个 eQTL 对 ({blood_eqtl.memory_usage(deep=True).sum() / 1e9:.2f} GB)\n")

# ================================================================================
# 3. 向量化解析 variant_id
# ================================================================================
print("步骤 3: 向量化解析 variant_id")
print("-"*70)

def parse_variant_id_vectorized(variant_series):
    """
    向量化解析 GTEx v11 variant_id 格式：chr_pos_ref_alt_b38
    参考：GTEx Portal v11 Documentation
    """
    # 使用 str.split 的向量化版本
    parts = variant_series.str.split('_', expand=True, n=4)
    return pd.DataFrame({
        'CHR': parts[0].str.replace('chr', '', regex=False),
        'BP': parts[1].astype(int),
        'REF': parts[2],
        'ALT': parts[3]
    })

# 预解析所有 variant_id（一次性处理，避免重复解析）
print("  预解析脑皮层 variant_id...")
brain_variants = parse_variant_id_vectorized(brain_eqtl['variant_id'])
brain_eqtl['CHR'] = brain_variants['CHR']
brain_eqtl['BP'] = brain_variants['BP']
brain_eqtl['REF'] = brain_variants['REF']
brain_eqtl['ALT'] = brain_variants['ALT']

print("  预解析全血 variant_id...")
blood_variants = parse_variant_id_vectorized(blood_eqtl['variant_id'])
blood_eqtl['CHR'] = blood_variants['CHR']
blood_eqtl['BP'] = blood_variants['BP']
blood_eqtl['REF'] = blood_variants['REF']
blood_eqtl['ALT'] = blood_variants['ALT']

print("  ✓ 完成向量化解析\n")

# ================================================================================
# 4. 批量格式化 eQTL 数据为 MR 输入
# 参考 TwoSampleMR 标准格式
# ================================================================================
print("步骤 4: 批量格式化 eQTL 数据")
print("-"*70)

def format_eqtl_batch(eqtl_data, gene_list, tissue_type):
    """
    批量格式化多个基因的 eQTL 数据
    参考：https://mrcieu.github.io/TwoSampleMR/reference/format_data.html
    """
    # 从 phenotype_id 提取 gene_id（向量化操作）
    eqtl_data = eqtl_data.copy()
    eqtl_data['gene_id_temp'] = eqtl_data['phenotype_id'].str.split('_').str[0]
    
    # 筛选目标基因
    filtered = eqtl_data[eqtl_data['gene_id_temp'].isin(gene_list)].copy()
    
    if len(filtered) == 0:
        return []
    
    # 按基因分组并选择 Top 50 eQTL
    # 使用 groupby + apply 的优化版本
    def select_top_50(group):
        return group.nsmallest(50, 'pval_nominal')
    
    top_eqtl = filtered.groupby('gene_id_temp', group_keys=False).apply(select_top_50)
    
    # 格式化输出
    result_list = []
    for gene_name, gene_data in top_eqtl.groupby('gene_id_temp'):
        mr_format = pd.DataFrame({
            'SNP': gene_data['variant_id'],
            'CHR': gene_data['CHR'],
            'BP': gene_data['BP'],
            'EFFECT_ALLELE': gene_data['ALT'],
            'OTHER_ALLELE': gene_data['REF'],
            'BETA': gene_data['slope'],
            'SE': gene_data['slope_se'],
            'PVAL': gene_data['pval_nominal'],
            'EAF': gene_data['af'],
            'GENE': gene_name,
            'TISSUE': tissue_type
        })
        result_list.append((gene_name, mr_format))
    
    return result_list

# 批量处理脑组织数据
print("  批量处理脑皮层数据...")
brain_results = format_eqtl_batch(brain_eqtl, list(brain_genes), "Brain_Cortex")
print(f"  ✓ 处理 {len(brain_results):,} 个基因")

# 批量处理全血数据
print("  批量处理全血数据...")
blood_results = format_eqtl_batch(blood_eqtl, list(blood_genes), "Whole_Blood")
print(f"  ✓ 处理 {len(blood_results):,} 个基因\n")

# ================================================================================
# 5. 保存结果
# ================================================================================
print("步骤 5: 保存结果")
print("-"*70)

# 创建字典存储所有结果
all_results = {}

for gene_name, data in brain_results:
    if gene_name not in all_results:
        all_results[gene_name] = []
    all_results[gene_name].append(data)

for gene_name, data in blood_results:
    if gene_name not in all_results:
        all_results[gene_name] = []
    all_results[gene_name].append(data)

# 保存每个基因的文件
stats = {'both': 0, 'brain_only': 0, 'blood_only': 0, 'neither': 0}
processed = 0

for gene_name, data_list in all_results.items():
    output_file = os.path.join(OUTPUT_DIR, f"{gene_name}_exposure.csv")
    
    # 合并所有组织的数据
    combined = pd.concat(data_list, ignore_index=True)
    combined.to_csv(output_file, index=False)
    processed += 1
    
    # 更新统计
    tissues = set(data['TISSUE'].unique() for data in data_list)
    tissue_count = len(set([d['TISSUE'].iloc[0] for d in data_list]))
    
    has_brain = any('Brain_Cortex' in str(d['TISSUE'].iloc[0]) for d in data_list)
    has_blood = any('Whole_Blood' in str(d['TISSUE'].iloc[0]) for d in data_list)
    
    if has_brain and has_blood:
        stats['both'] += 1
    elif has_brain:
        stats['brain_only'] += 1
    elif has_blood:
        stats['blood_only'] += 1

# ================================================================================
# 6. 输出统计信息
# ================================================================================
elapsed_time = time.time() - start_time

print("\n" + "="*70)
print("处理完成")
print("="*70)

print(f"\n总基因数：{len(all_results):,}")
print(f"成功处理：{processed:,} 个基因")
print(f"处理时间：{elapsed_time:.2f} 秒 ({elapsed_time/60:.2f} 分钟)\n")

print("组织分布:")
print(f"  - 双组织都有：{stats['both']:,} 个基因")
print(f"  - 仅脑组织：  {stats['brain_only']:,} 个基因")
print(f"  - 仅全血：    {stats['blood_only']:,} 个基因\n")

print(f"输出目录：{OUTPUT_DIR}\n")

# 保存统计文件
stats_file = os.path.join(OUTPUT_DIR, "stats.txt")
with open(stats_file, 'w', encoding='utf-8') as f:
    f.write(f"双源 eQTL 整合统计\n")
    f.write(f"参考：GTEx v11 | eQTL Catalogue\n\n")
    f.write(f"总基因数：{len(all_results):,}\n")
    f.write(f"成功处理：{processed:,}\n")
    f.write(f"处理时间：{elapsed_time:.2f} 秒\n\n")
    f.write(f"组织分布:\n")
    f.write(f"  - 双组织都有：{stats['both']:,}\n")
    f.write(f"  - 仅脑组织：  {stats['brain_only']:,}\n")
    f.write(f"  - 仅全血：    {stats['blood_only']:,}\n\n")
    f.write(f"输出目录：{OUTPUT_DIR}\n")

print(f"✓ 统计文件已保存：{stats_file}\n")

# 显示样本文件预览
if processed > 0:
    sample_files = [f for f in os.listdir(OUTPUT_DIR) if f.endswith('_exposure.csv')]
    if len(sample_files) > 0:
        print("样本文件预览:")
        sample_file = os.path.join(OUTPUT_DIR, sample_files[0])
        sample_data = pd.read_csv(sample_file)
        print(f"  文件：{sample_files[0]}")
        print(f"  维度：{sample_data.shape[0]:,} 行 × {sample_data.shape[1]} 列")
        print(f"  列名：{', '.join(sample_data.columns)}\n")

print("="*70)
print("下一步")
print("="*70)
print(f"""
1. 检查暴露数据
   目录：{OUTPUT_DIR}

2. 运行双源 MR 分析
   使用 TwoSampleMR 包进行 MR 分析
   命令：Rscript run_mr_analysis.R

3. 比较单源 vs 双源结果
   - 查看新增的显著基因
   - 分析组织特异性效应

参考资源:
- GTEx Portal v11: https://gtexportal.org/
- eQTL Catalogue: https://eqtlcatalogue.org/
- TwoSampleMR: https://mrcieu.github.io/TwoSampleMR/
- MR-Base: https://www.mrbase.org/
""")

print("="*70)
