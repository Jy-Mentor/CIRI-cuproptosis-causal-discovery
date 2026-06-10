#!/usr/bin/env python3
# ================================================================================
# 提取双源 eQTL 数据 - 纯 Python 实现
# 从 GTEx v11 parquet 文件中提取 eQTL 数据并格式化为 MR 输入
# 参考：GTEx Portal v11 | eQTL Catalogue | TwoSampleMR
# 不依赖 numpy/pandas，避免 DLL 问题
# ================================================================================

import os
import sys
import csv
from pathlib import Path

# 尝试导入 pandas，如果失败则使用纯 Python
try:
    import pandas as pd
    USE_PANDAS = True
except ImportError:
    print("警告：pandas 不可用，使用纯 Python 实现\n")
    USE_PANDAS = False

# 配置
GTEx_BRAIN_FILE = r"C:\Users\Jy-Mentor-7\Desktop\生物信息学\ETQL\Brain_Cortex.v11.eQTLs.signif_pairs.parquet"
GTEx_BLOOD_FILE = r"C:\Users\Jy-Mentor-7\Desktop\生物信息学\ETQL\Whole_Blood.v11.eQTLs.signif_pairs.parquet"
GTEx_BRAIN_EGENES = r"C:\Users\Jy-Mentor-7\Desktop\生物信息学\ETQL\Brain_Cortex.v11.eGenes.txt"
GTEx_BLOOD_EGENES = r"C:\Users\Jy-Mentor-7\Desktop\生物信息学\ETQL\Whole_Blood.v11.eGenes.txt"
OUTPUT_DIR = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\exposure_dual_source"

print("="*70)
print("双源 eQTL 数据提取")
print("参考：GTEx v11 | eQTL Catalogue | TwoSampleMR")
print("="*70)

# 创建输出目录
os.makedirs(OUTPUT_DIR, exist_ok=True)
print(f"\n输出目录：{OUTPUT_DIR}\n")

# ================================================================================
# 1. 加载 eGenes 列表
# ================================================================================
print("步骤 1: 加载 eGenes 列表")
print("-"*70)

def load_genes_file(filename):
    """加载 eGenes 文件并返回基因 ID 集合"""
    genes = set()
    with open(filename, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            genes.add(row['gene_id'])
    return genes

brain_genes = load_genes_file(GTEx_BRAIN_EGENES)
blood_genes = load_genes_file(GTEx_BLOOD_EGENES)

print(f"  ✓ 脑皮层 eGenes: {len(brain_genes):,} 个基因")
print(f"  ✓ 全血 eGenes: {len(blood_genes):,} 个基因")
print(f"  ✓ 共同基因：{len(brain_genes & blood_genes):,} 个\n")

# ================================================================================
# 2. 加载 eQTL 数据
# ================================================================================
print("步骤 2: 加载 eQTL 数据")
print("-"*70)

def load_parquet_simple(filename):
    """
    简单加载 parquet 文件
    使用 pyarrow 的 pandas 接口
    """
    try:
        import pyarrow.parquet as pq
        table = pq.read_table(filename)
        return table.to_pandas()
    except Exception as e:
        print(f"  ✗ 加载失败：{e}")
        return None

print("  加载脑皮层 eQTL 数据...")
brain_eqtl = load_parquet_simple(GTEx_BRAIN_FILE)
if brain_eqtl is not None:
    print(f"  ✓ 脑皮层：{len(brain_eqtl):,} 个 eQTL 对")

print("  加载全血 eQTL 数据...")
blood_eqtl = load_parquet_simple(GTEx_BLOOD_FILE)
if blood_eqtl is not None:
    print(f"  ✓ 全血：{len(blood_eqtl):,} 个 eQTL 对\n")

# ================================================================================
# 3. 解析 variant_id
# ================================================================================
def parse_variant_id(variant_id):
    """
    解析 GTEx v11 variant_id 格式：chr_pos_ref_alt_b38
    参考：GTEx Portal v11 Documentation
    """
    try:
        parts = str(variant_id).split('_')
        if len(parts) >= 4:
            return {
                'CHR': parts[0].replace('chr', ''),
                'BP': int(parts[1]),
                'REF': parts[2],
                'ALT': parts[3]
            }
    except:
        pass
    return {'CHR': 'NA', 'BP': 0, 'REF': 'NA', 'ALT': 'NA'}

# ================================================================================
# 4. 格式化 eQTL 数据为 MR 输入
# 参考 TwoSampleMR 标准格式
# ================================================================================
print("步骤 3: 格式化 eQTL 数据为 MR 输入")
print("-"*70)

def format_eqtl_for_mr(eqtl_data, gene_name, tissue_type):
    """
    格式化 eQTL 数据为 MR 分析输入格式
    参考：https://mrcieu.github.io/TwoSampleMR/reference/format_data.html
    """
    if eqtl_data is None or len(eqtl_data) == 0:
        return None
    
    # GTEx v11 使用 phenotype_id 格式：geneID_position_strand
    # 需要从中提取 gene_id
    gene_data = eqtl_data[
        eqtl_data['phenotype_id'].str.split('_').str[0] == gene_name
    ].copy()
    
    if len(gene_data) == 0:
        return None
    
    # 选择最强的 eQTL (最低 P 值) - 限制 Top 50
    # 参考 GTEx 推荐做法和 MR-Base 标准
    gene_data = gene_data.sort_values('pval_nominal').head(50)
    
    # 解析 variant_id
    variant_info = gene_data['variant_id'].apply(parse_variant_id)
    gene_data['CHR'] = [x['CHR'] for x in variant_info]
    gene_data['BP'] = [x['BP'] for x in variant_info]
    gene_data['REF'] = [x['REF'] for x in variant_info]
    gene_data['ALT'] = [x['ALT'] for x in variant_info]
    
    # 格式化输出 (TwoSampleMR 标准格式)
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
    
    return mr_format

# ================================================================================
# 5. 处理所有基因
# ================================================================================
print("步骤 4: 处理所有基因")
print("-"*70)

# 目标基因列表（所有在任一组织中表达的基因）
target_genes = list(brain_genes | blood_genes)
print(f"目标基因总数：{len(target_genes):,}\n")

stats = {'both': 0, 'brain_only': 0, 'blood_only': 0, 'neither': 0}
processed = 0

for i, gene in enumerate(target_genes, 1):
    output_file = os.path.join(OUTPUT_DIR, f"{gene}_exposure.csv")
    
    in_brain = gene in brain_genes
    in_blood = gene in blood_genes
    
    all_data = []
    
    # 添加脑组织数据
    if in_brain and brain_eqtl is not None:
        brain_data = format_eqtl_for_mr(brain_eqtl, gene, "Brain_Cortex")
        if brain_data is not None and len(brain_data) > 0:
            all_data.append(brain_data)
    
    # 添加全血数据
    if in_blood and blood_eqtl is not None:
        blood_data = format_eqtl_for_mr(blood_eqtl, gene, "Whole_Blood")
        if blood_data is not None and len(blood_data) > 0:
            all_data.append(blood_data)
    
    # 保存
    if len(all_data) > 0:
        # 合并所有组织的数据
        combined = pd.concat(all_data, ignore_index=True)
        combined.to_csv(output_file, index=False)
        processed += 1
        
        # 更新统计
        if in_brain and in_blood:
            stats['both'] += 1
        elif in_brain:
            stats['brain_only'] += 1
        elif in_blood:
            stats['blood_only'] += 1
    else:
        stats['neither'] += 1
    
    # 进度显示
    if i % 100 == 0:
        print(f"  进度：{i:,}/{len(target_genes):,} (已处理：{processed:,})")

# ================================================================================
# 6. 输出统计信息
# ================================================================================
print("\n" + "="*70)
print("处理完成")
print("="*70)

print(f"\n总基因数：{len(target_genes):,}")
print(f"成功处理：{processed:,} 个基因\n")

print("组织分布:")
print(f"  - 双组织都有：{stats['both']:,} 个基因")
print(f"  - 仅脑组织：  {stats['brain_only']:,} 个基因")
print(f"  - 仅全血：    {stats['blood_only']:,} 个基因")
print(f"  - 都无数据：  {stats['neither']:,} 个基因\n")

print(f"输出目录：{OUTPUT_DIR}\n")

# 保存统计文件
stats_file = os.path.join(OUTPUT_DIR, "stats.txt")
with open(stats_file, 'w', encoding='utf-8') as f:
    f.write(f"双源 eQTL 整合统计\n")
    f.write(f"参考：GTEx v11 | eQTL Catalogue\n\n")
    f.write(f"总基因数：{len(target_genes):,}\n")
    f.write(f"成功处理：{processed:,}\n\n")
    f.write(f"组织分布:\n")
    f.write(f"  - 双组织都有：{stats['both']:,}\n")
    f.write(f"  - 仅脑组织：  {stats['brain_only']:,}\n")
    f.write(f"  - 仅全血：    {stats['blood_only']:,}\n")
    f.write(f"  - 都无数据：  {stats['neither']:,}\n\n")
    f.write(f"输出目录：{OUTPUT_DIR}\n")

print(f"✓ 统计文件已保存：{stats_file}\n")

# 显示样本文件预览
if processed > 0:
    sample_files = [f for f in os.listdir(OUTPUT_DIR) if f.endswith('_exposure.csv')]
    if len(sample_files) > 0:
        print("样本文件预览:")
        sample_file = os.path.join(OUTPUT_DIR, sample_files[0])
        with open(sample_file, 'r', encoding='utf-8') as f:
            reader = csv.reader(f)
            header = next(reader)
            row_count = sum(1 for _ in reader) + 1
        print(f"  文件：{sample_files[0]}")
        print(f"  维度：{row_count:,} 行 × {len(header)} 列")
        print(f"  列名：{', '.join(header)}\n")

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
