#!/usr/bin/env python3
# ================================================================================
# 准备 GTEx v11 脑组织 eQTL 数据
# 处理本地 GTEx v11 parquet 文件并格式化为 MR 输入
# ================================================================================

import os
import sys
import pandas as pd
from pathlib import Path

# 配置 - 使用用户实际数据路径
GTEx_BRAIN_FILE = r"C:\Users\Jy-Mentor-7\Desktop\生物信息学\ETQL\Brain_Cortex.v11.eQTLs.signif_pairs.parquet"
GTEx_BLOOD_FILE = r"C:\Users\Jy-Mentor-7\Desktop\生物信息学\ETQL\Whole_Blood.v11.eQTLs.signif_pairs.parquet"
OUTPUT_DIR = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\exposure_gtex_brain"
GENE_LIST_FILE = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\gene_list_optimized.txt"

print("="*70)
print("准备 GTEx v11 脑组织 eQTL 数据")
print("="*70)

# 创建输出目录
os.makedirs(OUTPUT_DIR, exist_ok=True)
print(f"\n输出目录：{OUTPUT_DIR}")

# ================================================================================
# 加载基因列表
# ================================================================================

gene_list = []
if os.path.exists(GENE_LIST_FILE):
    with open(GENE_LIST_FILE, 'r', encoding='utf-8') as f:
        gene_list = [line.strip() for line in f if line.strip()]
    print(f"✓ 加载基因列表：{len(gene_list)} 个基因")
else:
    print(f"✗ 基因列表文件不存在：{GENE_LIST_FILE}")
    print("  将使用 eQTL 数据中的所有基因")

# ================================================================================
# 加载 GTEx 数据
# ================================================================================

print("\n" + "="*70)
print("加载 GTEx v11 数据")
print("="*70)

# 加载脑组织数据
if os.path.exists(GTEx_BRAIN_FILE):
    print(f"\n加载脑皮层 eQTL 数据...")
    try:
        brain_eqtl = pd.read_parquet(GTEx_BRAIN_FILE)
        print(f"  ✓ 脑皮层数据：{len(brain_eqtl):,} 个 eQTL 对")
        print(f"    列：{list(brain_eqtl.columns)}")
    except Exception as e:
        print(f"  ✗ 加载失败：{e}")
        brain_eqtl = pd.DataFrame()
else:
    print(f"  ✗ 文件不存在：{GTEx_BRAIN_FILE}")
    brain_eqtl = pd.DataFrame()

# 加载全血数据
if os.path.exists(GTEx_BLOOD_FILE):
    print(f"\n加载全血 eQTL 数据...")
    try:
        blood_eqtl = pd.read_parquet(GTEx_BLOOD_FILE)
        print(f"  ✓ 全血数据：{len(blood_eqtl):,} 个 eQTL 对")
    except Exception as e:
        print(f"  ✗ 加载失败：{e}")
        blood_eqtl = pd.DataFrame()
else:
    print(f"  ✗ 文件不存在：{GTEx_BLOOD_FILE}")
    blood_eqtl = pd.DataFrame()

# ================================================================================
# 格式化数据为 MR 输入
# ================================================================================

def format_eqtl_for_mr(eqtl_data, gene_name, tissue_type):
    """格式化 eQTL 数据为 MR 分析输入格式"""
    
    if len(eqtl_data) == 0:
        return pd.DataFrame()
    
    # 筛选该基因的 eQTL
    # phenotype_id 就是基因的 Ensembl ID (格式: ENSG000XXXXX.Y)
    gene_data = eqtl_data[eqtl_data['phenotype_id'] == gene_name].copy()
    
    if len(gene_data) == 0:
        return pd.DataFrame()
    
    # 选择最强的 eQTL (最低 P 值) - 限制 Top 50
    gene_data = gene_data.sort_values('pval_nominal').head(50)
    
    # 解析 variant_id 获取染色体和位置
    # variant_id 格式：chr_pos_ref_alt_b38
    def parse_variant_id(variant_id):
        try:
            parts = str(variant_id).split('_')
            if len(parts) >= 4:
                chr_ = parts[0].replace('chr', '')
                pos = int(parts[1])
                ref = parts[2]
                alt = parts[3]
                return chr_, pos, ref, alt
        except:
            pass
        return 'NA', 0, 'NA', 'NA'
    
    # 解析变异信息
    variant_info = gene_data['variant_id'].apply(parse_variant_id)
    gene_data['CHR'] = [x[0] for x in variant_info]
    gene_data['BP'] = [x[1] for x in variant_info]
    gene_data['REF'] = [x[2] for x in variant_info]
    gene_data['ALT'] = [x[3] for x in variant_info]
    
    # 格式化输出
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
        'TISSUE': tissue_type,
        'TSS_DISTANCE': gene_data['start_distance']
    })
    
    return mr_format

# ================================================================================
# 为每个基因创建暴露文件
# ================================================================================

print("\n" + "="*70)
print("格式化 eQTL 数据为 MR 输入")
print("="*70)

# 获取所有可用基因
if len(brain_eqtl) > 0:
    brain_genes = set(brain_eqtl['phenotype_id'].unique())
    print(f"\n脑皮层 eQTL 数据包含：{len(brain_genes):,} 个基因")
else:
    brain_genes = set()

if len(blood_eqtl) > 0:
    blood_genes = set(blood_eqtl['phenotype_id'].unique())
    print(f"全血 eQTL 数据包含：{len(blood_genes):,} 个基因")
else:
    blood_genes = set()

# 如果有基因列表，使用基因列表；否则使用所有 eQTL 基因
if len(gene_list) > 0:
    target_genes = gene_list
    print(f"\n目标基因：{len(target_genes)} 个（来自基因列表）")
else:
    # 使用所有在任一组织中表达的基因
    target_genes = list(brain_genes | blood_genes)
    print(f"\n目标基因：{len(target_genes)} 个（来自 eQTL 数据）")

# 统计
stats = {
    'brain_only': 0,
    'blood_only': 0,
    'both': 0,
    'neither': 0
}

# 处理每个基因
print(f"\n开始处理 {len(target_genes)} 个基因...")
processed_count = 0

for i, gene in enumerate(target_genes, 1):
    output_file = os.path.join(OUTPUT_DIR, f"{gene}_exposure.csv")
    
    # 检查基因是否在数据中
    in_brain = gene in brain_genes
    in_blood = gene in blood_genes
    
    all_data = []
    
    # 添加脑组织数据
    if in_brain:
        brain_data = format_eqtl_for_mr(brain_eqtl, gene, "Brain_Cortex")
        if len(brain_data) > 0:
            all_data.append(brain_data)
    
    # 添加全血数据
    if in_blood:
        blood_data = format_eqtl_for_mr(blood_eqtl, gene, "Whole_Blood")
        if len(blood_data) > 0:
            all_data.append(blood_data)
    
    # 保存
    if len(all_data) > 0:
        # 合并所有组织的数据
        combined = pd.concat(all_data, ignore_index=True)
        combined.to_csv(output_file, index=False)
        processed_count += 1
        
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
    if i % 50 == 0:
        print(f"  处理进度：{i}/{len(target_genes)} (已处理：{processed_count})")

# ================================================================================
# 输出统计信息
# ================================================================================

print("\n" + "="*70)
print("处理完成")
print("="*70)

print(f"\n总基因数：{len(target_genes)}")
print(f"成功处理：{processed_count} 个基因")
print(f"\n组织分布:")
print(f"  - 双组织都有：{stats['both']} 个基因")
print(f"  - 仅脑组织：  {stats['brain_only']} 个基因")
print(f"  - 仅全血：    {stats['blood_only']} 个基因")
print(f"  - 都无数据：  {stats['neither']} 个基因")
print(f"\n输出目录：{OUTPUT_DIR}")

# 创建统计文件
stats_file = os.path.join(OUTPUT_DIR, "processing_stats.txt")
with open(stats_file, 'w', encoding='utf-8') as f:
    f.write(f"GTEx v11 eQTL 数据处理统计\n")
    f.write(f"="*50 + "\n\n")
    f.write(f"总基因数：{len(target_genes)}\n")
    f.write(f"成功处理：{processed_count}\n\n")
    f.write(f"组织分布:\n")
    f.write(f"  - 双组织都有：{stats['both']}\n")
    f.write(f"  - 仅脑组织：  {stats['brain_only']}\n")
    f.write(f"  - 仅全血：    {stats['blood_only']}\n")
    f.write(f"  - 都无数据：  {stats['neither']}\n\n")
    f.write(f"输出目录：{OUTPUT_DIR}\n")

print(f"\n✓ 已保存统计信息：{stats_file}")

# 创建基因列表文件
gene_list_output = os.path.join(OUTPUT_DIR, "gene_list_processed.txt")
with open(gene_list_output, 'w', encoding='utf-8') as f:
    for gene in target_genes:
        f.write(f"{gene}\n")

print(f"✓ 已保存基因列表：{gene_list_output}")

print("\n" + "="*70)
print("下一步")
print("="*70)
print("""
1. 检查输出文件
   目录：C:\\Users\\Jy-Mentor-7\\Desktop\\你信我 这不是重蹈覆辙\\exposure_gtex_brain

2. 运行双源 MR 分析
   Rscript run_dual_source_mr.R

3. 比较单源 vs 双源结果
   - 查看新增的显著基因
   - 分析组织特异性效应
""")

print("\n" + "="*70)
