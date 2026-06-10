#!/usr/bin/env python3
# ================================================================================
# 双源 eQTL 匹配和整合
# 将脑组织和全血的 eQTL 数据进行匹配整合
# 参考：GTEx v11 | eQTL Catalogue | TwoSampleMR
# ================================================================================

import os
import sys
import csv
import time
from collections import defaultdict

print("="*70)
print("双源 eQTL 匹配和整合")
print("参考：GTEx v11 | eQTL Catalogue | TwoSampleMR")
print("="*70)

# 配置
BRAIN_CSV = r"C:\Users\Jy-Mentor-7\Desktop\生物信息学\ETQL\brain_eqtl.csv"
BLOOD_CSV = r"C:\Users\Jy-Mentor-7\Desktop\生物信息学\ETQL\blood_eqtl.csv"
OUTPUT_DIR = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\exposure_matched"
MATCHED_DIR = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\exposure_matched\matched_data"

os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(MATCHED_DIR, exist_ok=True)
print(f"\n输出目录：{OUTPUT_DIR}\n")

start_time = time.time()

# ================================================================================
# 1. 加载基因列表
# ================================================================================
print("步骤 1: 加载基因列表")
print("-"*70)

BRAIN_EGENES = r"C:\Users\Jy-Mentor-7\Desktop\生物信息学\ETQL\Brain_Cortex.v11.eGenes.txt"
BLOOD_EGENES = r"C:\Users\Jy-Mentor-7\Desktop\生物信息学\ETQL\Whole_Blood.v11.eGenes.txt"

def load_genes(filename):
    genes = set()
    with open(filename, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            genes.add(row['gene_id'])
    return genes

print("  加载脑皮层 eGenes...")
brain_genes = load_genes(BRAIN_EGENES)
print(f"  ✓ 脑皮层：{len(brain_genes):,} 基因")

print("  加载全血 eGenes...")
blood_genes = load_genes(BLOOD_EGENES)
print(f"  ✓ 全血：{len(blood_genes):,} 基因")

# 双源匹配
common_genes = brain_genes & blood_genes
brain_only = brain_genes - blood_genes
blood_only = blood_genes - brain_genes

print(f"\n双源匹配结果:")
print(f"  - 双组织都有：{len(common_genes):,} 基因")
print(f"  - 仅脑组织：{len(brain_only):,} 基因")
print(f"  - 仅全血：{len(blood_only):,} 基因")

target_genes = brain_genes | blood_genes
print(f"\n目标总数：{len(target_genes):,}\n")

# ================================================================================
# 2. 处理 eQTL 数据
# ================================================================================
print("步骤 2: 处理 eQTL 数据")
print("-"*70)

def process_and_match_eqtl(brain_file, blood_file, output_dir, matched_dir):
    """处理并匹配双源 eQTL 数据"""
    
    # 存储每个基因的数据
    gene_eqtl = defaultdict(lambda: {'brain': [], 'blood': []})
    stats = {'brain_processed': 0, 'blood_processed': 0, 'saved': 0}
    
    # 处理脑组织数据
    print("  处理脑皮层数据...")
    with open(brain_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            stats['brain_processed'] += 1
            gene_id = row['phenotype_id'].split('_')[0]
            
            if gene_id not in target_genes:
                continue
            
            # 解析 variant_id
            variant_parts = row['variant_id'].split('_')
            if len(variant_parts) >= 4:
                chr_ = variant_parts[0].replace('chr', '')
                bp = variant_parts[1]
                ref = variant_parts[2]
                alt = variant_parts[3]
            else:
                continue
            
            gene_eqtl[gene_id]['brain'].append({
                'SNP': row['variant_id'],
                'CHR': chr_,
                'BP': bp,
                'EFFECT_ALLELE': alt,
                'OTHER_ALLELE': ref,
                'BETA': row['slope'],
                'SE': row['slope_se'],
                'PVAL': float(row['pval_nominal']),
                'EAF': row['af'],
                'TISSUE': 'Brain_Cortex'
            })
            
            if stats['brain_processed'] % 500000 == 0:
                print(f"    已处理 {stats['brain_processed']:,} 行")
    
    print(f"  ✓ 脑皮层：{stats['brain_processed']:,} 行\n")
    
    # 处理全血数据
    print("  处理全血数据...")
    with open(blood_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            stats['blood_processed'] += 1
            gene_id = row['phenotype_id'].split('_')[0]
            
            if gene_id not in target_genes:
                continue
            
            # 解析 variant_id
            variant_parts = row['variant_id'].split('_')
            if len(variant_parts) >= 4:
                chr_ = variant_parts[0].replace('chr', '')
                bp = variant_parts[1]
                ref = variant_parts[2]
                alt = variant_parts[3]
            else:
                continue
            
            gene_eqtl[gene_id]['blood'].append({
                'SNP': row['variant_id'],
                'CHR': chr_,
                'BP': bp,
                'EFFECT_ALLELE': alt,
                'OTHER_ALLELE': ref,
                'BETA': row['slope'],
                'SE': row['slope_se'],
                'PVAL': float(row['pval_nominal']),
                'EAF': row['af'],
                'TISSUE': 'Whole_Blood'
            })
            
            if stats['blood_processed'] % 500000 == 0:
                print(f"    已处理 {stats['blood_processed']:,} 行")
    
    print(f"  ✓ 全血：{stats['blood_processed']:,} 行\n")
    
    # 保存匹配后的数据
    print("  保存匹配后的数据...")
    
    # 1. 双组织都有的基因 - 合并数据
    print(f"    处理双组织都有的 {len(common_genes):,} 个基因...")
    for gene_name in common_genes:
        brain_data = gene_eqtl[gene_name]['brain']
        blood_data = gene_eqtl[gene_name]['blood']
        
        if len(brain_data) > 0 or len(blood_data) > 0:
            # 合并两个组织的数据
            combined = brain_data + blood_data
            
            # 按 P 值排序，选择 Top 50
            combined.sort(key=lambda x: x['PVAL'])
            top_50 = combined[:50]
            
            # 保存
            output_file = os.path.join(matched_dir, f"{gene_name}_exposure.csv")
            with open(output_file, 'w', newline='', encoding='utf-8') as f:
                fieldnames = ['SNP', 'CHR', 'BP', 'EFFECT_ALLELE', 'OTHER_ALLELE', 'BETA', 'SE', 'PVAL', 'EAF', 'GENE', 'TISSUE']
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                for row in top_50:
                    row['GENE'] = gene_name
                    writer.writerow(row)
            
            stats['saved'] += 1
    
    # 2. 仅脑组织的基因
    print(f"    处理仅脑组织的 {len(brain_only):,} 个基因...")
    for gene_name in brain_only:
        brain_data = gene_eqtl[gene_name]['brain']
        
        if len(brain_data) > 0:
            brain_data.sort(key=lambda x: x['PVAL'])
            top_50 = brain_data[:50]
            
            output_file = os.path.join(matched_dir, f"{gene_name}_exposure.csv")
            with open(output_file, 'w', newline='', encoding='utf-8') as f:
                fieldnames = ['SNP', 'CHR', 'BP', 'EFFECT_ALLELE', 'OTHER_ALLELE', 'BETA', 'SE', 'PVAL', 'EAF', 'GENE', 'TISSUE']
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                for row in top_50:
                    row['GENE'] = gene_name
                    writer.writerow(row)
            
            stats['saved'] += 1
    
    # 3. 仅全血的基因
    print(f"    处理仅全血的 {len(blood_only):,} 个基因...")
    for gene_name in blood_only:
        blood_data = gene_eqtl[gene_name]['blood']
        
        if len(blood_data) > 0:
            blood_data.sort(key=lambda x: x['PVAL'])
            top_50 = blood_data[:50]
            
            output_file = os.path.join(matched_dir, f"{gene_name}_exposure.csv")
            with open(output_file, 'w', newline='', encoding='utf-8') as f:
                fieldnames = ['SNP', 'CHR', 'BP', 'EFFECT_ALLELE', 'OTHER_ALLELE', 'BETA', 'SE', 'PVAL', 'EAF', 'GENE', 'TISSUE']
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                for row in top_50:
                    row['GENE'] = gene_name
                    writer.writerow(row)
            
            stats['saved'] += 1
    
    return stats

# 运行匹配和整合
stats = process_and_match_eqtl(BRAIN_CSV, BLOOD_CSV, OUTPUT_DIR, MATCHED_DIR)

# ================================================================================
# 3. 统计结果
# ================================================================================
elapsed_time = time.time() - start_time

print("\n" + "="*70)
print("双源匹配完成")
print("="*70)

output_files = [f for f in os.listdir(MATCHED_DIR) if f.endswith('_exposure.csv')]

print(f"\n总基因数：{len(output_files):,}")
print(f"处理时间：{elapsed_time:.2f} 秒 ({elapsed_time/60:.2f} 分钟)")
print(f"处理速度：{(stats['brain_processed'] + stats['blood_processed']) / elapsed_time / 1000:.1f}K 行/秒\n")

print("双源匹配统计:")
print(f"  - 双组织都有：{len(common_genes):,} 基因 → 合并数据")
print(f"  - 仅脑组织：{len(brain_only):,} 基因 → 使用脑数据")
print(f"  - 仅全血：{len(blood_only):,} 基因 → 使用血数据\n")

# 保存统计
with open(os.path.join(OUTPUT_DIR, "matching_stats.txt"), 'w', encoding='utf-8') as f:
    f.write(f"双源 eQTL 匹配统计\n")
    f.write(f"参考：GTEx v11 | eQTL Catalogue\n\n")
    f.write(f"双源匹配结果:\n")
    f.write(f"  - 双组织都有：{len(common_genes):,} 基因\n")
    f.write(f"  - 仅脑组织：{len(brain_only):,} 基因\n")
    f.write(f"  - 仅全血：{len(blood_only):,} 基因\n\n")
    f.write(f"总基因数：{len(output_files):,}\n")
    f.write(f"处理时间：{elapsed_time:.2f} 秒\n\n")
    f.write(f"处理统计:\n")
    f.write(f"  - 脑组织处理：{stats['brain_processed']:,} 行\n")
    f.write(f"  - 全血处理：{stats['blood_processed']:,} 行\n")
    f.write(f"  - 保存基因：{stats['saved']:,} 个\n")

print(f"✓ 统计已保存\n")

# 样本预览
if len(output_files) > 0:
    # 检查双组织样本
    sample_file = os.path.join(MATCHED_DIR, output_files[0])
    with open(sample_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    
    # 统计组织分布
    tissues = defaultdict(int)
    for row in rows:
        tissues[row['TISSUE']] += 1
    
    print(f"样本预览 ({output_files[0]}):")
    print(f"  维度：{len(rows):,} × {len(rows[0])}")
    print(f"  组织分布:")
    for tissue, count in tissues.items():
        print(f"    - {tissue}: {count:,} SNP")
    print()

print("="*70)
print("完成！下一步：Rscript run_mr_analysis_matched.R")
print("="*70)
