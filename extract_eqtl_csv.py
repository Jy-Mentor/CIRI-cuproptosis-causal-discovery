#!/usr/bin/env python3
# ================================================================================
# 提取双源 eQTL 数据 - CSV 版本
# 使用纯 Python CSV 处理，避免 numpy/pandas 依赖
# 参考：GTEx v11 | eQTL Catalogue | TwoSampleMR
# ================================================================================

import os
import sys
import csv
import time
from collections import defaultdict

print("="*70)
print("双源 eQTL 数据提取 - CSV 版本")
print("参考：GTEx v11 | eQTL Catalogue | TwoSampleMR")
print("="*70)

# 配置
GTEx_BRAIN_FILE = r"C:\Users\Jy-Mentor-7\Desktop\生物信息学\ETQL\Brain_Cortex.v11.eQTLs.signif_pairs.parquet"
GTEx_BLOOD_FILE = r"C:\Users\Jy-Mentor-7\Desktop\生物信息学\ETQL\Whole_Blood.v11.eQTLs.signif_pairs.parquet"
GTEx_BRAIN_EGENES = r"C:\Users\Jy-Mentor-7\Desktop\生物信息学\ETQL\Brain_Cortex.v11.eGenes.txt"
GTEx_BLOOD_EGENES = r"C:\Users\Jy-Mentor-7\Desktop\生物信息学\ETQL\Whole_Blood.v11.eGenes.txt"
OUTPUT_DIR = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\exposure_dual_source"

os.makedirs(OUTPUT_DIR, exist_ok=True)
print(f"\n输出目录：{OUTPUT_DIR}\n")

start_time = time.time()

# ================================================================================
# 1. 加载基因列表
# ================================================================================
print("步骤 1: 加载基因列表")
print("-"*70)

def load_genes(filename):
    genes = set()
    with open(filename, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            genes.add(row['gene_id'])
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
# 2. 使用 pyarrow 读取 parquet
# ================================================================================
print("步骤 2: 读取 eQTL 数据")
print("-"*70)

try:
    import pyarrow.parquet as pq
    
    print("  加载脑皮层 eQTL...")
    brain_table = pq.read_table(GTEx_BRAIN_FILE)
    print(f"  ✓ 脑皮层：{brain_table.num_rows:,} 行")
    
    print("  加载全血 eQTL...")
    blood_table = pq.read_table(GTEx_BLOOD_FILE)
    print(f"  ✓ 全血：{blood_table.num_rows:,} 行\n")
    
    # 转换为 CSV 临时文件
    import tempfile
    
    print("  转换为临时 CSV 文件...")
    brain_csv = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.csv', encoding='utf-8')
    blood_csv = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.csv', encoding='utf-8')
    
    # 写入 CSV
    brain_df = brain_table.to_pandas()
    blood_df = blood_table.to_pandas()
    
    brain_df.to_csv(brain_csv.name, index=False)
    blood_df.to_csv(blood_csv.name, index=False)
    
    brain_csv.close()
    blood_csv.close()
    
    print(f"  ✓ 临时文件：{brain_csv.name}, {blood_csv.name}\n")
    
    USE_PARQUET = True
    
except Exception as e:
    print(f"  ✗ pyarrow 失败：{e}")
    print("  请手动转换 parquet 为 CSV 格式\n")
    sys.exit(1)

# ================================================================================
# 3. 处理 CSV 数据
# ================================================================================
print("步骤 3: 处理 eQTL 数据")
print("-"*70)

def process_csv_file(csv_file, target_genes, tissue_type, output_dir):
    """处理 CSV 文件"""
    print(f"  处理 {tissue_type}...")
    
    gene_data = defaultdict(list)
    stats = {'processed': 0, 'saved': 0}
    
    with open(csv_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        
        for row in reader:
            stats['processed'] += 1
            
            # 从 phenotype_id 提取 gene_id
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
            
            # 存储数据
            gene_data[gene_id].append({
                'SNP': row['variant_id'],
                'CHR': chr_,
                'BP': bp,
                'EFFECT_ALLELE': alt,
                'OTHER_ALLELE': ref,
                'BETA': row['slope'],
                'SE': row['slope_se'],
                'PVAL': float(row['pval_nominal']),
                'EAF': row['af'],
                'TISSUE': tissue_type
            })
            
            # 进度显示
            if stats['processed'] % 500000 == 0:
                print(f"    已处理 {stats['processed']:,} 行")
    
    # 保存每个基因的文件
    for gene_name, data_list in gene_data.items():
        # 按 P 值排序，选择 Top 50
        data_list.sort(key=lambda x: x['PVAL'])
        top_50 = data_list[:50]
        
        output_file = os.path.join(output_dir, f"{gene_name}_exposure.csv")
        
        # 如果文件已存在，追加数据
        file_exists = os.path.exists(output_file)
        
        with open(output_file, 'a', newline='', encoding='utf-8') as f:
            fieldnames = ['SNP', 'CHR', 'BP', 'EFFECT_ALLELE', 'OTHER_ALLELE', 'BETA', 'SE', 'PVAL', 'EAF', 'GENE', 'TISSUE']
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            
            if not file_exists:
                writer.writeheader()
                stats['saved'] += 1
            
            for row in top_50:
                row['GENE'] = gene_name
                writer.writerow(row)
        
        # 进度显示
        if stats['saved'] % 1000 == 0:
            print(f"    已保存 {stats['saved']:,} 个基因")
    
    return stats

# 处理脑组织数据
brain_stats = process_csv_file(brain_csv.name, brain_genes, "Brain_Cortex", OUTPUT_DIR)
print(f"  ✓ 脑皮层：处理 {brain_stats['processed']:,} 行，保存 {brain_stats['saved']:,} 基因\n")

# 处理全血数据
blood_stats = process_csv_file(blood_csv.name, blood_genes, "Whole_Blood", OUTPUT_DIR)
print(f"  ✓ 全血：处理 {blood_stats['processed']:,} 行，保存 {blood_stats['saved']:,} 基因\n")

# 清理临时文件
import os
try:
    os.unlink(brain_csv.name)
    os.unlink(blood_csv.name)
except:
    pass

# ================================================================================
# 4. 统计结果
# ================================================================================
elapsed_time = time.time() - start_time

print("="*70)
print("处理完成")
print("="*70)

output_files = [f for f in os.listdir(OUTPUT_DIR) if f.endswith('_exposure.csv')]

print(f"\n总基因数：{len(output_files):,}")
print(f"处理时间：{elapsed_time:.2f} 秒 ({elapsed_time/60:.2f} 分钟)")
print(f"处理速度：{(brain_stats['processed'] + blood_stats['processed']) / elapsed_time / 1000:.1f}K 行/秒\n")

# 统计组织分布
stats = {'both': 0, 'brain_only': 0, 'blood_only': 0}

for filename in output_files:
    filepath = os.path.join(OUTPUT_DIR, filename)
    with open(filepath, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        tissues = set(row['TISSUE'] for row in reader)
    
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

# 保存统计
with open(os.path.join(OUTPUT_DIR, "stats.txt"), 'w', encoding='utf-8') as f:
    f.write(f"双源 eQTL 整合统计\n")
    f.write(f"参考：GTEx v11 | eQTL Catalogue\n\n")
    f.write(f"总基因数：{len(output_files):,}\n")
    f.write(f"处理时间：{elapsed_time:.2f} 秒\n\n")
    f.write(f"组织分布:\n")
    f.write(f"  - 双组织都有：{stats['both']:,}\n")
    f.write(f"  - 仅脑组织：  {stats['brain_only']:,}\n")
    f.write(f"  - 仅全血：    {stats['blood_only']:,}\n")

print("✓ 统计已保存\n")

# 样本预览
if len(output_files) > 0:
    with open(os.path.join(OUTPUT_DIR, output_files[0]), 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    
    print(f"样本预览 ({output_files[0]}):")
    print(f"  维度：{len(rows):,} × {len(rows[0]) if rows else 0}")
    print(f"  列：{', '.join(rows[0].keys()) if rows else 'N/A'}\n")

print("="*70)
print("完成！")
print("="*70)
