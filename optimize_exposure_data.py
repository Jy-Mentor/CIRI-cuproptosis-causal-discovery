#!/usr/bin/env python3
"""
优化版暴露数据预处理脚本
目标：提高数据暴露率从 82% 到 95%+

优化策略:
1. 使用更宽松的 P 值阈值 (1e-5 代替 5e-8)
2. 优化 LD clump 参数 (r2 < 0.01 代替 0.001)
3. 分级阈值策略
4. 自动检测并使用最佳数据源
"""

import pandas as pd
import numpy as np
import os
import sys
from pathlib import Path
from datetime import datetime

# ==================== 配置参数 ====================

# 输入文件路径（按优先级排序）
INPUT_FILES = [
    r"D:\EQTL\clump\eQTLgen_allgene_p_1e-05_kb_10000_r2_0.01.xlsx",  # 推荐：P<1e-5
    r"D:\EQTL\clump\eQTLgen_allgene_p_5e-06_kb_10000_r2_0.01.xlsx",  # 备选：P<5e-6
    r"D:\EQTL\clump\eQTLgen_allgene_p_5e-08_kb_10000_r2_0.001.xlsx", # 严格：P<5e-8
]

# 输出目录
OUTPUT_DIR = r"D:\下载\MR_batch_results\exposure_optimized"
GENE_LIST_FILE = r"D:\下载\MR_batch_results\gene_list_optimized.txt"

# 候选基因列表（130 个）
CANDIDATE_GENES = [
    "ACAD11", "PTGR1", "ACADVL", "CPT2", "HSD17B4", "ACADM", "PDHX", "FABP4", "HIBADH", "CPT1A",
    "DLAT", "PPARG", "PDHB", "ALDH9A1", "SREBF1", "EPHX1", "ACTA2", "TIMP1", "TGFB1", "CCL2",
    "CTSS", "FNTA", "PTPRC", "EGFR", "ITGA1", "COL1A1", "ADRB1", "HTR2B", "HTR2C", "AKT1",
    "AIF1", "HMOX1", "NFKB1", "CTSC", "IL10RA", "CTSD", "C3", "GFAP", "STAT3", "ICAM1",
    "CCR5", "CTSB", "CASK", "TSPO", "IL6", "TNF", "CCND1", "NR3C1", "IGFBP2", "IRF1",
    "CDK4", "PRKCQ", "LEF1", "CTSK", "RHOC", "HPGDS", "FABP5", "PA2G4", "PTPN2", "SPHK1",
    "ZHX2", "IMPDH2", "HSPA5", "F3", "STAT5A", "CTSL", "ATP7A", "GPX4", "CASP8", "XRCC6",
    "STAT1", "ZEB1", "MTOR", "PTGS1", "PARP1", "MAPKAPK2", "PLA2G4A", "CXCR3", "MAOB", "XDH",
    "NFE2L2", "RELA", "PTPN6", "CHFR", "TOP2A", "LYN", "IKBKB", "HIF1A", "B2M", "CCNA2",
    "CAT", "CNDP2", "ATOX1", "ATP7B", "CP", "SLC31A1", "MKNK2", "MB", "SEC13", "OAZ1",
    "JAK1", "BRD3", "CITED2", "NR1H3", "CTSF", "MAN2B1", "CUL4B", "RBM39", "POLR2D", "DDC",
    "GCH1", "LIAS", "LIPT1", "PABPC1", "PTPRF", "S100A6", "PTPRJ", "PDCD6IP", "FABP2", "FDX1",
    "HSD17B10", "TBXAS1", "FLT4", "SAT1", "HBS1L", "SAT2", "KCNA5", "SCN9A", "NMT1", "NUDCD2",
    "PARP12", "TDP1", "PCTP", "STARD13", "PDCD6"
]

# 列名映射
COLUMN_MAPPING = {
    'SNP': ['SNP', 'SNPID', 'RSID', 'RSNUMBER'],
    'CHR': ['CHR', 'CHR.EXPOSURE', 'CHROMOSOME'],
    'BP': ['BP', 'POS.EXPOSURE', 'POSITION', 'BPOSITION'],
    'EFFECT_ALLELE': ['EFFECT_ALLELE', 'EA', 'A1', 'ALLELE1', 'EFFECT_ALLELE.EXPOSURE'],
    'OTHER_ALLELE': ['OTHER_ALLELE', 'OA', 'A2', 'ALLELE0', 'OTHER_ALLELE.EXPOSURE'],
    'BETA': ['BETA', 'BETA.EXPOSURE', 'B', 'Z'],
    'SE': ['SE', 'SE.EXPOSURE'],
    'PVAL': ['PVAL', 'PVAL.EXPOSURE', 'P', 'PVALUE', 'P-VALUE'],
    'EAF': ['EAF', 'EAF.EXPOSURE', 'FRQ', 'FREQ', 'MAF'],
    'GENE': ['GENE', 'GENENAME', 'SYMBOL']
}

# ==================== 工具函数 ====================

def find_best_input_file():
    """
    寻找最佳可用输入文件
    优先级：P<1e-5 > P<5e-6 > P<5e-8
    """
    print("=" * 60)
    print("步骤 1: 寻找最佳暴露数据源")
    print("=" * 60)
    
    for file_path in INPUT_FILES:
        if os.path.exists(file_path):
            file_size = os.path.getsize(file_path) / (1024 * 1024)  # MB
            print(f"✓ 找到文件：{file_path}")
            print(f"  文件大小：{file_size:.2f} MB")
            return file_path
    
    print("✗ 错误：未找到任何暴露数据文件")
    print(f" searched paths: {INPUT_FILES}")
    sys.exit(1)

def standardize_columns(df):
    """
    标准化列名
    """
    # 转换为大写
    df.columns = df.columns.str.upper()
    
    # 应用映射
    rename_dict = {}
    for target_col, source_cols in COLUMN_MAPPING.items():
        for source_col in source_cols:
            if source_col in df.columns:
                rename_dict[source_col] = target_col
                break
    
    if rename_dict:
        df = df.rename(columns=rename_dict)
        print(f"  列名转换：{len(rename_dict)} 列")
    
    return df

def check_data_quality(df):
    """
    检查数据质量
    """
    print("\n数据质量检查:")
    
    # 检查必需列
    required_cols = ['SNP', 'CHR', 'BP', 'EFFECT_ALLELE', 'OTHER_ALLELE', 'BETA', 'SE', 'PVAL', 'EAF', 'GENE']
    missing = [col for col in required_cols if col not in df.columns]
    
    if missing:
        print(f"  ✗ 缺少列：{missing}")
        return False
    
    print(f"  ✓ 所有必需列存在")
    
    # 检查缺失值
    for col in ['BETA', 'SE', 'PVAL', 'EAF']:
        n_na = df[col].isna().sum()
        if n_na > 0:
            print(f"  ⚠ {col} 列有 {n_na} 个缺失值")
    
    # 检查 P 值分布
    pval_min = df['PVAL'].min()
    pval_max = df['PVAL'].max()
    print(f"  ✓ P 值范围：{pval_min:.2e} - {pval_max:.2e}")
    
    # 检查等位基因频率
    eaf_min = df['EAF'].min()
    eaf_max = df['EAF'].max()
    print(f"  ✓ EAF 范围：{eaf_min:.3f} - {eaf_max:.3f}")
    
    return True

def process_with_optimized_parameters():
    """
    使用优化参数处理暴露数据
    """
    print("=" * 60)
    print("步骤 2: 使用优化参数处理暴露数据")
    print("=" * 60)
    
    # 1. 找到最佳输入文件
    input_file = find_best_input_file()
    
    # 2. 读取数据
    print(f"\n读取文件：{input_file}")
    try:
        df = pd.read_excel(input_file)
        print(f"  读取完成：{len(df):,} 行，{len(df.columns)} 列")
    except Exception as e:
        print(f"  ✗ 读取失败：{e}")
        sys.exit(1)
    
    # 3. 标准化列名
    df = standardize_columns(df)
    
    # 4. 数据质量检查
    if not check_data_quality(df):
        print("  ✗ 数据质量检查失败")
        sys.exit(1)
    
    # 5. 过滤候选基因
    print(f"\n过滤候选基因...")
    original_genes = df['GENE'].nunique()
    df = df[df['GENE'].isin(CANDIDATE_GENES)]
    filtered_genes = df['GENE'].nunique()
    print(f"  原始基因数：{original_genes}")
    print(f"  过滤后基因数：{filtered_genes}")
    print(f"  过滤后 SNP 数：{len(df):,}")
    
    # 6. 按基因分组并导出
    print(f"\n按基因分组并导出...")
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    gene_stats = []
    
    for gene in CANDIDATE_GENES:
        gene_df = df[df['GENE'] == gene].copy()
        
        if len(gene_df) == 0:
            print(f"  {gene:10s}: ✗ 无数据")
            gene_stats.append(f"{gene}\t0\tNO_DATA")
            continue
        
        # 添加 samplesize 列
        gene_df['samplesize'] = 31684
        
        # 选择输出列
        output_cols = ['SNP', 'CHR', 'BP', 'EFFECT_ALLELE', 'OTHER_ALLELE', 
                      'BETA', 'SE', 'PVAL', 'EAF', 'GENE', 'samplesize']
        available_cols = [col for col in output_cols if col in gene_df.columns]
        gene_df = gene_df[available_cols]
        
        # 保存
        output_file = os.path.join(OUTPUT_DIR, f"{gene}.exposure.csv")
        gene_df.to_csv(output_file, index=False)
        
        n_snps = len(gene_df)
        status = "OK"
        print(f"  {gene:10s}: ✓ {n_snps:3d} SNPs")
        gene_stats.append(f"{gene}\t{n_snps}\t{status}")
    
    # 7. 生成基因清单
    print(f"\n生成基因清单：{GENE_LIST_FILE}")
    with open(GENE_LIST_FILE, 'w') as f:
        f.write("GENE\tN_SNPS\tSTATUS\n")
        for stat in gene_stats:
            f.write(stat + "\n")
    
    # 8. 统计摘要
    n_with_data = sum(1 for stat in gene_stats if stat.endswith("OK"))
    n_without_data = sum(1 for stat in gene_stats if stat.endswith("NO_DATA"))
    total_snps = sum(int(stat.split('\t')[1]) for stat in gene_stats if stat.endswith("OK"))
    
    print(f"\n" + "=" * 60)
    print("优化处理完成 - 摘要统计")
    print("=" * 60)
    print(f"总候选基因数：     {len(CANDIDATE_GENES)}")
    print(f"有暴露数据基因：   {n_with_data} ({n_with_data/len(CANDIDATE_GENES)*100:.1f}%)")
    print(f"无暴露数据基因：   {n_without_data} ({n_without_data/len(CANDIDATE_GENES)*100:.1f}%)")
    print(f"总 SNP 数量：        {total_snps:,}")
    print(f"平均每个基因 SNP:  {total_snps/n_with_data:.2f}")
    print(f"\n输出目录：{OUTPUT_DIR}")
    print(f"基因清单：{GENE_LIST_FILE}")
    
    # 9. 对比优化效果
    print(f"\n" + "=" * 60)
    print("优化效果对比")
    print("=" * 60)
    
    # 原始数据（P<5e-8）
    original_exposure_rate = 107 / 130 * 100  # 82.3%
    
    # 优化后（P<1e-5）
    optimized_exposure_rate = n_with_data / len(CANDIDATE_GENES) * 100
    
    improvement = optimized_exposure_rate - original_exposure_rate
    
    print(f"优化前暴露率：     {original_exposure_rate:.1f}% (107/130)")
    print(f"优化后暴露率：     {optimized_exposure_rate:.1f}% ({n_with_data}/{len(CANDIDATE_GENES)})")
    print(f"暴露率提升：       +{improvement:.1f}%")
    print(f"新增基因数：       {n_with_data - 107}")
    
    if n_with_data - 107 > 0:
        print(f"\n✓ 优化成功！暴露率从 {original_exposure_rate:.1f}% 提升至 {optimized_exposure_rate:.1f}%")
    else:
        print(f"\n⚠ 优化效果不明显，请检查数据源")
    
    return n_with_data, len(CANDIDATE_GENES)

def compare_clump_parameters():
    """
    比较不同 clump 参数的效果
    """
    print("\n" + "=" * 60)
    print("LD Clump 参数对比分析")
    print("=" * 60)
    
    # 可用的 clump 文件
    clump_files = {
        "严格 (r2<0.001)": r"D:\EQTL\clump\eQTLgen_allgene_p_5e-08_kb_10000_r2_0.001.xlsx",
        "中等 (r2<0.01)": r"D:\EQTL\clump\eQTLgen_allgene_p_1e-05_kb_10000_r2_0.01.xlsx",
    }
    
    results = []
    
    for param_name, file_path in clump_files.items():
        if os.path.exists(file_path):
            df = pd.read_excel(file_path)
            df = df[df['GENE'].isin(CANDIDATE_GENES)]
            n_genes = df['GENE'].nunique()
            n_snps = len(df)
            avg_snps = n_snps / n_genes if n_genes > 0 else 0
            
            results.append({
                '参数': param_name,
                '基因数': n_genes,
                'SNP 总数': n_snps,
                '平均 SNP/基因': f"{avg_snps:.2f}"
            })
            
            print(f"\n{param_name}:")
            print(f"  基因数：{n_genes}")
            print(f"  SNP 总数：{n_snps:,}")
            print(f"  平均 SNP/基因：{avg_snps:.2f}")
    
    # 创建对比 DataFrame
    if results:
        对比_df = pd.DataFrame(results)
        print(f"\n参数对比表:")
        print(对比_df.to_string(index=False))
    
    return results

# ==================== 主程序 ====================

if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("暴露数据优化处理脚本")
    print("版本：2.0 (优化版)")
    print(f"时间：{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)
    
    try:
        # 1. 使用优化参数处理
        n_genes, total_genes = process_with_optimized_parameters()
        
        # 2. 比较 clump 参数
        compare_clump_parameters()
        
        # 3. 下一步建议
        print("\n" + "=" * 60)
        print("下一步操作建议")
        print("=" * 60)
        print("1. 检查优化后的暴露数据:")
        print(f"   目录：{OUTPUT_DIR}")
        print(f"   基因清单：{GENE_LIST_FILE}")
        print("\n2. 运行批量 MR 分析:")
        print(f'   Rscript mr_analysis_batch.R "{OUTPUT_DIR}" ./outcome "{GENE_LIST_FILE}"')
        print("\n3. 比较优化前后的结果:")
        print("   - 暴露率对比")
        print("   - 显著基因数量对比")
        print("   - F 统计量对比")
        
        print("\n✓ 优化处理完成！")
        
    except Exception as e:
        print(f"\n✗ 错误：{e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
