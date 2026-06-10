# -*- coding: utf-8 -*-
"""
GEO数据预处理 - 解析GSE16561和GSE37587并合并
使用pandas直接读取series matrix文件
"""
import os
import sys
import gzip
import pandas as pd
import numpy as np
from scipy import stats
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import RESULTS_DIR

# 数据路径 - 使用相对路径或环境变量，提高可移植性
DATA_DIR = os.environ.get('CUPROPTOSIS_DATA_DIR', 
                          str(Path.home() / 'Downloads'))

GSE16561_PATH = os.path.join(DATA_DIR, "GSE16561_series_matrix.txt.gz")
GSE37587_PATH = os.path.join(DATA_DIR, "GSE37587_series_matrix.txt.gz")
GPL6883_PATH = os.path.join(DATA_DIR, "GPL6883-11606.txt")
OUTPUT_DIR = os.path.join(RESULTS_DIR, "stage1_rma_degs")


def parse_geo_metadata(filepath):
    """解析GEO series matrix文件的元数据"""
    metadata = {}
    with gzip.open(filepath, 'rt', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line.startswith('!'):
                parts = line.split('\t')
                key = parts[0].lstrip('!')
                values = [v.strip('"') for v in parts[1:]]
                metadata[key] = values
            elif line.startswith('ID_REF') or line.startswith('"ID_REF"'):
                break
    return metadata


def read_geo_expression(filepath):
    """使用pandas直接读取GEO表达矩阵"""
    # 找到数据起始行
    skip_rows = 0
    with gzip.open(filepath, 'rt', encoding='utf-8') as f:
        for i, line in enumerate(f):
            if line.strip().startswith('ID_REF') or line.strip().startswith('"ID_REF"'):
                skip_rows = i
                break

    # 使用pandas读取
    df = pd.read_csv(filepath, compression='gzip', sep='\t',
                     skiprows=skip_rows, header=0, index_col=0)

    # 去除索引和列名的引号
    df.index = df.index.str.strip('"')
    df.columns = df.columns.str.strip('"')

    # 将null字符串转为NaN
    df = df.replace('null', np.nan)
    df = df.astype(float)

    return df


def load_gpl_annotation(gpl_path):
    """加载GPL平台注释文件"""
    print(f"加载GPL注释: {gpl_path}")

    with open(gpl_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    # 找到表头行
    header_idx = 0
    for i, line in enumerate(lines):
        if line.startswith('ID\t'):
            header_idx = i
            break

    header = lines[header_idx].strip().split('\t')
    print(f"  表头列: {header[:15]}")

    # 找到Symbol列索引
    symbol_idx = None
    for i, h in enumerate(header):
        if h == 'Symbol':
            symbol_idx = i
            break

    if symbol_idx is None:
        raise ValueError("未找到Symbol列")

    # 解析注释
    annot_data = []
    for line in lines[header_idx+1:]:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        parts = line.split('\t')
        if len(parts) > symbol_idx:
            probe_id = parts[0].strip('"')
            gene_symbol = parts[symbol_idx].strip('"')
            if gene_symbol and gene_symbol != '' and gene_symbol != '---':
                annot_data.append({'ProbeID': probe_id, 'GeneSymbol': gene_symbol})

    annot = pd.DataFrame(annot_data)
    print(f"  GPL注释: {len(annot)} 探针有GeneSymbol")
    return annot


def process_gse16561():
    """处理GSE16561: 39例卒中 vs 24例对照"""
    print("=" * 60)
    print("处理 GSE16561...")

    # 读取表达矩阵
    expr = read_geo_expression(GSE16561_PATH)

    # 读取元数据
    meta = parse_geo_metadata(GSE16561_PATH)
    sample_geo = meta.get('Sample_geo_accession', [])
    sample_titles = meta.get('Sample_title', [])
    sample_desc = meta.get('Sample_description', [])

    print(f"  样本数: {len(sample_geo)}")
    print(f"  表达矩阵: {expr.shape}")

    # 创建样本注释
    groups = []
    for desc in sample_desc:
        if 'Stroke' in str(desc):
            groups.append('Stroke')
        elif 'Control' in str(desc):
            groups.append('Control')
        else:
            groups.append('Unknown')

    sample_annot = pd.DataFrame({
        'SampleID': sample_geo,
        'Title': sample_titles,
        'Group': groups,
        'Dataset': 'GSE16561'
    })

    print(f"  分组: {sample_annot['Group'].value_counts().to_dict()}")

    return expr, sample_annot


def process_gse37587():
    """处理GSE37587: 34例卒中患者，Baseline vs Follow-Up"""
    print("=" * 60)
    print("处理 GSE37587...")

    # 读取表达矩阵
    expr = read_geo_expression(GSE37587_PATH)

    # 读取元数据
    meta = parse_geo_metadata(GSE37587_PATH)
    sample_geo = meta.get('Sample_geo_accession', [])
    sample_titles = meta.get('Sample_title', [])
    time_chars = meta.get('Sample_characteristics_ch1', [])

    # 提取时间点
    timepoints = []
    for tc in time_chars:
        if 'Baseline' in tc:
            timepoints.append('Baseline')
        elif 'Follow-Up' in tc:
            timepoints.append('FollowUp')
        else:
            timepoints.append('Unknown')

    print(f"  样本数: {len(sample_geo)}")
    print(f"  时间点: {pd.Series(timepoints).value_counts().to_dict()}")
    print(f"  表达矩阵: {expr.shape}")

    sample_annot = pd.DataFrame({
        'SampleID': sample_geo,
        'Title': sample_titles,
        'Group': timepoints,
        'Dataset': 'GSE37587'
    })

    return expr, sample_annot


def merge_and_normalize(expr1, annot1, expr2, annot2, gpl_annot):
    """合并数据集，映射到GeneSymbol，标准化"""
    print("=" * 60)
    print("合并数据集并映射到GeneSymbol...")

    # 取共同探针
    common_probes = expr1.index.intersection(expr2.index)
    print(f"  共同探针: {len(common_probes)}")

    expr1_common = expr1.loc[common_probes]
    expr2_common = expr2.loc[common_probes]

    # 合并
    expr_merged = pd.concat([expr1_common, expr2_common], axis=1)
    annot_merged = pd.concat([annot1, annot2], axis=0, ignore_index=True)

    # 处理null值（用样本均值填充）
    null_count = expr_merged.isna().sum().sum()
    if null_count > 0:
        print(f"  填充 {null_count} 个null值")
        expr_merged = expr_merged.apply(lambda x: x.fillna(x.mean()), axis=0)

    print(f"  合并后表达矩阵: {expr_merged.shape}")

    # 映射到GeneSymbol
    probe_to_gene = dict(zip(gpl_annot['ProbeID'], gpl_annot['GeneSymbol']))
    expr_merged['GeneSymbol'] = expr_merged.index.map(probe_to_gene)
    expr_mapped = expr_merged.dropna(subset=['GeneSymbol'])
    expr_mapped = expr_mapped[expr_mapped['GeneSymbol'] != '']
    expr_gene = expr_mapped.groupby('GeneSymbol').mean()

    print(f"  映射到GeneSymbol: {expr_gene.shape}")

    # 标准化 (z-score per gene)
    expr_norm = expr_gene.apply(lambda x: (x - x.mean()) / x.std() if x.std() > 0 else x, axis=1)

    return expr_norm, annot_merged


def run_ttest(expr_matrix, sample_annot):
    """运行t-test差异表达分析 with BH-FDR校正"""
    print("=" * 60)
    print("运行差异表达分析...")

    from statsmodels.stats.multitest import multipletests

    stroke_samples = sample_annot[sample_annot['Group'] == 'Stroke']['SampleID'].values
    control_samples = sample_annot[sample_annot['Group'] == 'Control']['SampleID'].values

    # 确保样本在表达矩阵中
    stroke_samples = [s for s in stroke_samples if s in expr_matrix.columns]
    control_samples = [s for s in control_samples if s in expr_matrix.columns]

    print(f"  Stroke样本: {len(stroke_samples)}")
    print(f"  Control样本: {len(control_samples)}")

    results = []
    for gene in expr_matrix.index:
        stroke_expr = expr_matrix.loc[gene, stroke_samples].values
        control_expr = expr_matrix.loc[gene, control_samples].values

        # 移除NaN
        stroke_expr = stroke_expr[~np.isnan(stroke_expr)]
        control_expr = control_expr[~np.isnan(control_expr)]

        if len(stroke_expr) < 2 or len(control_expr) < 2:
            continue

        t_stat, p_val = stats.ttest_ind(stroke_expr, control_expr)
        logfc = np.mean(stroke_expr) - np.mean(control_expr)

        results.append({
            'GeneSymbol': gene,
            'logFC': logfc,
            'AveExpr': np.mean(expr_matrix.loc[gene].values),
            't': t_stat,
            'P.Value': p_val,
            'adj.P.Val': np.nan,
            'B': 0
        })

    results_df = pd.DataFrame(results)

    # BH-FDR校正
    p_values = results_df['P.Value'].values
    p_values = np.nan_to_num(p_values, nan=1.0)
    _, fdr, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')
    results_df['adj.P.Val'] = fdr

    results_df = results_df.sort_values('P.Value')

    print(f"  DEGs: {len(results_df)}")
    sig = results_df[results_df['adj.P.Val'] < 0.05]
    print(f"  显著DEGs (FDR<0.05): {len(sig)}")
    if len(sig) > 0:
        print(f"  Top DEG: {sig.iloc[0]['GeneSymbol']} (logFC={sig.iloc[0]['logFC']:.3f}, FDR={sig.iloc[0]['adj.P.Val']:.2e})")

    return results_df


def main():
    """主函数"""
    print("=" * 60)
    print("GEO数据预处理")
    print("=" * 60)

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # 1. 加载GPL注释
    gpl_annot = load_gpl_annotation(GPL6883_PATH)

    # 2. 处理GSE16561
    expr1, annot1 = process_gse16561()

    # 3. 处理GSE37587
    expr2, annot2 = process_gse37587()

    # 4. 合并、映射、标准化
    expr_norm, annot_merged = merge_and_normalize(expr1, annot1, expr2, annot2, gpl_annot)

    # 5. 保存表达矩阵
    expr_file = os.path.join(OUTPUT_DIR, "expr_matrix.csv")
    expr_norm.to_csv(expr_file)
    print(f"  表达矩阵已保存: {expr_file}")

    # 6. 保存样本注释
    annot_file = os.path.join(OUTPUT_DIR, "sample_annotations.csv")
    annot_merged.to_csv(annot_file, index=False)
    print(f"  样本注释已保存: {annot_file}")

    # 7. 差异表达分析（仅GSE16561有对照组）
    degs = run_ttest(expr_norm, annot1)

    degs_file = os.path.join(OUTPUT_DIR, "limma_degs.csv")
    degs.to_csv(degs_file, index=False)
    print(f"  DEGs已保存: {degs_file}")

    # 8. 保存sample_groups.csv（兼容主控脚本）
    groups_df = annot_merged[['SampleID', 'Group']].rename(columns={'SampleID': 'Sample', 'Group': 'Group'})
    groups_file = os.path.join(OUTPUT_DIR, "sample_groups.csv")
    groups_df.to_csv(groups_file, index=False)
    print(f"  样本分组已保存: {groups_file}")

    print("=" * 60)
    print("数据预处理完成!")


if __name__ == "__main__":
    main()
