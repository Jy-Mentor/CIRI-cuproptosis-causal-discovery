#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GEO数据集差异表达分析管道
处理三个脑缺血再灌注损伤数据集：GSE61616 (7d, 大鼠芯片), GSE97537 (24H, 大鼠芯片), GSE104036 (多时序, 小鼠RNA-seq)

工作流程:
  1. 自动识别并解析Series Matrix或Count Matrix格式
  2. 芯片数据: 探针ID → Gene Symbol (GPL1355注释), 使用Welch t-test + BH校正
  3. RNA-seq数据: 使用PyDESeq2进行差异表达分析
  4. 保存各数据集DEGs、交集、并集
  5. 输出 disease_genes_for_GAT.csv

输出文件:
  - GSE61616_7d_DEGs.csv
  - GSE97537_24h_DEGs.csv
  - GSE104036_multitime_DEGs.csv
  - common_DEGs.csv
  - union_DEGs.csv
  - disease_genes_for_GAT.csv
"""

import pandas as pd
import numpy as np
import gzip
import os
import io
import re
import warnings
import logging
from scipy import stats
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings('ignore')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# ============================================================
# 路径配置
# ============================================================

BASE_DIRS = {
    'GSE61616': r'D:\反向网络药理学\L1 数据集\bulk\GSE61616（7d）',
    'GSE97537': r'D:\反向网络药理学\L1 数据集\bulk\GSE97537(24H)',
    'GSE104036': r'D:\反向网络药理学\L1 数据集\bulk\GSE104036（多时序）'
}
OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))
PLATFORM_FILE_GPL1355 = os.path.join(BASE_DIRS['GSE61616'], 'GPL1355-10794 (1).txt')

# ============================================================
# 样本分组映射（已验证自Series Matrix元数据）
# ============================================================

GSE61616_GROUPS = {
    'Sham':   ['GSM1509422', 'GSM1509423', 'GSM1509424', 'GSM1509425', 'GSM1509426'],
    'Model':  ['GSM1509427', 'GSM1509428', 'GSM1509429', 'GSM1509430', 'GSM1509431'],
}
GSE97537_GROUPS = {
    'Sham': ['GSM2571742', 'GSM2571743', 'GSM2571744', 'GSM2571745', 'GSM2571746'],
    'MCAO': ['GSM2571735', 'GSM2571736', 'GSM2571737', 'GSM2571738',
             'GSM2571739', 'GSM2571740', 'GSM2571741'],
}
GSE104036_COUNT_COLS = {
    'Sham': ['S1', 'S2', 'S3'],
    'Ipsilateral': ['I1_3hr', 'I2_3hr', 'I3_3hr', 'I1_6hr', 'I2_6hr', 'I3_6hr',
                    'I1_12hr', 'I2_12hr', 'I3_12hr', 'I1_24hr', 'I2_24hr', 'I3_24hr'],
}

DEG_CONFIG = {
    'GSE61616':  {'case': 'Model', 'control': 'Sham',  'label': 'GSE61616_7d'},
    'GSE97537':  {'case': 'MCAO',  'control': 'Sham',  'label': 'GSE97537_24h'},
    'GSE104036': {'case': 'Ipsilateral', 'control': 'Sham', 'label': 'GSE104036_multitime'},
}

# ============================================================
# 工具函数
# ============================================================

def find_file(dir_path, patterns):
    """在目录中查找匹配模式的文件"""
    for f in os.listdir(dir_path):
        for pat in patterns:
            if pat.lower() in f.lower():
                return os.path.join(dir_path, f)
    return None

def parse_series_matrix(filepath):
    """解析GEO Series Matrix文件，返回表达矩阵和列名"""
    open_func = gzip.open if filepath.endswith('.gz') else open
    with open_func(filepath, 'rt', encoding='latin-1') as f:
        content = f.read()

    meta = {}
    for line in content.splitlines():
        if line.startswith('!Sample_title'):
            parts = line.split('\t')
            meta['sample_titles'] = [p.strip('"') for p in parts[1:]]
        elif line.startswith('!Sample_geo_accession'):
            parts = line.split('\t')
            meta['sample_geo'] = [p.strip('"') for p in parts[1:]]

    data_start = content.find('!series_matrix_table_begin')
    data_end = content.find('!series_matrix_table_end')
    if data_start == -1 or data_end == -1:
        raise ValueError(f'无法在 {filepath} 中找到 series_matrix_table 标记')

    table_text = content[data_start:data_end]
    table_text = table_text.replace('!series_matrix_table_begin', '').strip()

    df = pd.read_csv(io.StringIO(table_text), sep='\t', quoting=1,
                     dtype=str, low_memory=False)
    if 'ID_REF' in df.columns:
        df = df.set_index('ID_REF')
    else:
        df = df.set_index(df.columns[0])

    df = df.apply(pd.to_numeric, errors='coerce')
    return df, meta


def parse_gpl1355_annotation(filepath):
    """解析GPL1355平台注释文件，返回探针→Gene Symbol映射"""
    mapping = {}
    with open(filepath, 'r', encoding='latin-1') as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            parts = line.strip().split('\t')
            if len(parts) < 11:
                continue
            probe_id = parts[0].strip()
            gene_symbol = parts[10].strip()  # Gene Symbol列
            if gene_symbol and gene_symbol != '---':
                mapping[probe_id] = gene_symbol.split('///')[0].strip()
    return mapping


def collapse_probes_to_genes(expr_df, probe_to_gene):
    """将探针级别的表达矩阵折叠为基因级别（取每个基因平均表达最高的探针）"""
    mapped_probes = set(probe_to_gene.keys()) & set(expr_df.index)
    expr_mapped = expr_df.loc[list(mapped_probes)]
    probe_to_gene_sub = {p: probe_to_gene[p] for p in mapped_probes}

    gene_rows = []
    for gene in set(probe_to_gene_sub.values()):
        probes = [p for p in mapped_probes if probe_to_gene_sub[p] == gene]
        if len(probes) == 1:
            gene_rows.append((gene, expr_mapped.loc[probes[0]]))
        else:
            sub = expr_mapped.loc[probes]
            mean_expr = sub.mean(axis=1)
            best_probe = mean_expr.idxmax()
            gene_rows.append((gene, expr_mapped.loc[best_probe]))

    result = pd.DataFrame([row[1] for row in gene_rows], index=[r[0] for r in gene_rows])
    return result


def deg_microarray_t_test(expr_df, case_samples, control_samples):
    """对芯片表达矩阵执行Welch t-test + Benjamini-Hochberg校正"""
    results = []
    case = expr_df[case_samples].values
    control = expr_df[control_samples].values

    for i, gene in enumerate(expr_df.index):
        c = control[i, :].astype(float)
        t = case[i, :].astype(float)
        if np.all(np.isnan(c)) or np.all(np.isnan(t)):
            continue
        c = c[~np.isnan(c)]
        t = t[~np.isnan(t)]
        if len(c) < 2 or len(t) < 2:
            continue
        log2fc = np.mean(t) - np.mean(c)
        stat, pval = stats.ttest_ind(t, c, equal_var=False)
        results.append({
            'gene_symbol': gene,
            'log2FoldChange': log2fc,
            'stat': stat,
            'pvalue': pval
        })

    res_df = pd.DataFrame(results)
    if res_df.empty:
        return res_df

    reject, padj, _, _ = multipletests(res_df['pvalue'], method='fdr_bh')
    res_df['padj'] = padj
    res_df = res_df.sort_values('pvalue')
    return res_df


def deg_rnaseq_pydeseq2(counts_df, metadata, case_label, control_label):
    """对RNA-seq计数矩阵执行PyDESeq2差异表达分析"""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    samples = metadata.index.tolist()
    common_genes = counts_df.index.tolist()

    counts_sub = counts_df.loc[common_genes, samples].astype(int)
    counts_sub = counts_sub[~(counts_sub == 0).all(axis=1)]

    dds = DeseqDataSet(
        counts=counts_sub.T,
        metadata=metadata,
        design='~condition',
    )
    dds.deseq2()

    stat_res = DeseqStats(dds, contrast=['condition', case_label, control_label])
    stat_res.summary()

    result = stat_res.results_df.copy()
    result = result.reset_index()
    first_col = result.columns[0]
    if first_col != 'gene_symbol':
        result = result.rename(columns={first_col: 'gene_symbol'})
    result = result[['gene_symbol', 'log2FoldChange', 'pvalue', 'padj']].dropna()
    return result


# ============================================================
# 数据集处理函数
# ============================================================

def process_gse61616():
    """处理GSE61616数据集：大鼠芯片，Model (7d) vs Sham"""
    logger.info('=' * 60)
    logger.info('开始处理 GSE61616 (7d 大鼠芯片，Model vs Sham)')
    dir_path = BASE_DIRS['GSE61616']

    sm_file = find_file(dir_path, ['series_matrix.txt'])
    if not sm_file:
        raise FileNotFoundError(f'未在 {dir_path} 中找到series matrix文件')
    logger.info(f'发现Series Matrix文件: {sm_file}')

    expr_df, meta = parse_series_matrix(sm_file)
    logger.info(f'表达矩阵维度: {expr_df.shape}')

    avail_samples = set(expr_df.columns)
    sham_sams = [s for s in GSE61616_GROUPS['Sham'] if s in avail_samples]
    model_sams = [s for s in GSE61616_GROUPS['Model'] if s in avail_samples]
    logger.info(f'Sham样本数: {len(sham_sams)}, Model样本数: {len(model_sams)}')

    if not os.path.exists(PLATFORM_FILE_GPL1355):
        alt_path = os.path.join(dir_path, 'GPL1355-10794 (1).txt')
        plat_file = alt_path if os.path.exists(alt_path) else None
    else:
        plat_file = PLATFORM_FILE_GPL1355

    if plat_file:
        logger.info(f'读取GPL1355平台注释: {plat_file}')
        probe_to_gene = parse_gpl1355_annotation(plat_file)
        logger.info(f'注释映射数: {len(probe_to_gene)}')
        expr_gene = collapse_probes_to_genes(expr_df, probe_to_gene)
        logger.info(f'折叠为基因后维度: {expr_gene.shape}')
    else:
        logger.warning('未找到GPL1355注释文件，直接使用探针ID')
        expr_gene = expr_df

    result = deg_microarray_t_test(expr_gene, model_sams, sham_sams)
    return result


def process_gse97537():
    """处理GSE97537数据集：大鼠芯片，MCAO (24H) vs Sham"""
    logger.info('=' * 60)
    logger.info('开始处理 GSE97537 (24H 大鼠芯片，MCAO vs Sham)')
    dir_path = BASE_DIRS['GSE97537']

    sm_file = find_file(dir_path, ['series_matrix.txt'])
    if not sm_file:
        raise FileNotFoundError(f'未在 {dir_path} 中找到series matrix文件')
    logger.info(f'发现Series Matrix文件: {sm_file}')

    expr_df, meta = parse_series_matrix(sm_file)
    logger.info(f'表达矩阵维度: {expr_df.shape}')

    avail_samples = set(expr_df.columns)
    sham_sams = [s for s in GSE97537_GROUPS['Sham'] if s in avail_samples]
    mcao_sams = [s for s in GSE97537_GROUPS['MCAO'] if s in avail_samples]
    logger.info(f'Sham样本数: {len(sham_sams)}, MCAO样本数: {len(mcao_sams)}')

    plat_dir = BASE_DIRS['GSE97537']
    plat_file = find_file(plat_dir, ['GPL1355'])
    if not plat_file:
        plat_file = PLATFORM_FILE_GPL1355 if os.path.exists(PLATFORM_FILE_GPL1355) else None

    if plat_file:
        logger.info(f'读取GPL1355平台注释: {plat_file}')
        probe_to_gene = parse_gpl1355_annotation(plat_file)
        logger.info(f'注释映射数: {len(probe_to_gene)}')
        expr_gene = collapse_probes_to_genes(expr_df, probe_to_gene)
        logger.info(f'折叠为基因后维度: {expr_gene.shape}')
    else:
        logger.warning('未找到GPL1355注释文件，直接使用探针ID')
        expr_gene = expr_df

    result = deg_microarray_t_test(expr_gene, mcao_sams, sham_sams)
    return result


def process_gse104036():
    """处理GSE104036数据集：小鼠RNA-seq，Ipsilateral (多时序) vs Sham"""
    logger.info('=' * 60)
    logger.info('开始处理 GSE104036 (多时序 小鼠RNA-seq，Ipsilateral vs Sham)')
    dir_path = BASE_DIRS['GSE104036']

    count_file = find_file(dir_path, ['counts.txt'])
    if not count_file:
        raise FileNotFoundError(f'未在 {dir_path} 中找到count matrix文件')
    logger.info(f'发现Count Matrix文件: {count_file}')

    open_func = gzip.open if count_file.endswith('.gz') else open
    with open_func(count_file, 'rt') as f:
        header = f.readline().strip().split('\t')
        n_cols = len(header)
    logger.info(f'Count矩阵列数(含gene): {n_cols}')

    df = pd.read_csv(count_file, sep='\t', dtype=str, low_memory=False)
    df = df.set_index(df.columns[0])
    df = df.apply(pd.to_numeric, errors='coerce').fillna(0).astype(int)

    dup_genes = df.index[df.index.duplicated()]
    if len(dup_genes) > 0:
        logger.info(f'发现 {len(dup_genes)} 个重复基因名，合并计数')
        df = df.groupby(df.index).sum()

    logger.info(f'Count矩阵原始维度: {df.shape}')

    avail_cols = set(df.columns)
    sham_cols = [c for c in GSE104036_COUNT_COLS['Sham'] if c in avail_cols]
    ipsi_cols = [c for c in GSE104036_COUNT_COLS['Ipsilateral'] if c in avail_cols]
    logger.info(f'Sham样本列: {sham_cols}')
    logger.info(f'Ipsilateral样本列: {ipsi_cols}')

    if len(ipsi_cols) < 2 or len(sham_cols) < 2:
        raise ValueError(f'样本列不足: sham {len(sham_cols)}, ipsi {len(ipsi_cols)}')

    study_cols = sham_cols + ipsi_cols
    metadata = pd.DataFrame({
        'condition': [sham_cols[0].replace('S1', 'Sham')] * len(sham_cols) +
                     [ipsi_cols[0].split('_')[0] + '_Ipsilateral'] * len(ipsi_cols)
    }, index=sham_cols + ipsi_cols)

    metadata.loc[sham_cols, 'condition'] = 'Sham'
    metadata.loc[ipsi_cols, 'condition'] = 'Ipsilateral'

    result = deg_rnaseq_pydeseq2(df, metadata, 'Ipsilateral', 'Sham')
    return result


def filter_and_save_deg(deg_df, output_path, lfc_thresh=1.0, p_thresh=0.05):
    """筛选差异基因并保存"""
    required = {'gene_symbol', 'log2FoldChange', 'padj'}
    if not required.issubset(deg_df.columns):
        missing = required - set(deg_df.columns)
        logger.warning(f'缺少列 {missing}，跳过保存 {output_path}')
        return None

    deg_df['abs_log2FC'] = deg_df['log2FoldChange'].abs()
    filtered = deg_df[
        (deg_df['abs_log2FC'] > lfc_thresh) &
        (deg_df['padj'] < p_thresh)
    ].copy()

    filtered = filtered.sort_values('padj')
    out = filtered[['gene_symbol', 'log2FoldChange', 'pvalue', 'padj']]
    out.to_csv(output_path, index=False, encoding='utf-8-sig')
    logger.info(f'保存 {len(out)} 个差异基因 -> {output_path}')
    return out


# ============================================================
# 主流程
# ============================================================

def main():
    logger.info('=' * 60)
    logger.info('GEO差异表达分析管道启动')
    logger.info(f'输出目录: {OUTPUT_DIR}')
    logger.info('=' * 60)

    deg_files = {
        'GSE61616_7d_DEGs.csv': 'GSE61616',
        'GSE97537_24h_DEGs.csv': 'GSE97537',
        'GSE104036_multitime_DEGs.csv': 'GSE104036',
    }

    all_deg_sets = {}

    # 处理GSE61616
    result_61616 = process_gse61616()
    out = filter_and_save_deg(result_61616, os.path.join(OUTPUT_DIR, 'GSE61616_7d_DEGs.csv'))
    if out is not None:
        all_deg_sets['GSE61616'] = set(out['gene_symbol'].dropna().unique())
    else:
        all_deg_sets['GSE61616'] = set()
    logger.info(f'GSE61616 差异基因数: {len(all_deg_sets["GSE61616"])}')

    # 处理GSE97537
    result_97537 = process_gse97537()
    out = filter_and_save_deg(result_97537, os.path.join(OUTPUT_DIR, 'GSE97537_24h_DEGs.csv'))
    if out is not None:
        all_deg_sets['GSE97537'] = set(out['gene_symbol'].dropna().unique())
    else:
        all_deg_sets['GSE97537'] = set()
    logger.info(f'GSE97537 差异基因数: {len(all_deg_sets["GSE97537"])}')

    # 处理GSE104036
    result_104036 = process_gse104036()
    out = filter_and_save_deg(result_104036, os.path.join(OUTPUT_DIR, 'GSE104036_multitime_DEGs.csv'))
    if out is not None:
        all_deg_sets['GSE104036'] = set(out['gene_symbol'].dropna().unique())
    else:
        all_deg_sets['GSE104036'] = set()
    logger.info(f'GSE104036 差异基因数: {len(all_deg_sets["GSE104036"])}')

    # 计算交集和并集
    logger.info('=' * 60)
    logger.info('计算三个数据集交集和并集')

    all_symbols = [all_deg_sets[k] for k in ['GSE61616', 'GSE97537', 'GSE104036']]

    common_symbols = set.intersection(*all_symbols) if all(all_symbols) else set()
    union_symbols = set.union(*all_symbols)

    logger.info(f'三个数据集交集差异基因数: {len(common_symbols)}')
    logger.info(f'三个数据集并集差异基因数: {len(union_symbols)}')

    # 保存交集
    common_df = pd.DataFrame({'gene_symbol': sorted(common_symbols)})
    common_df.to_csv(os.path.join(OUTPUT_DIR, 'common_DEGs.csv'), index=False, encoding='utf-8-sig')
    logger.info(f'保存交集 -> common_DEGs.csv ({len(common_df)} 个基因)')

    # 保存并集
    union_df = pd.DataFrame({'gene_symbol': sorted(union_symbols)})
    union_df.to_csv(os.path.join(OUTPUT_DIR, 'union_DEGs.csv'), index=False, encoding='utf-8-sig')
    logger.info(f'保存并集 -> union_DEGs.csv ({len(union_df)} 个基因)')

    # 保存disease_genes_for_GAT.csv（包含交集和并集）
    combined = pd.concat([
        pd.DataFrame({'gene_symbol': sorted(common_symbols), 'source': 'common'}),
        pd.DataFrame({'gene_symbol': sorted(union_symbols - common_symbols), 'source': 'union_only'})
    ], ignore_index=True)
    combined.to_csv(os.path.join(OUTPUT_DIR, 'disease_genes_for_GAT.csv'), index=False, encoding='utf-8-sig')
    logger.info(f'保存GAT输入 -> disease_genes_for_GAT.csv ({len(combined)} 个基因)')

    # 输出摘要
    logger.info('=' * 60)
    logger.info('分析完成！摘要:')
    logger.info(f'  GSE61616 (7d, 大鼠芯片): {len(all_deg_sets["GSE61616"])} DEGs')
    logger.info(f'  GSE97537 (24H, 大鼠芯片): {len(all_deg_sets["GSE97537"])} DEGs')
    logger.info(f'  GSE104036 (多时序, 小鼠RNA-seq): {len(all_deg_sets["GSE104036"])} DEGs')
    logger.info(f'  交集 (common): {len(common_symbols)} 个基因')
    logger.info(f'  并集 (union): {len(union_symbols)} 个基因')
    logger.info(f'  GAT输入总基因: {len(combined)} 个')
    logger.info('=' * 60)


if __name__ == '__main__':
    main()