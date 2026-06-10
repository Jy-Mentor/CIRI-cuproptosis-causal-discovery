#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
5个GEO数据集：分时序DEG分析 + 交集/并集 / 转人类靶点 → 汇总Excel

数据集：
  GSE16561  (人, Illumina GPL6883): 39 Stroke vs 24 Control
  GSE37587  (人, Illumina GPL6883): 34 Baseline vs 34 Follow-Up (24-48h, 配对)
  GSE61616  (大鼠, Affy GPL1355):   5 Sham / 5 Model / 5 XST (7d)
  GSE97537  (大鼠, Affy GPL1355):   5 Sham / 7 MCAO (24h)
  GSE104036 (小鼠, RNA-seq):        3hr / 6hr / 12hr / 24hr Ipsi vs Sham

输出: deg_5datasets_summary.xlsx (多Sheet)
"""

import pandas as pd
import numpy as np
import gzip
import os
import io
import warnings
import logging
from scipy import stats
from scripts.geo_data_processor import GEODataProcessor

warnings.filterwarnings('ignore')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# ============================================================
# 路径配置
# ============================================================
BASE_DIRS = {
    'GSE16561': r'D:\反向网络药理学\L1 数据集\bulk\GSE16561',
    'GSE37587': r'D:\反向网络药理学\L1 数据集\bulk\GSE37587',
    'GSE61616': r'D:\反向网络药理学\L1 数据集\bulk\GSE61616（7d）',
    'GSE97537': r'D:\反向网络药理学\L1 数据集\bulk\GSE97537(24H)',
    'GSE104036': r'D:\反向网络药理学\L1 数据集\bulk\GSE104036（多时序）',
}
OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))
GPL6883_ANNOT = os.path.join(OUTPUT_DIR, 'GPL6883.annot.gz')
PLATFORM_FILE_GPL1355 = os.path.join(BASE_DIRS['GSE61616'], 'GPL1355-10794 (1).txt')

# ============================================================
# 样本分组映射
# ============================================================
GSE61616_GROUPS = {
    'Sham':  ['GSM1509422', 'GSM1509423', 'GSM1509424', 'GSM1509425', 'GSM1509426'],
    'Model': ['GSM1509427', 'GSM1509428', 'GSM1509429', 'GSM1509430', 'GSM1509431'],
}
GSE97537_GROUPS = {
    'Sham': ['GSM2571742', 'GSM2571743', 'GSM2571744', 'GSM2571745', 'GSM2571746'],
    'MCAO': ['GSM2571735', 'GSM2571736', 'GSM2571737', 'GSM2571738',
             'GSM2571739', 'GSM2571740', 'GSM2571741'],
}
GSE104036_TIME_MAP = {
    'Sham': ['S1', 'S2', 'S3'],
    '3hr':  ['I1_3hr', 'I2_3hr', 'I3_3hr'],
    '6hr':  ['I1_6hr', 'I2_6hr', 'I3_6hr'],
    '12hr': ['I1_12hr', 'I2_12hr', 'I3_12hr'],
    '24hr': ['I1_24hr', 'I2_24hr', 'I3_24hr'],
}

DEG_THRESHOLDS = {'lfc': 1.0, 'padj': 0.05}

# ============================================================
# Illumina GPL6883 探针注释
# ============================================================
def parse_gpl6883_annotation(annot_path):
    probe_map = {}
    with gzip.open(annot_path, 'rt', encoding='latin-1') as f:
        in_table = False
        header = None
        for line in f:
            l = line.strip()
            if l == '!platform_table_begin':
                in_table = True
                header = f.readline().strip().split('\t')
                if 'Gene symbol' in header:
                    gs_idx = header.index('Gene symbol')
                elif 'Symbol' in header:
                    gs_idx = header.index('Symbol')
                else:
                    gs_idx = 2
                continue
            if not in_table:
                continue
            if l == '':
                continue
            fields = l.split('\t')
            if len(fields) <= max(gs_idx, 0):
                continue
            probe_id = fields[0].strip('"').strip()
            gene_symbol = fields[gs_idx].strip('"').strip()
            if gene_symbol and gene_symbol != '':
                probe_map[probe_id] = gene_symbol
    logger.info(f'  GPL6883: {len(probe_map)} 探针注释加载')
    return probe_map

def collapse_illumina_probes(expr_df, probe_map):
    expr_df = expr_df.copy()
    mapped = expr_df[expr_df.index.isin(probe_map.keys())]
    if mapped.empty:
        return expr_df
    gene_series = pd.Series(mapped.index.map(probe_map), index=mapped.index)
    gene_symbols = gene_series.fillna(gene_series.index.to_series())
    mapped.index = gene_symbols
    mapped = mapped.groupby(mapped.index).max()
    return mapped

# ============================================================
# 通用工具
# ============================================================
def find_file(dir_path, keywords):
    if not os.path.isdir(dir_path):
        return None
    for f in os.listdir(dir_path):
        fl = f.lower()
        if all(k.lower() in fl for k in keywords):
            return os.path.join(dir_path, f)
    return None

def deg_t_test(expr_df, case_cols, control_cols):
    results = []
    case_data = expr_df[case_cols].apply(pd.to_numeric, errors='coerce').values
    ctrl_data = expr_df[control_cols].apply(pd.to_numeric, errors='coerce').values
    for i, gene in enumerate(expr_df.index):
        c = case_data[i, :]
        t = ctrl_data[i, :]
        c = c[~np.isnan(c)]
        t = t[~np.isnan(t)]
        if len(c) < 2 or len(t) < 2:
            results.append({'gene_symbol': gene, 'log2FoldChange': np.nan,
                            'pvalue': 1.0, 'padj': 1.0})
            continue
        mean_c = np.mean(c)
        mean_t = np.mean(t)
        log2fc = mean_c - mean_t
        try:
            t_stat, p_val = stats.ttest_ind(c, t, equal_var=False)
        except:
            p_val = 1.0
        results.append({'gene_symbol': gene, 'log2FoldChange': log2fc,
                        'pvalue': p_val, 'padj': np.nan})
    df = pd.DataFrame(results)
    pvals = df['pvalue'].values
    n = len(pvals)
    ranks = np.zeros(n)
    sort_idx = np.argsort(pvals)
    ranks[sort_idx] = np.arange(1, n + 1)
    padj = pvals * n / ranks
    padj = np.minimum(padj, 1.0)
    for i in range(n - 1, 0, -1):
        if padj[i] > padj[i - 1]:
            padj[i] = padj[i - 1]
    df['padj'] = padj
    return df

def deg_paired_t_test(expr_df, case_cols, control_cols):
    """配对t检验"""
    results = []
    case_data = expr_df[case_cols].apply(pd.to_numeric, errors='coerce').values
    ctrl_data = expr_df[control_cols].apply(pd.to_numeric, errors='coerce').values
    for i, gene in enumerate(expr_df.index):
        c = case_data[i, :]
        t = ctrl_data[i, :]
        valid = ~np.isnan(c) & ~np.isnan(t)
        if valid.sum() < 3:
            results.append({'gene_symbol': gene, 'log2FoldChange': np.nan,
                            'pvalue': 1.0, 'padj': 1.0})
            continue
        diff = c[valid] - t[valid]
        log2fc = np.mean(diff)
        try:
            t_stat, p_val = stats.ttest_1samp(diff, 0)
        except:
            p_val = 1.0
        results.append({'gene_symbol': gene, 'log2FoldChange': log2fc,
                        'pvalue': p_val, 'padj': np.nan})
    df = pd.DataFrame(results)
    pvals = df['pvalue'].values
    n = len(pvals)
    ranks = np.zeros(n)
    sort_idx = np.argsort(pvals)
    ranks[sort_idx] = np.arange(1, n + 1)
    padj = pvals * n / ranks
    padj = np.minimum(padj, 1.0)
    for i in range(n - 1, 0, -1):
        if padj[i] > padj[i - 1]:
            padj[i] = padj[i - 1]
    df['padj'] = padj
    return df

def filter_deg(deg_df):
    if deg_df is None or deg_df.empty:
        return set()
    deg_df = deg_df.copy()
    deg_df['abs_log2FC'] = deg_df['log2FoldChange'].abs()
    filtered = deg_df[(deg_df['abs_log2FC'] > DEG_THRESHOLDS['lfc']) &
                      (deg_df['padj'] < DEG_THRESHOLDS['padj'])]
    return set(filtered['gene_symbol'].dropna().unique())

# ============================================================
# 数据集处理：GSE16561 (人, Stroke vs Control)
# ============================================================
def process_gse16561():
    dir_path = BASE_DIRS['GSE16561']
    sm_file = find_file(dir_path, ['series_matrix'])
    logger.info(f'[GSE16561] 加载: {os.path.basename(sm_file)}')
    expr_df = GEODataProcessor.parse_series_matrix(sm_file)
    logger.info(f'  矩阵: {expr_df.shape}')

    with gzip.open(sm_file, 'rt', encoding='latin-1') as f:
        lines = f.readlines()
    desc_line = None
    sample_ids_line = None
    for l in lines:
        if l.startswith('!Sample_description'):
            desc_line = l.strip().split('\t')
        if l.startswith('!Sample_geo_accession'):
            sample_ids_line = l.strip().split('\t')

    stroke_cols = []
    control_cols = []
    for i, gsm_id in enumerate(sample_ids_line[1:], 1):
        gsm = gsm_id.strip('"').strip()
        desc = desc_line[i].strip('"').strip()
        if 'Stroke' in desc:
            stroke_cols.append(gsm)
        elif 'Control' in desc or 'control' in desc.lower():
            control_cols.append(gsm)

    avail = set(expr_df.columns)
    stroke_cols = [c for c in stroke_cols if c in avail]
    control_cols = [c for c in control_cols if c in avail]
    logger.info(f'  Stroke={len(stroke_cols)}, Control={len(control_cols)}')

    # GPL6883 注释
    probe_map = parse_gpl6883_annotation(GPL6883_ANNOT)
    expr_df.columns = [c.strip('"').strip() for c in expr_df.columns]
    expr_df.index = [str(idx).strip('"').strip() for idx in expr_df.index]
    expr_gene = collapse_illumina_probes(expr_df, probe_map)
    logger.info(f'  基因水平: {expr_gene.shape}')

    deg = deg_t_test(expr_gene, stroke_cols, control_cols)
    genes = filter_deg(deg)
    logger.info(f'  DEGs: {len(genes)}')
    return deg, genes

# ============================================================
# 数据集处理：GSE37587 (人, Follow-Up vs Baseline, 配对)
# ============================================================
def process_gse37587():
    dir_path = BASE_DIRS['GSE37587']
    sm_file = find_file(dir_path, ['series_matrix'])
    logger.info(f'[GSE37587] 加载: {os.path.basename(sm_file)}')
    expr_df = GEODataProcessor.parse_series_matrix(sm_file)
    logger.info(f'  矩阵: {expr_df.shape}')

    with gzip.open(sm_file, 'rt', encoding='latin-1') as f:
        lines = f.readlines()
    desc_line = None
    sample_ids_line = None
    char_lines = []
    for l in lines:
        if l.startswith('!Sample_description'):
            desc_line = l.strip().split('\t')
        if l.startswith('!Sample_geo_accession'):
            sample_ids_line = l.strip().split('\t')
        if l.startswith('!Sample_characteristics'):
            char_lines.append(l.strip().split('\t'))

    # 获取 patient number
    patient_col = None
    for cl in char_lines:
        for val in cl[1:]:
            if 'patient number' in val.lower():
                patient_col = cl
                break
        if patient_col is not None:
            break

    # 解析: Baseline vs Follow-Up
    baseline_cols = []
    followup_cols = []
    patient_map = {}
    for i, gsm_id in enumerate(sample_ids_line[1:], 1):
        gsm = gsm_id.strip('"').strip()
        desc = desc_line[i].strip('"').strip()
        if 'Baseline' in desc:
            baseline_cols.append(gsm)
        elif 'Follow-Up' in desc or 'Follow' in desc:
            followup_cols.append(gsm)
        if patient_col:
            pn = patient_col[i].strip('"').strip()
            if 'patient number:' in pn.lower():
                pn = pn.split(':', 1)[-1].strip()
            patient_map[gsm] = pn

    avail = set(expr_df.columns)
    baseline_cols = [c for c in baseline_cols if c in avail]
    followup_cols = [c for c in followup_cols if c in avail]

    # 用病人ID进行配对
    baseline_by_pat = {}
    for c in baseline_cols:
        p = patient_map.get(c, c)
        baseline_by_pat[p] = c
    followup_by_pat = {}
    for c in followup_cols:
        p = patient_map.get(c, c)
        followup_by_pat[p] = c
    common_pats = set(baseline_by_pat.keys()) & set(followup_by_pat.keys())
    paired_baseline = [baseline_by_pat[p] for p in sorted(common_pats)]
    paired_followup = [followup_by_pat[p] for p in sorted(common_pats)]
    logger.info(f'  Baseline={len(baseline_cols)}, Follow-Up={len(followup_cols)}, Paired={len(paired_baseline)}')

    # GPL6883 注释
    probe_map = parse_gpl6883_annotation(GPL6883_ANNOT)
    expr_df.columns = [c.strip('"').strip() for c in expr_df.columns]
    expr_df.index = [str(idx).strip('"').strip() for idx in expr_df.index]
    expr_gene = collapse_illumina_probes(expr_df, probe_map)
    logger.info(f'  基因水平: {expr_gene.shape}')

    deg = deg_paired_t_test(expr_gene, paired_followup, paired_baseline)
    genes = filter_deg(deg)
    logger.info(f'  DEGs: {len(genes)}')
    return deg, genes

# ============================================================
# 数据集处理：GSE61616 & GSE97537 (复用现有逻辑)
# ============================================================
def process_microarray(dataset_key):
    dir_path = BASE_DIRS[dataset_key]
    sm_file = find_file(dir_path, ['series_matrix.txt'])
    logger.info(f'[{dataset_key}] 加载: {os.path.basename(sm_file)}')
    expr_df = GEODataProcessor.parse_series_matrix(sm_file)
    logger.info(f'  矩阵: {expr_df.shape}')

    groups = GSE61616_GROUPS if dataset_key == 'GSE61616' else GSE97537_GROUPS
    ctrl, case = ('Sham', 'Model') if dataset_key == 'GSE61616' else ('Sham', 'MCAO')
    avail = set(expr_df.columns)
    ctrl_cols = [s for s in groups[ctrl] if s in avail]
    case_cols = [s for s in groups[case] if s in avail]

    plat_file = find_file(dir_path, ['GPL1355'])
    if not plat_file and os.path.exists(PLATFORM_FILE_GPL1355):
        plat_file = PLATFORM_FILE_GPL1355
    if plat_file:
        probe_map = GEODataProcessor.parse_gpl1355_annotation(plat_file)
        expr_gene = GEODataProcessor.collapse_probes_to_genes(expr_df, probe_map)
    else:
        expr_gene = expr_df

    deg = GEODataProcessor.deg_microarray_t_test(expr_gene, case_cols, ctrl_cols)
    genes = filter_deg(deg)
    logger.info(f'  DEGs: {len(genes)}')
    return deg, genes

# ============================================================
# 数据集处理：GSE104036 (小鼠, RNA-seq, 分时序)
# ============================================================
def load_gse104036_full_counts():
    dir_path = BASE_DIRS['GSE104036']
    count_file = find_file(dir_path, ['counts.txt'])
    df = pd.read_csv(count_file, sep='\t', dtype=str, low_memory=False)
    df = df.set_index(df.columns[0])
    df = df.apply(pd.to_numeric, errors='coerce').fillna(0).astype(int)
    dup = df.index[df.index.duplicated()]
    if len(dup) > 0:
        df = df.groupby(df.index).sum()
    return df

def process_gse104036_temporal(counts_df):
    sham_cols = GSE104036_TIME_MAP['Sham']
    time_labels = ['3hr', '6hr', '12hr', '24hr']
    results = {}
    all_time_ipsi = []
    for t in time_labels:
        ipsi_cols = GSE104036_TIME_MAP[t]
        all_time_ipsi.extend(ipsi_cols)
        sample_cols = sham_cols + ipsi_cols
        meta = pd.DataFrame(index=sample_cols)
        meta['condition'] = ['Sham'] * len(sham_cols) + [f'Ipsi_{t}'] * len(ipsi_cols)
        logger.info(f'  [GSE104036_{t}]')
        deg = GEODataProcessor.deg_rnaseq_pydeseq2(counts_df, meta, f'Ipsi_{t}', 'Sham')
        genes = filter_deg(deg)
        results[t] = {'deg': deg, 'genes': genes, 'n': len(genes)}

    sample_cols = sham_cols + all_time_ipsi
    meta_all = pd.DataFrame(index=sample_cols)
    meta_all['condition'] = ['Sham'] * len(sham_cols) + ['Ipsilateral'] * len(all_time_ipsi)
    logger.info(f'  [GSE104036_all]')
    deg_all = GEODataProcessor.deg_rnaseq_pydeseq2(counts_df, meta_all, 'Ipsilateral', 'Sham')
    genes_all = filter_deg(deg_all)
    results['all'] = {'deg': deg_all, 'genes': genes_all, 'n': len(genes_all)}
    return results

# ============================================================
# 交集/并集组合计算
# ============================================================
def compute_all_combinations(all_gene_sets):
    """all_gene_sets: {name: set_of_genes}"""
    combos = {}

    for name, s in all_gene_sets.items():
        combos[name] = s

    ds_names = list(all_gene_sets.keys())

    # 两两交叉
    for i, n1 in enumerate(ds_names):
        for n2 in ds_names[i+1:]:
            combos[f'{n1}_&_{n2}_intersect'] = all_gene_sets[n1] & all_gene_sets[n2]
            combos[f'{n1}_&_{n2}_union'] = all_gene_sets[n1] | all_gene_sets[n2]

    # 全部交集
    combos['ALL_5_intersect'] = set.intersection(*all_gene_sets.values()) if all_gene_sets else set()
    # 全部并集
    combos['ALL_5_union'] = set.union(*all_gene_sets.values())

    # GSE104036 时序内部
    time_keys = [k for k in ds_names if 'GSE104036_' in k and k != 'GSE104036_all']
    if time_keys:
        combos['104036_4time_common'] = set.intersection(*[all_gene_sets[k] for k in time_keys])
        combos['104036_4time_union'] = set.union(*[all_gene_sets[k] for k in time_keys])
        for t in time_keys:
            others = set.union(*[all_gene_sets[o] for o in time_keys if o != t])
            label = t.replace('GSE104036_', '')
            combos[f'104036_unique_{label}'] = all_gene_sets[t] - others

    # 三数据集交集 (人 + 大鼠 + 小鼠 × 各时序)
    human_keys = [k for k in ds_names if k in ('GSE16561', 'GSE37587')]
    rat_keys = [k for k in ds_names if k in ('GSE61616', 'GSE97537')]
    mouse_keys = [k for k in ds_names if 'GSE104036' in k]

    for mk in mouse_keys:
        t_label = mk.replace('GSE104036_', '')
        for hk in human_keys:
            for rk in rat_keys:
                combos[f'3way_{hk}_{rk}_{mk}_intersect'] = (
                    all_gene_sets.get(hk, set()) & all_gene_sets.get(rk, set()) & all_gene_sets.get(mk, set())
                )

    return combos

# ============================================================
# Excel输出
# ============================================================
def write_excel(deg_results, combos, output_path):
    logger.info(f'写入Excel: {output_path}')
    writer = pd.ExcelWriter(output_path, engine='openpyxl')

    # --- Summary ---
    summary_rows = []
    for ds_name in ['GSE16561', 'GSE37587', 'GSE61616', 'GSE97537',
                     'GSE104036_all', 'GSE104036_3hr', 'GSE104036_6hr',
                     'GSE104036_12hr', 'GSE104036_24hr']:
        if ds_name in combos:
            summary_rows.append((ds_name, len(combos[ds_name])))

    summary_rows.append(('--- 5数据集全局 ---', ''))
    summary_rows.append(('5数据集交集', len(combos.get('ALL_5_intersect', set()))))
    summary_rows.append(('5数据集并集', len(combos.get('ALL_5_union', set()))))

    summary_rows.append(('--- GSE104036 时序内部 ---', ''))
    if '104036_4time_common' in combos:
        summary_rows.append(('4时间点交集', len(combos['104036_4time_common'])))
        summary_rows.append(('4时间点并集', len(combos['104036_4time_union'])))
        for t in ['3hr', '6hr', '12hr', '24hr']:
            key = f'104036_unique_{t}'
            if key in combos:
                summary_rows.append((f'仅{t}期', len(combos[key])))

    summary_rows.append(('--- 两两交集 ---', ''))
    for name, s in combos.items():
        if '_&_' in name and 'intersect' in name:
            summary_rows.append((name.replace('_&_', ' ∩ ').replace('_intersect', ''),
                                 len(s)))

    summary_df = pd.DataFrame(summary_rows, columns=['组合', '基因数'])
    summary_df.to_excel(writer, sheet_name='Summary', index=False)

    # --- 各数据集DEG详情 ---
    for ds_name, deg_df in deg_results.items():
        deg_out = deg_df.copy()
        deg_out['significant'] = ((deg_out['log2FoldChange'].abs() > DEG_THRESHOLDS['lfc']) &
                                   (deg_out['padj'] < DEG_THRESHOLDS['padj']))
        sheet_name = ds_name.replace('GSE104036_', '104036_')[:31]
        deg_out.to_excel(writer, sheet_name=sheet_name, index=False)

    # --- 组合基因列表 ---
    for name, s in combos.items():
        if name in deg_results:
            continue
        genes = sorted(s)
        df = pd.DataFrame({'gene_symbol': genes})
        sheet_name = name.replace('_&_', '-').replace('_intersect', '-交')
        sheet_name = sheet_name.replace('_union', '-并')[:31]
        df.to_excel(writer, sheet_name=sheet_name, index=False)

    writer.close()
    logger.info('Excel写入完成')

# ============================================================
# 主流程
# ============================================================
def main():
    logger.info('=' * 70)
    logger.info('5个GEO数据集 D差异表达 + 分时序整合启动')
    logger.info('=' * 70)

    deg_results = {}
    all_gene_sets = {}

    # === 1: GSE16561 (人, Stroke vs Control) ===
    logger.info('--- [1/5] GSE16561 ---')
    deg_16561, genes_16561 = process_gse16561()
    deg_results['GSE16561'] = deg_16561
    all_gene_sets['GSE16561'] = genes_16561

    # === 2: GSE37587 (人, 配对 Follow-Up vs Baseline) ===
    logger.info('--- [2/5] GSE37587 ---')
    deg_37587, genes_37587 = process_gse37587()
    deg_results['GSE37587'] = deg_37587
    all_gene_sets['GSE37587'] = genes_37587

    # === 3: GSE61616 (大鼠) ===
    logger.info('--- [3/5] GSE61616 ---')
    deg_61616, genes_61616 = process_microarray('GSE61616')
    deg_results['GSE61616'] = deg_61616
    all_gene_sets['GSE61616'] = genes_61616

    # === 4: GSE97537 (大鼠) ===
    logger.info('--- [4/5] GSE97537 ---')
    deg_97537, genes_97537 = process_microarray('GSE97537')
    deg_results['GSE97537'] = deg_97537
    all_gene_sets['GSE97537'] = genes_97537

    # === 5: GSE104036 (小鼠, 分时序) ===
    logger.info('--- [5/5] GSE104036 ---')
    counts_df = load_gse104036_full_counts()
    deg_104036 = process_gse104036_temporal(counts_df)
    deg_results['GSE104036_all'] = deg_104036['all']['deg']
    for t in ['3hr', '6hr', '12hr', '24hr', 'all']:
        deg_results[f'GSE104036_{t}'] = deg_104036[t]['deg']
        all_gene_sets[f'GSE104036_{t}'] = deg_104036[t]['genes']

    # === 计算所有组合 ===
    logger.info('--- 计算组合 ---')
    combos = compute_all_combinations(all_gene_sets)
    for name, s in combos.items():
        logger.info(f'  {name}: {len(s)}')

    # === 输出Excel ===
    output_path = os.path.join(OUTPUT_DIR, 'deg_5datasets_summary.xlsx')
    write_excel(deg_results, combos, output_path)
    logger.info(f'输出: {output_path}')

    # === 最终摘要 ===
    logger.info('=' * 70)
    logger.info('DEG结果摘要:')
    for k, v in all_gene_sets.items():
        logger.info(f'  {k}: {len(v)}')
    logger.info(f'  5数据集交集: {len(combos.get("ALL_5_intersect", set()))}')
    logger.info(f'  5数据集并集: {len(combos.get("ALL_5_union", set()))}')
    logger.info('=' * 70)


if __name__ == '__main__':
    main()