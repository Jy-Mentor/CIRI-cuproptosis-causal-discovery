#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GEO差异表达分时序整合管道
对GSE104036按3hr/6hr/12hr/24hr分时序分析，与GSE61616/GSE97537交叉整合
输出汇总Excel（多Sheet）

输出: deg_temporal_summary.xlsx
"""

import pandas as pd
import numpy as np
import gzip
import os
import io
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
    'GSE104036': r'D:\反向网络药理学\L1 数据集\bulk\GSE104036（多时序）',
}
OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))
PLATFORM_FILE_GPL1355 = os.path.join(BASE_DIRS['GSE61616'], 'GPL1355-10794 (1).txt')

# ============================================================
# 样本分组映射（已验证）
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
# 解析工具（复用已有逻辑）
# ============================================================

def find_file(dir_path, patterns):
    for f in os.listdir(dir_path):
        for pat in patterns:
            if pat.lower() in f.lower():
                return os.path.join(dir_path, f)
    return None

def parse_series_matrix(filepath):
    open_func = gzip.open if filepath.endswith('.gz') else open
    with open_func(filepath, 'rt', encoding='latin-1') as f:
        content = f.read()
    data_start = content.find('!series_matrix_table_begin')
    data_end = content.find('!series_matrix_table_end')
    table_text = content[data_start:data_end].replace('!series_matrix_table_begin', '').strip()
    df = pd.read_csv(io.StringIO(table_text), sep='\t', quoting=1, dtype=str, low_memory=False)
    if 'ID_REF' in df.columns:
        df = df.set_index('ID_REF')
    else:
        df = df.set_index(df.columns[0])
    df = df.apply(pd.to_numeric, errors='coerce')
    return df

def parse_gpl1355_annotation(filepath):
    mapping = {}
    with open(filepath, 'r', encoding='latin-1') as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            parts = line.strip().split('\t')
            if len(parts) < 11:
                continue
            probe_id = parts[0].strip()
            gene_symbol = parts[10].strip()
            if gene_symbol and gene_symbol != '---':
                mapping[probe_id] = gene_symbol.split('///')[0].strip()
    return mapping

def collapse_probes_to_genes(expr_df, probe_to_gene):
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
            best_probe = sub.mean(axis=1).idxmax()
            gene_rows.append((gene, expr_mapped.loc[best_probe]))
    return pd.DataFrame([r[1] for r in gene_rows], index=[r[0] for r in gene_rows])

def deg_microarray_t_test(expr_df, case_samples, control_samples):
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
        results.append({'gene_symbol': gene, 'log2FoldChange': log2fc,
                        'stat': stat, 'pvalue': pval})
    res_df = pd.DataFrame(results)
    if res_df.empty:
        return res_df
    _, padj, _, _ = multipletests(res_df['pvalue'], method='fdr_bh')
    res_df['padj'] = padj
    return res_df.sort_values('pvalue')

def deg_rnaseq_pydeseq2(counts_df, metadata, case_label, control_label):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    samples = metadata.index.tolist()
    counts_sub = counts_df.loc[:, samples].astype(int)
    counts_sub = counts_sub.loc[~(counts_sub == 0).all(axis=1)]
    dds = DeseqDataSet(counts=counts_sub.T, metadata=metadata, design='~condition')
    dds.deseq2()
    stat_res = DeseqStats(dds, contrast=['condition', case_label, control_label])
    stat_res.summary()
    result = stat_res.results_df.copy().reset_index()
    first_col = result.columns[0]
    if first_col != 'gene_symbol':
        result = result.rename(columns={first_col: 'gene_symbol'})
    return result[['gene_symbol', 'log2FoldChange', 'pvalue', 'padj']].dropna()

def filter_deg(deg_df):
    if deg_df is None or deg_df.empty:
        return set()
    deg_df = deg_df.copy()
    deg_df['abs_log2FC'] = deg_df['log2FoldChange'].abs()
    filtered = deg_df[(deg_df['abs_log2FC'] > DEG_THRESHOLDS['lfc']) &
                      (deg_df['padj'] < DEG_THRESHOLDS['padj'])]
    return set(filtered['gene_symbol'].dropna().unique())

# ============================================================
# 数据集处理
# ============================================================

def load_gse104036_full_counts():
    """加载GSE104036完整count矩阵，去重"""
    dir_path = BASE_DIRS['GSE104036']
    count_file = find_file(dir_path, ['counts.txt'])
    df = pd.read_csv(count_file, sep='\t', dtype=str, low_memory=False)
    df = df.set_index(df.columns[0])
    df = df.apply(pd.to_numeric, errors='coerce').fillna(0).astype(int)
    dup = df.index[df.index.duplicated()]
    if len(dup) > 0:
        df = df.groupby(df.index).sum()
    return df

def process_microarray(dataset_key):
    """统一处理芯片数据（GSE61616或GSE97537）"""
    dir_path = BASE_DIRS[dataset_key]
    sm_file = find_file(dir_path, ['series_matrix.txt'])
    logger.info(f'[{dataset_key}] 加载 Series Matrix: {os.path.basename(sm_file)}')
    expr_df = parse_series_matrix(sm_file)
    logger.info(f'  表达矩阵: {expr_df.shape}')

    groups = GSE61616_GROUPS if dataset_key == 'GSE61616' else GSE97537_GROUPS
    control_lbl, case_lbl = ('Sham', 'Model') if dataset_key == 'GSE61616' else ('Sham', 'MCAO')

    avail = set(expr_df.columns)
    control_sams = [s for s in groups[control_lbl] if s in avail]
    case_sams = [s for s in groups[case_lbl] if s in avail]
    logger.info(f'  对照{len(control_sams)}个, 处理{len(case_sams)}个')

    plat_file = find_file(dir_path, ['GPL1355'])
    if not plat_file and os.path.exists(PLATFORM_FILE_GPL1355):
        plat_file = PLATFORM_FILE_GPL1355

    if plat_file:
        probe_map = parse_gpl1355_annotation(plat_file)
        expr_gene = collapse_probes_to_genes(expr_df, probe_map)
        logger.info(f'  探针→基因: {expr_gene.shape}')
    else:
        expr_gene = expr_df

    deg = deg_microarray_t_test(expr_gene, case_sams, control_sams)
    genes = filter_deg(deg)
    logger.info(f'  DEGs: {len(genes)}')
    return deg, genes

def process_gse104036_temporal(counts_df):
    """分时序处理GSE104036，返回每个时间点的DEG结果"""
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

        logger.info(f'  [GSE104036_{t}] Sham={len(sham_cols)}, Ipsi={len(ipsi_cols)}')
        deg = deg_rnaseq_pydeseq2(counts_df, meta, f'Ipsi_{t}', 'Sham')
        genes = filter_deg(deg)
        results[t] = {'deg': deg, 'genes': genes, 'n': len(genes)}
        logger.info(f'  DEGs @{t}: {len(genes)}')

    sample_cols = sham_cols + all_time_ipsi
    meta_all = pd.DataFrame(index=sample_cols)
    meta_all['condition'] = ['Sham'] * len(sham_cols) + ['Ipsilateral'] * len(all_time_ipsi)
    logger.info(f'  [GSE104036_all] Sham={len(sham_cols)}, Ipsi={len(all_time_ipsi)}')
    deg_all = deg_rnaseq_pydeseq2(counts_df, meta_all, 'Ipsilateral', 'Sham')
    genes_all = filter_deg(deg_all)
    results['all'] = {'deg': deg_all, 'genes': genes_all, 'n': len(genes_all)}
    logger.info(f'  DEGs @all: {len(genes_all)}')

    return results

# ============================================================
# 交集/并集组合计算
# ============================================================

def compute_combinations(deg_61616, deg_97537, deg_104036_time):
    """计算所有组合的交集/并集"""
    combos = {}

    sets = {}
    sets['GSE61616'] = deg_61616
    sets['GSE97537'] = deg_97537
    for t in ['all', '3hr', '6hr', '12hr', '24hr']:
        sets[f'GSE104036_{t}'] = deg_104036_time[t]['genes']

    # --- 单数据集 ---
    for name, s in sets.items():
        combos[name] = s

    # --- GSE61616 vs GSE97537 ---
    combos['Intersect_61616_97537'] = sets['GSE61616'] & sets['GSE97537']
    combos['Union_61616_97537'] = sets['GSE61616'] | sets['GSE97537']

    # --- GSE61616 vs GSE104036 (各时序) ---
    for t in ['all', '3hr', '6hr', '12hr', '24hr']:
        combos[f'Intersect_61616_104036_{t}'] = sets['GSE61616'] & sets[f'GSE104036_{t}']
        combos[f'Union_61616_104036_{t}'] = sets['GSE61616'] | sets[f'GSE104036_{t}']

    # --- GSE97537 vs GSE104036 (各时序) ---
    for t in ['all', '3hr', '6hr', '12hr', '24hr']:
        combos[f'Intersect_97537_104036_{t}'] = sets['GSE97537'] & sets[f'GSE104036_{t}']
        combos[f'Union_97537_104036_{t}'] = sets['GSE97537'] | sets[f'GSE104036_{t}']

    # --- 三数据集交集 (各时序) ---
    for t in ['all', '3hr', '6hr', '12hr', '24hr']:
        combos[f'Intersect_3datasets_{t}'] = (
            sets['GSE61616'] & sets['GSE97537'] & sets[f'GSE104036_{t}']
        )

    # --- 三数据集并集 (各时序) ---
    for t in ['all', '3hr', '6hr', '12hr', '24hr']:
        combos[f'Union_3datasets_{t}'] = (
            sets['GSE61616'] | sets['GSE97537'] | sets[f'GSE104036_{t}']
        )

    # --- GSE104036 时序内部组合 ---
    time_sets = [sets[f'GSE104036_{t}'] for t in ['3hr', '6hr', '12hr', '24hr']]
    combos['GSE104036_common_4time'] = set.intersection(*time_sets) if all(time_sets) else set()
    combos['GSE104036_union_4time'] = set.union(*time_sets)

    time_3way_sets = [sets[f'GSE104036_{t}'] for t in ['3hr', '6hr', '12hr']]
    combos['GSE104036_common_3hr_6hr_12hr'] = set.intersection(*time_3way_sets) if all(time_3way_sets) else set()

    time_3way_sets = [sets[f'GSE104036_{t}'] for t in ['6hr', '12hr', '24hr']]
    combos['GSE104036_common_6hr_12hr_24hr'] = set.intersection(*time_3way_sets) if all(time_3way_sets) else set()

    # --- GSE104036 时序特有基因 (仅在该时间点出现) ---
    for t in ['3hr', '6hr', '12hr', '24hr']:
        others = set.union(*[sets[f'GSE104036_{o}'] for o in ['3hr', '6hr', '12hr', '24hr'] if o != t])
        combos[f'GSE104036_unique_{t}'] = sets[f'GSE104036_{t}'] - others

    return combos

# ============================================================
# Excel输出
# ============================================================

def write_excel(deg_results, combos, output_path):
    """将DEG详情和组合结果写入Excel"""
    logger.info(f'写入Excel: {output_path}')
    writer = pd.ExcelWriter(output_path, engine='openpyxl')
    deg_104036_data = deg_results['GSE104036']

    # --- Sheet 1: Summary ---
    summary_rows = []
    summary_rows.append(('GSE61616 (7d, 大鼠芯片, Model vs Sham)', len(combos['GSE61616'])))

    summary_rows.append(('GSE97537 (24H, 大鼠芯片, MCAO vs Sham)', len(combos['GSE97537'])))

    for t in ['all', '3hr', '6hr', '12hr', '24hr']:
        label = {'all': '全部时序', '3hr': '3hr', '6hr': '6hr', '12hr': '12hr', '24hr': '24hr'}[t]
        summary_rows.append((f'GSE104036_{label} (Ipsilateral vs Sham)', len(combos[f'GSE104036_{t}'])))

    summary_rows.append(('--- GSE104036 时序内部 ---', ''))
    summary_rows.append(('GSE104036 4时间点交集', len(combos['GSE104036_common_4time'])))
    summary_rows.append(('GSE104036 4时间点并集', len(combos['GSE104036_union_4time'])))
    summary_rows.append(('GSE104036 3hr+6hr+12hr交集', len(combos['GSE104036_common_3hr_6hr_12hr'])))
    summary_rows.append(('GSE104036 6hr+12hr+24hr交集', len(combos['GSE104036_common_6hr_12hr_24hr'])))
    for t in ['3hr', '6hr', '12hr', '24hr']:
        summary_rows.append((f'GSE104036 特有_{t}', len(combos[f'GSE104036_unique_{t}'])))

    summary_rows.append(('--- 两数据集交叉 ---', ''))
    summary_rows.append(('GSE61616 ∩ GSE97537', len(combos['Intersect_61616_97537'])))
    for t in ['all', '3hr', '6hr', '12hr', '24hr']:
        summary_rows.append((f'GSE61616 ∩ GSE104036_{t}', len(combos[f'Intersect_61616_104036_{t}'])))
        summary_rows.append((f'GSE97537 ∩ GSE104036_{t}', len(combos[f'Intersect_97537_104036_{t}'])))

    summary_rows.append(('--- 三数据集交叉 ---', ''))
    for t in ['all', '3hr', '6hr', '12hr', '24hr']:
        summary_rows.append((f'3数据集交集_{t}', len(combos[f'Intersect_3datasets_{t}'])))
        summary_rows.append((f'3数据集并集_{t}', len(combos[f'Union_3datasets_{t}'])))

    summary_df = pd.DataFrame(summary_rows, columns=['组合', '基因数'])
    summary_df.to_excel(writer, sheet_name='Summary', index=False)

    # --- Sheet 2-6: GSE104036 分时序 DEG 详情 ---
    for t in ['all', '3hr', '6hr', '12hr', '24hr']:
        label = {'all': 'all', '3hr': '3hr', '6hr': '6hr', '12hr': '12hr', '24hr': '24hr'}[t]
        deg_df = deg_104036_data[t]['deg'].copy()
        if 'abs_log2FC' in deg_df.columns:
            deg_df = deg_df.drop(columns=['abs_log2FC'])
        deg_df['significant'] = ((deg_df['log2FoldChange'].abs() > DEG_THRESHOLDS['lfc']) &
                                 (deg_df['padj'] < DEG_THRESHOLDS['padj']))
        sheet_name = f'GSE104036_{label}_DEGs'[:31]
        deg_df.to_excel(writer, sheet_name=sheet_name, index=False)

    # --- Sheet 7: GSE61616 DEG 详情 ---
    deg_61616_df = deg_results['GSE61616'].copy()
    deg_61616_df['significant'] = ((deg_61616_df['log2FoldChange'].abs() > DEG_THRESHOLDS['lfc']) &
                                    (deg_61616_df['padj'] < DEG_THRESHOLDS['padj']))
    deg_61616_df.to_excel(writer, sheet_name='GSE61616_DEGs', index=False)

    # --- Sheet 8: GSE97537 DEG 详情 ---
    deg_97537_df = deg_results['GSE97537'].copy()
    deg_97537_df['significant'] = ((deg_97537_df['log2FoldChange'].abs() > DEG_THRESHOLDS['lfc']) &
                                    (deg_97537_df['padj'] < DEG_THRESHOLDS['padj']))
    deg_97537_df.to_excel(writer, sheet_name='GSE97537_DEGs', index=False)

    # --- 组合输出（仅基因列表） ---
    combo_sheet_config = [
        ('Intersect_61616_97537', '二交_61616_97537'),
        ('GSE104036_common_4time', '104036_4时间点交集'),
        ('GSE104036_common_3hr_6hr_12hr', '104036_3hr_6hr_12hr交集'),
        ('GSE104036_common_6hr_12hr_24hr', '104036_6hr_12hr_24hr交集'),
        ('GSE104036_union_4time', '104036_4时间点并集'),
    ]
    for t in ['all', '3hr', '6hr', '12hr', '24hr']:
        combo_sheet_config.append((f'Intersect_3datasets_{t}', f'三交集_{t}'))
        combo_sheet_config.append((f'Union_3datasets_{t}', f'三并集_{t}'))
        for ds in ['61616', '97537']:
            combo_sheet_config.append((f'Intersect_{ds}_104036_{t}', f'{ds}∩104036_{t}'))

    for key, sheet_name in combo_sheet_config:
        genes = sorted(combos[key])
        df = pd.DataFrame({'gene_symbol': genes})
        df.to_excel(writer, sheet_name=sheet_name[:31], index=False)

    # --- 时序特有基因 ---
    for t in ['3hr', '6hr', '12hr', '24hr']:
        genes = sorted(combos[f'GSE104036_unique_{t}'])
        df = pd.DataFrame({'gene_symbol': genes})
        ws_name = f'仅_{t}期'[:31]
        df.to_excel(writer, sheet_name=ws_name, index=False)

    writer.close()
    logger.info('Excel写入完成')

# ============================================================
# 主流程
# ============================================================

def main():
    logger.info('=' * 60)
    logger.info('分时序差异基因整合管道启动')
    logger.info('=' * 60)

    deg_results = {}
    combos = {}

    # Step 1: GSE61616
    logger.info('-' * 50)
    deg_61616, genes_61616 = process_microarray('GSE61616')
    deg_results['GSE61616'] = deg_61616

    # Step 2: GSE97537
    logger.info('-' * 50)
    deg_97537, genes_97537 = process_microarray('GSE97537')
    deg_results['GSE97537'] = deg_97537

    # Step 3: GSE104036 分时序
    logger.info('-' * 50)
    logger.info('[GSE104036] 加载Count矩阵...')
    counts_df = load_gse104036_full_counts()
    deg_104036_time = process_gse104036_temporal(counts_df)
    deg_results['GSE104036'] = deg_104036_time

    # Step 4: 计算所有组合
    logger.info('-' * 50)
    logger.info('计算交集/并集组合...')
    combos = compute_combinations(genes_61616, genes_97537, deg_104036_time)

    for name, s in combos.items():
        logger.info(f'  {name}: {len(s)} 个基因')

    # Step 5: 输出Excel
    logger.info('-' * 50)
    output_path = os.path.join(OUTPUT_DIR, 'deg_temporal_summary.xlsx')
    write_excel(deg_results, combos, output_path)
    logger.info(f'Excel文件: {output_path}')

    # Final summary
    logger.info('=' * 60)
    logger.info('完成！关键结果:')
    logger.info(f'  GSE61616: {len(genes_61616)} DEGs')
    logger.info(f'  GSE97537:  {len(genes_97537)} DEGs')
    for t in ['all', '3hr', '6hr', '12hr', '24hr']:
        logger.info(f'  GSE104036_{t}: {len(deg_104036_time[t]["genes"])} DEGs')
    logger.info(f'  3数据集交集(all): {len(combos["Intersect_3datasets_all"])}')
    logger.info(f'  GSE104036 4时间点交集: {len(combos["GSE104036_common_4time"])}')
    logger.info(f'  GSE104036 4时间点并集: {len(combos["GSE104036_union_4time"])}')
    logger.info('=' * 60)


if __name__ == '__main__':
    main()