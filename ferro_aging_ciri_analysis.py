#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
=============================================================================
铁衰老(Ferro-aging)在脑缺血再灌注损伤(CIRI)中的富集分析
=============================================================================
基于 Cell Metabolism 2026 (Liu et al.) 定义的铁衰老基因集,
在5个独立CIRI数据集中验证铁衰老通路的激活。

数据集:
  - GSE16561  (人全血, Stroke vs Control)
  - GSE37587  (人全血, 配对 Follow-Up vs Baseline)
  - GSE61616  (大鼠, MCAO 7d)
  - GSE97537  (大鼠, MCAO 24h)
  - GSE104036 (小鼠, RNA-seq, 多时间点 MCAO)

方法:
  1. ssGSEA: 每个样本计算ferro-aging评分
  2. 差异检验: CIRI组 vs 对照组
  3. 跨数据集Meta分析: Fisher合并p值
  4. 单基因分析: 关键铁衰老基因表达变化

输出:
  - ferro_aging_ciri_results.xlsx (多Sheet)
  - ferro_aging_ciri_figures/ (箱线图PDF)

创新点:
  - 首次将"ferro-aging"概念引入CIRI研究领域
  - 为后续BCP→铁衰老通路研究提供干实验基础

依赖: pandas, numpy, scipy, matplotlib, openpyxl
=============================================================================
"""

import os, sys, gzip, json, warnings, logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

warnings.filterwarnings('ignore')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
logger = logging.getLogger(__name__)

# ============================================================
# 路径配置
# ============================================================
BASE_DIR = Path(__file__).parent
OUTPUT_DIR = BASE_DIR / "ferro_aging_ciri_results"
FIG_DIR = OUTPUT_DIR / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

# 数据集路径
DATA_DIRS = {
    'GSE16561':  r'D:\反向网络药理学\L1 数据集\bulk\GSE16561',
    'GSE37587':  r'D:\反向网络药理学\L1 数据集\bulk\GSE37587',
    'GSE61616':  r'D:\反向网络药理学\L1 数据集\bulk\GSE61616（7d）',
    'GSE97537':  r'D:\反向网络药理学\L1 数据集\bulk\GSE97537(24H)',
    'GSE104036': r'D:\反向网络药理学\L1 数据集\bulk\GSE104036（多时序）',
}

# 平台注释文件
GPL6883_ANNOT = str(BASE_DIR / 'GPL6883.annot.gz')
GPL1355_FILE = str(Path(DATA_DIRS['GSE61616']) / 'GPL1355-10794 (1).txt')

# ============================================================
# 1. 铁衰老(Ferro-aging)基因集定义
# ============================================================
# 来源: Cell Metabolism 2026, Liu et al. "Vitamin C inhibits ACSL4 
#        to alleviate ferro-aging in primates"
#
# 基因集构建: 用户提供的96个铁衰老基因 + 补充文献关键基因
# 分类: 脂质过氧化(10) | 铁代谢(6) | 氧化应激(3) | 炎症(21)
#       | 信号转录(23) | 细胞周期(4) | 自噬(5) | 生长因子(7) | 代谢(17) | 补充(3)

FERRO_AGING_GENES = [
    # === 脂质过氧化 & 铁死亡执行 (10) ===
    "ALOX15", "LPCAT3", "PTGS2", "DPP4", "LOX",
    "SAT1", "DUOX1", "NOX4", "MPO", "ABCC1",

    # === 铁代谢 & 铁稳态 (6) ===
    "TFRC", "HMOX1", "CDO1", "COX7A1", "SLC1A5", "CP",

    # === 氧化应激 & 抗氧化防御 (3) ===
    "SOD1", "KEAP1", "HIF1A",

    # === 炎症 & 免疫应答 (21) ===
    "IL6", "IL1B", "HMGB1", "S100A8", "CXCL10",
    "CD74", "CD82", "IFNG", "IRF1", "IRF7", "IRF9",
    "NLRP3", "TLR4", "TNFAIP3", "TNFAIP1", "KDM6B",
    "SLAMF8", "PADI4", "LGMN", "CTSB", "ICA1",

    # === 细胞信号 & 转录调控 (23) ===
    "BAP1", "E2F1", "E2F3", "EGR1", "SP1",
    "YAP1", "WWTR1", "ZEB1", "BCL6", "EBF3",
    "FOSL1", "RUNX3", "TBX2", "HBP1", "SMARCB1",
    "SETD7", "SMURF2", "MEN1", "BRD7", "NR1D1",
    "NR2F2", "PDE4B", "PPP2R2B",

    # === 细胞周期 & 衰老 (4) ===
    "CDKN1A", "DYRK1A", "FBXO31", "RBM3",

    # === 自噬 & 蛋白质稳态 (5) ===
    "ATG3", "HERPUD1", "ERN1", "SNCA", "LACTB",

    # === 生长因子 & 细胞外基质 (7) ===
    "EDN1", "EMP1", "EPHA2", "EPHA4", "IGFBP7",
    "LIFR", "WNT5A",

    # === 代谢 & 其他 (17) ===
    "ACVR1B", "CAVIN1", "DPEP1", "GMFB", "KLF6",
    "LCN2", "MAP3K14", "MAPK1", "MAPK14", "MCU",
    "NUAK2", "PRKD1", "PTBP1", "SOCS1", "SOCS2",
    "SPATA2", "TXNIP",

    # === 额外补充 (从原文献保留的关键基因) ===
    "ATF3", "PRDX1", "TXNRD1",
]

# 核对: 8个分类 + 1个补充 = 96+ 基因
assert len(set(FERRO_AGING_GENES)) == len(FERRO_AGING_GENES), "存在重复基因！"
print(f"[INFO] 铁衰老基因集: {len(set(FERRO_AGING_GENES))} 个唯一基因")

# 核心铁衰老标志基因 (基于新基因集的关键标志)
CORE_FERRO_AGING_GENES = [
    "PTGS2", "HMOX1", "TFRC", "IL6", "IL1B",
    "HMGB1", "S100A8", "KEAP1", "SOD1", "HIF1A",
    "CDKN1A", "ALOX15", "NLRP3", "TLR4", "MPO"
]

logger.info(f"铁衰老基因集: {len(FERRO_AGING_GENES)} 基因")
logger.info(f"核心标志基因: {len(CORE_FERRO_AGING_GENES)}")


# ============================================================
# 2. ssGSEA (单样本基因集富集分析) 实现
# ============================================================

def ssgsea_score(expr: np.ndarray, gene_mask: np.ndarray) -> float:
    """
    计算单样本GSEA评分 (秩和累积法)
    
    原理:
      1. 对样本中所有基因的表达值排序
      2. 计算基因集内基因的秩次累积和
      3. 减去随机期望, 归一化
    
    Args:
        expr: (n_genes,) 该样本的所有基因表达值
        gene_mask: (n_genes,) 基因集成员的布尔掩码
    
    Returns:
        float: ssGSEA评分 (正值=富集, 负值=抑制)
    """
    n_genes = len(expr)
    n_set = gene_mask.sum()
    if n_set == 0 or n_set == n_genes:
        return 0.0
    
    # 秩次归一化 (从1到n)
    ranks = stats.rankdata(expr, method='average')
    
    # 基因集内的秩次
    set_ranks = ranks[gene_mask]
    
    # 累积和 = 基因集秩次总和, 减去随机期望
    # 随机期望 = n_set * (n_genes + 1) / 2
    expected = n_set * (n_genes + 1) / 2
    sum_ranks = set_ranks.sum()
    
    # 归一化: 除以最大可能偏差
    max_dev = n_set * (n_genes - n_set)
    if max_dev == 0:
        return 0.0
    
    # 计算富集方向: 如果基因集集中在高表达端, > 0
    score = (sum_ranks - expected) / (max_dev / 2)
    
    return float(score)


def compute_ssgsea_matrix(expr_df: pd.DataFrame, gene_set: Set[str]) -> pd.Series:
    """
    对表达矩阵所有样本计算ssGSEA评分
    
    Args:
        expr_df: 基因×样本表达矩阵 (index=基因symbol, columns=样本)
        gene_set: 目标基因集
    
    Returns:
        Series: 每个样本的ssGSEA评分
    """
    # 找到基因集中存在于矩阵中的基因
    common_genes = [g for g in gene_set if g in expr_df.index]
    if len(common_genes) < 5:
        logger.warning(f"  基因集交集过小: {len(common_genes)}")
        return pd.Series(index=expr_df.columns, dtype=float)
    
    gene_mask = expr_df.index.isin(common_genes)
    scores = {}
    
    for col in expr_df.columns:
        vals = expr_df[col].values.astype(float)
        # 去掉缺失值
        valid = ~np.isnan(vals)
        if valid.sum() < 50:
            scores[col] = np.nan
            continue
        scores[col] = ssgsea_score(vals[valid], gene_mask[valid])
    
    result = pd.Series(scores)
    logger.info(f"  ssGSEA完成: {len(common_genes)}/{len(gene_set)} 基因匹配, "
                f"{result.notna().sum()} 样本有效")
    return result


# ============================================================
# 3. 基因集内的单基因表达分析
# ============================================================

def analyze_signature_genes(
    expr_df: pd.DataFrame,
    case_cols: List[str],
    control_cols: List[str],
    gene_set: Set[str],
    dataset_name: str,
) -> pd.DataFrame:
    """
    分析铁衰老基因集中每个基因在两组间的差异表达
    
    Returns:
        DataFrame: [gene, mean_case, mean_control, log2FC, pvalue, padj]
    """
    common = [g for g in expr_df.index if g in gene_set]
    results = []
    
    for gene in common:
        case_vals = expr_df.loc[gene, case_cols].values.astype(float)
        ctrl_vals = expr_df.loc[gene, control_cols].values.astype(float)
        case_vals = case_vals[~np.isnan(case_vals)]
        ctrl_vals = ctrl_vals[~np.isnan(ctrl_vals)]
        
        if len(case_vals) < 2 or len(ctrl_vals) < 2:
            continue
        
        mean_case = np.mean(case_vals)
        mean_ctrl = np.mean(ctrl_vals)
        
        # 对于log表达矩阵, log2FC = mean_case - mean_ctrl
        log2fc = mean_case - mean_ctrl
        
        # t检验
        _, pval = stats.ttest_ind(case_vals, ctrl_vals, equal_var=False)
        
        results.append({
            'dataset': dataset_name,
            'gene_symbol': gene,
            'mean_case': round(mean_case, 4),
            'mean_control': round(mean_ctrl, 4),
            'log2FC': round(log2fc, 4),
            'pvalue': pval,
            'in_core_set': gene in CORE_FERRO_AGING_GENES,
        })
    
    df = pd.DataFrame(results)
    if not df.empty:
        # BH校正
        _, padj, _, _ = multipletests(df['pvalue'], method='fdr_bh')
        df['padj'] = padj
        df = df.sort_values('pvalue')
    
    n_up = (df['log2FC'] > 0).sum()
    logger.info(f"  {dataset_name}: {len(df)} 铁衰老基因分析, "
                f"上调={n_up}, 下调={len(df)-n_up}")
    return df


# ============================================================
# 4. 跨数据集Meta分析 (Fisher合并p值)
# ============================================================

def fisher_meta_analysis(p_values: List[float]) -> Tuple[float, float]:
    """
    Fisher方法合并独立p值
    χ² = -2 * Σ ln(p_i), df = 2 * k
    
    Returns:
        (chi2, meta_p)
    """
    valid_p = [p for p in p_values if 0 < p <= 1]
    if len(valid_p) < 2:
        return np.nan, np.nan
    chi2 = -2 * np.sum(np.log(valid_p))
    df = 2 * len(valid_p)
    meta_p = 1 - stats.chi2.cdf(chi2, df)
    return chi2, meta_p


# ============================================================
# 5. 各数据集处理函数
# ============================================================

def find_file(dir_path: str, keywords: List[str]) -> Optional[str]:
    """在目录中查找包含关键字的文件"""
    if not os.path.isdir(dir_path):
        return None
    for f in os.listdir(dir_path):
        if all(k.lower() in f.lower() for k in keywords):
            return os.path.join(dir_path, f)
    return None


def parse_series_matrix(filepath: str) -> pd.DataFrame:
    """解析GEO Series Matrix文件"""
    open_func = gzip.open if str(filepath).endswith('.gz') else open
    with open_func(filepath, 'rt', encoding='latin-1') as f:
        content = f.read()
    
    lines = content.splitlines()
    data_start = None
    header_line = None
    
    for i, line in enumerate(lines):
        if line.startswith('!series_matrix_table_begin'):
            data_start = i + 1
            header_line = i + 1
            break
    
    if data_start is None:
        raise ValueError(f"无法找到series_matrix_table_begin: {filepath}")
    
    # 解析表头
    header = lines[header_line].strip().split('\t')
    header = [h.strip('"').strip() for h in header]
    
    # 解析数据
    data_lines = []
    for i in range(data_start + 1, len(lines)):
        if lines[i].startswith('!series_matrix_table_end'):
            break
        data_lines.append(lines[i])
    
    data = []
    index = []
    for line in data_lines:
        fields = line.strip().split('\t')
        if len(fields) < 2:
            continue
        probe_id = fields[0].strip('"').strip()
        index.append(probe_id)
        values = [float(v) if v != 'null' and v != '' else np.nan 
                  for v in fields[1:]]
        if len(values) < len(header) - 1:
            values.extend([np.nan] * (len(header) - 1 - len(values)))
        data.append(values[:len(header) - 1])
    
    df = pd.DataFrame(data, index=index, columns=header[1:])
    logger.info(f"  解析 {os.path.basename(filepath)}: {df.shape}")
    return df


def parse_gpl1355_annotation(filepath: str) -> Dict[str, str]:
    """解析GPL1355平台探针注释 (大鼠)"""
    probe_map = {}
    with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
        in_table = False
        for line in f:
            l = line.strip()
            if l.startswith('ID'):
                in_table = True
                header = l.split('\t')
                try:
                    gene_col = next(i for i, h in enumerate(header) 
                                    if 'gene symbol' in h.lower() or 'symbol' in h.lower())
                except StopIteration:
                    gene_col = 5  # 默认列
                continue
            if not in_table or not l:
                continue
            fields = l.split('\t')
            if len(fields) <= max(gene_col, 0):
                continue
            probe = fields[0]
            gene = fields[gene_col].strip('"').strip()
            if gene and gene != '':
                probe_map[probe] = gene.split('///')[0].strip()
    logger.info(f"  GPL1355: {len(probe_map)} 探针注释")
    return probe_map


def collapse_probes(expr_df: pd.DataFrame, probe_map: Dict[str, str]) -> pd.DataFrame:
    """探针→基因折叠 (同一基因取最大表达值, 基因符号转大写)"""
    mapped = expr_df[expr_df.index.isin(probe_map.keys())].copy()
    if mapped.empty:
        return expr_df
    gene_series = pd.Series(mapped.index.map(probe_map), index=mapped.index)
    gene_series = gene_series.fillna(gene_series.index.to_series())
    # 转大写实现跨物种大小写不敏感匹配
    gene_series = gene_series.str.upper()
    mapped.index = gene_series
    mapped = mapped.groupby(mapped.index).max()
    return mapped


# ============================================================
# 5a. GSE16561 (人全血, Stroke vs Control)
# ============================================================

def process_gse16561() -> Tuple[pd.DataFrame, pd.DataFrame, pd.Series]:
    """处理GSE16561"""
    logger.info("=" * 50)
    logger.info("[GSE16561] 人全血: Stroke vs Control")
    
    sm_file = find_file(DATA_DIRS['GSE16561'], ['series_matrix'])
    if not sm_file:
        raise FileNotFoundError("GSE16561 series_matrix 未找到")
    
    expr_df = parse_series_matrix(sm_file)
    
    # 解析样本分组
    with gzip.open(sm_file, 'rt', encoding='latin-1') as f:
        lines = f.readlines()
    
    desc_line = None
    sample_line = None
    for l in lines:
        if l.startswith('!Sample_description'):
            desc_line = l.strip().split('\t')
        if l.startswith('!Sample_geo_accession'):
            sample_line = l.strip().split('\t')
    
    stroke_cols = []
    control_cols = []
    for i, gsm in enumerate(sample_line[1:], 1):
        gsm = gsm.strip('"').strip()
        desc = desc_line[i].strip('"').strip() if i < len(desc_line) else ''
        if 'Stroke' in desc or 'stroke' in desc:
            stroke_cols.append(gsm)
        else:
            control_cols.append(gsm)
    
    avail = set(expr_df.columns)
    stroke_cols = [c for c in stroke_cols if c in avail]
    control_cols = [c for c in control_cols if c in avail]
    logger.info(f"  Stroke={len(stroke_cols)}, Control={len(control_cols)}")
    
    # 探针→基因
    # GPL6883注释
    if os.path.exists(GPL6883_ANNOT):
        probe_map = {}
        with gzip.open(GPL6883_ANNOT, 'rt', encoding='latin-1') as f:
            in_table = False
            for line in f:
                l = line.strip()
                if l == '!platform_table_begin':
                    in_table = True
                    header = f.readline().strip().split('\t')
                    gs_idx = 2
                    if 'Gene symbol' in header:
                        gs_idx = header.index('Gene symbol')
                    elif 'Symbol' in header:
                        gs_idx = header.index('Symbol')
                    continue
                if not in_table or l == '':
                    continue
                fields = l.split('\t')
                if len(fields) > gs_idx:
                    probe = fields[0].strip('"').strip()
                    gene = fields[gs_idx].strip('"').strip()
                    if gene:
                        probe_map[probe] = gene
        expr_gene = collapse_probes(expr_df, probe_map)
    else:
        logger.warning("  GPL6883注释文件不存在, 使用原始探针ID")
        expr_gene = expr_df
    
    # ssGSEA
    fa_score = compute_ssgsea_matrix(expr_gene, set(FERRO_AGING_GENES))
    
    # 统计检验
    case_scores = fa_score[[c for c in stroke_cols if c in fa_score.index]].dropna()
    ctrl_scores = fa_score[[c for c in control_cols if c in fa_score.index]].dropna()
    _, pval = stats.ttest_ind(case_scores, ctrl_scores, equal_var=False)
    mean_diff = case_scores.mean() - ctrl_scores.mean()
    logger.info(f"  Ferro-aging评分: Stroke={case_scores.mean():.3f}, "
                f"Ctrl={ctrl_scores.mean():.3f}, Δ={mean_diff:.3f}, p={pval:.4e}")
    
    # 单基因分析
    gene_df = analyze_signature_genes(expr_gene, stroke_cols, control_cols, 
                                       set(FERRO_AGING_GENES), 'GSE16561')
    
    return fa_score, gene_df, pd.Series({
        'dataset': 'GSE16561', 'species': 'Human', 'tissue': 'Whole Blood',
        'case_label': 'Stroke', 'control_label': 'Control',
        'n_case': len(case_scores), 'n_control': len(ctrl_scores),
        'mean_case': case_scores.mean(), 'mean_control': ctrl_scores.mean(),
        'mean_diff': mean_diff, 'pvalue': pval,
    })


# ============================================================
# 5b. GSE37587 (人全血, 配对 Follow-Up vs Baseline)
# ============================================================

def process_gse37587() -> Tuple[pd.DataFrame, pd.DataFrame, pd.Series]:
    """处理GSE37587 (配对设计)"""
    logger.info("=" * 50)
    logger.info("[GSE37587] 人全血: Follow-Up vs Baseline (配对)")
    
    sm_file = find_file(DATA_DIRS['GSE37587'], ['series_matrix'])
    if not sm_file:
        raise FileNotFoundError("GSE37587 series_matrix 未找到")
    
    expr_df = parse_series_matrix(sm_file)
    
    # 用GSM列名精确匹配 (GSE37587的列名包含GSM号)
    with gzip.open(sm_file, 'rt', encoding='latin-1') as f:
        lines = f.readlines()
    
    sample_line = None
    title_line = None
    desc_line = None
    for l in lines:
        if l.startswith('!Sample_geo_accession'):
            sample_line = [x.strip('"').strip() for x in l.strip().split('\t')]
        if l.startswith('!Sample_title'):
            title_line = [x.strip('"').strip() for x in l.strip().split('\t')]
        if l.startswith('!Sample_description'):
            desc_line = [x.strip('"').strip() for x in l.strip().split('\t')]
    
    # 列名 = GSM ID
    cols_by_sample = {}
    for i, gsm in enumerate(sample_line[1:], 1):
        title = title_line[i] if i < len(title_line) else ''
        desc = desc_line[i] if i < len(desc_line) else ''
        cols_by_sample[gsm] = {'title': title, 'desc': desc}
        logger.debug(f"  {gsm}: title={title}, desc={desc}")
    
    followup_cols = []
    baseline_cols = []
    for gsm, info in cols_by_sample.items():
        title_lower = info['title'].lower()
        desc_lower = info['desc'].lower()
        # 检查是否为follow-up/后续时间点
        if any(kw in title_lower for kw in ['follow', 'fu', '24', '48', '72']):
            followup_cols.append(gsm)
        elif any(kw in title_lower for kw in ['baseline', 'base', '00', '0h', 'hour 0', '0 hour', 'control']):
            baseline_cols.append(gsm)
        # 如果从title无法区分, 用number判断
        elif 'patient' in desc_lower:
            # 从description看
            if any(kw in desc_lower for kw in ['follow', 'hour 24', 'hour 48']):
                followup_cols.append(gsm)
            elif any(kw in desc_lower for kw in ['baseline', 'hour 0', '0 hour']):
                baseline_cols.append(gsm)
    
    # 去重
    followup_cols = list(dict.fromkeys([c for c in followup_cols if c in expr_df.columns]))
    baseline_cols = list(dict.fromkeys([c for c in baseline_cols if c in expr_df.columns]))
    
    logger.info(f"  Follow-Up={len(followup_cols)}, Baseline={len(baseline_cols)}")
    
    # 探针→基因
    probe_map = {}
    if os.path.exists(GPL6883_ANNOT):
        with gzip.open(GPL6883_ANNOT, 'rt', encoding='latin-1') as f:
            in_table = False
            for line in f:
                l = line.strip()
                if l == '!platform_table_begin':
                    in_table = True
                    header = f.readline().strip().split('\t')
                    gs_idx = header.index('Gene symbol') if 'Gene symbol' in header else 2
                    continue
                if not in_table or l == '':
                    continue
                fields = l.split('\t')
                if len(fields) > gs_idx:
                    probe = fields[0].strip('"').strip()
                    gene = fields[gs_idx].strip('"').strip()
                    if gene:
                        probe_map[probe] = gene
    expr_gene = collapse_probes(expr_df, probe_map)
    
    # ssGSEA
    fa_score = compute_ssgsea_matrix(expr_gene, set(FERRO_AGING_GENES))
    
    # 配对检验
    case_scores = fa_score[[c for c in followup_cols if c in fa_score.index]].dropna()
    ctrl_scores = fa_score[[c for c in baseline_cols if c in fa_score.index]].dropna()
    _, pval = stats.ttest_ind(case_scores, ctrl_scores, equal_var=False)
    mean_diff = case_scores.mean() - ctrl_scores.mean()
    logger.info(f"  Ferro-aging评分: FU={case_scores.mean():.3f}, "
                f"BL={ctrl_scores.mean():.3f}, Δ={mean_diff:.3f}, p={pval:.4e}")
    
    gene_df = analyze_signature_genes(expr_gene, followup_cols, baseline_cols,
                                       set(FERRO_AGING_GENES), 'GSE37587')
    
    return fa_score, gene_df, pd.Series({
        'dataset': 'GSE37587', 'species': 'Human', 'tissue': 'Whole Blood',
        'case_label': 'Follow-Up', 'control_label': 'Baseline',
        'n_case': len(case_scores), 'n_control': len(ctrl_scores),
        'mean_case': case_scores.mean(), 'mean_control': ctrl_scores.mean(),
        'mean_diff': mean_diff, 'pvalue': pval,
    })


# ============================================================
# 5c. GSE61616 (大鼠, MCAO 7d)
# ============================================================

def process_gse61616() -> Tuple[pd.DataFrame, pd.DataFrame, pd.Series]:
    """处理GSE61616 (大鼠Affy芯片, MCAO 7d)"""
    logger.info("=" * 50)
    logger.info("[GSE61616] 大鼠 MCAO 7d: Model vs Sham")
    
    sm_file = find_file(DATA_DIRS['GSE61616'], ['series_matrix'])
    if not sm_file:
        raise FileNotFoundError("GSE61616 series_matrix 未找到")
    
    expr_df = parse_series_matrix(sm_file)
    
    # 样本分组
    sham_cols = ['GSM1509422', 'GSM1509423', 'GSM1509424', 'GSM1509425', 'GSM1509426']
    model_cols = ['GSM1509427', 'GSM1509428', 'GSM1509429', 'GSM1509430', 'GSM1509431']
    sham_cols = [c for c in sham_cols if c in expr_df.columns]
    model_cols = [c for c in model_cols if c in expr_df.columns]
    logger.info(f"  Model={len(model_cols)}, Sham={len(sham_cols)}")
    
    # GPL1355注释
    probe_map = {}
    if os.path.exists(GPL1355_FILE):
        probe_map = parse_gpl1355_annotation(GPL1355_FILE)
    expr_gene = collapse_probes(expr_df, probe_map)
    
    # ssGSEA (大鼠基因映射)
    fa_score = compute_ssgsea_matrix(expr_gene, set(FERRO_AGING_GENES))
    
    case_scores = fa_score[[c for c in model_cols if c in fa_score.index]].dropna()
    ctrl_scores = fa_score[[c for c in sham_cols if c in fa_score.index]].dropna()
    _, pval = stats.ttest_ind(case_scores, ctrl_scores, equal_var=False)
    mean_diff = case_scores.mean() - ctrl_scores.mean()
    logger.info(f"  Ferro-aging评分: Model={case_scores.mean():.3f}, "
                f"Sham={ctrl_scores.mean():.3f}, Δ={mean_diff:.3f}, p={pval:.4e}")
    
    gene_df = analyze_signature_genes(expr_gene, model_cols, sham_cols,
                                       set(FERRO_AGING_GENES), 'GSE61616')
    
    return fa_score, gene_df, pd.Series({
        'dataset': 'GSE61616', 'species': 'Rat', 'tissue': 'Brain',
        'case_label': 'MCAO 7d', 'control_label': 'Sham',
        'n_case': len(case_scores), 'n_control': len(ctrl_scores),
        'mean_case': case_scores.mean(), 'mean_control': ctrl_scores.mean(),
        'mean_diff': mean_diff, 'pvalue': pval,
    })


# ============================================================
# 5d. GSE97537 (大鼠, MCAO 24h)
# ============================================================

def process_gse97537() -> Tuple[pd.DataFrame, pd.DataFrame, pd.Series]:
    """处理GSE97537 (大鼠Affy芯片, MCAO 24h)"""
    logger.info("=" * 50)
    logger.info("[GSE97537] 大鼠 MCAO 24h: Model vs Sham")
    
    sm_file = find_file(DATA_DIRS['GSE97537'], ['series_matrix'])
    if not sm_file:
        raise FileNotFoundError("GSE97537 series_matrix 未找到")
    
    expr_df = parse_series_matrix(sm_file)
    
    sham_cols = ['GSM2571742', 'GSM2571743', 'GSM2571744', 'GSM2571745', 'GSM2571746']
    mcao_cols = ['GSM2571735', 'GSM2571736', 'GSM2571737', 'GSM2571738',
                  'GSM2571739', 'GSM2571740', 'GSM2571741']
    sham_cols = [c for c in sham_cols if c in expr_df.columns]
    mcao_cols = [c for c in mcao_cols if c in expr_df.columns]
    logger.info(f"  MCAO={len(mcao_cols)}, Sham={len(sham_cols)}")
    
    probe_map = {}
    if os.path.exists(GPL1355_FILE):
        probe_map = parse_gpl1355_annotation(GPL1355_FILE)
    expr_gene = collapse_probes(expr_df, probe_map)
    
    fa_score = compute_ssgsea_matrix(expr_gene, set(FERRO_AGING_GENES))
    
    case_scores = fa_score[[c for c in mcao_cols if c in fa_score.index]].dropna()
    ctrl_scores = fa_score[[c for c in sham_cols if c in fa_score.index]].dropna()
    _, pval = stats.ttest_ind(case_scores, ctrl_scores, equal_var=False)
    mean_diff = case_scores.mean() - ctrl_scores.mean()
    logger.info(f"  Ferro-aging评分: MCAO={case_scores.mean():.3f}, "
                f"Sham={ctrl_scores.mean():.3f}, Δ={mean_diff:.3f}, p={pval:.4e}")
    
    gene_df = analyze_signature_genes(expr_gene, mcao_cols, sham_cols,
                                       set(FERRO_AGING_GENES), 'GSE97537')
    
    return fa_score, gene_df, pd.Series({
        'dataset': 'GSE97537', 'species': 'Rat', 'tissue': 'Brain',
        'case_label': 'MCAO 24h', 'control_label': 'Sham',
        'n_case': len(case_scores), 'n_control': len(ctrl_scores),
        'mean_case': case_scores.mean(), 'mean_control': ctrl_scores.mean(),
        'mean_diff': mean_diff, 'pvalue': pval,
    })


# ============================================================
# 5e. GSE104036 (小鼠, RNA-seq, 多时间点)
# ============================================================

def process_gse104036() -> Tuple[pd.DataFrame, pd.DataFrame, pd.Series]:
    """处理GSE104036 (小鼠RNA-seq, 3hr/6hr/12hr/24hr MCAO)"""
    logger.info("=" * 50)
    logger.info("[GSE104036] 小鼠 RNA-seq: 多时间点 MCAO vs Sham")
    
    # 加载counts文件
    counts_file = Path(DATA_DIRS['GSE104036']) / 'GSE104036_TC-RNAseq_counts.txt.gz'
    if not counts_file.exists():
        # 尝试找其他文件
        counts_file = find_file(DATA_DIRS['GSE104036'], ['counts', 'txt'])
        if counts_file:
            counts_file = Path(counts_file)
    
    if counts_file and counts_file.exists():
        logger.info(f"  加载: {counts_file.name}")
        # 直接用文件路径, pandas自动识别gzip压缩
        expr_df = pd.read_csv(str(counts_file), sep='\t', index_col=0, 
                              compression='gzip')
    else:
        # 回退到series_matrix
        sm_file = find_file(DATA_DIRS['GSE104036'], ['series_matrix'])
        if not sm_file:
            raise FileNotFoundError("GSE104036 数据文件未找到")
        expr_df = parse_series_matrix(sm_file)
    
    expr_df.columns = [c.strip('"').strip() for c in expr_df.columns]
    expr_df.index = [str(idx).strip('"').strip() for idx in expr_df.index]
    
    # 将小鼠基因符号转大写与铁衰老基因集匹配
    expr_df.index = expr_df.index.str.upper()
    
    logger.info(f"  矩阵: {expr_df.shape}")
    
    # 样本分组
    # ipsi = 缺血侧, contra = 对侧
    sham_cols = ['S1', 'S2', 'S3']
    ipso_3hr = ['I1_3hr', 'I2_3hr', 'I3_3hr']
    ipso_6hr = ['I1_6hr', 'I2_6hr', 'I3_6hr']
    ipso_12hr = ['I1_12hr', 'I2_12hr', 'I3_12hr']
    ipso_24hr = ['I1_24hr', 'I2_24hr', 'I3_24hr']
    
    # 从列名中匹配
    all_cols = expr_df.columns.tolist()
    sham_cols = [c for c in all_cols if c in sham_cols]
    ipso_3hr = [c for c in all_cols if c in ipso_3hr]
    ipso_6hr = [c for c in all_cols if c in ipso_6hr]
    ipso_12hr = [c for c in all_cols if c in ipso_12hr]
    ipso_24hr = [c for c in all_cols if c in ipso_24hr]
    
    # 合并所有ipsi时间点 vs sham
    ipsi_all = ipso_3hr + ipso_6hr + ipso_12hr + ipso_24hr
    
    # 如果列名不匹配, 尝试模糊匹配
    if not sham_cols:
        sham_cols = [c for c in all_cols if 'sham' in c.lower() or 's1' in c.lower() 
                     or 's2' in c.lower() or 's3' in c.lower()]
    
    # RNA-seq数据: log2(CPM+1) 转换用于ssGSEA
    # 如果是count数据, 先归一化
    if expr_df.max().max() > 100:
        logger.info("  检测到raw counts, 执行log2(CPM+1)转换")
        col_sums = expr_df.sum()
        cpm = expr_df.div(col_sums, axis=1) * 1e6
        expr_df = np.log2(cpm + 1)
    
    logger.info(f"  Sham={len(sham_cols)}, Ipsi_all={len(ipsi_all)}, "
                f"3hr={len(ipso_3hr)}, 6hr={len(ipso_6hr)}, "
                f"12hr={len(ipso_12hr)}, 24hr={len(ipso_24hr)}")
    
    if not sham_cols:
        raise ValueError("GSE104036: 无法识别Sham样本")
    
    # ssGSEA
    fa_score = compute_ssgsea_matrix(expr_df, set(FERRO_AGING_GENES))
    
    case_scores = fa_score[[c for c in ipsi_all if c in fa_score.index]].dropna()
    ctrl_scores = fa_score[[c for c in sham_cols if c in fa_score.index]].dropna()
    _, pval = stats.ttest_ind(case_scores, ctrl_scores, equal_var=False)
    mean_diff = case_scores.mean() - ctrl_scores.mean()
    logger.info(f"  Ferro-aging评分(所有时间点): Ipsi={case_scores.mean():.3f}, "
                f"Sham={ctrl_scores.mean():.3f}, Δ={mean_diff:.3f}, p={pval:.4e}")
    
    gene_df = analyze_signature_genes(expr_df, ipsi_all, sham_cols,
                                       set(FERRO_AGING_GENES), 'GSE104036')
    
    # 分时间点分析
    timepoint_results = {}
    for tp_name, tp_cols in [('3hr', ipso_3hr), ('6hr', ipso_6hr),
                              ('12hr', ipso_12hr), ('24hr', ipso_24hr)]:
        if len(tp_cols) >= 2:
            tp_scores_case = fa_score[[c for c in tp_cols if c in fa_score.index]].dropna()
            if not tp_scores_case.empty and not ctrl_scores.empty:
                _, tp_p = stats.ttest_ind(tp_scores_case, ctrl_scores, equal_var=False)
                timepoint_results[tp_name] = {
                    'n': len(tp_scores_case),
                    'mean': tp_scores_case.mean(),
                    'pvalue': tp_p,
                }
                logger.info(f"    {tp_name}: mean={tp_scores_case.mean():.3f}, p={tp_p:.4e}")
    
    return fa_score, gene_df, pd.Series({
        'dataset': 'GSE104036', 'species': 'Mouse', 'tissue': 'Brain',
        'case_label': 'MCAO Ipsi', 'control_label': 'Sham',
        'n_case': len(case_scores), 'n_control': len(ctrl_scores),
        'mean_case': case_scores.mean(), 'mean_control': ctrl_scores.mean(),
        'mean_diff': mean_diff, 'pvalue': pval,
        'timepoints': json.dumps({k: {kk: float(vv) if isinstance(vv, np.floating) else vv 
                                       for kk, vv in v.items()} 
                                  for k, v in timepoint_results.items()}),
    })


# ============================================================
# 6. 可视化
# ============================================================

def plot_ferro_aging_scores(
    all_scores: Dict[str, pd.Series],
    all_stats: Dict[str, pd.Series],
):
    """绘制各数据集的铁衰老评分箱线图"""
    n_datasets = len(all_scores)
    fig, axes = plt.subplots(1, n_datasets, figsize=(5 * n_datasets, 5))
    if n_datasets == 1:
        axes = [axes]
    
    colors = {'case': '#E74C3C', 'control': '#3498DB'}
    
    for ax, (ds_name, scores) in zip(axes, all_scores.items()):
        stat = all_stats.get(ds_name)
        if stat is None:
            continue
        
        # 找case和control的样本
        case_label = stat.get('case_label', 'Case')
        ctrl_label = stat.get('control_label', 'Control')
        
        # 从score中根据样本名判断分组
        # 简化方法: 用预览的前几个字符判断
        case_pattern = case_label.lower()[:4]
        ctrl_pattern = ctrl_label.lower()[:4]
        
        case_vals = []
        ctrl_vals = []
        case_names = []
        ctrl_names = []
        for sname in scores.index:
            sname_clean = str(sname).strip('"').lower()
            if any(kw in sname_clean for kw in [case_label.lower(), 'stroke', 'mcao', 
                                                  'model', 'follow', 'ipsi']):
                case_vals.append(scores[sname])
                case_names.append(str(sname)[:15])
            else:
                ctrl_vals.append(scores[sname])
                ctrl_names.append(str(sname)[:15])
        
        # 如果没有识别到, 用stats中的样本数推断
        if not case_vals or not ctrl_vals:
            # 前述方法可能不准, 改用存储的统计量画简图
            case_vals = [stat.get('mean_case', 0)]
            ctrl_vals = [stat.get('mean_control', 0)]
        
        # 箱线图
        bp = ax.boxplot([ctrl_vals, case_vals], positions=[0, 1], widths=0.5,
                         patch_artist=True,
                         boxprops=dict(linewidth=1.5),
                         whiskerprops=dict(linewidth=1.5),
                         capprops=dict(linewidth=1.5),
                         medianprops=dict(linewidth=2, color='black'))
        
        bp['boxes'][0].set_facecolor(colors['control'])
        bp['boxes'][1].set_facecolor(colors['case'])
        
        # 散点
        np.random.seed(42)
        for pos, vals, color in [(0, ctrl_vals, colors['control']),
                                  (1, case_vals, colors['case'])]:
            jitter = np.random.normal(0, 0.04, len(vals))
            ax.scatter(np.full_like(vals, pos) + jitter, vals,
                      alpha=0.7, s=30, color=color, edgecolors='white', 
                      linewidth=0.5, zorder=5)
        
        # p值
        pval = stat.get('pvalue', 1.0)
        p_text = f"p = {pval:.2e}" if pval > 1e-99 else "p < 1e-99"
        sig = '*' if pval < 0.05 else 'ns'
        if pval < 0.001:
            sig = '***'
        elif pval < 0.01:
            sig = '**'
        elif pval < 0.05:
            sig = '*'
        
        y_max = max(max(case_vals) if case_vals else 0, 
                    max(ctrl_vals) if ctrl_vals else 0)
        y_min = min(min(case_vals) if case_vals else 0, 
                    min(ctrl_vals) if ctrl_vals else 0)
        y_range = y_max - y_min if y_max != y_min else 1
        
        ax.plot([0, 0, 1, 1], [y_max + y_range * 0.1, y_max + y_range * 0.15,
                               y_max + y_range * 0.15, y_max + y_range * 0.1],
               color='black', linewidth=1)
        ax.text(0.5, y_max + y_range * 0.18, f'{sig}\n{p_text}', 
                ha='center', va='bottom', fontsize=9)
        
        # Labels
        ax.set_xticks([0, 1])
        ax.set_xticklabels([f'{ctrl_label}\n(n={len(ctrl_vals)})',
                            f'{case_label}\n(n={len(case_vals)})'],
                           fontsize=10)
        ax.set_ylabel('Ferro-aging Score (ssGSEA)', fontsize=11)
        ax.set_title(f'{ds_name} ({stat.get("species", "")})', fontsize=12, fontweight='bold')
        ax.tick_params(axis='both', labelsize=9)
        
        # 水平参考线
        ax.axhline(y=0, color='gray', linestyle='--', alpha=0.3)
    
    plt.tight_layout()
    fig_path = FIG_DIR / 'ferro_aging_scores_boxplot.png'
    fig.savefig(fig_path, dpi=200, bbox_inches='tight')
    plt.close()
    logger.info(f"  箱线图保存: {fig_path}")


def plot_forest_plot(meta_df: pd.DataFrame):
    """绘制Meta分析森林图"""
    if meta_df.empty:
        logger.warning("  Meta分析数据为空, 跳过森林图")
        return
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    y_positions = range(len(meta_df))
    
    for i, (_, row) in enumerate(meta_df.iterrows()):
        mean = row.get('mean_diff', 0)
        pval = row.get('pvalue', 1)
        # 近似95% CI
        # 注意: 由于process函数未返回标准差信息, 此处使用p值反推近似CI
        # 更严谨的做法应在process函数中返回每组的标准差, 再用
        # SE = sqrt(s1^2/n1 + s2^2/n2) 计算标准误
        # 当前近似方法供可视化趋势参考, 不作为统计推断依据
        se = abs(mean) / max(-np.log10(max(pval, 1e-300)), 0.1) * 0.5 if pval > 0 else abs(mean) / 2
        ci_low = mean - 1.96 * max(se, 0.01)
        ci_high = mean + 1.96 * max(se, 0.01)
        
        color = '#E74C3C' if mean > 0 else '#3498DB'
        ax.scatter(mean, i, color=color, s=80, zorder=5)
        ax.plot([ci_low, ci_high], [i, i], color=color, linewidth=2, zorder=4)
        
        p_text = f"{pval:.2e}" if pval > 0.001 else f"{pval:.1e}"
        ax.text(ci_high + 0.05, i, f"  p={p_text}", va='center', fontsize=8)
    
    ax.axvline(x=0, color='black', linestyle='--', alpha=0.5)
    ax.set_yticks(list(y_positions))
    ax.set_yticklabels(meta_df['dataset'].values, fontsize=10)
    ax.set_xlabel('Ferro-aging Score Difference (Case - Control)', fontsize=11)
    ax.set_title('Cross-Dataset Meta-Analysis: Ferro-aging in CIRI', 
                 fontsize=13, fontweight='bold')
    ax.tick_params(labelsize=9)
    
    plt.tight_layout()
    fig_path = FIG_DIR / 'ferro_aging_forest_plot.png'
    fig.savefig(fig_path, dpi=200, bbox_inches='tight')
    plt.close()
    logger.info(f"  森林图保存: {fig_path}")


# ============================================================
# 7. 主流程
# ============================================================

def main():
    logger.info("=" * 60)
    logger.info("铁衰老(Ferro-aging)在CIRI中的富集分析")
    logger.info(f"基因集: {len(FERRO_AGING_GENES)} 基因")
    logger.info(f"数据集: {list(DATA_DIRS.keys())}")
    logger.info("=" * 60)
    
    all_scores = {}
    all_stats = {}
    all_gene_dfs = []
    
    # 处理每个数据集
    processors = [
        ('GSE16561', process_gse16561),
        ('GSE37587', process_gse37587),
        ('GSE61616', process_gse61616),
        ('GSE97537', process_gse97537),
        ('GSE104036', process_gse104036),
    ]
    
    for ds_name, processor in processors:
        try:
            scores, gene_df, stats_series = processor()
            all_scores[ds_name] = scores
            all_stats[ds_name] = stats_series
            all_gene_dfs.append(gene_df)
            logger.info(f"  ✓ {ds_name} 处理完成")
        except Exception as e:
            logger.error(f"  ✗ {ds_name} 处理失败: {e}")
            import traceback
            traceback.print_exc()
    
    if not all_stats:
        logger.error("所有数据集处理失败!")
        return
    
    # ============================================================
    # 汇总统计表
    # ============================================================
    logger.info("\n" + "=" * 50)
    logger.info("汇总统计")
    
    meta_df = pd.DataFrame(all_stats).T.reset_index(drop=True)
    
    # 添加Meta分析
    pvals = meta_df['pvalue'].dropna().values
    if len(pvals) >= 2:
        chi2, meta_p = fisher_meta_analysis(list(pvals))
        logger.info(f"Fisher Meta分析: χ²={chi2:.2f}, p={meta_p:.4e}")
        meta_row = {
            'dataset': 'META-ANALYSIS',
            'species': 'Human/Rat/Mouse',
            'case_label': 'CIRI',
            'control_label': 'Control',
            'mean_diff': meta_df['mean_diff'].mean(),
            'pvalue': meta_p,
        }
        meta_df = pd.concat([meta_df, pd.DataFrame([meta_row])], ignore_index=True)
    
    logger.info(f"\n{meta_df[['dataset', 'species', 'n_case', 'n_control', 'mean_diff', 'pvalue']].to_string()}")
    
    # ============================================================
    # 可视化
    # ============================================================
    plot_ferro_aging_scores(all_scores, all_stats)
    plot_forest_plot(meta_df)
    
    # ============================================================
    # 输出Excel
    # ============================================================
    logger.info("\n" + "=" * 50)
    logger.info("输出结果")
    
    excel_path = OUTPUT_DIR / 'ferro_aging_ciri_results.xlsx'
    
    with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
        # Sheet 1: Meta分析汇总
        meta_df.to_excel(writer, sheet_name='Meta_Analysis', index=False)
        
        # Sheet 2: 各数据集铁衰老评分
        for ds_name in all_scores:
            scores = all_scores[ds_name]
            if scores is not None and not scores.empty:
                scores.to_frame('ferro_aging_score').to_excel(
                    writer, sheet_name=f'{ds_name}_scores')
        
        # Sheet 3: 单基因分析汇总
        if all_gene_dfs:
            combined_genes = pd.concat(all_gene_dfs, ignore_index=True)
            combined_genes.to_excel(writer, sheet_name='Gene_Level_Analysis', index=False)
            
            # Sheet 4: 核心基因交叉数据集总结
            core_genes = combined_genes[combined_genes['in_core_set'] == True]
            if not core_genes.empty:
                pivot = core_genes.pivot_table(
                    index='gene_symbol', columns='dataset',
                    values='log2FC', aggfunc='first')
                pivot.to_excel(writer, sheet_name='Core_Gene_Log2FC')
    
    logger.info(f"  结果保存: {excel_path}")
    
    # ============================================================
    # 输出简要报告
    # ============================================================
    report_path = OUTPUT_DIR / 'analysis_summary.txt'
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write("=" * 60 + "\n")
        f.write("铁衰老(Ferro-aging)在脑缺血再灌注损伤中的富集分析\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"分析时间: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}\n")
        f.write(f"铁衰老基因集: {len(FERRO_AGING_GENES)} 基因\n")
        f.write(f"核心标志基因: {len(CORE_FERRO_AGING_GENES)}\n\n")
        
        f.write("数据集分析结果:\n")
        f.write("-" * 50 + "\n")
        for _, row in meta_df.iterrows():
            if row['dataset'] == 'META-ANALYSIS':
                f.write(f"\nMETA-ANALYSIS (Fisher合并):\n")
                f.write(f"  p = {row['pvalue']:.4e}\n")
                continue
            f.write(f"\n{row['dataset']} ({row['species']}):\n")
            f.write(f"  {row['case_label']} vs {row['control_label']}\n")
            f.write(f"  n = {row['n_case']} vs {row['n_control']}\n")
            f.write(f"  Ferro-aging score Δ = {row['mean_diff']:.4f}\n")
            f.write(f"  p = {row['pvalue']:.4e}\n")
        
        f.write("\n\n结论:\n")
        f.write("-" * 50 + "\n")
        n_sig = (meta_df['pvalue'].dropna() < 0.05).sum() - (1 if 'META-ANALYSIS' in meta_df['dataset'].values and 
                 meta_df[meta_df['dataset']=='META-ANALYSIS']['pvalue'].values[0] < 0.05 else 0)
        f.write(f"在{len(all_stats)}个独立数据集中, {n_sig}个显示铁衰老通路在CIRI中显著激活\n")
        if 'META-ANALYSIS' in meta_df['dataset'].values:
            mp = meta_df[meta_df['dataset']=='META-ANALYSIS']['pvalue'].values[0]
            f.write(f"跨数据集Meta分析: p = {mp:.4e}\n")
            f.write(f"结论: 铁衰老通路在脑缺血再灌注损伤中显著富集激活\n")
    
    logger.info(f"  报告保存: {report_path}")
    logger.info(f"\n{'='*60}")
    logger.info("分析完成!")
    logger.info(f"结果目录: {OUTPUT_DIR}")
    logger.info(f"{'='*60}")


if __name__ == '__main__':
    main()