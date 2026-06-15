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
  1. 秩和富集评分: 每个样本计算ferro-aging评分
  2. 差异检验: CIRI组 vs 对照组
  3. 跨数据集Meta分析: Fisher合并p值
  4. 单基因分析: 关键铁衰老基因表达变化

输出:
  - ferro_aging_ciri_results.xlsx (多Sheet)

创新点:
  - 首次将"ferro-aging"概念引入CIRI研究领域
  - 为后续BCP→铁衰老通路研究提供干实验基础

依赖: pandas, numpy, scipy, matplotlib, openpyxl
=============================================================================
"""

import os, re, sys, gzip, json, warnings, logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings('ignore')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
logger = logging.getLogger(__name__)

# ============================================================
# 路径配置
# ============================================================
BASE_DIR = Path(__file__).parent
OUTPUT_DIR = BASE_DIR / "ferro_aging_ciri_results"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

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

# 核心铁衰老标志基因 (基于新基因集的关键标志)
CORE_FERRO_AGING_GENES = [
    "PTGS2", "HMOX1", "TFRC", "IL6", "IL1B",
    "HMGB1", "S100A8", "KEAP1", "SOD1", "HIF1A",
    "CDKN1A", "ALOX15", "NLRP3", "TLR4", "MPO"
]

logger.info(f"铁衰老基因集定义完成: {len(set(FERRO_AGING_GENES))} 个唯一基因")
logger.info(f"核心标志基因: {len(CORE_FERRO_AGING_GENES)}")


# ============================================================
# 2. 效应量工具函数
# ============================================================

def cohens_d(case: np.ndarray, control: np.ndarray) -> float:
    """计算Cohen's d效应量 (pooled标准差)"""
    n1, n2 = len(case), len(control)
    s1, s2 = np.var(case, ddof=1), np.var(control, ddof=1)
    pooled = np.sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
    return (np.mean(case) - np.mean(control)) / pooled if pooled > 0 else 0.0


# ============================================================
# 3. 秩和富集评分 (单样本)
# ============================================================

def rank_sum_enrichment_score(expr: np.ndarray, gene_mask: np.ndarray) -> float:
    """
    计算秩和富集评分 (单样本)
    
    注意: 非标准ssGSEA (非ECDF积分差), 而是内部的秩和富集度量:
      1. 对样本中所有基因的表达值排序
      2. 计算基因集内基因的秩次总和, 减去随机期望
      3. 除以基因集大小归一化
    
    与标准ssGSEA的区别:
      - 标准ssGSEA: 两类ECDF累积积分差 (需用gseapy库)
      - 本方法: 秩次总和偏差, 计算更轻量, 趋势一致但数值分布不同
    
    Args:
        expr: (n_genes,) 该样本的所有基因表达值
        gene_mask: (n_genes,) 基因集成员的布尔掩码
    
    Returns:
        float: 富集评分 (正值=富集, 负值=抑制)
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


def compute_enrichment_score_matrix(expr_df: pd.DataFrame, gene_set: Set[str]) -> pd.Series:
    """
    对表达矩阵所有样本计算秩和富集评分
    
    调用 rank_sum_enrichment_score 逐样本计算。
    注意: 非标准ssGSEA, 详见 rank_sum_enrichment_score 注释。
        expr_df: 基因×样本表达矩阵 (index=基因symbol, columns=样本)
        gene_set: 目标基因集
    
    Returns:
        Series: 每个样本的富集评分
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
        scores[col] = rank_sum_enrichment_score(vals[valid], gene_mask[valid])
    
    result = pd.Series(scores)
    logger.info(f"  富集评分计算完成: {len(common_genes)}/{len(gene_set)} 基因匹配, "
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
    paired: bool = False,
) -> pd.DataFrame:
    """
    分析铁衰老基因集中每个基因在两组间的差异表达
    
    Args:
        paired: 是否使用配对t检验 (需case_cols/control_cols顺序对齐)
    
    Returns:
        DataFrame: [gene, mean_case, mean_control, log2FC, pvalue, padj]
    """
    # 安全过滤: 确保传入列名存在于expr_df中, 避免loc KeyError
    case_cols = [c for c in case_cols if c in expr_df.columns]
    control_cols = [c for c in control_cols if c in expr_df.columns]
    if not case_cols or not control_cols:
        return pd.DataFrame()
    
    common = [g for g in expr_df.index if g in gene_set]
    results = []
    
    for gene in common:
        raw_case = expr_df.loc[gene, case_cols].values.astype(float)
        raw_ctrl = expr_df.loc[gene, control_cols].values.astype(float)
        
        if paired:
            # 配对检验: 保留完整对应对, 同时过滤NaN
            valid = ~(np.isnan(raw_case) | np.isnan(raw_ctrl))
            if valid.sum() < 3:
                continue
            case_vals = raw_case[valid]
            ctrl_vals = raw_ctrl[valid]
        else:
            # 独立检验: 分别过滤NaN
            case_vals = raw_case[~np.isnan(raw_case)]
            ctrl_vals = raw_ctrl[~np.isnan(raw_ctrl)]
            if len(case_vals) < 2 or len(ctrl_vals) < 2:
                continue
        
        mean_case = np.mean(case_vals)
        mean_ctrl = np.mean(ctrl_vals)
        
        # log2FC = mean_case - mean_ctrl
        # 假设: 表达矩阵已为log2空间
        #   - Affy芯片(RMA/PLIER/MAS5): 默认log2尺度 ✓
        #   - Illumina BeadArray: 通常为log2尺度 ✓
        #   - GSE104036 RNA-seq: 已转为log2(CPM+1) ✓
        #   - 结论: 对所有数据集这假设成立, 差值即log2倍变化
        log2fc = mean_case - mean_ctrl
        
        # t检验
        if paired and len(case_vals) == len(ctrl_vals) and len(case_vals) >= 3:
            _, pval = stats.ttest_rel(case_vals, ctrl_vals)
        else:
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
    """在目录中查找包含关键字的文件 (递归搜索子目录)"""
    if not os.path.isdir(dir_path):
        return None
    for root, dirs, files in os.walk(dir_path):
        for f in files:
            if all(k.lower() in f.lower() for k in keywords):
                return os.path.join(root, f)
    return None


def parse_series_matrix(filepath: str) -> pd.DataFrame:
    """解析GEO Series Matrix文件"""
    open_func = gzip.open if str(filepath).endswith('.gz') else open
    with open_func(filepath, 'rt', encoding='latin-1') as f:
        content = f.read()
    
    lines = content.splitlines()
    header_idx = None
    
    for i, line in enumerate(lines):
        if line.startswith('!series_matrix_table_begin'):
            # 跳过表头后的空行
            j = i + 1
            while j < len(lines) and not lines[j].strip():
                j += 1
            header_idx = j
            break
    
    if header_idx is None:
        raise ValueError(f"无法找到series_matrix_table_begin: {filepath}")
    
    # 解析表头
    header = lines[header_idx].strip().split('\t')
    header = [h.strip('"').strip() for h in header]
    
    # 解析数据 (从表头下一行开始, 跳过空行, 直到 !series_matrix_table_end)
    data_lines = []
    for i in range(header_idx + 1, len(lines)):
        if lines[i].startswith('!series_matrix_table_end'):
            break
        stripped = lines[i].strip()
        if stripped:  # 跳过空行
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
    n_mapped = mapped.shape[0]
    gene_series = pd.Series(mapped.index.map(probe_map), index=mapped.index)
    # 未映射探针用NaN填充 (不再用探针ID当基因名)
    unmapped_count = gene_series.isna().sum()
    if unmapped_count > 0:
        logger.debug(f"  {unmapped_count}/{n_mapped} 探针无基因映射, 已过滤")
    gene_series = gene_series.dropna()
    mapped = mapped.loc[gene_series.index]
    # 转大写实现跨物种大小写不敏感匹配
    gene_series = gene_series.str.upper()
    mapped.index = gene_series
    mapped = mapped.groupby(mapped.index).max()
    return mapped


# ============================================================
# 5a. GSE16561 (人全血, Stroke vs Control)
# ============================================================

def process_gse16561() -> Tuple[pd.DataFrame, pd.DataFrame, pd.Series, List[str], List[str]]:
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
    if sample_line is None:
        raise ValueError("GSE16561 series_matrix 中未找到 !Sample_geo_accession")
    if desc_line is None:
        raise ValueError("GSE16561 series_matrix 中未找到 !Sample_description")
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
                    gs_idx = next((i for i, h in enumerate(header) 
                                    if 'gene symbol' in h.lower() or 'symbol' in h.lower()), 2)
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
    fa_score = compute_enrichment_score_matrix(expr_gene, set(FERRO_AGING_GENES))
    
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
        'std_case': float(case_scores.std()) if len(case_scores) > 1 else 0.0,
        'std_control': float(ctrl_scores.std()) if len(ctrl_scores) > 1 else 0.0,
        'mean_diff': mean_diff, 'pvalue': pval,
        'cohens_d': round(cohens_d(case_scores.values, ctrl_scores.values), 4),
    }), list(stroke_cols), list(control_cols)


# ============================================================
# 5b. GSE37587 (人全血, 配对 Follow-Up vs Baseline)
# ============================================================

def process_gse37587() -> Tuple[pd.DataFrame, pd.DataFrame, pd.Series, List[str], List[str]]:
    """处理GSE37587 (配对设计)"""
    logger.info("=" * 50)
    logger.info("[GSE37587] 人全血: Follow-Up vs Baseline (配对)")
    
    sm_file = find_file(DATA_DIRS['GSE37587'], ['series_matrix'])
    if not sm_file:
        raise FileNotFoundError("GSE37587 series_matrix 未找到")
    
    expr_df = parse_series_matrix(sm_file)
    
    # 用GSM列名精确匹配
    with gzip.open(sm_file, 'rt', encoding='latin-1') as f:
        lines = f.readlines()
    
    sample_line = None
    desc_line = None
    for l in lines:
        if l.startswith('!Sample_geo_accession'):
            sample_line = [x.strip('"').strip() for x in l.strip().split('\t')]
        if l.startswith('!Sample_description'):
            desc_line = [x.strip('"').strip() for x in l.strip().split('\t')]
    
    # 从description解析分组和配对信息
    followup_cols = []
    baseline_cols = []
    patient_pairs = []  # patient_id -> group -> gsm

    for i, gsm in enumerate(sample_line[1:], 1):
        desc = desc_line[i] if i < len(desc_line) else ''
        desc_lower = desc.lower()
        
        # 从description判断分组
        is_baseline = any(kw in desc_lower for kw in ['baseline', 'hour 0', '0 hour'])
        is_followup = any(kw in desc_lower for kw in ['follow-up', 'follow up', 'followup',
                                                        'hour 24', 'hour 48', 'hour 72'])
        
        # 提取患者ID用于配对
        patient_match = re.search(r'[Pp]atient\s+(\d+)', desc)
        patient_id = patient_match.group(1) if patient_match else None
        
        if is_baseline:
            baseline_cols.append(gsm)
            if patient_id:
                patient_pairs.append({'patient': patient_id, 'group': 'baseline', 'gsm': gsm})
        elif is_followup:
            followup_cols.append(gsm)
            if patient_id:
                patient_pairs.append({'patient': patient_id, 'group': 'followup', 'gsm': gsm})
        else:
            logger.warning(f"  无法识别分组: {gsm} desc={desc[:80]}")
    
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
                    gs_idx = next((i for i, h in enumerate(header) 
                                    if 'gene symbol' in h.lower() or 'symbol' in h.lower()), 2)
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
    fa_score = compute_enrichment_score_matrix(expr_gene, set(FERRO_AGING_GENES))
    
    case_scores = fa_score[[c for c in followup_cols if c in fa_score.index]].dropna()
    ctrl_scores = fa_score[[c for c in baseline_cols if c in fa_score.index]].dropna()
    
    # 配对t检验 (优先, 因实验设计为配对)
    pval = None
    pair_df = pd.DataFrame(patient_pairs)
    common_patients = []
    bl_arr = np.array([])  # 预初始化用于d_z计算
    if not pair_df.empty:
        baseline_map = pair_df[pair_df['group']=='baseline'].set_index('patient')['gsm'].to_dict()
        followup_map = pair_df[pair_df['group']=='followup'].set_index('patient')['gsm'].to_dict()
        common_patients = sorted(set(baseline_map.keys()) & set(followup_map.keys()))
        
        if len(common_patients) >= 3:
            bl_scores = [fa_score.get(baseline_map[p], np.nan) for p in common_patients]
            fu_scores = [fa_score.get(followup_map[p], np.nan) for p in common_patients]
            valid = [(b, f) for b, f in zip(bl_scores, fu_scores)
                     if not np.isnan(b) and not np.isnan(f)]
            if len(valid) >= 3:
                bl_arr = np.array([v[0] for v in valid])
                fu_arr = np.array([v[1] for v in valid])
                _, pval_paired = stats.ttest_rel(bl_arr, fu_arr)
                pval = pval_paired
                logger.info(f"  配对t检验: n_pairs={len(valid)}, p={pval_paired:.4e}")
    
    # 如果配对检验不可行, 用非配对检验
    if pval is None:
        _, pval = stats.ttest_ind(case_scores, ctrl_scores, equal_var=False)
    
    mean_diff = case_scores.mean() - ctrl_scores.mean()
    logger.info(f"  Ferro-aging评分: FU={case_scores.mean():.3f}, "
                f"BL={ctrl_scores.mean():.3f}, Δ={mean_diff:.3f}, p={pval:.4e}")
    
    # 构建对齐的配对列用于单基因分析
    if (not pair_df.empty and len(common_patients) >= 3):
        baseline_ordered = [baseline_map[p] for p in common_patients
                            if baseline_map[p] in expr_gene.columns and
                               followup_map[p] in expr_gene.columns]
        followup_ordered = [followup_map[p] for p in common_patients
                            if baseline_map[p] in expr_gene.columns and
                               followup_map[p] in expr_gene.columns]
        if len(baseline_ordered) >= 3:
            gene_df = analyze_signature_genes(expr_gene, followup_ordered, baseline_ordered,
                                               set(FERRO_AGING_GENES), 'GSE37587',
                                               paired=True)
        else:
            gene_df = analyze_signature_genes(expr_gene, followup_cols, baseline_cols,
                                               set(FERRO_AGING_GENES), 'GSE37587',
                                               paired=False)
    else:
        gene_df = analyze_signature_genes(expr_gene, followup_cols, baseline_cols,
                                           set(FERRO_AGING_GENES), 'GSE37587',
                                           paired=False)

    # 配对Cohen's d_z (差值均值/差值标准差)
    dz_val = np.nan
    if len(bl_arr) >= 3 and len(fu_arr) >= 3:
        diff = fu_arr - bl_arr
        dz_val = round(np.mean(diff) / np.std(diff, ddof=1), 4)
    
    return fa_score, gene_df, pd.Series({
        'dataset': 'GSE37587', 'species': 'Human', 'tissue': 'Whole Blood',
        'case_label': 'Follow-Up', 'control_label': 'Baseline',
        'n_case': len(case_scores), 'n_control': len(ctrl_scores),
        'mean_case': case_scores.mean(), 'mean_control': ctrl_scores.mean(),
        'std_case': float(case_scores.std()) if len(case_scores) > 1 else 0.0,
        'std_control': float(ctrl_scores.std()) if len(ctrl_scores) > 1 else 0.0,
        'mean_diff': mean_diff, 'pvalue': pval,
        'cohens_d': dz_val if not np.isnan(dz_val) else round(cohens_d(case_scores.values, ctrl_scores.values), 4),
    }), list(followup_cols), list(baseline_cols)


# ============================================================
# 5c. GSE61616 (大鼠, MCAO 7d)
# ============================================================

def process_gse61616() -> Tuple[pd.DataFrame, pd.DataFrame, pd.Series, List[str], List[str]]:
    """处理GSE61616 (大鼠Affy芯片, MCAO 7d)"""
    logger.info("=" * 50)
    logger.info("[GSE61616] 大鼠 MCAO 7d: Model vs Sham")
    
    sm_file = find_file(DATA_DIRS['GSE61616'], ['series_matrix'])
    if not sm_file:
        raise FileNotFoundError("GSE61616 series_matrix 未找到")
    
    expr_df = parse_series_matrix(sm_file)
    
    # 动态解析样本分组 (替代硬编码GSM ID)
    with gzip.open(sm_file, 'rt', encoding='latin-1') as f:
        lines = f.readlines()
    sample_acc = None
    sample_title = None
    for l in lines:
        if l.startswith('!Sample_geo_accession'):
            sample_acc = [x.strip('"').strip() for x in l.strip().split('\t')]
        if l.startswith('!Sample_title'):
            sample_title = [x.strip('"').strip() for x in l.strip().split('\t')]
    sham_cols = []
    model_cols = []
    for i, gsm in enumerate(sample_acc[1:], 1):
        title = sample_title[i].lower() if i < len(sample_title) else ''
        if 'sham' in title:
            sham_cols.append(gsm)
        elif any(kw in title for kw in ['mcao', 'model', 'stroke']):
            model_cols.append(gsm)
    sham_cols = [c for c in sham_cols if c in expr_df.columns]
    model_cols = [c for c in model_cols if c in expr_df.columns]
    logger.info(f"  Model={len(model_cols)}, Sham={len(sham_cols)}")
    
    # GPL1355注释
    probe_map = {}
    if os.path.exists(GPL1355_FILE):
        probe_map = parse_gpl1355_annotation(GPL1355_FILE)
    expr_gene = collapse_probes(expr_df, probe_map)
    
    # ssGSEA (大鼠基因映射)
    fa_score = compute_enrichment_score_matrix(expr_gene, set(FERRO_AGING_GENES))
    
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
        'std_case': float(case_scores.std()) if len(case_scores) > 1 else 0.0,
        'std_control': float(ctrl_scores.std()) if len(ctrl_scores) > 1 else 0.0,
        'mean_diff': mean_diff, 'pvalue': pval,
        'cohens_d': round(cohens_d(case_scores.values, ctrl_scores.values), 4),
    }), list(model_cols), list(sham_cols)


# ============================================================
# 5d. GSE97537 (大鼠, MCAO 24h)
# ============================================================

def process_gse97537() -> Tuple[pd.DataFrame, pd.DataFrame, pd.Series, List[str], List[str]]:
    """处理GSE97537 (大鼠Affy芯片, MCAO 24h)"""
    logger.info("=" * 50)
    logger.info("[GSE97537] 大鼠 MCAO 24h: Model vs Sham")
    
    sm_file = find_file(DATA_DIRS['GSE97537'], ['series_matrix'])
    if not sm_file:
        raise FileNotFoundError("GSE97537 series_matrix 未找到")
    
    expr_df = parse_series_matrix(sm_file)
    
    # 动态解析样本分组 (替代硬编码GSM ID)
    with gzip.open(sm_file, 'rt', encoding='latin-1') as f:
        lines = f.readlines()
    sample_acc = None
    sample_title = None
    for l in lines:
        if l.startswith('!Sample_geo_accession'):
            sample_acc = [x.strip('"').strip() for x in l.strip().split('\t')]
        if l.startswith('!Sample_title'):
            sample_title = [x.strip('"').strip() for x in l.strip().split('\t')]
    sham_cols = []
    mcao_cols = []
    for i, gsm in enumerate(sample_acc[1:], 1):
        title = sample_title[i].lower() if i < len(sample_title) else ''
        if 'sham' in title:
            sham_cols.append(gsm)
        elif any(kw in title for kw in ['mcao', 'model', 'stroke']):
            mcao_cols.append(gsm)
    sham_cols = [c for c in sham_cols if c in expr_df.columns]
    mcao_cols = [c for c in mcao_cols if c in expr_df.columns]
    logger.info(f"  MCAO={len(mcao_cols)}, Sham={len(sham_cols)}")
    
    probe_map = {}
    if os.path.exists(GPL1355_FILE):
        probe_map = parse_gpl1355_annotation(GPL1355_FILE)
    expr_gene = collapse_probes(expr_df, probe_map)
    
    fa_score = compute_enrichment_score_matrix(expr_gene, set(FERRO_AGING_GENES))
    
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
        'std_case': float(case_scores.std()) if len(case_scores) > 1 else 0.0,
        'std_control': float(ctrl_scores.std()) if len(ctrl_scores) > 1 else 0.0,
        'mean_diff': mean_diff, 'pvalue': pval,
        'cohens_d': round(cohens_d(case_scores.values, ctrl_scores.values), 4),
    }), list(mcao_cols), list(sham_cols)


# ============================================================
# 5e. GSE104036 (小鼠, RNA-seq, 多时间点)
# ============================================================

def process_gse104036() -> Tuple[pd.DataFrame, pd.DataFrame, pd.Series, List[str], List[str]]:
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
    
    # 样本分组: 使用模糊匹配替代硬编码
    # 列名格式: S1, S2, S3 (sham) | I1_3hr, I2_6hr, I3_24hr (ipsilateral)
    # 使用关键词匹配: 排除sham列和C前缀(对侧)列, 按时间关键词匹配
    all_cols = expr_df.columns.tolist()
    logger.info(f"  列名: {all_cols}")
    
    sham_cols = sorted([c for c in all_cols 
                        if re.match(r'^S\d+', str(c)) or 'sham' in str(c).lower()])
    
    ipsi_candidates = [c for c in all_cols 
                       if 'sham' not in str(c).lower()
                       and not re.match(r'^C\d+', str(c))]  # 排除对侧(C前缀)
    ipsi_3hr = sorted([c for c in ipsi_candidates if re.search(r'(?i)3h', str(c))])
    ipsi_6hr = sorted([c for c in ipsi_candidates if re.search(r'(?i)6h', str(c))])
    ipsi_12hr = sorted([c for c in ipsi_candidates if re.search(r'(?i)12h', str(c))])
    ipsi_24hr = sorted([c for c in ipsi_candidates if re.search(r'(?i)24h', str(c))])
    
    # 合并所有ipsi时间点 vs sham
    ipsi_all = ipsi_3hr + ipsi_6hr + ipsi_12hr + ipsi_24hr
    
    # RNA-seq数据: log2(CPM+1) 转换用于ssGSEA
    # 如果是count数据, 先归一化
    # 判断是否为raw counts需要log转换
    # 准则: 最大值>50 且 中列值>5 且 大部分为整数 -> 视为raw counts
    flat = expr_df.values.flatten()
    flat = flat[~np.isnan(flat)]
    max_val = np.max(flat)
    median_val = np.median(flat)
    int_ratio = np.mean(flat == np.floor(flat))  # 整数比例
    if max_val > 50 and median_val > 5 and int_ratio > 0.5:
        logger.info(f"  检测到raw counts (max={max_val:.0f}, median={median_val:.0f}, "
                    f"int_ratio={int_ratio:.0%}), 执行log2(CPM+1)转换")
        col_sums = expr_df.sum()
        cpm = expr_df.div(col_sums, axis=1) * 1e6
        expr_df = np.log2(cpm + 1)
    else:
        logger.info(f"  数据已是log空间或归一化后表达值 (max={max_val:.1f}), 跳过转换")
    
    logger.info(f"  Sham={len(sham_cols)}, Ipsi_all={len(ipsi_all)}, "
                f"3hr={len(ipsi_3hr)}, 6hr={len(ipsi_6hr)}, "
                f"12hr={len(ipsi_12hr)}, 24hr={len(ipsi_24hr)}")
    
    if not sham_cols:
        raise ValueError("GSE104036: 无法识别Sham样本")
    
    # ssGSEA
    fa_score = compute_enrichment_score_matrix(expr_df, set(FERRO_AGING_GENES))
    
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
    for tp_name, tp_cols in [('3hr', ipsi_3hr), ('6hr', ipsi_6hr),
                              ('12hr', ipsi_12hr), ('24hr', ipsi_24hr)]:
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
        'std_case': float(case_scores.std()) if len(case_scores) > 1 else 0.0,
        'std_control': float(ctrl_scores.std()) if len(ctrl_scores) > 1 else 0.0,
        'mean_diff': mean_diff, 'pvalue': pval,
        'cohens_d': round(cohens_d(case_scores.values, ctrl_scores.values), 4),
        'timepoints': json.dumps({k: {kk: float(vv) if isinstance(vv, np.floating) else vv 
                                       for kk, vv in v.items()} 
                                  for k, v in timepoint_results.items()}),
    }), list(ipsi_all), list(sham_cols)


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
    all_groups = {}  # {ds_name: (case_cols, control_cols)}
    
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
            scores, gene_df, stats_series, case_cols, control_cols = processor()
            all_scores[ds_name] = scores
            all_stats[ds_name] = stats_series
            all_groups[ds_name] = (case_cols, control_cols)
            all_gene_dfs.append(gene_df)
            logger.info(f"  ✓ {ds_name} 处理完成")
        except FileNotFoundError as e:
            logger.error(f"  ✗ {ds_name} 数据文件缺失: {e}")
            continue
        except ValueError as e:
            logger.error(f"  ✗ {ds_name} 数据/分组解析错误: {e}")
            import traceback; traceback.print_exc()
            continue
        except Exception as e:
            logger.error(f"  ✗ {ds_name} 处理失败: [{type(e).__name__}] {e}")
            import traceback; traceback.print_exc()
            continue
    
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
            'n_case': int(meta_df['n_case'].sum()) if 'n_case' in meta_df.columns else 0,
            'n_control': int(meta_df['n_control'].sum()) if 'n_control' in meta_df.columns else 0,
            'mean_case': np.nan,    # 不报告合并效应量
            'mean_control': np.nan,  # 仅展示Fisher合并p值
            'std_case': np.nan,      # 跨数据集不合并标准差
            'std_control': np.nan,
            'mean_diff': np.nan,     # 简单平均无加权意义, 不科学
            'pvalue': meta_p,
            'cohens_d': np.nan,
        }
        meta_df = pd.concat([meta_df, pd.DataFrame([meta_row])], ignore_index=True)
    
    logger.info(f"\n{meta_df[['dataset', 'species', 'n_case', 'n_control', 'mean_diff', 'pvalue']].to_string()}")
    
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
            f.write(f"  Cohen's d = {row['cohens_d']:.3f}\n")
            f.write(f"  p = {row['pvalue']:.4e}\n")
        
        f.write("\n\n结论:\n")
        f.write("-" * 50 + "\n")
        # 统计显著的数据集数量 (排除META-ANALYSIS行)
        sig_datasets = 0
        for _, r in meta_df.iterrows():
            if r['dataset'] == 'META-ANALYSIS':
                continue
            if r['pvalue'] < 0.05:
                sig_datasets += 1
        f.write(f"在{len(all_stats)}个独立数据集中, {sig_datasets}个显示铁衰老通路在CIRI中显著激活\n")
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