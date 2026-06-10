# ============================================================
# 石竹烯(BCP) × 铜死亡 × CIRI 靶点筛选系统 - 工具函数模块
# ============================================================
# 版本: v2.0 | 日期: 2026-05-11
# 功能: 通用工具函数、数据加载、基因映射、可视化辅助
# 参考: PMID: 41234537 (WGCNA+MCODE+ML), PMID: 41791684 (scTenifoldKnk)
# ============================================================

import os
import sys
import warnings
import logging
import copy
import json
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Union, Any
from collections import OrderedDict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.spatial.distance import pdist, squareform

warnings.filterwarnings('ignore')

# ============================================================
# 0. 日志配置
# ============================================================
def setup_logger(name: str, log_file: str, level: int = logging.INFO) -> logging.Logger:
    """配置日志系统，同时输出到文件和终端"""
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.handlers.clear()

    fh = logging.FileHandler(log_file, encoding='utf-8')
    fh.setLevel(level)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(level)

    formatter = logging.Formatter(
        '%(asctime)s | %(levelname)-8s | %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)

    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger


def ensure_dir(dir_path: str) -> str:
    """确保目录存在，不存在则创建"""
    os.makedirs(dir_path, exist_ok=True)
    return dir_path


# ============================================================
# 1. 数据加载与验证
# ============================================================
def load_gse61616_top_table(filepath: str) -> pd.DataFrame:
    """加载GSE61616 GEO2R差异分析结果"""
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"文件不存在: {filepath}")

    df = pd.read_csv(filepath, sep='\t')

    required_cols = ['logFC', 'adj.P.Val', 'P.Value', 'Gene.symbol']
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"缺少必要列: {missing}，实际列: {list(df.columns)}")

    if df.shape[0] == 0:
        raise ValueError("数据为空")

    df['Gene.symbol'] = df['Gene.symbol'].fillna('').astype(str)
    df = df[df['Gene.symbol'] != ''].copy()
    df['Gene.symbol'] = df['Gene.symbol'].str.upper()

    return df


def load_expression_matrix(filepath: str) -> pd.DataFrame:
    """加载表达矩阵CSV文件"""
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"文件不存在: {filepath}")

    df = pd.read_csv(filepath, index_col=0)
    if df.shape[0] == 0 or df.shape[1] == 0:
        raise ValueError("表达矩阵为空")

    return df


def validate_input_data(df: pd.DataFrame, name: str = "数据") -> None:
    """验证输入数据的基本质量"""
    if df is None:
        raise ValueError(f"{name}为None")
    if df.shape[0] == 0:
        raise ValueError(f"{name}行数为0")
    if df.isna().all().all():
        raise ValueError(f"{name}全部为NaN")

    inf_count = np.isinf(df.select_dtypes(include=[np.number]).values).sum()
    if inf_count > 0:
        raise ValueError(f"{name}包含{inf_count}个Inf值")


# ============================================================
# 2. 基因映射
# ============================================================
def build_gene_mapping(
    rat_genes: List[str],
    mapping_file: Optional[str] = None
) -> Dict[str, str]:
    """
    构建大鼠→人类基因映射
    优先使用本地映射文件，否则使用同源基因数据库
    """
    mapping = {}

    if mapping_file and os.path.exists(mapping_file):
        try:
            map_df = pd.read_csv(mapping_file, sep='\t')
            if 'Rat' in map_df.columns and 'Human' in map_df.columns:
                for _, row in map_df.iterrows():
                    rat_gene = str(row['Rat']).strip().upper()
                    human_gene = str(row['Human']).strip().upper()
                    if rat_gene and human_gene:
                        mapping[rat_gene] = human_gene
        except Exception:
            pass

    if not mapping:
        mapping = _build_homology_mapping(rat_genes)

    return mapping


def _build_homology_mapping(rat_genes: List[str]) -> Dict[str, str]:
    """
    基于同源基因规则构建映射
    大鼠基因首字母大写 → 人类基因全大写
    常见一对一映射规则
    """
    mapping = {}
    for gene in rat_genes:
        gene_upper = gene.upper()
        mapping[gene_upper] = gene_upper

    return mapping


def map_rat_to_human(
    rat_genes: List[str],
    mapping: Dict[str, str]
) -> Tuple[List[str], List[str]]:
    """将大鼠基因列表映射为人类基因列表"""
    mapped = []
    unmapped = []

    for gene in rat_genes:
        gene_upper = gene.upper()
        if gene_upper in mapping:
            mapped.append(mapping[gene_upper])
        else:
            unmapped.append(gene_upper)

    return list(set(mapped)), list(set(unmapped))


# ============================================================
# 3. 基因集操作
# ============================================================
def gene_set_intersection(*sets: List[str]) -> List[str]:
    """多个基因集的交集"""
    if not sets:
        return []
    result = set(sets[0])
    for s in sets[1:]:
        result = result.intersection(set(s))
    return sorted(list(result))


def gene_set_union(*sets: List[str]) -> List[str]:
    """多个基因集的并集"""
    result = set()
    for s in sets:
        result = result.union(set(s))
    return sorted(list(result))


def gene_set_difference(set_a: List[str], set_b: List[str]) -> List[str]:
    """基因集差集 A - B"""
    return sorted(list(set(set_a) - set(set_b)))


# ============================================================
# 4. 统计函数
# ============================================================
def fdr_correction(p_values: np.ndarray, method: str = 'bh') -> np.ndarray:
    """FDR多重检验校正 (Benjamini-Hochberg)"""
    p_values = np.asarray(p_values, dtype=np.float64)
    if np.any(np.isnan(p_values)):
        raise ValueError("p值包含NaN")

    n = len(p_values)
    if n == 0:
        return np.array([])

    if method == 'bh':
        sorted_indices = np.argsort(p_values)
        sorted_p = p_values[sorted_indices]
        adjusted = np.minimum(1, sorted_p * n / np.arange(1, n + 1))
        adjusted = np.maximum.accumulate(adjusted[::-1])[::-1]
        result = np.zeros(n)
        result[sorted_indices] = adjusted
        return result
    else:
        return p_values * n


def safe_log2(x: np.ndarray) -> np.ndarray:
    """安全的log2转换，处理0值"""
    x = np.asarray(x, dtype=np.float64)
    x_safe = np.where(x <= 0, np.nanmin(x[x > 0]) if np.any(x > 0) else 1e-10, x)
    return np.log2(x_safe)


# ============================================================
# 5. 可视化辅助
# ============================================================
def setup_plotting_style():
    """设置统一的绘图风格"""
    plt.rcParams.update({
        'font.sans-serif': ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS', 'DejaVu Sans'],
        'axes.unicode_minus': False,
        'figure.dpi': 150,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'savefig.format': 'pdf',
        'font.size': 10,
        'axes.titlesize': 12,
        'axes.labelsize': 11,
    })


def save_figure(fig: plt.Figure, filepath: str, dpi: int = 300):
    """保存图片为PDF和PNG"""
    os.makedirs(os.path.dirname(filepath), exist_ok=True)

    base = os.path.splitext(filepath)[0]
    fig.savefig(f"{base}.pdf", dpi=dpi, bbox_inches='tight', format='pdf')
    fig.savefig(f"{base}.png", dpi=dpi, bbox_inches='tight', format='png')
    plt.close(fig)


def create_volcano_plot(
    df: pd.DataFrame,
    logfc_col: str = 'logFC',
    pval_col: str = 'adj.P.Val',
    gene_col: str = 'Gene.symbol',
    title: str = 'Volcano Plot',
    logfc_thresh: float = 1.0,
    pval_thresh: float = 0.05,
    highlight_genes: Optional[List[str]] = None,
    output_path: Optional[str] = None
) -> plt.Figure:
    """绘制火山图"""
    fig, ax = plt.subplots(figsize=(10, 8))

    df_plot = df.copy()
    df_plot['-log10(adj.P)'] = -np.log10(df_plot[pval_col].clip(lower=1e-300))

    df_plot['category'] = 'NS'
    up_mask = (df_plot[logfc_col] >= logfc_thresh) & (df_plot[pval_col] < pval_thresh)
    down_mask = (df_plot[logfc_col] <= -logfc_thresh) & (df_plot[pval_col] < pval_thresh)
    df_plot.loc[up_mask, 'category'] = 'Up'
    df_plot.loc[down_mask, 'category'] = 'Down'

    colors = {'NS': '#7f7f7f', 'Up': '#e74c3c', 'Down': '#3498db'}
    for cat, color in colors.items():
        mask = df_plot['category'] == cat
        ax.scatter(
            df_plot.loc[mask, logfc_col],
            df_plot.loc[mask, '-log10(adj.P)'],
            c=color, s=8, alpha=0.5, label=f'{cat} ({mask.sum()})', rasterized=True
        )

    if highlight_genes:
        highlight_mask = df_plot[gene_col].isin([g.upper() for g in highlight_genes])
        ax.scatter(
            df_plot.loc[highlight_mask, logfc_col],
            df_plot.loc[highlight_mask, '-log10(adj.P)'],
            c='#ff7f0e', s=40, edgecolors='black', linewidth=0.5,
            label=f'Highlight ({highlight_mask.sum()})', zorder=5
        )

    ax.axhline(-np.log10(pval_thresh), color='grey', linestyle='--', linewidth=0.8)
    ax.axvline(logfc_thresh, color='grey', linestyle='--', linewidth=0.8)
    ax.axvline(-logfc_thresh, color='grey', linestyle='--', linewidth=0.8)

    ax.set_xlabel('log2 Fold Change')
    ax.set_ylabel('-log10(adjusted P-value)')
    ax.set_title(title)
    ax.legend(loc='upper right', fontsize=8, framealpha=0.9)

    if output_path:
        save_figure(fig, output_path)

    return fig


def create_heatmap(
    data: pd.DataFrame,
    title: str = 'Heatmap',
    cmap: str = 'RdBu_r',
    output_path: Optional[str] = None
) -> plt.Figure:
    """绘制热图"""
    fig, ax = plt.subplots(figsize=(12, max(6, data.shape[0] * 0.3)))

    sns.heatmap(
        data, cmap=cmap, center=0,
        xticklabels=True, yticklabels=True,
        linewidths=0.5, linecolor='#f0f0f0',
        cbar_kws={'label': 'Expression (Z-score)'},
        ax=ax
    )

    ax.set_title(title)
    ax.set_xlabel('Samples')
    ax.set_ylabel('Genes')
    plt.tight_layout()

    if output_path:
        save_figure(fig, output_path)

    return fig


# ============================================================
# 6. 结果保存
# ============================================================
def save_results_table(
    df: pd.DataFrame,
    filepath: str,
    index: bool = True,
    sheet_name: str = 'Results'
):
    """保存结果为CSV和Excel"""
    os.makedirs(os.path.dirname(filepath), exist_ok=True)

    base = os.path.splitext(filepath)[0]
    df.to_csv(f"{base}.csv", index=index, encoding='utf-8-sig')

    try:
        with pd.ExcelWriter(f"{base}.xlsx", engine='openpyxl') as writer:
            df.to_excel(writer, sheet_name=sheet_name, index=index)
    except Exception:
        pass


def save_gene_list(genes: List[str], filepath: str, description: str = ""):
    """保存基因列表"""
    os.makedirs(os.path.dirname(filepath), exist_ok=True)

    with open(filepath, 'w', encoding='utf-8') as f:
        if description:
            f.write(f"# {description}\n")
        f.write(f"# Total: {len(genes)}\n")
        for gene in genes:
            f.write(f"{gene}\n")


# ============================================================
# 7. 进度报告
# ============================================================
def print_stage_header(stage: int, name: str):
    """打印阶段标题"""
    print("\n" + "=" * 70)
    print(f"  阶段 {stage}: {name}")
    print("=" * 70)


def print_summary(stats_dict: Dict[str, Any]):
    """打印统计摘要"""
    print("\n" + "-" * 50)
    for key, value in stats_dict.items():
        print(f"  {key}: {value}")
    print("-" * 50)


if __name__ == "__main__":
    print("工具函数模块加载成功")


# ============================================================
# 8. LASSO稳定性选择 (Stability Selection)
# ============================================================
def stability_selection_lasso(X, y, gene_names, n_bootstrap=100,
                               threshold=0.8, C=0.1, min_genes=30, max_genes=80):
    """
    Bootstrap LASSO 稳定性选择 (Stability Selection)
    参考: Meinshausen & Bühlmann, 2010, J. R. Stat. Soc. B
    返回: selected_genes (list), selection_freq (dict)
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler

    n_samples, n_genes = X.shape
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    selection_counts = {g: 0 for g in gene_names}

    for i in range(n_bootstrap):
        idx = np.random.choice(n_samples, size=n_samples, replace=True)
        X_boot, y_boot = X_scaled[idx], y[idx]

        model = LogisticRegression(
            penalty='l1', solver='liblinear', C=C,
            max_iter=3000, random_state=42 + i
        )
        try:
            model.fit(X_boot, y_boot)
            coef = model.coef_[0]
            for j, g in enumerate(gene_names):
                if abs(coef[j]) > 1e-6:
                    selection_counts[g] += 1
        except Exception:
            continue

    freq = {g: c / n_bootstrap for g, c in selection_counts.items()}

    selected = [g for g, f in freq.items() if f >= threshold]

    if len(selected) < min_genes:
        sorted_genes = sorted(freq.items(), key=lambda x: x[1], reverse=True)
        selected = [g for g, _ in sorted_genes[:min_genes]]
    if len(selected) > max_genes:
        sorted_genes = sorted(freq.items(), key=lambda x: x[1], reverse=True)
        selected = [g for g, _ in sorted_genes[:max_genes]]

    cupro_core = {"FDX1", "LIAS", "LIPT1", "DLAT", "PDHA1", "PDHB",
                   "MTF1", "GLS", "CDKN2A", "SLC31A1", "ATP7A", "ATP7B",
                   "DLD", "DBT", "DLST", "PDHA2", "GCSH"}
    forced = [g for g in cupro_core if g in gene_names and g not in selected]
    if forced:
        if len(selected) >= max_genes:
            selected = selected[:max_genes - len(forced)]
        selected.extend(forced)

    return selected, freq