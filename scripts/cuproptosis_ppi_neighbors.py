# -*- coding: utf-8 -*-
"""
============================================================================
模块6: PPI扩展网络的差异表达分析
============================================================================
目的: 分析铜死亡基因的PPI邻居基因在CIRI中的差异表达
输入: PPI网络文件 (STRING) + DEGs文件 + 铜死亡基因集
输出: 邻居基因差异表达统计 + 网络可视化 + 功能富集

参考文献:
- Szklarczyk D, et al. Nucleic Acids Res 2021 (PMID:33237321) - STRING数据库
- Tsvetkov P, et al. Science 2022 (PMID:35298263) - 铜死亡核心基因
- Yang S, et al. Sci Rep 2024 (PMID:38956251) - FDX1在CIRI中的作用

GitHub最佳实践参考:
- Changwuuu/Cuproptosis-pancancer (PPI网络分析)
============================================================================
"""

import os
import sys
import logging
import warnings
from typing import Dict, List, Optional, Tuple, Union, Set
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import fisher_exact, hypergeom

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import (
    BASE_DIR, RESULTS_DIR, DATA_DIR,
    CUPROPTOSIS_GENES, CUPROPTOSIS_RELATED,
    FIG_FORMAT, FIG_DPI
)
from scripts.utils import setup_logger, ensure_dir

warnings.filterwarnings('ignore')

STAGE_DIR = os.path.join(RESULTS_DIR, "cuproptosis_ppi")
ensure_dir(STAGE_DIR)

logger = setup_logger("cuproptosis_ppi", os.path.join(STAGE_DIR, "cuproptosis_ppi.log"))


# ============================================================
# 铜死亡核心基因
# ============================================================
CUPROPTOSIS_CORE = CUPROPTOSIS_GENES


# ============================================================
# PPI网络加载与验证
# ============================================================
def load_ppi_network(ppi_file: Optional[str] = None,
                     min_score: float = 0.4) -> Optional[pd.DataFrame]:
    """
    加载PPI网络数据
    
    Parameters:
    -----------
    ppi_file : str, optional
        PPI文件路径 (TSV格式)
    min_score : float
        最小置信度分数 (STRING: 0-1000)
    
    Returns:
    --------
    pd.DataFrame : PPI网络数据框
    """
    if ppi_file is None:
        # 尝试从stage5查找
        possible_paths = [
            os.path.join(RESULTS_DIR, "stage5_string_ppi", "string_ppi.tsv"),
            os.path.join(RESULTS_DIR, "stage5_string_ppi", "string_interactions.tsv"),
            os.path.join(RESULTS_DIR, "ppi_network", "ppi_edges.csv"),
        ]
        for path in possible_paths:
            if os.path.exists(path):
                ppi_file = path
                logger.info(f"找到PPI文件: {ppi_file}")
                break
    
    if ppi_file is None or not os.path.exists(ppi_file):
        logger.error("PPI网络文件不存在")
        return None
    
    # 读取PPI文件
    try:
        ppi = pd.read_csv(ppi_file, sep='\t')
    except Exception:
        try:
            ppi = pd.read_csv(ppi_file)
        except Exception as e:
            logger.error(f"读取PPI文件失败: {e}")
            return None
    
    # 检查必要列
    required_cols = ['preferredName_A', 'preferredName_B']
    alt_cols = {
        'preferredName_A': ['protein1', 'gene_a', 'Gene_A', 'node1', 'source'],
        'preferredName_B': ['protein2', 'gene_b', 'Gene_B', 'node2', 'target'],
        'score': ['score', 'combined_score', 'confidence', 'weight']
    }
    
    for req, alts in alt_cols.items():
        if req not in ppi.columns:
            for alt in alts:
                if alt in ppi.columns:
                    ppi[req] = ppi[alt]
                    logger.info(f"  列映射: {alt} -> {req}")
                    break
    
    if 'preferredName_A' not in ppi.columns or 'preferredName_B' not in ppi.columns:
        logger.error(f"PPI文件缺少必要列，可用列: {list(ppi.columns)}")
        return None
    
    # 标准化基因名
    ppi['preferredName_A'] = ppi['preferredName_A'].astype(str).str.upper()
    ppi['preferredName_B'] = ppi['preferredName_B'].astype(str).str.upper()
    
    # 过滤低置信度互作
    if 'score' in ppi.columns:
        ppi = ppi[ppi['score'] >= min_score]
        logger.info(f"  按score>={min_score}过滤后: {len(ppi)} 条互作")
    
    # 去除自环
    ppi = ppi[ppi['preferredName_A'] != ppi['preferredName_B']].copy()
    
    logger.info(f"加载PPI网络: {len(ppi)} 条互作")
    return ppi


# ============================================================
# 邻居基因提取
# ============================================================
def extract_ppi_neighbors(ppi: pd.DataFrame,
                          seed_genes: List[str],
                          order: int = 1) -> Dict[str, Set[str]]:
    """
    提取PPI邻居基因
    
    Parameters:
    -----------
    ppi : pd.DataFrame
        PPI网络数据
    seed_genes : List[str]
        种子基因列表 (铜死亡基因)
    order : int
        邻居阶数 (1=直接邻居, 2=二阶邻居)
    
    Returns:
    --------
    Dict[str, Set[str]] : 每个种子基因的邻居集合
    """
    seed_set = set(g.upper() for g in seed_genes)
    
    # 构建邻接表
    adjacency = defaultdict(set)
    for _, row in ppi.iterrows():
        a = row['preferredName_A']
        b = row['preferredName_B']
        adjacency[a].add(b)
        adjacency[b].add(a)
    
    neighbors = {}
    
    for gene in seed_set:
        if gene not in adjacency:
            logger.warning(f"  {gene}: 不在PPI网络中")
            continue
        
        # 一阶邻居
        first_order = adjacency[gene] - seed_set  # 排除其他种子基因
        
        if order == 1:
            neighbors[gene] = first_order
        else:
            # 二阶邻居
            second_order = set()
            for neighbor in first_order:
                if neighbor in adjacency:
                    second_order.update(adjacency[neighbor])
            second_order = second_order - seed_set - first_order - {gene}
            neighbors[gene] = first_order | second_order
        
        logger.info(f"  {gene}: {len(neighbors[gene])} 个{order}阶邻居")
    
    return neighbors


def get_all_neighbors(neighbors_dict: Dict[str, Set[str]]) -> Set[str]:
    """获取所有邻居基因的并集"""
    all_neighbors = set()
    for gene, neighs in neighbors_dict.items():
        all_neighbors.update(neighs)
    return all_neighbors


# ============================================================
# 差异表达分析
# ============================================================
def analyze_neighbor_degs(ppi: pd.DataFrame,
                          degs_file: str,
                          cupro_genes: Optional[List[str]] = None,
                          neighbor_order: int = 1) -> Dict:
    """
    分析铜死亡PPI邻居基因的差异表达情况
    
    Parameters:
    -----------
    ppi : pd.DataFrame
        PPI网络
    degs_file : str
        DEGs文件路径
    cupro_genes : List[str], optional
        铜死亡基因列表
    neighbor_order : int
        邻居阶数
    
    Returns:
    --------
    Dict : 分析结果字典
    """
    cupro_genes = cupro_genes or CUPROPTOSIS_CORE
    
    # 1. 提取邻居
    logger.info("[1/3] 提取PPI邻居基因...")
    neighbors_dict = extract_ppi_neighbors(ppi, cupro_genes, neighbor_order)
    all_neighbors = get_all_neighbors(neighbors_dict)
    
    logger.info(f"总邻居基因数: {len(all_neighbors)}")
    
    # 2. 加载DEGs
    logger.info("[2/3] 加载DEGs...")
    if not os.path.exists(degs_file):
        logger.error(f"DEGs文件不存在: {degs_file}")
        return {}
    
    degs = pd.read_csv(degs_file)
    
    # 标准化列名
    col_mapping = {
        'GeneSymbol': ['GeneSymbol', 'Gene.symbol', 'gene_symbol', 'Gene', 'SYMBOL', 'gene'],
        'logFC': ['logFC', 'log2FoldChange', 'log2fc'],
        'adj.P.Val': ['adj.P.Val', 'padj', 'fdr', 'adj_pval'],
        'P.Value': ['P.Value', 'pvalue', 'p_value']
    }
    
    for std, alts in col_mapping.items():
        if std not in degs.columns:
            for alt in alts:
                if alt in degs.columns:
                    degs[std] = degs[alt]
                    break
    
    if 'GeneSymbol' not in degs.columns:
        logger.error(f"DEGs文件缺少基因名列，可用列: {list(degs.columns)}")
        return {}
    
    degs['GeneSymbol'] = degs['GeneSymbol'].astype(str).str.upper()
    
    # 3. 分析邻居基因的差异表达
    logger.info("[3/3] 分析邻居基因差异表达...")
    
    # 邻居基因在DEGs中
    deg_neighbors = degs[degs['GeneSymbol'].isin(all_neighbors)]
    
    # 显著差异的邻居
    if 'adj.P.Val' in degs.columns:
        sig_neighbors = deg_neighbors[deg_neighbors['adj.P.Val'] < 0.05]
    elif 'P.Value' in degs.columns:
        sig_neighbors = deg_neighbors[deg_neighbors['P.Value'] < 0.05]
    else:
        sig_neighbors = deg_neighbors
    
    # 上调和下调
    if 'logFC' in degs.columns:
        up_neighbors = sig_neighbors[sig_neighbors['logFC'] > 0]
        down_neighbors = sig_neighbors[sig_neighbors['logFC'] < 0]
    else:
        up_neighbors = pd.DataFrame()
        down_neighbors = pd.DataFrame()
    
    logger.info(f"PPI邻居基因差异表达分析:")
    logger.info(f"  邻居基因总数: {len(all_neighbors)}")
    logger.info(f"  在DEGs中的邻居: {len(deg_neighbors)}")
    logger.info(f"  显著差异邻居: {len(sig_neighbors)}")
    if len(sig_neighbors) > 0:
        logger.info(f"    上调: {len(up_neighbors)}")
        logger.info(f"    下调: {len(down_neighbors)}")
    
    # 按铜死亡基因统计
    per_gene_stats = []
    for cupro_gene, neighs in neighbors_dict.items():
        neigh_degs = degs[degs['GeneSymbol'].isin(neighs)]
        
        if 'adj.P.Val' in degs.columns and 'logFC' in degs.columns:
            sig = neigh_degs[neigh_degs['adj.P.Val'] < 0.05]
            up = sig[sig['logFC'] > 0]
            down = sig[sig['logFC'] < 0]
        else:
            sig = neigh_degs
            up = pd.DataFrame()
            down = pd.DataFrame()
        
        per_gene_stats.append({
            'Cuproptosis_Gene': cupro_gene,
            'N_Neighbors': len(neighs),
            'N_in_DEGs': len(neigh_degs),
            'N_Significant': len(sig),
            'N_Up': len(up),
            'N_Down': len(down),
            'Top_Up': ', '.join(up.nlargest(3, 'logFC')['GeneSymbol'].tolist()) if len(up) > 0 else '',
            'Top_Down': ', '.join(down.nsmallest(3, 'logFC')['GeneSymbol'].tolist()) if len(down) > 0 else ''
        })
    
    per_gene_df = pd.DataFrame(per_gene_stats)
    
    results = {
        'total_neighbors': len(all_neighbors),
        'deg_neighbors': len(deg_neighbors),
        'sig_neighbors': len(sig_neighbors),
        'up_neighbors': len(up_neighbors),
        'down_neighbors': len(down_neighbors),
        'sig_neighbor_genes': sig_neighbors['GeneSymbol'].tolist() if 'GeneSymbol' in sig_neighbors.columns else [],
        'per_gene_stats': per_gene_df,
        'sig_neighbors_df': sig_neighbors,
        'neighbors_dict': neighbors_dict
    }
    
    return results


# ============================================================
# 邻居基因功能富集
# ============================================================
def enrich_neighbor_functions(sig_neighbor_genes: List[str],
                               background_genes: Optional[int] = None) -> Optional[pd.DataFrame]:
    """
    对显著差异的邻居基因进行功能富集 (简化版)
    
    Parameters:
    -----------
    sig_neighbor_genes : List[str]
        显著差异的邻居基因列表
    background_genes : int, optional
        背景基因数
    
    Returns:
    --------
    pd.DataFrame : 富集结果
    """
    if len(sig_neighbor_genes) < 5:
        logger.warning("显著差异邻居基因不足5个，跳过富集")
        return None
    
    try:
        import gseapy as gp
        
        enr = gp.enrichr(
            gene_list=sig_neighbor_genes,
            gene_sets=['GO_Biological_Process_2023', 'KEGG_2021_Human'],
            organism='human',
            outdir=None,
            no_plot=True,
            cutoff=0.05
        )
        
        if enr.results is not None and len(enr.results) > 0:
            results = enr.results.copy()
            results['-log10(adjP)'] = -np.log10(results['Adjusted P-value'].clip(lower=1e-300))
            logger.info(f"邻居基因富集: {len(results)} 条")
            return results
        
    except ImportError:
        logger.warning("gseapy未安装，跳过富集分析")
    except Exception as e:
        logger.error(f"富集分析失败: {e}")
    
    return results


def identify_hub_genes(ppi_df: pd.DataFrame,
                       degs: pd.DataFrame,
                       top_n: int = 20,
                       output_dir: str = None) -> pd.DataFrame:
    """
    v13 FIX: 鉴定PPI网络hub基因 (文献标准: PMID:39391037, 41194198)
    
    使用度中心性和差异表达显著性综合评分
    
    Parameters:
    -----------
    ppi_df : pd.DataFrame
        PPI互作数据
    degs : pd.DataFrame
        差异表达结果
    top_n : int
        返回前N个hub基因
    output_dir : str
        输出目录
    
    Returns:
    --------
    pd.DataFrame : hub基因列表
    """
    if output_dir is None:
        output_dir = STAGE_DIR
    
    logger.info("[额外] Hub基因鉴定...")
    
    # 1. 计算每个基因的度中心性
    col1 = 'preferredName_A' if 'preferredName_A' in ppi_df.columns else 'node1'
    col2 = 'preferredName_B' if 'preferredName_B' in ppi_df.columns else 'node2'
    
    gene_degree = ppi_df[col1].value_counts() + ppi_df[col2].value_counts()
    gene_degree = gene_degree[~gene_degree.index.duplicated(keep='first')]
    
    # 2. 合并差异表达信息
    hub_data = []
    for gene, degree in gene_degree.items():
        if gene in degs.index:
            row = degs.loc[gene]
            logFC = row.get('logFC', 0)
            adj_P_Val = row.get('adj.P.Val', 1)
            avg_logFC = abs(logFC)
            neg_logFDR = -np.log10(max(adj_P_Val, 1e-300))
            
            hub_data.append({
                'Gene': gene,
                'Degree': degree,
                'logFC': logFC,
                'avg_logFC': avg_logFC,
                'adj_P_Val': adj_P_Val,
                'neg_logFDR': neg_logFDR,
                'Direction': 'Up' if logFC > 0 else 'Down'
            })
    
    if not hub_data:
        logger.warning("无hub基因数据")
        return pd.DataFrame()
    
    hub_df = pd.DataFrame(hub_data)
    
    # 3. 计算综合评分: HubScore = 归一化度 × 0.4 + 归一化-logFDR × 0.3 + 归一化|logFC| × 0.3
    hub_df['Degree_norm'] = hub_df['Degree'] / hub_df['Degree'].max()
    hub_df['neg_logFDR_norm'] = hub_df['neg_logFDR'] / hub_df['neg_logFDR'].max()
    hub_df['avg_logFC_norm'] = hub_df['avg_logFC'] / hub_df['avg_logFC'].max()
    
    hub_df['HubScore'] = (
        0.4 * hub_df['Degree_norm'] +
        0.3 * hub_df['neg_logFDR_norm'] +
        0.3 * hub_df['avg_logFC_norm']
    )
    
    hub_df = hub_df.sort_values('HubScore', ascending=False).head(top_n)
    
    # 保存结果
    hub_df.to_csv(os.path.join(output_dir, "ppi_hub_genes.csv"), index=False)
    
    logger.info(f"  Top {top_n} Hub基因:")
    for i, (_, row) in enumerate(hub_df.iterrows()):
        logger.info(f"    {i+1:2d}. {row['Gene']:10s}: Degree={int(row['Degree']):4d}, "
                   f"logFC={row['logFC']:+.2f}, FDR={row['adj_P_Val']:.2e}, "
                   f"HubScore={row['HubScore']:.3f}")
    
    # 绘制hub基因图
    _plot_hub_genes(hub_df, output_dir)
    
    return hub_df


# ============================================================
# 可视化
# ============================================================
def plot_neighbor_stats(per_gene_df: pd.DataFrame, output_dir: str):
    """绘制每个铜死亡基因的邻居统计图"""
    if per_gene_df.empty:
        return
    
    fig, axes = plt.subplots(1, 2, figsize=(14, max(5, len(per_gene_df) * 0.4)))
    
    # 左图: 邻居基因数
    y_pos = range(len(per_gene_df))
    axes[0].barh(y_pos, per_gene_df['N_Neighbors'], color='#3498DB', alpha=0.8, label='Total')
    axes[0].barh(y_pos, per_gene_df['N_Significant'], color='#E74C3C', alpha=0.8, label='Significant DEGs')
    axes[0].set_yticks(y_pos)
    axes[0].set_yticklabels(per_gene_df['Cuproptosis_Gene'], fontsize=9)
    axes[0].set_xlabel('Number of Neighbors', fontsize=11)
    axes[0].set_title('PPI Neighbors per Cuproptosis Gene', fontsize=12, fontweight='bold')
    axes[0].legend()
    
    # 右图: 上调vs下调
    width = 0.35
    axes[1].barh([y - width/2 for y in y_pos], per_gene_df['N_Up'], 
                 width, color='#E74C3C', alpha=0.8, label='Up-regulated')
    axes[1].barh([y + width/2 for y in y_pos], per_gene_df['N_Down'],
                 width, color='#3498DB', alpha=0.8, label='Down-regulated')
    axes[1].set_yticks(y_pos)
    axes[1].set_yticklabels(per_gene_df['Cuproptosis_Gene'], fontsize=9)
    axes[1].set_xlabel('Number of DEGs', fontsize=11)
    axes[1].set_title('Differential Expression of Neighbors', fontsize=12, fontweight='bold')
    axes[1].legend()
    
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, f"ppi_neighbor_stats.{FIG_FORMAT}"),
                dpi=FIG_DPI, bbox_inches='tight')
    plt.close()
    logger.info(f"  邻居统计图已保存")


def plot_neighbor_volcano(degs: pd.DataFrame,
                           neighbor_genes: Set[str],
                           cupro_genes: List[str],
                           output_dir: str):
    """绘制邻居基因的火山图"""
    if 'logFC' not in degs.columns or 'adj.P.Val' not in degs.columns:
        logger.warning("DEGs缺少必要列，跳过火山图")
        return
    
    degs = degs.copy()
    degs['-log10(adjP)'] = -np.log10(degs['adj.P.Val'].clip(lower=1e-300))
    
    # 分类
    degs['category'] = 'Other'
    degs.loc[degs['GeneSymbol'].isin(cupro_genes), 'category'] = 'Cuproptosis'
    degs.loc[degs['GeneSymbol'].isin(neighbor_genes), 'category'] = 'Neighbor'
    
    sig_mask = degs['adj.P.Val'] < 0.05
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # 绘制背景
    other = degs[degs['category'] == 'Other']
    ax.scatter(other['logFC'], other['-log10(adjP)'], c='grey', alpha=0.3, s=5, label='Other')
    
    # 绘制邻居
    neighbors = degs[(degs['category'] == 'Neighbor') & sig_mask]
    ax.scatter(neighbors['logFC'], neighbors['-log10(adjP)'], 
               c='#3498DB', alpha=0.7, s=20, label=f'Neighbor DEGs ({len(neighbors)})')
    
    # 绘制铜死亡基因
    cupro = degs[degs['category'] == 'Cuproptosis']
    ax.scatter(cupro['logFC'], cupro['-log10(adjP)'],
               c='#E74C3C', alpha=0.9, s=50, edgecolors='black', linewidth=0.5,
               marker='D', label=f'Cuproptosis ({len(cupro)})', zorder=5)
    
    # 标注铜死亡基因
    for _, row in cupro.iterrows():
        ax.annotate(row['GeneSymbol'], (row['logFC'], row['-log10(adjP)']),
                   fontsize=8, xytext=(5, 5), textcoords='offset points')
    
    ax.axhline(y=-np.log10(0.05), color='grey', linestyle='--', linewidth=0.5)
    ax.axvline(x=0, color='grey', linestyle='-', linewidth=0.5)
    
    ax.set_xlabel('log2 Fold Change', fontsize=11)
    ax.set_ylabel('-log10(Adjusted P-value)', fontsize=11)
    ax.set_title('PPI Neighbor Genes in CIRI DEGs', fontsize=13, fontweight='bold')
    ax.legend(fontsize=9, loc='upper right')
    
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, f"ppi_neighbor_volcano.{FIG_FORMAT}"),
                dpi=FIG_DPI, bbox_inches='tight')
    plt.close()
    logger.info(f"  邻居火山图已保存")


def _plot_hub_genes(hub_df: pd.DataFrame, output_dir: str):
    """绘制Hub基因散点图 (Degree vs -logFDR, 大小=HubScore)"""
    if hub_df.empty:
        return
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    colors = ['#E74C3C' if d == 'Up' else '#3498DB' for d in hub_df['Direction']]
    
    scatter = ax.scatter(hub_df['Degree'], hub_df['neg_logFDR'], 
                        s=hub_df['HubScore'] * 500, c=colors, alpha=0.7, edgecolors='white', linewidth=0.5)
    
    for i, (_, row) in enumerate(hub_df.iterrows()):
        ax.annotate(row['Gene'], (row['Degree'], row['neg_logFDR']),
                    fontsize=8, ha='center', va='bottom',
                    fontweight='bold' if row['HubScore'] > 0.5 else 'normal')
    
    ax.set_xlabel('Degree Centrality', fontsize=11)
    ax.set_ylabel('-log10(FDR)', fontsize=11)
    ax.set_title('PPI Hub Gene Identification\n(Size = HubScore, Color = Direction)',
                 fontsize=12, fontweight='bold')
    
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#E74C3C', alpha=0.7, label='Up-regulated'),
        Patch(facecolor='#3498DB', alpha=0.7, label='Down-regulated')
    ]
    ax.legend(handles=legend_elements, fontsize=9)
    
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, f"ppi_hub_genes.{FIG_FORMAT}"),
                dpi=FIG_DPI, bbox_inches='tight')
    plt.close()
    logger.info(f"  Hub基因散点图已保存")


# ============================================================
# 主函数
# ============================================================
def main():
    """主函数 - PPI邻居差异表达分析"""
    logger.info("=" * 60)
    logger.info("模块6: PPI扩展网络的差异表达分析")
    logger.info("=" * 60)
    
    # 1. 加载PPI网络
    logger.info("[1/4] 加载PPI网络...")
    ppi = load_ppi_network()
    if ppi is None:
        logger.error("无法加载PPI网络，退出")
        return None
    
    # 2. 加载DEGs
    logger.info("[2/4] 加载DEGs...")
    degs_file = os.path.join(RESULTS_DIR, "stage1_rma_degs", "limma_degs.csv")
    
    if not os.path.exists(degs_file):
        alt_paths = [
            os.path.join(RESULTS_DIR, "stage3_enrichment", "human_degs.csv"),
            os.path.join(RESULTS_DIR, "deg_analysis", "degs.csv"),
        ]
        for alt in alt_paths:
            if os.path.exists(alt):
                degs_file = alt
                logger.info(f"使用替代DEGs文件: {degs_file}")
                break
    
    if not os.path.exists(degs_file):
        logger.error(f"DEGs文件不存在")
        return None
    
    # 3. 分析邻居差异表达
    logger.info("[3/4] 分析PPI邻居差异表达...")
    results = analyze_neighbor_degs(ppi, degs_file, CUPROPTOSIS_CORE, neighbor_order=1)
    
    if not results:
        logger.error("邻居分析失败")
        return None
    
    # 保存结果
    if 'per_gene_stats' in results:
        results['per_gene_stats'].to_csv(
            os.path.join(STAGE_DIR, "ppi_neighbor_per_gene_stats.csv"), index=False)
        plot_neighbor_stats(results['per_gene_stats'], STAGE_DIR)
    
    if 'sig_neighbors_df' in results and not results['sig_neighbors_df'].empty:
        results['sig_neighbors_df'].to_csv(
            os.path.join(STAGE_DIR, "sig_neighbor_degs.csv"), index=False)
        
        # 火山图
        degs = pd.read_csv(degs_file)
        if 'GeneSymbol' not in degs.columns:
            for col in ['Gene.symbol', 'gene_symbol', 'Gene']:
                if col in degs.columns:
                    degs['GeneSymbol'] = degs[col]
                    break
        degs['GeneSymbol'] = degs['GeneSymbol'].astype(str).str.upper()
        
        all_neighbors = get_all_neighbors(results['neighbors_dict'])
        plot_neighbor_volcano(degs, all_neighbors, CUPROPTOSIS_CORE, STAGE_DIR)
    
    # 4. Hub基因鉴定与功能富集
    logger.info("[4/5] 鉴定PPI网络Hub基因...")
    degs_for_hub = pd.read_csv(degs_file)
    if 'GeneSymbol' not in degs_for_hub.columns:
        for col in ['Gene.symbol', 'gene_symbol', 'Gene']:
            if col in degs_for_hub.columns:
                degs_for_hub['GeneSymbol'] = degs_for_hub[col]
                break
    degs_for_hub['GeneSymbol'] = degs_for_hub['GeneSymbol'].astype(str).str.upper()
    degs_for_hub.set_index('GeneSymbol', inplace=True)
    
    hub_genes = identify_hub_genes(ppi, degs_for_hub, top_n=20, output_dir=STAGE_DIR)
    
    logger.info("[5/5] 邻居基因功能富集...")
    if results.get('sig_neighbor_genes'):
        enrichment = enrich_neighbor_functions(results['sig_neighbor_genes'])
        if enrichment is not None:
            enrichment.to_csv(os.path.join(STAGE_DIR, "neighbor_enrichment.csv"), index=False)
            
            # 绘制富集图
            top_enr = enrichment.head(15)
            fig, ax = plt.subplots(figsize=(10, max(4, len(top_enr) * 0.35)))
            
            colors = ['#E74C3C' if 'GO' in gs else '#3498DB' for gs in top_enr['Gene_set']]
            ax.barh(range(len(top_enr)), top_enr['-log10(adjP)'], color=colors, alpha=0.8)
            ax.set_yticks(range(len(top_enr)))
            ax.set_yticklabels([t[:60] for t in top_enr['Term']], fontsize=9)
            ax.set_xlabel('-log10(Adjusted P-value)', fontsize=11)
            ax.set_title('Neighbor Genes Functional Enrichment', fontsize=12, fontweight='bold')
            
            plt.tight_layout()
            fig.savefig(os.path.join(STAGE_DIR, f"neighbor_enrichment.{FIG_FORMAT}"),
                        dpi=FIG_DPI, bbox_inches='tight')
            plt.close()
    
    # 输出摘要
    logger.info("\n" + "=" * 60)
    logger.info("模块6完成!")
    logger.info(f"  总邻居基因: {results['total_neighbors']}")
    logger.info(f"  DEGs中邻居: {results['deg_neighbors']}")
    logger.info(f"  显著差异邻居: {results['sig_neighbors']}")
    logger.info(f"    上调: {results['up_neighbors']}")
    logger.info(f"    下调: {results['down_neighbors']}")
    logger.info("=" * 60)
    
    return results


if __name__ == "__main__":
    main()
