# -*- coding: utf-8 -*-
"""
============================================================================
模块4: WGCNA模块与铜死亡基因的关联分析
============================================================================
目的: 寻找与铜死亡基因共表达的基因模块
输入: WGCNA模块分配文件 + 铜死亡基因集
输出: 模块富集统计 + 可视化 + 共表达网络摘要

参考文献:
- Langfelder P & Horvath S, BMC Bioinformatics 2008 (PMID:19114008) - WGCNA
- Tsvetkov P, et al. Science 2022 (PMID:35298263) - 铜死亡基因集
- Chen S, et al. Cancer Cell Int 2025 (PMID:41194198) - WGCNA+铜死亡

GitHub最佳实践参考:
- Changwuuu/Cuproptosis-pancancer (WGCNA模块分析)
============================================================================
"""

import os
import sys
import logging
import warnings
from typing import Dict, List, Optional, Tuple, Union
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

STAGE_DIR = os.path.join(RESULTS_DIR, "cuproptosis_wgcna")
ensure_dir(STAGE_DIR)

logger = setup_logger("cuproptosis_wgcna", os.path.join(STAGE_DIR, "cuproptosis_wgcna.log"))


# ============================================================
# 铜死亡基因集
# ============================================================
CUPROPTOSIS_CORE = CUPROPTOSIS_GENES

CUPROPTOSIS_EXTENDED = CUPROPTOSIS_CORE + [
    "ATOX1", "COX17", "CCS", "COX11", "SCO1", "SCO2",
    "STEAP1", "STEAP2", "STEAP3", "STEAP4",
    "CP", "COMMD1", "MT1A", "MT2A"
]


# ============================================================
# WGCNA结果加载与验证
# ============================================================
def load_wgcna_modules(modules_file: Optional[str] = None) -> Optional[pd.DataFrame]:
    """
    加载WGCNA模块分配结果
    
    Parameters:
    -----------
    modules_file : str, optional
        WGCNA模块分配文件路径
    
    Returns:
    --------
    pd.DataFrame : 模块分配数据框
    """
    if modules_file is None:
        # 尝试从stage4查找
        possible_paths = [
            os.path.join(RESULTS_DIR, "stage4_seed_wgcna", "wgcna_modules.csv"),
            os.path.join(RESULTS_DIR, "wgcna", "module_colors.csv"),
            os.path.join(RESULTS_DIR, "stage4_seed_wgcna", "module_assignments.csv"),
        ]
        for path in possible_paths:
            if os.path.exists(path):
                modules_file = path
                logger.info(f"找到WGCNA模块文件: {modules_file}")
                break
    
    if modules_file is None or not os.path.exists(modules_file):
        logger.error("WGCNA模块文件不存在")
        return None
    
    # v9 FIX: 优先使用已映射的wgcna_modules_mapped.csv
    mapped_file = os.path.join(RESULTS_DIR, "stage4_seed_wgcna", "wgcna_modules_mapped.csv")
    if os.path.exists(mapped_file):
        logger.info("  使用已映射的wgcna_modules_mapped.csv")
        modules = pd.read_csv(mapped_file)
    else:
        modules = pd.read_csv(modules_file)
        logger.info(f"  原始列: {list(modules.columns)}")
        
        # 处理ProbeID格式，映射到GeneSymbol
        if 'ProbeID' in modules.columns and 'Module' in modules.columns:
            logger.info("  检测到ProbeID格式，映射到GeneSymbol...")
            pg_file = os.path.join(RESULTS_DIR, "stage4_seed_wgcna", "wgcna_probe_gene_map.csv")
            if os.path.exists(pg_file):
                pg = pd.read_csv(pg_file)
                probe_to_gene = dict(zip(pg['ProbeID'], pg['GeneSymbol']))
                modules['Gene'] = modules['ProbeID'].map(probe_to_gene)
                modules = modules.dropna(subset=['Gene'])
                modules['Gene'] = modules['Gene'].astype(str).str.upper().str.strip()
                modules = modules[modules['Gene'] != ''].copy()
                logger.info(f"  映射后: {len(modules)} 基因")
            else:
                logger.error("  wgcna_probe_gene_map.csv不存在，无法映射")
                return None
        else:
            # 检查必要列
            if 'Module' not in modules.columns:
                col_map = {
                    'Module': ['module', 'ModuleColor', 'moduleColor', 'color'],
                    'Gene': ['gene', 'GeneSymbol', 'gene_symbol', 'ProbeID', 'probe_id']
                }
                for std, alts in col_map.items():
                    if std not in modules.columns:
                        for alt in alts:
                            if alt in modules.columns:
                                modules[std] = modules[alt]
                                logger.info(f"  列映射: {alt} -> {std}")
                                break
    
    if 'Module' not in modules.columns or 'Gene' not in modules.columns:
        logger.error(f"WGCNA模块文件缺少必要列，可用列: {list(modules.columns)}")
        return None
    
    modules['Gene'] = modules['Gene'].astype(str).str.upper().str.strip()
    modules = modules[modules['Gene'] != ''].copy()
    # 去重: 同一基因取第一个模块
    modules = modules.drop_duplicates(subset=['Gene'], keep='first')
    
    logger.info(f"加载WGCNA模块: {len(modules)} 基因, {modules['Module'].nunique()} 模块")
    return modules


def load_module_eigengenes(eigengenes_file: Optional[str] = None) -> Optional[pd.DataFrame]:
    """加载模块特征基因(ME)"""
    if eigengenes_file is None:
        possible_paths = [
            os.path.join(RESULTS_DIR, "stage4_seed_wgcna", "module_eigengenes.csv"),
            os.path.join(RESULTS_DIR, "wgcna", "MEs.csv"),
        ]
        for path in possible_paths:
            if os.path.exists(path):
                eigengenes_file = path
                break
    
    if eigengenes_file and os.path.exists(eigengenes_file):
        me = pd.read_csv(eigengenes_file, index_col=0)
        logger.info(f"加载模块特征基因: {me.shape}")
        return me
    
    return None


# ============================================================
# 模块富集分析
# ============================================================
def analyze_cuproptosis_modules(modules: pd.DataFrame,
                                 cupro_genes: Optional[List[str]] = None,
                                 background_genes: Optional[int] = None) -> pd.DataFrame:
    """
    分析WGCNA模块中铜死亡基因的富集情况
    
    Parameters:
    -----------
    modules : pd.DataFrame
        WGCNA模块分配数据框
    cupro_genes : List[str], optional
        铜死亡基因集
    background_genes : int, optional
        背景基因总数 (用于超几何检验)
    
    Returns:
    --------
    pd.DataFrame : 富集结果
    """
    cupro_genes = cupro_genes or CUPROPTOSIS_CORE
    cupro_set = set(g.upper() for g in cupro_genes)
    
    if background_genes is None:
        background_genes = modules['Gene'].nunique()
    
    n_cupro_total = len(cupro_set)
    
    results = []
    module_colors = modules['Module'].unique()
    
    for module in module_colors:
        if str(module).lower() == 'grey':
            continue  # 跳过未分配模块
        
        module_genes = set(modules[modules['Module'] == module]['Gene'].str.upper())
        n_module = len(module_genes)
        
        cupro_in_module = cupro_set & module_genes
        n_cupro_in_module = len(cupro_in_module)
        
        if n_cupro_in_module > 0:
            # 超几何检验
            # 总体: background_genes
            # 铜死亡基因: n_cupro_total
            # 模块基因: n_module
            # 模块中铜死亡基因: n_cupro_in_module
            
            pvalue = hypergeom.sf(n_cupro_in_module - 1, background_genes, 
                                  n_cupro_total, n_module)
            
            # 富集比
            enrichment_ratio = (n_cupro_in_module / n_module) / (n_cupro_total / background_genes) \
                               if n_module > 0 and n_cupro_total > 0 else 0
            
            results.append({
                'Module': module,
                'Module_Size': n_module,
                'N_Cuproptosis': n_cupro_in_module,
                'Cuproptosis_Genes': ','.join(sorted(cupro_in_module)),
                'Coverage': n_cupro_in_module / n_cupro_total,
                'Enrichment_Ratio': enrichment_ratio,
                'P_value': pvalue,
                'Significant': pvalue < 0.05
            })
    
    if not results:
        logger.warning("未找到富集铜死亡基因的模块")
        return pd.DataFrame()
    
    results_df = pd.DataFrame(results).sort_values('P_value')
    
    # v12 FIX: 添加BH-FDR多重检验校正 (检查清单2.4要求)
    from statsmodels.stats.multitest import multipletests
    _, fdr_values, _, _ = multipletests(results_df['P_value'].values, method='fdr_bh')
    results_df['FDR'] = fdr_values
    results_df['FDR_Significant'] = results_df['FDR'] < 0.05
    
    logger.info("铜死亡基因WGCNA模块富集分析:")
    for _, row in results_df.iterrows():
        sig = "***" if row['P_value'] < 0.001 else "**" if row['P_value'] < 0.01 else "*" if row['P_value'] < 0.05 else ""
        fdr_sig = "✓" if row.get('FDR_Significant', False) else "✗"
        logger.info(f"  {row['Module']:15s}: {row['N_Cuproptosis']:2d}/{row['Module_Size']:4d} 基因, "
                   f"ratio={row['Enrichment_Ratio']:.2f}, p={row['P_value']:.4e} {sig}, "
                   f"FDR={row.get('FDR', 1.0):.4e} {fdr_sig}")
    
    return results_df


def analyze_module_trait_correlation(modules: pd.DataFrame,
                                      eigengenes: pd.DataFrame,
                                      trait_data: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """
    分析模块特征基因与性状(如condition)的相关性
    
    Parameters:
    -----------
    modules : pd.DataFrame
        WGCNA模块分配
    eigengenes : pd.DataFrame
        模块特征基因
    trait_data : pd.DataFrame, optional
        性状数据 (样本×性状)
    
    Returns:
    --------
    pd.DataFrame : 相关性结果
    """
    if trait_data is None:
        # 尝试从样本注释推断
        annot_file = os.path.join(RESULTS_DIR, "stage1_rma_degs", "sample_annotations.csv")
        if os.path.exists(annot_file):
            annot = pd.read_csv(annot_file)
            if 'SampleID' in annot.columns and 'Group' in annot.columns:
                trait_data = pd.get_dummies(annot.set_index('SampleID')['Group'])
                logger.info(f"从样本注释构建性状数据: {trait_data.shape}")
    
    if trait_data is None:
        logger.warning("无性状数据，跳过模块-性状相关性分析")
        return pd.DataFrame()
    
    # 对齐样本
    common_samples = eigengenes.index.intersection(trait_data.index)
    if len(common_samples) == 0:
        logger.warning("模块特征基因与性状数据样本不匹配")
        return pd.DataFrame()
    
    me_aligned = eigengenes.loc[common_samples]
    trait_aligned = trait_data.loc[common_samples]
    
    from scipy.stats import pearsonr
    
    correlations = []
    for module_col in me_aligned.columns:
        for trait_col in trait_aligned.columns:
            me_values = me_aligned[module_col].values
            trait_values = trait_aligned[trait_col].values
            
            if np.std(me_values) > 0 and np.std(trait_values) > 0:
                corr, pval = pearsonr(me_values, trait_values)
                correlations.append({
                    'Module': module_col,
                    'Trait': trait_col,
                    'Pearson_R': corr,
                    'P_value': pval,
                    'Significant': pval < 0.05
                })
    
    corr_df = pd.DataFrame(correlations)
    if not corr_df.empty:
        # v13 FIX: 添加BH-FDR多重检验校正 (文献标准: PMID:39853851, 39513039)
        from statsmodels.stats.multitest import multipletests
        _, fdr_values, _, _ = multipletests(corr_df['P_value'].values, method='fdr_bh')
        corr_df['FDR'] = fdr_values
        corr_df['FDR_Significant'] = corr_df['FDR'] < 0.05
        
        corr_df = corr_df.sort_values('P_value')
        logger.info(f"模块-性状相关性: {len(corr_df)} 对")
        
        sig_count = corr_df['FDR_Significant'].sum()
        if sig_count > 0:
            logger.info(f"  FDR显著相关: {sig_count} 对")
            for _, row in corr_df[corr_df['FDR_Significant']].head(10).iterrows():
                logger.info(f"    {row['Module']} x {row['Trait']}: R={row['Pearson_R']:+.3f}, FDR={row['FDR']:.2e}")
    
    return corr_df


# ============================================================
# 可视化
# ============================================================
def plot_module_enrichment(enrichment_df: pd.DataFrame, output_dir: str):
    """绘制模块富集结果图"""
    if enrichment_df.empty:
        logger.warning("富集结果为空，跳过可视化")
        return
    
    fig, axes = plt.subplots(1, 2, figsize=(14, max(5, len(enrichment_df) * 0.5)))
    
    # 左图: 富集比
    colors = ['#E74C3C' if sig else '#95A5A6' for sig in enrichment_df['Significant']]
    axes[0].barh(range(len(enrichment_df)), enrichment_df['Enrichment_Ratio'], 
                 color=colors, alpha=0.8, edgecolor='black', linewidth=0.5)
    axes[0].set_yticks(range(len(enrichment_df)))
    axes[0].set_yticklabels(enrichment_df['Module'], fontsize=10)
    axes[0].set_xlabel('Enrichment Ratio', fontsize=11)
    axes[0].set_title('Cuproptosis Gene Enrichment by Module', fontsize=12, fontweight='bold')
    axes[0].axvline(x=1, color='grey', linestyle='--', linewidth=0.8)
    
    # 右图: -log10(p)
    neg_log_p = -np.log10(enrichment_df['P_value'].clip(lower=1e-300))
    axes[1].barh(range(len(enrichment_df)), neg_log_p,
                 color=colors, alpha=0.8, edgecolor='black', linewidth=0.5)
    axes[1].set_yticks(range(len(enrichment_df)))
    axes[1].set_yticklabels(enrichment_df['Module'], fontsize=10)
    axes[1].set_xlabel('-log10(P-value)', fontsize=11)
    axes[1].set_title('Statistical Significance', fontsize=12, fontweight='bold')
    axes[1].axvline(x=-np.log10(0.05), color='grey', linestyle='--', linewidth=0.8)
    
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, f"wgcna_module_enrichment.{FIG_FORMAT}"),
                dpi=FIG_DPI, bbox_inches='tight')
    plt.close()
    logger.info(f"  模块富集图已保存")


def plot_module_trait_heatmap(corr_df: pd.DataFrame, output_dir: str):
    """绘制模块-性状相关性热图"""
    if corr_df.empty:
        return
    
    pivot = corr_df.pivot(index='Module', columns='Trait', values='Pearson_R')
    
    fig, ax = plt.subplots(figsize=(max(6, pivot.shape[1] * 1.5), max(5, pivot.shape[0] * 0.5)))
    
    sns.heatmap(pivot, annot=True, fmt='.3f', cmap='RdBu_r', center=0,
                vmin=-1, vmax=1, linewidths=0.5, ax=ax)
    ax.set_title('Module-Trait Correlation', fontsize=12, fontweight='bold')
    ax.set_xlabel('Trait', fontsize=11)
    ax.set_ylabel('Module', fontsize=11)
    
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, f"wgcna_module_trait_corr.{FIG_FORMAT}"),
                dpi=FIG_DPI, bbox_inches='tight')
    plt.close()
    logger.info(f"  模块-性状热图已保存")


def plot_cuproptosis_module_network(modules: pd.DataFrame,
                                     enrichment_df: pd.DataFrame,
                                     output_dir: str):
    """
    绘制铜死亡基因在模块中的分布网络图 (简化版)
    """
    if enrichment_df.empty:
        return
    
    # 获取显著富集的模块
    sig_modules = enrichment_df[enrichment_df['Significant']]['Module'].tolist()
    
    if not sig_modules:
        logger.warning("无显著富集模块，跳过网络图")
        return
    
    # 统计每个模块中的铜死亡基因
    cupro_set = set(g.upper() for g in CUPROPTOSIS_CORE)
    
    module_cupro_counts = []
    for module in sig_modules:
        module_genes = set(modules[modules['Module'] == module]['Gene'].str.upper())
        cupro_in_module = cupro_set & module_genes
        for gene in cupro_in_module:
            module_cupro_counts.append({'Module': module, 'Gene': gene})
    
    if not module_cupro_counts:
        return
    
    network_df = pd.DataFrame(module_cupro_counts)
    
    # 创建简单的矩阵可视化
    pivot = pd.crosstab(network_df['Gene'], network_df['Module'])
    
    fig, ax = plt.subplots(figsize=(max(6, len(sig_modules) * 1.2), max(5, len(pivot) * 0.5)))
    
    sns.heatmap(pivot, annot=True, fmt='d', cmap='YlOrRd', 
                cbar_kws={'label': 'Presence'}, linewidths=0.5, ax=ax)
    ax.set_title('Cuproptosis Genes in Significant WGCNA Modules', 
                 fontsize=12, fontweight='bold')
    ax.set_xlabel('WGCNA Module', fontsize=11)
    ax.set_ylabel('Cuproptosis Gene', fontsize=11)
    
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, f"cuproptosis_module_network.{FIG_FORMAT}"),
                dpi=FIG_DPI, bbox_inches='tight')
    plt.close()
    logger.info(f"  网络图已保存")


# ============================================================
# 主函数
# ============================================================
def main():
    """主函数 - WGCNA模块与铜死亡基因关联分析"""
    logger.info("=" * 60)
    logger.info("模块4: WGCNA模块与铜死亡基因的关联分析")
    logger.info("=" * 60)
    
    # 1. 加载WGCNA模块
    logger.info("[1/4] 加载WGCNA模块...")
    modules = load_wgcna_modules()
    if modules is None:
        logger.error("无法加载WGCNA模块，退出")
        return None
    
    # 2. 模块富集分析
    logger.info("[2/4] 铜死亡基因模块富集分析...")
    enrichment_df = analyze_cuproptosis_modules(modules, CUPROPTOSIS_CORE)
    
    if not enrichment_df.empty:
        enrichment_file = os.path.join(STAGE_DIR, "cuproptosis_module_enrichment.csv")
        enrichment_df.to_csv(enrichment_file, index=False)
        logger.info(f"  富集结果已保存: {enrichment_file}")
        
        # 可视化
        plot_module_enrichment(enrichment_df, STAGE_DIR)
        plot_cuproptosis_module_network(modules, enrichment_df, STAGE_DIR)
    
    # 3. 扩展基因集分析
    logger.info("[3/4] 扩展铜死亡基因集分析...")
    enrichment_ext = analyze_cuproptosis_modules(modules, CUPROPTOSIS_EXTENDED)
    if not enrichment_ext.empty:
        enrichment_ext.to_csv(os.path.join(STAGE_DIR, "cuproptosis_module_enrichment_extended.csv"), 
                              index=False)
    
    # 4. 模块-性状相关性
    logger.info("[4/4] 模块-性状相关性分析...")
    eigengenes = load_module_eigengenes()
    if eigengenes is not None:
        trait_corr = analyze_module_trait_correlation(modules, eigengenes)
        if not trait_corr.empty:
            trait_corr.to_csv(os.path.join(STAGE_DIR, "module_trait_correlation.csv"), index=False)
            plot_module_trait_heatmap(trait_corr, STAGE_DIR)
    
    # 输出摘要
    logger.info("\n" + "=" * 60)
    logger.info("模块4完成!")
    if not enrichment_df.empty:
        sig_modules = enrichment_df[enrichment_df['Significant']]
        logger.info(f"  显著富集模块: {len(sig_modules)}")
        for _, row in sig_modules.iterrows():
            logger.info(f"    {row['Module']}: {row['N_Cuproptosis']} 铜死亡基因, p={row['P_value']:.4e}")
    else:
        logger.info("  无显著富集模块")
    logger.info("=" * 60)
    
    return enrichment_df


if __name__ == "__main__":
    main()
