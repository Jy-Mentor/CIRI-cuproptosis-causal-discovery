# -*- coding: utf-8 -*-
"""
============================================================================
模块2: 铜死亡基因集GSEA富集分析
============================================================================
目的: 检测铜死亡基因在CIRI差异表达排序中的非随机分布
输入: DEGs文件 (含logFC和P值) + 铜死亡基因集
输出: GSEA富集结果 + 富集图 + 统计摘要

参考文献:
- Subramanian A, et al. PNAS 2005 (PMID:16199517) - 原始GSEA算法
- Mootha VK, et al. Nat Genet 2003 (PMID:12808457) - 基因集富集概念
- Tsvetkov P, et al. Science 2022 (PMID:35298263) - 铜死亡基因集

GitHub最佳实践参考:
- Changwuuu/Cuproptosis-pancancer (GSEA富集分析)
- 16622911388/cuproptosis (prerank GSEA实现)
============================================================================
"""

import os
import sys
import logging
import warnings
from typing import Dict, List, Optional, Tuple, Union
from pathlib import Path
from collections import OrderedDict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import (
    BASE_DIR, RESULTS_DIR, DATA_DIR,
    CUPROPTOSIS_GENES, CUPROPTOSIS_RELATED,
    FIG_FORMAT, FIG_DPI
)
from scripts.utils import setup_logger, ensure_dir, validate_input_data

warnings.filterwarnings('ignore')

STAGE_DIR = os.path.join(RESULTS_DIR, "cuproptosis_gsea")
ensure_dir(STAGE_DIR)

logger = setup_logger("cuproptosis_gsea", os.path.join(STAGE_DIR, "cuproptosis_gsea.log"))


# ============================================================
# 铜死亡基因集定义
# 基于: Tsvetkov P, et al. Science 2022 (PMID:35298263)
#       Liu H, et al. 2022 (PMID:36119826) - 12核心基因
# ============================================================
CUPROPTOSIS_GENE_SETS = {
    # 12个核心基因 (Liu H 2022 权威定义)
    'Cuproptosis_Core_12': [
        "FDX1", "LIAS", "LIPT1", "DLD", "DLAT", "PDHA1", "PDHB",
        "MTF1", "GLS", "CDKN2A", "SLC31A1", "ATP7B"
    ],
    # 17个扩展核心基因
    'Cuproptosis_Core': [
        "FDX1", "LIAS", "LIPT1", "DLAT", "PDHA1", "PDHB",
        "MTF1", "GLS", "CDKN2A", "SLC31A1", "ATP7A", "ATP7B",
        "DLD", "DBT", "DLST", "PDHA2", "GCSH"
    ],
    # 扩展基因集 (含铜代谢相关)
    'Cuproptosis_Extended': [
        "FDX1", "LIAS", "LIPT1", "DLAT", "PDHA1", "PDHB",
        "MTF1", "GLS", "CDKN2A", "SLC31A1", "ATP7A", "ATP7B",
        "DLD", "DBT", "DLST", "PDHA2", "GCSH",
        "ATOX1", "COX17", "CCS", "COX11", "SCO1", "SCO2",
        "STEAP1", "STEAP2", "STEAP3", "STEAP4",
        "CP", "COMMD1", "MT1A", "MT2A", "HSPA4"
    ],
    # 脂酰化相关基因 (铜死亡直接靶点)
    'Cuproptosis_Lipoylation': [
        "FDX1", "LIAS", "LIPT1", "DLAT", "PDHA1", "PDHB",
        "DLD", "DBT", "DLST", "PDHA2", "GCSH"
    ],
    # 铜转运蛋白
    'Cuproptosis_Copper_Transport': [
        "SLC31A1", "ATP7A", "ATP7B", "ATOX1", "COX17", "CCS"
    ],
    # 促铜死亡基因 (7个)
    'Cuproptosis_Pro': [
        "FDX1", "LIAS", "LIPT1", "DLD", "DLAT", "PDHA1", "PDHB"
    ],
    # 抗铜死亡基因 (3个)
    'Cuproptosis_Anti': [
        "MTF1", "GLS", "CDKN2A"
    ]
}


# ============================================================
# 输入验证
# ============================================================
def validate_degs_file(degs_file: str) -> pd.DataFrame:
    """验证DEGs文件格式和内容"""
    if not os.path.exists(degs_file):
        raise FileNotFoundError(f"DEGs文件不存在: {degs_file}")
    
    degs = pd.read_csv(degs_file)
    
    if degs.shape[0] == 0:
        raise ValueError("DEGs文件为空")
    
    # 检查必要列
    required_cols = ['logFC', 'P.Value']
    missing = [c for c in required_cols if c not in degs.columns]
    if missing:
        # 尝试常见变体
        col_mapping = {
            'logFC': ['logFC', 'log2FoldChange', 'Log2FoldChange', 'log2fc'],
            'P.Value': ['P.Value', 'pvalue', 'PValue', 'p_value', 'adj.P.Val', 'padj'],
            'GeneSymbol': ['GeneSymbol', 'Gene.symbol', 'gene_symbol', 'Gene', 'SYMBOL']
        }
        
        for req, alts in col_mapping.items():
            if req in missing:
                for alt in alts:
                    if alt in degs.columns:
                        degs[req] = degs[alt]
                        missing.remove(req)
                        logger.info(f"  列映射: {alt} -> {req}")
                        break
        
        if missing:
            raise ValueError(f"DEGs文件缺少必要列: {missing}，可用列: {list(degs.columns)}")
    
    # 确保有基因名列
    if 'GeneSymbol' not in degs.columns:
        if 'Gene.symbol' in degs.columns:
            degs['GeneSymbol'] = degs['Gene.symbol']
        elif 'Gene' in degs.columns:
            degs['GeneSymbol'] = degs['Gene']
        else:
            raise ValueError("DEGs文件缺少基因名列")
    
    degs['GeneSymbol'] = degs['GeneSymbol'].fillna('').astype(str).str.upper()
    degs = degs[degs['GeneSymbol'] != ''].copy()
    
    # 检查数值范围
    if degs['logFC'].isna().all():
        raise ValueError("logFC列全部为NaN")
    
    logger.info(f"DEGs文件验证通过: {len(degs)} 行, 列: {list(degs.columns)}")
    return degs


# ============================================================
# GSEA分析核心函数
# ============================================================
def create_ranked_gene_list(degs: pd.DataFrame, 
                            ranking_method: str = 'logfc') -> pd.Series:
    """
    创建用于GSEA的排序基因列表
    
    Parameters:
    -----------
    degs : pd.DataFrame
        DEGs数据框，需包含GeneSymbol, logFC, P.Value列
    ranking_method : str
        'logfc' | 'signed_pvalue' | 'combined'
    
    Returns:
    --------
    pd.Series : 基因→排序值的Series，按值降序排列
    """
    degs = degs.copy()
    
    if ranking_method == 'logfc':
        # 仅使用logFC排序
        ranked = degs.set_index('GeneSymbol')['logFC']
        
    elif ranking_method == 'signed_pvalue':
        # 使用 -log10(p) * sign(logFC) 排序
        degs['signed_p'] = -np.log10(degs['P.Value'].clip(lower=1e-300)) * np.sign(degs['logFC'])
        ranked = degs.set_index('GeneSymbol')['signed_p']
        
    elif ranking_method == 'combined':
        # 组合排序: logFC加权p值
        degs['ranking'] = degs['logFC'] * (-np.log10(degs['P.Value'].clip(lower=1e-300)))
        ranked = degs.set_index('GeneSymbol')['ranking']
        
    else:
        raise ValueError(f"未知排序方法: {ranking_method}")
    
    # 去重: 保留绝对值最大的
    # v9 FIX: 处理重复基因，保留绝对logFC最大的
    ranked_df = ranked.reset_index()
    ranked_df.columns = ['GeneSymbol', 'score']
    ranked_df = ranked_df.loc[ranked_df['score'].abs().groupby(ranked_df['GeneSymbol']).idxmax()]
    ranked = ranked_df.set_index('GeneSymbol')['score'].sort_values(ascending=False)
    
    logger.info(f"排序基因列表: {len(ranked)} 基因, 方法: {ranking_method}")
    logger.info(f"  Top 5: {ranked.head(5).index.tolist()}")
    logger.info(f"  Bottom 5: {ranked.tail(5).index.tolist()}")
    
    return ranked


def run_cuproptosis_gsea(degs_file: str,
                          output_dir: Optional[str] = None,
                          gene_sets: Optional[Dict[str, List[str]]] = None,
                          ranking_method: str = 'logfc',
                          permutation_num: int = 1000) -> Optional[Dict]:
    """
    对铜死亡基因集进行GSEA分析
    
    Parameters:
    -----------
    degs_file : str
        DEGs CSV文件路径
    output_dir : str, optional
        输出目录
    gene_sets : Dict[str, List[str]], optional
        自定义基因集字典
    ranking_method : str
        基因排序方法
    permutation_num : int
        置换检验次数
    
    Returns:
    --------
    Dict : GSEA结果字典
    """
    if output_dir is None:
        output_dir = STAGE_DIR
    ensure_dir(output_dir)
    
    gene_sets = gene_sets or CUPROPTOSIS_GENE_SETS
    
    # 1. 验证并加载DEGs
    logger.info(f"[1/4] 加载DEGs: {degs_file}")
    degs = validate_degs_file(degs_file)
    
    # 2. 创建排序基因列表
    logger.info(f"[2/4] 创建排序基因列表 (方法: {ranking_method})...")
    ranked_genes = create_ranked_gene_list(degs, ranking_method)
    
    # 3. 运行GSEA
    logger.info(f"[3/4] 运行GSEA分析...")
    
    results = {}
    
    try:
        import gseapy as gp
        logger.info("  使用gseapy进行prerank GSEA")
        
        for set_name, gene_list in gene_sets.items():
            logger.info(f"  分析基因集: {set_name} ({len(gene_list)} 基因)")
            
            # 检查基因集覆盖度
            genes_in_ranking = [g for g in gene_list if g in ranked_genes.index]
            coverage = len(genes_in_ranking) / len(gene_list)
            logger.info(f"    基因集覆盖度: {len(genes_in_ranking)}/{len(gene_list)} ({coverage:.1%})")
            
            if len(genes_in_ranking) < 5:
                logger.warning(f"    基因集 {set_name} 覆盖度不足，跳过")
                continue
            
            try:
                pre_res = gp.prerank(
                    rnk=ranked_genes,
                    gene_sets={set_name: gene_list},
                    outdir=os.path.join(output_dir, f"gsea_{set_name.lower()}"),
                    min_size=3,
                    max_size=500,
                    permutation_num=permutation_num,
                    seed=42,
                    no_plot=False
                )
                
                # v10 FIX: 直接从prerank报告CSV读取结果，而非res2d
                report_file = os.path.join(output_dir, f"gsea_{set_name.lower()}", "gseapy.gene_set.prerank.report.csv")
                
                if os.path.exists(report_file):
                    res_df = pd.read_csv(report_file)
                    
                    if len(res_df) > 0:
                        row = res_df.iloc[0]
                        
                        # 从Term列获取基因集名，或直接使用set_name
                        nes_val = float(row.get('NES', row.get('nes', None)))
                        pval_val = float(row.get('NOM p-val', row.get('nom p-val', None)))
                        fdr_val = float(row.get('FDR q-val', row.get('fdr q-val', None)))
                        tag_pct = row.get('Tag %', None)

                        results[set_name] = {
                            'result': res_df,
                            'nes': nes_val,
                            'pval': pval_val,
                            'fdr': fdr_val,
                            'coverage': coverage,
                            'tag_pct': tag_pct,
                            'genes_in_ranking': genes_in_ranking
                        }

                        logger.info(f"    NES: {nes_val:.4f}, p={pval_val:.4e}, FDR={fdr_val:.4e}")
                    else:
                        logger.warning(f"    基因集 {set_name} GSEA报告为空")
                elif hasattr(pre_res, 'res2d') and pre_res.res2d is not None:
                    res_df = pre_res.res2d.copy()
                    nes_col = 'NES' if 'NES' in res_df.columns else 'es'
                    pval_col = 'NOM p-val' if 'NOM p-val' in res_df.columns else 'pval'
                    fdr_col = 'FDR q-val' if 'FDR q-val' in res_df.columns else 'fdr'

                    nes_val = res_df[nes_col].iloc[0] if nes_col in res_df.columns else None
                    pval_val = res_df[pval_col].iloc[0] if pval_col in res_df.columns else None
                    fdr_val = res_df[fdr_col].iloc[0] if fdr_col in res_df.columns else None

                    try:
                        nes_val = float(nes_val) if nes_val is not None else None
                    except (ValueError, TypeError):
                        nes_val = None
                    try:
                        pval_val = float(pval_val) if pval_val is not None else None
                    except (ValueError, TypeError):
                        pval_val = None
                    try:
                        fdr_val = float(fdr_val) if fdr_val is not None else None
                    except (ValueError, TypeError):
                        fdr_val = None

                    results[set_name] = {
                        'result': res_df,
                        'nes': nes_val,
                        'pval': pval_val,
                        'fdr': fdr_val,
                        'coverage': coverage,
                        'genes_in_ranking': genes_in_ranking
                    }

                    logger.info(f"    NES: {nes_val:.4f}" if nes_val is not None else "    NES: None")
                    logger.info(f"    p-value: {pval_val:.4e}" if pval_val is not None else "    p-value: None")
                    logger.info(f"    FDR: {fdr_val:.4e}" if fdr_val is not None else "    FDR: None")
                else:
                    logger.warning(f"    基因集 {set_name} GSEA无结果")
                    
            except Exception as e:
                logger.error(f"    基因集 {set_name} GSEA失败: {e}")
                continue
        
    except ImportError:
        logger.warning("gseapy未安装，使用简化版GSEA近似计算")
        results = _run_simplified_gsea(ranked_genes, gene_sets)
    
    # 4. 汇总结果
    logger.info("[4/4] 汇总结果...")
    
    if results:
        summary_data = []
        for set_name, res in results.items():
            pval = res.get('pval')
            fdr = res.get('fdr')
            nes = res.get('nes')
            tag_pct = res.get('tag_pct', res.get('coverage'))

            summary_data.append({
                'Gene_Set': set_name,
                'NES': nes,
                'P_value': pval,
                'FDR': fdr,
                'Coverage': res.get('coverage'),
                'Tag_%': tag_pct,
                'Significant': pval < 0.05 if pval is not None else False
            })
        
        summary_df = pd.DataFrame(summary_data)
        summary_file = os.path.join(output_dir, "cuproptosis_gsea_summary.csv")
        summary_df.to_csv(summary_file, index=False)
        logger.info(f"  GSEA汇总已保存: {summary_file}")
        
        # 可视化
        plot_gsea_summary(summary_df, output_dir)
        
        # 保存详细结果
        for set_name, res in results.items():
            if 'result' in res and res['result'] is not None:
                detail_file = os.path.join(output_dir, f"gsea_detail_{set_name.lower()}.csv")
                res['result'].to_csv(detail_file, index=False)
    else:
        logger.warning("GSEA无显著结果")
    
    return results


def _run_simplified_gsea(ranked_genes: pd.Series,
                         gene_sets: Dict[str, List[str]]) -> Dict:
    """
    简化版GSEA近似计算 (gseapy不可用时使用)
    基于Kolmogorov-Smirnov-like统计量
    """
    results = {}
    n = len(ranked_genes)
    
    for set_name, gene_list in gene_sets.items():
        genes_in_set = [g for g in gene_list if g in ranked_genes.index]
        if len(genes_in_set) < 5:
            continue
        
        # 获取基因集在排序列表中的位置
        positions = [ranked_genes.index.get_loc(g) for g in genes_in_set]
        positions = np.array(positions)
        
        # 计算ES (简化版)
        positions_sorted = np.sort(positions)
        n_hit = len(positions_sorted)
        n_miss = n - n_hit
        
        # 计算running sum
        running_sum = np.zeros(n + 1)
        hit_indices = set(positions_sorted)
        
        for i in range(n):
            if i in hit_indices:
                running_sum[i + 1] = running_sum[i] + 1.0 / n_hit
            else:
                running_sum[i + 1] = running_sum[i] - 1.0 / n_miss
        
        es = running_sum.max() if abs(running_sum.max()) > abs(running_sum.min()) else running_sum.min()
        
        # 近似NES (简化)
        nes = es / np.sqrt(n_hit / n) if n_hit > 0 else 0
        
        results[set_name] = {
            'nes': nes,
            'pval': None,  # 简化版不计算p值
            'fdr': None,
            'coverage': len(genes_in_set) / len(gene_list),
            'genes_in_ranking': genes_in_set
        }
        
        logger.info(f"  {set_name}: ES={es:.4f}, approx NES={nes:.4f}")
    
    return results


# ============================================================
# 可视化
# ============================================================
def plot_gsea_summary(summary_df: pd.DataFrame, output_dir: str):
    if summary_df.empty or summary_df['NES'].isna().all():
        logger.warning("GSEA摘要数据为空，跳过可视化")
        return
    
    summary_df = summary_df.copy()
    summary_df['NES'] = pd.to_numeric(summary_df['NES'], errors='coerce')
    summary_df = summary_df.dropna(subset=['NES'])
    
    fig, ax = plt.subplots(figsize=(10, max(4, len(summary_df) * 0.5)))
    
    colors = ['#E74C3C' if nes > 0 else '#3498DB' for nes in summary_df['NES']]
    
    y_pos = range(len(summary_df))
    bars = ax.barh(y_pos, summary_df['NES'], color=colors, alpha=0.8, edgecolor='black', linewidth=0.5)
    
    for i, (_, row) in enumerate(summary_df.iterrows()):
        sig_text = ""
        try:
            pval = float(row.get('P_value', None)) if row.get('P_value') is not None else None
            if pval is not None:
                if pval < 0.001:
                    sig_text = "***"
                elif pval < 0.01:
                    sig_text = "**"
                elif pval < 0.05:
                    sig_text = "*"
        except (ValueError, TypeError):
            pass
        
        ax.text(row['NES'] + 0.05 if row['NES'] > 0 else row['NES'] - 0.05, 
                i, sig_text, ha='left' if row['NES'] > 0 else 'right', 
                va='center', fontsize=10, fontweight='bold')
    
    ax.set_yticks(y_pos)
    ax.set_yticklabels(summary_df['Gene_Set'], fontsize=10)
    ax.set_xlabel('Normalized Enrichment Score (NES)', fontsize=11)
    ax.set_title('Cuproptosis Gene Sets GSEA Results', fontsize=13, fontweight='bold')
    ax.axvline(x=0, color='black', linestyle='-', linewidth=0.8)
    
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#E74C3C', alpha=0.8, label='Enriched in MCAO'),
        Patch(facecolor='#3498DB', alpha=0.8, label='Enriched in Sham')
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=9)
    
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, f"cuproptosis_gsea_summary.{FIG_FORMAT}"),
                dpi=FIG_DPI, bbox_inches='tight')
    plt.close()
    logger.info(f"  GSEA摘要图已保存")


# ============================================================
# 主函数
# ============================================================
def main():
    """主函数 - 从现有管线数据运行GSEA分析"""
    logger.info("=" * 60)
    logger.info("模块2: 铜死亡基因集GSEA富集分析")
    logger.info("=" * 60)
    
    # 尝试从现有管线加载DEGs
    degs_file = os.path.join(RESULTS_DIR, "stage1_rma_degs", "limma_degs.csv")
    
    if not os.path.exists(degs_file):
        # 尝试其他可能的位置
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
        logger.error(f"DEGs文件不存在: {degs_file}")
        logger.info("请确保阶段1 (stage1_rma_degs) 已成功运行")
        return None
    
    # 运行GSEA
    results = run_cuproptosis_gsea(degs_file)
    
    if results:
        logger.info("\n" + "=" * 60)
        logger.info("模块2完成!")
        logger.info(f"  分析基因集数: {len(results)}")
        for set_name, res in results.items():
            nes = res.get('nes', 'N/A')
            pval = res.get('pval', 'N/A')
            logger.info(f"  {set_name}: NES={nes}, p={pval}")
        logger.info("=" * 60)
    
    return results


if __name__ == "__main__":
    main()
