# -*- coding: utf-8 -*-
"""
阶段3: DEG功能富集分析 (GO/KEGG/GSEA)
输入: 阶段1的limma_degs.csv + 表达矩阵
输出: GO/KEGG富集结果 + GSEA结果 + 可视化
"""

import os
import sys
import json
import logging
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import (
    BASE_DIR, RESULTS_DIR, DATA_DIR,
    RAT_MOUSE_HUMAN_MAP, CUPROPTOSIS_GENES, CUPROPTOSIS_RELATED,
    BCP_TARGETS, FIG_FORMAT, FIG_DPI
)
from scripts.utils import setup_logger, ensure_dir

STAGE_DIR = os.path.join(RESULTS_DIR, "stage3_enrichment")
ensure_dir(STAGE_DIR)

logger = setup_logger("stage3", os.path.join(STAGE_DIR, "stage3.log"))


def load_rat_human_mapping(map_file):
    """加载大鼠-人类基因映射库"""
    logger.info(f"加载映射库: {map_file}")
    
    if not os.path.exists(map_file):
        raise FileNotFoundError(f"映射库不存在: {map_file}")
    
    mapping = {}
    with open(map_file, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                rat_symbol = parts[0].strip().upper()
                human_symbol = parts[3].strip().upper()
                if rat_symbol and human_symbol and human_symbol != '':
                    mapping[rat_symbol] = human_symbol
    
    logger.info(f"  加载 {len(mapping)} 个大鼠-人类映射对")
    return mapping


def load_degs_and_map(degs_file, mapping):
    """加载DEGs并映射到人类同源基因"""
    logger.info(f"加载DEGs: {degs_file}")
    
    degs = pd.read_csv(degs_file)
    logger.info(f"  总探针: {len(degs)}")
    
    degs['GeneSymbol'] = degs['GeneSymbol'].fillna('').astype(str).str.upper()
    degs = degs[degs['GeneSymbol'] != ''].copy()
    
    sig_degs = degs[(degs['adj.P.Val'] < 0.05) & (abs(degs['logFC']) > 1)].copy()
    logger.info(f"  显著DEGs (大鼠): {len(sig_degs)}")
    
    sig_degs['HumanSymbol'] = sig_degs['GeneSymbol'].map(mapping)
    
    unmapped = sig_degs[sig_degs['HumanSymbol'].isna()]['GeneSymbol'].unique()
    if len(unmapped) > 0:
        unmapped_file = os.path.join(STAGE_DIR, "unmapped_rat_genes.txt")
        with open(unmapped_file, 'w', encoding='utf-8') as f:
            for g in sorted(unmapped):
                f.write(f"{g}\n")
        logger.info(f"  未映射大鼠基因: {len(unmapped)} → {unmapped_file}")
    
    mapped_degs = sig_degs.dropna(subset=['HumanSymbol']).copy()
    
    gene_degs = mapped_degs.groupby('HumanSymbol').agg({
        'logFC': 'mean',
        'P.Value': 'min',
        'adj.P.Val': 'min'
    }).reset_index()
    gene_degs.columns = ['Gene', 'logFC', 'PValue', 'adjPVal']
    
    logger.info(f"  映射后人类DEGs: {len(gene_degs)}")
    logger.info(f"    上调: {(gene_degs['logFC'] > 0).sum()}")
    logger.info(f"    下调: {(gene_degs['logFC'] < 0).sum()}")
    
    return gene_degs, sig_degs


def run_go_enrichment(gene_list, gene_degs):
    """运行GO富集分析"""
    logger.info("运行GO富集分析...")
    
    try:
        import gseapy as gp
        
        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets=['GO_Biological_Process_2023',
                       'GO_Cellular_Component_2023',
                       'GO_Molecular_Function_2023'],
            organism='Human',
            outdir=None,
            no_plot=True,
            cutoff=0.05
        )
        
        if enr.results is None or len(enr.results) == 0:
            logger.warning("  GO富集无显著结果")
            return None
        
        results = enr.results.copy()
        results['-log10(adjP)'] = -np.log10(results['Adjusted P-value'].clip(lower=1e-300))
        
        go_file = os.path.join(STAGE_DIR, "go_enrichment.csv")
        results.to_csv(go_file, index=False)
        logger.info(f"  GO富集结果: {len(results)} 条 → {go_file}")
        
        for go_type in ['GO_Biological_Process_2023', 'GO_Cellular_Component_2023',
                         'GO_Molecular_Function_2023']:
            subset = results[results['Gene_set'] == go_type].head(10)
            if len(subset) > 0:
                plot_go_dotplot(subset, go_type, gene_degs)
        
        return results
        
    except ImportError:
        logger.warning("  gseapy未安装，跳过GO富集")
        return None
    except Exception as e:
        logger.error(f"  GO富集失败: {e}")
        return None


def plot_go_dotplot(go_df, go_type, gene_degs):
    """绘制GO富集点图"""
    fig, ax = plt.subplots(figsize=(10, max(4, len(go_df) * 0.35)))
    
    go_df = go_df.sort_values('-log10(adjP)', ascending=True)
    
    up_genes = set(gene_degs[gene_degs['logFC'] > 0]['Gene'])
    down_genes = set(gene_degs[gene_degs['logFC'] < 0]['Gene'])
    
    colors = []
    for _, row in go_df.iterrows():
        genes_in_term = set(row['Genes'].split(';'))
        up_count = len(genes_in_term & up_genes)
        down_count = len(genes_in_term & down_genes)
        if up_count > down_count:
            colors.append('#E74C3C')
        else:
            colors.append('#3498DB')
    
    terms = [t[:60] for t in go_df['Term']]
    ax.scatter(
        go_df['-log10(adjP)'], range(len(go_df)),
        s=go_df['Overlap'].apply(lambda x: int(x.split('/')[0]) if '/' in str(x) else 1) * 15,
        c=colors, alpha=0.8, edgecolors='black', linewidth=0.5
    )
    
    ax.set_yticks(range(len(go_df)))
    ax.set_yticklabels(terms, fontsize=9)
    ax.set_xlabel('-log10(Adjusted P-value)', fontsize=11)
    
    go_name = go_type.replace('GO_', '').replace('_2023', '').replace('_', ' ')
    ax.set_title(f'GO {go_name} Enrichment', fontsize=13, fontweight='bold')
    
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#E74C3C',
               markersize=10, label='Up-regulated'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#3498DB',
               markersize=10, label='Down-regulated')
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=8)
    
    plt.tight_layout()
    fig.savefig(os.path.join(STAGE_DIR, f"go_{go_name.lower().replace(' ', '_')}_dotplot.{FIG_FORMAT}"),
                dpi=FIG_DPI, bbox_inches='tight')
    plt.close()
    logger.info(f"    GO点图已保存: go_{go_name.lower().replace(' ', '_')}_dotplot.{FIG_FORMAT}")


def run_kegg_enrichment(gene_list, gene_degs):
    """运行KEGG通路富集分析"""
    logger.info("运行KEGG通路富集分析...")
    
    try:
        import gseapy as gp
        
        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets=['KEGG_2021_Human'],
            organism='Human',
            outdir=None,
            no_plot=True,
            cutoff=0.05
        )
        
        if enr.results is None or len(enr.results) == 0:
            logger.warning("  KEGG富集无显著结果")
            return None
        
        results = enr.results.copy()
        results['-log10(adjP)'] = -np.log10(results['Adjusted P-value'].clip(lower=1e-300))
        
        kegg_file = os.path.join(STAGE_DIR, "kegg_enrichment.csv")
        results.to_csv(kegg_file, index=False)
        logger.info(f"  KEGG富集结果: {len(results)} 条 → {kegg_file}")
        
        top_kegg = results.head(20)
        plot_kegg_bar(top_kegg, gene_degs)
        
        return results
        
    except ImportError:
        logger.warning("  gseapy未安装，跳过KEGG富集")
        return None
    except Exception as e:
        logger.error(f"  KEGG富集失败: {e}")
        return None


def plot_kegg_bar(kegg_df, gene_degs):
    """绘制KEGG富集条形图"""
    fig, ax = plt.subplots(figsize=(10, max(4, len(kegg_df) * 0.35)))
    
    kegg_df = kegg_df.sort_values('-log10(adjP)', ascending=True)
    
    up_genes = set(gene_degs[gene_degs['logFC'] > 0]['Gene'])
    down_genes = set(gene_degs[gene_degs['logFC'] < 0]['Gene'])
    
    colors = []
    for _, row in kegg_df.iterrows():
        genes_in_term = set(row['Genes'].split(';'))
        up_count = len(genes_in_term & up_genes)
        down_count = len(genes_in_term & down_genes)
        if up_count > down_count:
            colors.append('#E74C3C')
        else:
            colors.append('#3498DB')
    
    terms = [t[:60] for t in kegg_df['Term']]
    ax.barh(range(len(kegg_df)), kegg_df['-log10(adjP)'], color=colors, alpha=0.8)
    
    ax.set_yticks(range(len(kegg_df)))
    ax.set_yticklabels(terms, fontsize=9)
    ax.set_xlabel('-log10(Adjusted P-value)', fontsize=11)
    ax.set_title('KEGG Pathway Enrichment', fontsize=13, fontweight='bold')
    
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color='#E74C3C', lw=4, label='Up-regulated'),
        Line2D([0], [0], color='#3498DB', lw=4, label='Down-regulated')
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=8)
    
    plt.tight_layout()
    fig.savefig(os.path.join(STAGE_DIR, f"kegg_barplot.{FIG_FORMAT}"),
                dpi=FIG_DPI, bbox_inches='tight')
    plt.close()
    logger.info("    KEGG条形图已保存")


def run_gsea(expr_matrix_file, degs_file, mapping):
    """运行GSEA分析"""
    logger.info("运行GSEA分析...")
    
    try:
        import gseapy as gp
        
        expr = pd.read_csv(expr_matrix_file, index_col=0)
        logger.info(f"  表达矩阵: {expr.shape[0]} 探针 x {expr.shape[1]} 样本")
        
        degs = pd.read_csv(degs_file)
        degs['GeneSymbol'] = degs['GeneSymbol'].fillna('').astype(str).str.upper()
        degs = degs[degs['GeneSymbol'] != ''].copy()
        degs['HumanSymbol'] = degs['GeneSymbol'].map(mapping)
        degs = degs.dropna(subset=['HumanSymbol'])
        
        gene_degs = degs.groupby('HumanSymbol').agg({
            'logFC': 'mean',
            'P.Value': 'min'
        }).reset_index()
        gene_degs.columns = ['Gene', 'logFC', 'PValue']
        
        ranked_genes = gene_degs.set_index('Gene')['logFC'].sort_values(ascending=False)
        
        logger.info(f"  排序基因数: {len(ranked_genes)}")
        
        pre_res = gp.prerank(
            rnk=ranked_genes,
            gene_sets=['KEGG_2021_Human', 'GO_Biological_Process_2023'],
            outdir=os.path.join(STAGE_DIR, "gsea_output"),
            seed=123,
            permutation_num=1000,
            min_size=10,
            max_size=500,
            no_plot=True
        )
        
        if pre_res.results is not None and len(pre_res.results) > 0:
            gsea_results = pre_res.results.copy()
            gsea_file = os.path.join(STAGE_DIR, "gsea_results.csv")
            gsea_results.to_csv(gsea_file, index=False)
            
            sig_gsea = gsea_results[gsea_results['FDR q-val'] < 0.25]
            logger.info(f"  GSEA结果: {len(gsea_results)} 条, 显著(FDR<0.25): {len(sig_gsea)} 条")
            
            if len(sig_gsea) > 0:
                plot_gsea_summary(sig_gsea)
        else:
            logger.warning("  GSEA无结果")
        
        return pre_res
        
    except ImportError:
        logger.warning("  gseapy未安装，跳过GSEA")
        return None
    except Exception as e:
        logger.error(f"  GSEA失败: {e}")
        import traceback
        traceback.print_exc()
        return None


def plot_gsea_summary(gsea_df):
    """绘制GSEA结果摘要图"""
    gsea_df = gsea_df.copy()
    gsea_df['abs_NES'] = gsea_df['NES'].abs()
    gsea_df = gsea_df.sort_values('abs_NES', ascending=False).head(20)
    
    fig, ax = plt.subplots(figsize=(10, max(4, len(gsea_df) * 0.35)))
    
    colors = ['#E74C3C' if x > 0 else '#3498DB' for x in gsea_df['NES']]
    terms = [t[:60] for t in gsea_df['Term']]
    
    ax.barh(range(len(gsea_df)), gsea_df['NES'], color=colors, alpha=0.8)
    ax.set_yticks(range(len(gsea_df)))
    ax.set_yticklabels(terms, fontsize=9)
    ax.set_xlabel('Normalized Enrichment Score (NES)', fontsize=11)
    ax.set_title('GSEA - Top Enriched Pathways', fontsize=13, fontweight='bold')
    ax.axvline(x=0, color='black', linestyle='-', linewidth=0.5)
    
    plt.tight_layout()
    fig.savefig(os.path.join(STAGE_DIR, f"gsea_summary.{FIG_FORMAT}"),
                dpi=FIG_DPI, bbox_inches='tight')
    plt.close()
    logger.info("    GSEA摘要图已保存")


def analyze_cuproptosis_enrichment(gene_degs):
    """分析铜死亡相关基因在DEGs中的富集情况"""
    logger.info("分析铜死亡基因富集...")
    
    deg_genes = set(gene_degs['Gene'])
    
    cupro_human = set(CUPROPTOSIS_GENES)
    cupro_related_human = set(CUPROPTOSIS_RELATED)
    
    cupro_in_degs = cupro_human & deg_genes
    cupro_related_in_degs = cupro_related_human & deg_genes
    
    logger.info(f"  铜死亡核心基因在DEGs中: {len(cupro_in_degs)}/{len(cupro_human)}")
    if cupro_in_degs:
        for g in sorted(cupro_in_degs):
            row = gene_degs[gene_degs['Gene'] == g].iloc[0]
            logger.info(f"    {g}: logFC={row['logFC']:.3f}, adjP={row['adjPVal']:.2e}")
    
    logger.info(f"  铜死亡相关基因在DEGs中: {len(cupro_related_in_degs)}/{len(cupro_related_human)}")
    
    cupro_info = {
        'core_in_degs': sorted(cupro_in_degs),
        'related_in_degs': sorted(cupro_related_in_degs),
        'core_total': len(cupro_human),
        'related_total': len(cupro_related_human)
    }
    
    with open(os.path.join(STAGE_DIR, "cuproptosis_enrichment.json"), 'w', encoding='utf-8') as f:
        json.dump(cupro_info, f, indent=2, ensure_ascii=False)
    
    return cupro_info


def plot_volcano(gene_degs):
    """绘制火山图"""
    logger.info("绘制火山图...")
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    gene_degs = gene_degs.copy()
    gene_degs['-log10(adjP)'] = -np.log10(gene_degs['adjPVal'].clip(lower=1e-300))
    
    sig_up = gene_degs[(gene_degs['adjPVal'] < 0.05) & (gene_degs['logFC'] > 1)]
    sig_down = gene_degs[(gene_degs['adjPVal'] < 0.05) & (gene_degs['logFC'] < -1)]
    ns = gene_degs[(gene_degs['adjPVal'] >= 0.05) | (abs(gene_degs['logFC']) <= 1)]
    
    ax.scatter(ns['logFC'], ns['-log10(adjP)'], c='grey', alpha=0.3, s=5, label='NS')
    ax.scatter(sig_up['logFC'], sig_up['-log10(adjP)'], c='#E74C3C', alpha=0.5, s=10, label=f'Up ({len(sig_up)})')
    ax.scatter(sig_down['logFC'], sig_down['-log10(adjP)'], c='#3498DB', alpha=0.5, s=10, label=f'Down ({len(sig_down)})')
    
    cupro_human = set(CUPROPTOSIS_GENES)
    cupro_in_data = gene_degs[gene_degs['Gene'].isin(cupro_human)]
    if len(cupro_in_data) > 0:
        ax.scatter(cupro_in_data['logFC'], cupro_in_data['-log10(adjP)'],
                   c='#F39C12', s=30, edgecolors='black', linewidth=0.5,
                   marker='D', label=f'Cuproptosis ({len(cupro_in_data)})', zorder=5)
        for _, row in cupro_in_data.iterrows():
            ax.annotate(row['Gene'], (row['logFC'], row['-log10(adjP)']),
                       fontsize=7, xytext=(5, 5), textcoords='offset points')
    
    ax.axhline(y=-np.log10(0.05), color='grey', linestyle='--', linewidth=0.5)
    ax.axvline(x=1, color='grey', linestyle='--', linewidth=0.5)
    ax.axvline(x=-1, color='grey', linestyle='--', linewidth=0.5)
    
    ax.set_xlabel('log2 Fold Change', fontsize=11)
    ax.set_ylabel('-log10(Adjusted P-value)', fontsize=11)
    ax.set_title('GSE61616: Sham vs MCAO (Human Orthologs)', fontsize=13, fontweight='bold')
    ax.legend(fontsize=8, loc='upper right')
    
    plt.tight_layout()
    fig.savefig(os.path.join(STAGE_DIR, f"volcano_plot.{FIG_FORMAT}"),
                dpi=FIG_DPI, bbox_inches='tight')
    plt.close()
    logger.info("  火山图已保存")


def main():
    logger.info("=" * 60)
    logger.info("阶段3: DEG功能富集分析 (GO/KEGG/GSEA)")
    logger.info("=" * 60)
    
    logger.info("[1/6] 加载大鼠-人类基因映射库...")
    mapping = load_rat_human_mapping(RAT_MOUSE_HUMAN_MAP)
    
    logger.info("[2/6] 加载DEGs并映射到人类同源基因...")
    degs_file = os.path.join(RESULTS_DIR, "stage1_rma_degs", "limma_degs.csv")
    gene_degs, raw_degs = load_degs_and_map(degs_file, mapping)
    
    gene_degs.to_csv(os.path.join(STAGE_DIR, "human_degs.csv"), index=False)
    
    logger.info("[3/6] 绘制火山图...")
    plot_volcano(gene_degs)
    
    logger.info("[4/6] GO/KEGG富集分析...")
    deg_gene_list = gene_degs['Gene'].tolist()
    
    go_results = run_go_enrichment(deg_gene_list, gene_degs)
    kegg_results = run_kegg_enrichment(deg_gene_list, gene_degs)
    
    logger.info("[5/6] GSEA分析...")
    expr_file = os.path.join(RESULTS_DIR, "stage1_rma_degs", "expr_matrix.csv")
    gsea_results = run_gsea(expr_file, degs_file, mapping)
    
    logger.info("[6/6] 铜死亡基因富集分析...")
    cupro_info = analyze_cuproptosis_enrichment(gene_degs)
    
    logger.info("\n" + "=" * 60)
    logger.info("阶段3完成! 摘要:")
    logger.info(f"  人类DEGs: {len(gene_degs)} (上调: {(gene_degs['logFC']>0).sum()}, 下调: {(gene_degs['logFC']<0).sum()})")
    logger.info(f"  铜死亡核心基因在DEGs: {cupro_info['core_in_degs']}")
    logger.info(f"  铜死亡相关基因在DEGs: {cupro_info['related_in_degs']}")
    if go_results is not None:
        logger.info(f"  GO富集: {len(go_results)} 条")
    if kegg_results is not None:
        logger.info(f"  KEGG富集: {len(kegg_results)} 条")
    logger.info("=" * 60)
    
    return gene_degs, go_results, kegg_results, cupro_info


if __name__ == "__main__":
    main()