#!/usr/bin/env python3
"""
GO/KEGG Enrichment Analysis for Top 10 Hub Genes
Based on cytoNCA + MCODE + MCC Combined Ranking
Using gseapy
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import warnings

warnings.filterwarnings('ignore')
matplotlib.use('Agg')

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.unicode_minus'] = False

def smart_install(package):
    import subprocess
    import sys
    try:
        __import__(package)
    except ImportError:
        subprocess.check_call([sys.executable, "-m", "pip", "install", package, "-q"])

for pkg in ["seaborn", "openpyxl", "scipy", "gseapy"]:
    smart_install(pkg)

def get_top10_genes():
    """获取基于cytoNCA+MCODE+MCC综合排名的前10基因"""
    xls = pd.ExcelFile('cytoNCA_Complete_Results.xlsx')
    df_complete = pd.read_excel(xls, sheet_name='Complete_85_Genes')
    df_hub = pd.read_excel(xls, sheet_name='Hub_Genes_8')
    df_mcode = pd.read_excel(xls, sheet_name='MCODE_Cluster_13')
    df_median = pd.read_excel(xls, sheet_name='Median_Filtered_14')

    hub_genes = set(df_hub['Gene'].tolist())
    mcode_genes = set(df_mcode['Gene'].tolist())
    median_genes = set(df_median['Gene'].tolist())

    df = df_complete.copy()

    df['cytoNCA_Score'] = df[['BC', 'CC', 'DC', 'EC', 'LAC', 'NC', 'SC', 'IC']].mean(axis=1)

    df['Hub_Flag'] = df['Gene'].isin(hub_genes).astype(int)
    df['MCODE_Flag'] = df['Gene'].isin(mcode_genes).astype(int)
    df['Median_Flag'] = df['Gene'].isin(median_genes).astype(int)

    df['Combined_Score'] = (
        df['cytoNCA_Score'] * 0.4 +
        df['MCC'] * 0.4 +
        df['Hub_Flag'] * 0.1 +
        df['MCODE_Flag'] * 0.05 +
        df['Median_Flag'] * 0.05
    )

    df_sorted = df.sort_values('Combined_Score', ascending=False)

    top10 = df_sorted.head(10)['Gene'].tolist()

    print("Top 10 genes by cytoNCA + MCODE + MCC combined ranking:")
    for i, gene in enumerate(top10, 1):
        score = df_sorted[df_sorted['Gene']==gene]['Combined_Score'].values[0]
        mcc = df_sorted[df_sorted['Gene']==gene]['MCC'].values[0]
        dc = df_sorted[df_sorted['Gene']==gene]['DC'].values[0]
        print(f"  {i}. {gene}: Combined={score:.4f}, MCC={mcc:.4f}, DC={dc:.4f}")

    return top10, df_sorted

def simplify_term(term):
    """简化GO/KEGG term名称，去除编号，保留核心含义"""
    import re
    term = re.sub(r'\s*\(GO:\d+\)', '', term)
    term = re.sub(r'\s*\(hsa\d+\)', '', term)
    term = re.sub(r'\s*\(Rattus norvegicus\)', '', term)
    term = re.sub(r'\s*\(Human\)', '', term)
    term = re.sub(r'\s*\(mouse\)', '', term)
    term = re.sub(r'\s*\(Drosophila melanogaster\)', '', term)
    term = re.sub(r'\s*\(Caenorhabditis elegans\)', '', term)
    term = re.sub(r'\s*\(Saccharomyces cerevisiae\)', '', term)
    term = re.sub(r'_Homo_sapiens_Gn', '', term)
    term = re.sub(r'_Mus_musculus_Gn', '', term)
    term = re.sub(r'_Rattus_norvegicus_Gn', '', term)
    term = re.sub(r'\s*\(GR.*?\)', '', term)
    term = re.sub(r'\s*WP\d+', '', term)
    term = re.sub(r'\s*\d+\.\d+', '', term)
    term = term.strip()
    term = ' '.join(term.split())

    kegg_shorten = {
        'AGE-RAGE signaling pathway in diabetic complications': 'AGE-RAGE signaling in diabetic complications',
        'AGE-RAGE signaling pathway in diabetic': 'AGE-RAGE signaling in diabetic',
        'Toll-like receptor signaling pathway': 'TLR signaling pathway',
        'NF-kappa B signaling pathway': 'NF-κB signaling pathway',
        'NOD-like receptor signaling pathway': 'NLR signaling pathway',
        'TNF signaling pathway': 'TNF signaling pathway',
        'IL-17 signaling pathway': 'IL-17 signaling pathway',
        'C-type lectin receptor signaling pathway': 'C-type lectin receptor signaling',
        'RIG-I-like receptor signaling pathway': 'RIG-I-like receptor signaling',
        'Cytosolic DNA-sensing pathway': 'Cytosolic DNA-sensing pathway',
        'Ferroptosis': 'Ferroptosis',
        'Mitophagy': 'Mitophagy',
    }

    for long_name, short_name in kegg_shorten.items():
        if long_name in term:
            term = short_name
            break

    if len(term) > 55:
        words = term.split()
        result = []
        current_len = 0
        for i, word in enumerate(words):
            if current_len + len(word) + (1 if result else 0) > 52:
                break
            result.append(word)
            current_len += len(word) + (1 if result else 0)
        term = ' '.join(result)
        if len(result) < len(words) and result:
            if not result[-1].endswith((':', '-', '/')):
                term = term.rstrip(',;:') + '...'
    return term

def go_kegg_analysis(top10_genes):
    """使用gseapy进行GO/KEGG富集分析"""
    import gseapy as gp

    print(f"\nPerforming GO/KEGG enrichment analysis for {len(top10_genes)} genes...")
    print(f"Organism: Human (human)")
    print(f"Genes: {top10_genes}")

    gene_sets = ['GO_Biological_Process_2021', 'GO_Molecular_Function_2021',
                 'GO_Cellular_Component_2021', 'KEGG_2021_Human']

    results = {}

    for gs in gene_sets:
        try:
            print(f"\n  Analyzing {gs}...")
            enr = gp.enrichr(
                gene_list=top10_genes,
                gene_sets=gs,
                organism='Human',
                outdir=None,
                no_plot=True
            )
            if enr is not None and len(enr.results) > 0:
                results[gs] = enr.results
                print(f"    Found {len(enr.results)} terms")
            else:
                results[gs] = pd.DataFrame()
                print(f"    No significant results")
        except Exception as e:
            print(f"    Error: {e}")
            results[gs] = pd.DataFrame()

    return results

def plot_enrichment_results(results, top10_genes):
    """绘制GO/KEGG富集分析结果 - 含Hub基因通路映射热图"""
    import seaborn as sns
    from matplotlib.colors import LinearSegmentedColormap

    fig = plt.figure(figsize=(22, 16))

    gs = fig.add_gridspec(3, 3, width_ratios=[1.0, 1.0, 0.8], height_ratios=[1, 1, 1.2],
                         wspace=0.35, hspace=0.45,
                         left=0.05, right=0.95, top=0.92, bottom=0.06)

    ax_go_bp = fig.add_subplot(gs[0, 0])
    ax_go_mf = fig.add_subplot(gs[1, 0])
    ax_go_cc = fig.add_subplot(gs[2, 0])
    ax_kegg = fig.add_subplot(gs[:2, 1])
    ax_heatmap = fig.add_subplot(gs[2, 1:])
    ax_gene_list = fig.add_subplot(gs[0:2, 2])

    go_colors = ['#1565C0', '#1976D2', '#1E88E5']
    kegg_color = '#C62828'

    def plot_go_enrichment(ax, df, title, bar_color, n_show=10):
        if len(df) == 0:
            ax.text(0.5, 0.5, 'No significant results', ha='center', va='center',
                   fontsize=11, transform=ax.transAxes)
            ax.set_title(title, fontsize=12, fontweight='bold', color=bar_color)
            return

        df = df.sort_values('Adjusted P-value').head(n_show).copy()
        sig_df = df[df['Adjusted P-value'] < 0.05]
        n_terms = len(df)
        n_sig = len(sig_df)

        y_pos = np.arange(n_terms)
        bar_values = -np.log10(df['Adjusted P-value'] + 1e-10)

        colors = [bar_color if p < 0.05 else '#BDBDBD' for p in df['Adjusted P-value']]
        bars = ax.barh(y_pos, bar_values, height=0.7,
                      color=colors, alpha=0.85, edgecolor='white', linewidth=0.8)

        ax.set_yticks(y_pos)
        ax.set_yticklabels(df['Term_Simple'].tolist(), fontsize=9)
        ax.set_xlabel('−lg(P-value)', fontsize=10, fontweight='bold')
        ax.set_title(title, fontsize=12, fontweight='bold', pad=8, color=bar_color)

        p05_val = -np.log10(0.05)
        max_x = max(bar_values) * 1.35
        ax.set_xlim(0, max_x)
        ax.set_ylim(-0.5, n_terms - 0.5)

        ax.axvline(p05_val, color='#E74C3C', linestyle='--', linewidth=1.2, alpha=0.8)
        ax.text(p05_val + 0.03, n_terms - 0.2, 'P=0.05', ha='left', va='top',
               fontsize=7, color='#E74C3C', fontweight='bold')

        ax.grid(axis='x', alpha=0.2, linestyle='-', linewidth=0.5)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.invert_yaxis()

        for bar, row in zip(bars, df.itertuples()):
            overlap_parts = row.Overlap.split('/')
            gene_count = overlap_parts[0] if len(overlap_parts) >= 1 else 'N/A'
            width = bar.get_width()
            ax.text(width + 0.02, bar.get_y() + bar.get_height()/2,
                   f'{gene_count}', va='center', ha='left', fontsize=7, color='#333333')

        ax.text(0.97, 0.02, f'Sig: {n_sig}/{n_terms}',
               transform=ax.transAxes, fontsize=7, va='bottom', ha='right',
               bbox=dict(boxstyle='round,pad=0.25', facecolor='white', alpha=0.95,
                        edgecolor=bar_color, linewidth=1),
               color=bar_color, fontweight='bold')

    go_configs = [
        ('GO_Biological_Process_2021', 'GO Biological Process', ax_go_bp, go_colors[0]),
        ('GO_Molecular_Function_2021', 'GO Molecular Function', ax_go_mf, go_colors[1]),
        ('GO_Cellular_Component_2021', 'GO Cellular Component', ax_go_cc, go_colors[2])
    ]

    for key, title, ax, bar_color in go_configs:
        if key in results and len(results[key]) > 0:
            df = results[key].copy()
            df['Term_Simple'] = df['Term'].apply(simplify_term)
            plot_go_enrichment(ax, df, title, bar_color, n_show=10)
        else:
            ax.text(0.5, 0.5, 'No significant results', ha='center', va='center',
                   fontsize=11, transform=ax.transAxes)
            ax.set_title(title, fontsize=12, fontweight='bold', color=bar_color)

    if 'KEGG_2021_Human' in results and len(results['KEGG_2021_Human']) > 0:
        df = results['KEGG_2021_Human'].copy()
        df['Term_Simple'] = df['Term'].apply(simplify_term)

        df = df.sort_values('Adjusted P-value').head(12).copy()
        sig_df = df[df['Adjusted P-value'] < 0.05]
        n_terms = len(df)
        n_sig = len(sig_df)

        y_pos = np.arange(n_terms)
        bar_values = -np.log10(df['Adjusted P-value'] + 1e-10)

        colors = [kegg_color if p < 0.05 else '#BDBDBD' for p in df['Adjusted P-value']]
        bars = ax_kegg.barh(y_pos, bar_values, height=0.65,
                           color=colors, alpha=0.85, edgecolor='white', linewidth=0.8)

        ax_kegg.set_yticks(y_pos)
        ax_kegg.set_yticklabels(df['Term_Simple'].tolist(), fontsize=10)
        ax_kegg.set_xlabel('−lg(P-value)', fontsize=11, fontweight='bold')
        ax_kegg.set_title('KEGG Pathway', fontsize=14, fontweight='bold', pad=10, color=kegg_color)

        p05_val = -np.log10(0.05)
        max_x = max(bar_values) * 1.3
        ax_kegg.set_xlim(0, max_x)
        ax_kegg.set_ylim(-0.5, n_terms - 0.5)

        ax_kegg.axvline(p05_val, color='#1565C0', linestyle='--', linewidth=1.3, alpha=0.8)
        ax_kegg.text(p05_val + 0.03, n_terms - 0.2, 'P=0.05', ha='left', va='top',
                    fontsize=8, color='#1565C0', fontweight='bold')

        ax_kegg.grid(axis='x', alpha=0.2, linestyle='-', linewidth=0.5)
        ax_kegg.spines['top'].set_visible(False)
        ax_kegg.spines['right'].set_visible(False)
        ax_kegg.invert_yaxis()

        for bar, row in zip(bars, df.itertuples()):
            overlap_parts = row.Overlap.split('/')
            gene_count = overlap_parts[0] if len(overlap_parts) >= 1 else 'N/A'
            width = bar.get_width()
            ax_kegg.text(width + 0.05, bar.get_y() + bar.get_height()/2,
                        f'{gene_count}', va='center', ha='left', fontsize=8, color='#333333')

        ax_kegg.text(0.97, 0.02, f'Sig: {n_sig}/{n_terms}',
                    transform=ax_kegg.transAxes, fontsize=9, va='bottom', ha='right',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.95,
                             edgecolor=kegg_color, linewidth=1.2),
                    color=kegg_color, fontweight='bold')

        pathway_genes = {}
        for _, row in df.iterrows():
            genes_in_pathway = set(str(row['Genes']).split(';')) & set(top10_genes)
            if genes_in_pathway:
                pathway_genes[row['Term_Simple']] = list(genes_in_pathway)

        ax_gene_list.text(0.5, 0.95, 'Hub Genes in\nTop Pathways', ha='center', va='top',
                         fontsize=11, fontweight='bold', transform=ax_gene_list.transAxes)
        ax_gene_list.axis('off')

        y_offset = 0.88
        for pathway, genes in pathway_genes.items():
            if y_offset < 0.3:
                break
            ax_gene_list.text(0.05, y_offset, f'• {pathway[:35]}',
                            fontsize=8, fontweight='bold', va='top', transform=ax_gene_list.transAxes,
                            wrap=True)
            y_offset -= 0.06
            for gene in genes[:5]:
                ax_gene_list.text(0.12, y_offset, f'  → {gene}',
                                fontsize=8, color='#1565C0', va='top', transform=ax_gene_list.transAxes)
                y_offset -= 0.05
            y_offset -= 0.02

        kegg_terms_for_heatmap = list(pathway_genes.keys())[:10]
        if kegg_terms_for_heatmap:
            heatmap_data = np.zeros((len(top10_genes), len(kegg_terms_for_heatmap)))

            for j, term in enumerate(kegg_terms_for_heatmap):
                for i, gene in enumerate(top10_genes):
                    if gene in pathway_genes.get(term, []):
                        heatmap_data[i, j] = 1

            im = ax_heatmap.imshow(heatmap_data, cmap='Blues', aspect='auto', vmin=0, vmax=1)

            ax_heatmap.set_xticks(np.arange(len(kegg_terms_for_heatmap)))
            short_labels = [t[:18] + '...' if len(t) > 18 else t for t in kegg_terms_for_heatmap]
            ax_heatmap.set_xticklabels(short_labels, rotation=45, ha='right', fontsize=8)
            ax_heatmap.set_yticks(np.arange(len(top10_genes)))
            ax_heatmap.set_yticklabels(top10_genes, fontsize=9, fontweight='bold')
            ax_heatmap.set_xlabel('Pathways', fontsize=10, fontweight='bold')
            ax_heatmap.set_title('Hub Genes-Pathway Membership', fontsize=12, fontweight='bold', pad=8, color='#1565C0')

            for i in range(len(top10_genes)):
                for j in range(len(kegg_terms_for_heatmap)):
                    if heatmap_data[i, j] == 1:
                        ax_heatmap.text(j, i, '●', ha='center', va='center', fontsize=8, color='white')

            cbar = plt.colorbar(im, ax=ax_heatmap, shrink=0.6, pad=0.02)
            cbar.set_ticks([0, 1])
            cbar.set_ticklabels(['No', 'Yes'])
        else:
            ax_heatmap.text(0.5, 0.5, 'No gene-pathway mapping available', ha='center', va='center',
                           fontsize=11, transform=ax_heatmap.transAxes)
            ax_heatmap.axis('off')
    else:
        ax_kegg.text(0.5, 0.5, 'No significant enrichment', ha='center', va='center',
                    fontsize=12, transform=ax_kegg.transAxes)
        ax_kegg.set_title('KEGG Pathway', fontsize=14, fontweight='bold', color=kegg_color)

    gene_str = ', '.join(top10_genes)
    fig.text(0.5, 0.01, f'Top 10 Hub Genes: {gene_str}',
            ha='center', fontsize=10, style='italic', color='#444444')

    plt.suptitle('GO/KEGG Enrichment Analysis | BCP Copper Death Hub Genes',
                fontsize=16, fontweight='bold', y=0.98)

    plt.savefig('GO_KEGG_Enrichment_v6.pdf', format='pdf', dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig('GO_KEGG_Enrichment_v6.png', format='png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    print("\nGO/KEGG enrichment figures saved!")

def save_enrichment_results(results, top10_genes):
    """保存富集分析结果到Excel"""
    with pd.ExcelWriter('GO_KEGG_Enrichment_Results.xlsx') as writer:
        for key in ['GO_Biological_Process_2021', 'GO_Molecular_Function_2021',
                    'GO_Cellular_Component_2021', 'KEGG_2021_Human']:
            if key in results and len(results[key]) > 0:
                df = results[key].copy()
                df['Term_Simple'] = df['Term'].apply(simplify_term)
                df.to_excel(writer, sheet_name=key.replace('_2021', ''), index=False)

        pd.DataFrame({'Top10_Genes': top10_genes}).to_excel(writer, sheet_name='Input_Genes', index=False)

    print("GO/KEGG enrichment results saved to GO_KEGG_Enrichment_Results.xlsx")

def main():
    print("=" * 70)
    print("GO/KEGG Enrichment Analysis")
    print("Based on cytoNCA + MCODE + MCC Combined Ranking")
    print("=" * 70)

    top10_genes, df_sorted = get_top10_genes()

    try:
        results = go_kegg_analysis(top10_genes)
        plot_enrichment_results(results, top10_genes)
        save_enrichment_results(results, top10_genes)
    except Exception as e:
        print(f"\nError in enrichment analysis: {e}")
        import traceback
        traceback.print_exc()

    print("\n" + "=" * 70)
    print("Analysis Complete!")
    print("=" * 70)

if __name__ == "__main__":
    main()
