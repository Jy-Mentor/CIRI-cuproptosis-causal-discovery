#!/usr/bin/env python3
"""
cytoNCA + MCC + MCODE Analysis - 4 Visualization Charts (v9)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
from scipy import stats
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

for pkg in ["seaborn", "openpyxl", "scipy"]:
    smart_install(pkg)

def avoid_overlap(labels_info, offset, existing_positions, min_dist=0.02):
    """检测offset是否会与其他标注重叠，重叠则调整"""
    x, y = labels_info['x'] + offset[0], labels_info['y'] + offset[1]

    for pos in existing_positions:
        dist = np.sqrt((x - pos[0])**2 + (y - pos[1])**2)
        if dist < min_dist:
            return False
    return True

def compute_smart_offset(x, y, existing_positions, base_offset=0.03):
    """计算智能偏移量，避免重叠"""
    angles = np.linspace(0, 2*np.pi, 36, endpoint=False)

    for dist in [base_offset, base_offset*0.8, base_offset*0.6, base_offset*0.4]:
        for angle in angles:
            offset = (dist * np.cos(angle), dist * np.sin(angle))
            if avoid_overlap({'x': x, 'y': y}, offset, existing_positions):
                return offset

    return (base_offset * 0.3, base_offset * 0.3)

def main():
    xls = pd.ExcelFile('cytoNCA_Complete_Results.xlsx')
    df_complete = pd.read_excel(xls, sheet_name='Complete_85_Genes')
    df_hub = pd.read_excel(xls, sheet_name='Hub_Genes_8')
    df_mcode = pd.read_excel(xls, sheet_name='MCODE_Cluster_13')
    df_median = pd.read_excel(xls, sheet_name='Median_Filtered_14')

    print(f"Complete: {df_complete.shape}, Hub: {df_hub.shape}, MCODE: {df_mcode.shape}, Median: {df_median.shape}")

    sns.set_style("white")
    fig = plt.figure(figsize=(22, 20))

    ax1 = fig.add_subplot(2, 2, 1)
    hub_genes = df_hub['Gene'].tolist()
    top8_complete = df_complete[df_complete['Gene'].isin(hub_genes)].copy()
    top8_complete = top8_complete.sort_values('MCC', ascending=False)
    hub_genes_sorted = top8_complete['Gene'].tolist()

    metrics_display = ['BC', 'CC', 'DC', 'EC', 'NC', 'IC', 'MCC']
    x_pos = np.arange(len(hub_genes_sorted))
    width = 0.115

    colors = ['#E74C3C', '#3498DB', '#27AE60', '#8E44AD', '#F39C12', '#E67E22', '#2C3E50']

    for i, metric in enumerate(metrics_display):
        values = [top8_complete[top8_complete['Gene']==g][metric].values[0] for g in hub_genes_sorted]
        ax1.bar(x_pos + i*width, values, width, label=metric, color=colors[i], alpha=0.85, edgecolor='white', linewidth=0.5)

    ax1.set_xlabel('Hub Genes (sorted by MCC descending)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Centrality Score', fontsize=12, fontweight='bold')
    ax1.set_title('Fig1: Top 8 Hub Genes Multi-dimensional Topological Parameters', fontsize=13, fontweight='bold', pad=15)
    ax1.set_xticks(x_pos + width * 3)
    ax1.set_xticklabels([f'{g}\n(MCC={top8_complete[top8_complete["Gene"]==g]["MCC"].values[0]:.3f})' for g in hub_genes_sorted],
                       rotation=45, ha='right', fontsize=9)

    legend_labels = ['BC: Betweenness', 'CC: Closeness', 'DC: Degree', 'EC: Eigenvector', 'NC: Network', 'IC: Information', 'MCC: Maximal Clique']
    ax1.legend(legend_labels, loc='upper center', bbox_to_anchor=(0.5, -0.12), framealpha=0.95, fontsize=9, ncol=4)

    ax1.set_ylim(0, 1.15)
    ax1.grid(axis='y', alpha=0.3, linestyle='--', linewidth=0.5)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    ax2 = fig.add_subplot(2, 2, 2)
    mcc_values = df_complete['MCC'].values
    n_total = len(mcc_values)

    bins = np.linspace(0, 1, 26)
    n, bins_out, patches = ax2.hist(mcc_values, bins=bins, color='#3498DB', alpha=0.7, edgecolor='white', linewidth=0.8)

    kde_x = np.linspace(0, 1, 200)
    kde = stats.gaussian_kde(mcc_values)
    ax2.plot(kde_x, kde(kde_x) * n_total * (bins[1] - bins[0]), 'r-', linewidth=2.5, label='Density Curve (KDE)')

    hub_threshold = df_hub['MCC'].min()
    ax2.axvline(hub_threshold, color='#E74C3C', linestyle='--', linewidth=2.5, label=f'Hub Threshold (MCC={hub_threshold:.3f})')

    ax2.text(0.98, 0.95, f'n = {n_total} genes\nTop 8 Hub Genes Selected\nMCC Threshold = {hub_threshold:.3f}',
            transform=ax2.transAxes, fontsize=10, va='top', ha='right',
            bbox=dict(boxstyle='round', facecolor='lightyellow', edgecolor='gray', alpha=0.9))

    top3 = df_complete.nlargest(3, 'MCC')
    for idx, (_, row) in enumerate(top3.iterrows()):
        ax2.annotate(row['Gene'], xy=(row['MCC'], n.max() * 0.75 - idx * 3),
                    fontsize=11, fontweight='bold', color='#C0392B',
                    bbox=dict(boxstyle='round,pad=0.2', facecolor='#FDEBD0', edgecolor='#E74C3B', alpha=0.9))

    ax2.text(0.02, 0.98, 'Top 3 Hub Genes\nby MCC Score:', transform=ax2.transAxes,
            fontsize=9, va='top', ha='left', fontweight='bold', color='#2C3E50')

    ax2.set_xlabel('MCC Score', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Frequency', fontsize=12, fontweight='bold')
    ax2.set_title('Fig2: MCC Score Distribution (n=85 genes)', fontsize=13, fontweight='bold', pad=15)
    ax2.legend(loc='upper right', fontsize=9)
    ax2.grid(axis='y', alpha=0.3, linestyle='--', linewidth=0.5)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    ax3 = fig.add_subplot(2, 2, 3)
    mcode_genes = df_mcode['Gene'].tolist()
    mcode_full = df_complete[df_complete['Gene'].isin(mcode_genes)].copy()
    mcode_full = mcode_full.sort_values('MCC', ascending=False)

    heatmap_metrics = ['MCC', 'DC', 'CC', 'BC', 'EC', 'NC']
    heatmap_data = mcode_full[heatmap_metrics].values.astype(float)

    for j in range(len(heatmap_metrics)):
        col_min = heatmap_data[:, j].min()
        col_max = heatmap_data[:, j].max()
        if col_max > col_min:
            heatmap_data[:, j] = (heatmap_data[:, j] - col_min) / (col_max - col_min)

    cmap = sns.color_palette("coolwarm", as_cmap=True)
    im = ax3.imshow(heatmap_data, cmap=cmap, aspect='auto', vmin=0, vmax=1)

    ax3.set_xticks(np.arange(len(heatmap_metrics)))
    ax3.set_yticks(np.arange(len(mcode_full)))
    metric_labels = ['MCC', 'DC', 'CC', 'BC', 'EC', 'NC']
    ax3.set_xticklabels(metric_labels, fontsize=11, fontweight='bold')
    ax3.set_yticklabels(mcode_full['Gene'].tolist(), fontsize=11)

    for i in range(len(mcode_full)):
        for j in range(len(heatmap_metrics)):
            val = heatmap_data[i, j]
            text_color = 'white' if val > 0.5 else 'black'
            ax3.text(j, i, f'{val:.2f}', ha='center', va='center',
                    color=text_color, fontsize=10, fontweight='bold')

    ax3.set_title('Fig3: MCODE Cluster Genes (13 genes)\nParameters: nodeScoreCut=0.2, K-core=2 | Sorted by MCC',
                 fontsize=13, fontweight='bold', pad=15)

    cbar = plt.colorbar(im, ax=ax3, shrink=0.8)
    cbar.set_label('Normalized Score (0-1)', fontsize=10)

    ax4 = fig.add_subplot(2, 2, 4)

    bc_median = df_complete['BC'].median()
    cc_median = df_complete['CC'].median()

    dc_normalized = (df_complete['DC'] - df_complete['DC'].min()) / (df_complete['DC'].max() - df_complete['DC'].min())
    point_sizes = 30 + dc_normalized * 120

    scatter = ax4.scatter(df_complete['BC'], df_complete['CC'],
                          c=df_complete['MCC'], s=point_sizes, cmap='Spectral_r', alpha=0.7,
                          edgecolors='white', linewidth=0.5, vmin=0, vmax=1)

    z = np.polyfit(df_complete['BC'], df_complete['CC'], 1)
    p = np.poly1d(z)
    x_line = np.linspace(df_complete['BC'].min(), df_complete['BC'].max(), 100)
    r_val = np.corrcoef(df_complete['BC'], df_complete['CC'])[0,1]
    ax4.plot(x_line, p(x_line), "k--", alpha=0.7, linewidth=2, label=f'Linear Fit (r={r_val:.3f})')

    ax4.axvline(bc_median, color='gray', linestyle=':', alpha=0.6, linewidth=1.5)
    ax4.axhline(cc_median, color='gray', linestyle=':', alpha=0.6, linewidth=1.5)

    hub_sorted = df_hub.sort_values('MCC', ascending=False).reset_index(drop=True)

    existing_positions = []
    hub_coords = []

    for _, row in hub_sorted.iterrows():
        x, y = row['BC'], row['CC']
        offset = compute_smart_offset(x, y, existing_positions, base_offset=0.035)
        existing_positions.append((x + offset[0], y + offset[1]))
        hub_coords.append((x, y, offset[0], offset[1]))

    for x, y, off_x, off_y in hub_coords:
        ax4.annotate(hub_sorted[hub_sorted['BC']==x]['Gene'].values[0],
                    xy=(x, y),
                    xytext=(x + off_x, y + off_y),
                    fontsize=9, fontweight='bold', color='#1A1A1A',
                    bbox=dict(boxstyle='round,pad=0.12', facecolor='#FFFFCC', edgecolor='#E74C3C', alpha=0.95),
                    arrowprops=dict(arrowstyle='->', color='#E74C3C', lw=0.7))

    ax4.set_xlabel('Betweenness Centrality (BC)', fontsize=12, fontweight='bold')
    ax4.set_ylabel('Closeness Centrality (CC)', fontsize=12, fontweight='bold')
    ax4.set_title('Fig4: Network Centrality Scatter (BC vs CC)\n(Point Size=DC, Color=MCC, Dashed=Linear Fit, Dotted=Medians)',
                 fontsize=12, fontweight='bold', pad=15)

    ax4.text(0.02, 0.98, f'BC Median={bc_median:.3f}\nCC Median={cc_median:.3f}',
            transform=ax4.transAxes, fontsize=9, va='top', ha='left',
            bbox=dict(boxstyle='round', facecolor='white', edgecolor='gray', alpha=0.8))

    cbar = plt.colorbar(scatter, ax=ax4, shrink=0.8)
    cbar.set_label('MCC Score', fontsize=10)
    ax4.legend(loc='lower right', fontsize=9)
    ax4.grid(alpha=0.2, linestyle='--', linewidth=0.5)
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)

    plt.tight_layout(pad=4.0)
    plt.savefig('cytoNCA_Analysis_4Charts.pdf', format='pdf', dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.savefig('cytoNCA_Analysis_4Charts.png', format='png', dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close()

    print("\n4 visualization charts (v9) generated!")
    print("  cytoNCA_Analysis_4Charts.pdf")
    print("  cytoNCA_Analysis_4Charts.png")

    print("\n" + "="*70)
    print("Summary Statistics:")
    print("="*70)
    print(f"Total genes: {n_total}")
    print(f"Hub threshold MCC: {hub_threshold:.3f}")
    print(f"Hub genes: {hub_genes_sorted}")
    print(f"BC-CC correlation: r={r_val:.3f}")
    print(f"MCODE cluster: {len(mcode_genes)} genes")

if __name__ == "__main__":
    main()
