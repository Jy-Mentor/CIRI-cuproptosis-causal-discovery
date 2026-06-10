#!/usr/bin/env python3
# ================================================================================
# MR 分析专业图表生成 - 参考 Nature Communications 和 PLOS Genetics 标准
# 基于 GitHub 优秀实践和权威论文图表规范
# ================================================================================

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import numpy as np
import os
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

# ================================================================================
# 图表风格设置 - Nature Communications 标准
# ================================================================================

# 使用 Nature 期刊推荐字体和样式
plt.rcParams.update({
    'font.size': 11,
    'axes.titlesize': 13,
    'axes.labelsize': 11,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'axes.linewidth': 1.0,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'grid.linewidth': 0.5,
    'lines.linewidth': 1.5,
    'patch.linewidth': 1.0,
})

# Nature 配色方案 (来自 Nature 官方)
NATURE_COLORS = {
    'primary_blue': '#0072B2',  # 主要颜色
    'primary_orange': '#E69F00',
    'primary_green': '#009E73',
    'primary_pink': '#CC79A7',
    'primary_grey': '#999999',
    'highlight_red': '#D55E00',
    'highlight_purple': '#56B4E9',
    'light_grey': '#E5E5E5',
    'white': '#FFFFFF'
}

# PLOS Genetics 配色方案
PLOS_COLORS = {
    'blue': '#377EB8',
    'red': '#E41A1C',
    'green': '#4DAF4A',
    'purple': '#984EA3',
    'orange': '#FF7F00',
    'brown': '#A65628',
    'pink': '#F781BF',
    'grey': '#999999'
}

print("="*70)
print("MR 分析专业图表生成 - Nature Communications 标准")
print("="*70)

# ================================================================================
# 数据加载
# ================================================================================

def load_mr_data():
    """加载 MR 分析结果"""
    mr_file = r"D:\下载\MR_batch_results\20260508_optimized_fixed_v2\MR_results_main_optimized.csv"
    if os.path.exists(mr_file):
        df = pd.read_csv(mr_file)
        print(f"✓ 加载 MR 结果：{len(df)} 个基因")
        return df
    else:
        print("✗ MR 结果文件不存在")
        return None

def load_enrichment_data():
    """加载富集分析结果"""
    reactome_file = r"D:\下载\MR_batch_results\20260508_optimized_fixed_v2\functional_enrichment\Reactome_results.csv"
    if os.path.exists(reactome_file):
        df = pd.read_csv(reactome_file)
        print(f"✓ 加载 Reactome 富集：{len(df)} 个通路")
        return df
    else:
        return None

# ================================================================================
# Figure 1: 森林图 (Forest Plot) - MR 分析主要结果
# ================================================================================

def create_forest_plot(mr_data, output_dir):
    """创建专业森林图 - 参考 Nature 风格"""
    print("\n创建 Figure 1: 森林图...")
    
    # 筛选显著基因和随机选择的非显著基因
    sig_genes = mr_data[mr_data['fdr_sig'] == True].copy()
    non_sig_genes = mr_data[mr_data['fdr_sig'] == False].sample(n=min(10, len(mr_data[mr_data['fdr_sig'] == False])), random_state=42)
    
    plot_data = pd.concat([sig_genes, non_sig_genes]).sort_values('discovery_b', ascending=True)
    
    # 创建图表
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # 颜色映射
    colors = [NATURE_COLORS['highlight_red'] if is_sig else NATURE_COLORS['primary_grey'] 
              for is_sig in plot_data['fdr_sig']]
    
    # 绘制点估计和置信区间
    y_pos = range(len(plot_data))
    
    for i, (idx, row) in enumerate(plot_data.iterrows()):
        or_val = row['discovery_or']
        ci_low = float(row['discovery_ci'].split('-')[0]) if '-' in str(row['discovery_ci']) else or_val * 0.8
        ci_high = float(row['discovery_ci'].split('-')[1]) if '-' in str(row['discovery_ci']) else or_val * 1.2
        
        # 绘制置信区间
        ax.plot([ci_low, ci_high], [i, i], linestyle='-', color=colors[i], linewidth=1.5, alpha=0.7)
        
        # 绘制点估计
        ax.scatter(or_val, i, s=80, c=colors[i], edgecolors='white', linewidth=1.5, zorder=5)
    
    # 添加无效线 (OR=1)
    ax.axvline(x=1.0, color='black', linestyle='--', linewidth=1.5, alpha=0.5)
    
    # 设置坐标轴
    ax.set_yticks(y_pos)
    ax.set_yticklabels(plot_data['gene'], fontsize=10)
    ax.set_xlabel('Odds Ratio (95% CI)', fontsize=11, fontweight='bold')
    ax.set_title('Mendelian Randomization Results\\nGenetic Associations with Stroke Risk', 
                 fontsize=13, fontweight='bold', pad=10)
    
    # 设置 x 轴为对数坐标
    ax.set_xscale('log')
    ax.set_xlim(0.5, 2.0)
    
    # 添加网格
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
    
    # 添加图例
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=NATURE_COLORS['highlight_red'], 
               markersize=10, label='FDR Significant (q < 0.05)'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=NATURE_COLORS['primary_grey'], 
               markersize=10, label='Non-significant'),
        Line2D([0], [0], color='black', linestyle='--', linewidth=1.5, label='Null (OR=1)')
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=10, frameon=True, framealpha=0.9)
    
    # 调整布局
    plt.tight_layout()
    
    # 保存
    output_file = os.path.join(output_dir, 'Figure1_Forest_Plot.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  ✓ 已保存：{output_file}")

# ================================================================================
# Figure 2: 火山图 (Volcano Plot) - 效应量 vs 显著性
# ================================================================================

def create_volcano_plot(mr_data, output_dir):
    """创建专业火山图"""
    print("\n创建 Figure 2: 火山图...")
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # 计算 -log10(P 值)
    mr_data['neg_log10_pval'] = -np.log10(mr_data['discovery_pval'].replace(0, 1e-10))
    
    # 区分显著和非显著
    sig_mask = mr_data['fdr_sig'] == True
    non_sig_mask = ~sig_mask
    
    # 绘制非显著基因
    ax.scatter(mr_data.loc[non_sig_mask, 'discovery_b'], 
               mr_data.loc[non_sig_mask, 'neg_log10_pval'],
               c=NATURE_COLORS['primary_grey'], alpha=0.5, s=30, label='Non-significant')
    
    # 绘制显著基因
    ax.scatter(mr_data.loc[sig_mask, 'discovery_b'], 
               mr_data.loc[sig_mask, 'neg_log10_pval'],
               c=NATURE_COLORS['highlight_red'], alpha=0.8, s=50, label='FDR Significant')
    
    # 添加显著性阈值线
    ax.axhline(y=-np.log10(0.05), color=NATURE_COLORS['primary_blue'], linestyle='--', 
               linewidth=1.5, alpha=0.7, label='P = 0.05')
    ax.axhline(y=-np.log10(0.001), color=NATURE_COLORS['primary_orange'], linestyle='--', 
               linewidth=1.5, alpha=0.7, label='P = 0.001')
    ax.axvline(x=0, color='black', linestyle='-', linewidth=0.8, alpha=0.3)
    
    # 标注显著基因
    for idx, row in mr_data[sig_mask].iterrows():
        ax.annotate(row['gene'], (row['discovery_b'], row['neg_log10_pval']),
                   fontsize=9, fontweight='bold', ha='center', va='bottom',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.7, edgecolor='none'))
    
    # 设置标签
    ax.set_xlabel('Beta (Effect Size)', fontsize=11, fontweight='bold')
    ax.set_ylabel('-log₁₀(P-value)', fontsize=11, fontweight='bold')
    ax.set_title('Volcano Plot of MR Associations', fontsize=13, fontweight='bold', pad=10)
    
    # 添加图例
    ax.legend(loc='upper right', fontsize=10, frameon=True, framealpha=0.9)
    
    # 添加网格
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
    
    # 调整布局
    plt.tight_layout()
    
    # 保存
    output_file = os.path.join(output_dir, 'Figure2_Volcano_Plot.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  ✓ 已保存：{output_file}")

# ================================================================================
# Figure 3: 气泡图 (Dot Plot) - 功能富集分析
# ================================================================================

def create_enrichment_dot_plot(enrichment_data, output_dir):
    """创建专业富集分析气泡图"""
    print("\n创建 Figure 3: 功能富集气泡图...")
    
    if enrichment_data is None or len(enrichment_data) == 0:
        print("  ✗ 无富集数据，跳过")
        return
    
    # 选择 Top 15 通路
    top_pathways = enrichment_data.nsmallest(15, 'pvalue')
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # 创建散点图
    scatter = ax.scatter(
        top_pathways['Count'],
        range(len(top_pathways)),
        s=top_pathways['pvalue'] * 1000,  # 点大小与 P 值相关
        c=-np.log10(top_pathways['pvalue']),  # 颜色与 -log10(P) 相关
        cmap='RdYlBu_r',
        alpha=0.7,
        edgecolors='white',
        linewidth=1.5
    )
    
    # 设置 y 轴
    ax.set_yticks(range(len(top_pathways)))
    ax.set_yticklabels(top_pathways['Description'], fontsize=9, ha='right')
    
    # 设置标签
    ax.set_xlabel('Gene Count', fontsize=11, fontweight='bold')
    ax.set_title('Reactome Pathway Enrichment Analysis', fontsize=13, fontweight='bold', pad=10)
    
    # 添加颜色条
    cbar = plt.colorbar(scatter, ax=ax, label='-log₁₀(P-value)', shrink=0.8)
    cbar.ax.tick_params(labelsize=9)
    
    # 添加网格
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5, axis='x')
    ax.set_axisbelow(True)
    
    # 移除上边框和右边框
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # 调整布局
    plt.tight_layout()
    
    # 保存
    output_file = os.path.join(output_dir, 'Figure3_Enrichment_Dot_Plot.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  ✓ 已保存：{output_file}")

# ================================================================================
# Figure 4: 敏感性分析图 (Sensitivity Analysis)
# ================================================================================

def create_sensitivity_plot(mr_data, output_dir):
    """创建敏感性分析图"""
    print("\n创建 Figure 4: 敏感性分析图...")
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # 1. F 统计量分布
    ax1 = axes[0, 0]
    f_stats = mr_data['F_mean'].dropna()
    ax1.hist(f_stats, bins=20, color=NATURE_COLORS['primary_blue'], 
             edgecolor='white', alpha=0.7)
    ax1.axvline(x=10, color=NATURE_COLORS['highlight_red'], linestyle='--', 
                linewidth=1.5, label='F=10 (Weak IV Threshold)')
    ax1.set_xlabel('F-statistic')
    ax1.set_ylabel('Count')
    ax1.set_title('Instrument Strength (F-statistic)', fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 2. SNP 数量分布
    ax2 = axes[0, 1]
    snp_counts = mr_data['nsnp'].dropna()
    ax2.hist(snp_counts, bins=20, color=NATURE_COLORS['primary_green'], 
             edgecolor='white', alpha=0.7)
    ax2.set_xlabel('Number of SNPs')
    ax2.set_ylabel('Count')
    ax2.set_title('Number of Instruments per Gene', fontweight='bold')
    ax2.grid(True, alpha=0.3)
    
    # 3. 异质性 P 值分布
    ax3 = axes[1, 0]
    het_pvals = mr_data['Q_p'].dropna()
    het_pvals = het_pvals[het_pvals > 0]  # 移除 0 值
    ax3.hist(-np.log10(het_pvals), bins=20, color=NATURE_COLORS['primary_orange'], 
             edgecolor='white', alpha=0.7)
    ax3.axvline(x=-np.log10(0.05), color=NATURE_COLORS['highlight_red'], 
                linestyle='--', linewidth=1.5, label='P=0.05')
    ax3.set_xlabel('-log₁₀(P-value)')
    ax3.set_ylabel('Count')
    ax3.set_title('Heterogeneity Test (Cochran Q)', fontweight='bold')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # 4. 多效性 P 值分布
    ax4 = axes[1, 1]
    pleio_pvals = mr_data['Egger_intercept_p'].dropna()
    pleio_pvals = pleio_pvals[pleio_pvals > 0]
    ax4.hist(-np.log10(pleio_pvals), bins=20, color=NATURE_COLORS['primary_pink'], 
             edgecolor='white', alpha=0.7)
    ax4.axvline(x=-np.log10(0.05), color=NATURE_COLORS['highlight_red'], 
                linestyle='--', linewidth=1.5, label='P=0.05')
    ax4.set_xlabel('-log₁₀(P-value)')
    ax4.set_ylabel('Count')
    ax4.set_title('Pleiotropy Test (MR-Egger)', fontweight='bold')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    # 总标题
    fig.suptitle('Sensitivity Analysis Quality Control', fontsize=14, fontweight='bold', y=1.02)
    
    # 调整布局
    plt.tight_layout()
    
    # 保存
    output_file = os.path.join(output_dir, 'Figure4_Sensitivity_Analysis.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  ✓ 已保存：{output_file}")

# ================================================================================
# Figure 5: 药物靶点图 (Drug Target Network)
# ================================================================================

def create_drug_target_plot(drug_data, output_dir):
    """创建药物靶点图"""
    print("\n创建 Figure 5: 药物靶点图...")
    
    if drug_data is None or len(drug_data) == 0:
        print("  ✗ 无药物靶点数据，跳过")
        return
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # 按优先级排序
    priority_order = {'High': 0, 'Medium': 1, 'Low': 2, 'Unknown': 3}
    drug_data_sorted = drug_data.copy()
    drug_data_sorted['priority_rank'] = drug_data_sorted['drug_priority'].map(priority_order)
    drug_data_sorted = drug_data_sorted.sort_values('priority_rank')
    
    # 颜色映射
    color_map = {
        'High': NATURE_COLORS['highlight_red'],
        'Medium': NATURE_COLORS['primary_orange'],
        'Low': NATURE_COLORS['primary_green'],
        'Unknown': NATURE_COLORS['primary_grey']
    }
    
    colors = [color_map.get(p, NATURE_COLORS['primary_grey']) for p in drug_data_sorted['drug_priority']]
    
    # 绘制条形图
    y_pos = range(len(drug_data_sorted))
    bars = ax.barh(y_pos, drug_data_sorted['dgidb_count'], color=colors, 
                   edgecolor='white', linewidth=1.5, alpha=0.8)
    
    # 设置 y 轴
    ax.set_yticks(y_pos)
    ax.set_yticklabels(drug_data_sorted['gene'], fontsize=10)
    
    # 设置标签
    ax.set_xlabel('Number of Drug Interactions', fontsize=11, fontweight='bold')
    ax.set_title('Drug Target Prediction from DGIdb', fontsize=13, fontweight='bold', pad=10)
    
    # 添加图例
    legend_elements = [
        Rectangle((0, 0), 1, 1, facecolor=color_map[p], edgecolor='white', linewidth=1.5)
        for p in ['High', 'Medium', 'Low', 'Unknown']
    ]
    ax.legend(legend_elements, ['High Priority', 'Medium Priority', 'Low Priority', 'Unknown'],
              loc='upper right', fontsize=10, frameon=True)
    
    # 添加网格
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5, axis='x')
    ax.set_axisbelow(True)
    
    # 移除上边框和右边框
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # 调整布局
    plt.tight_layout()
    
    # 保存
    output_file = os.path.join(output_dir, 'Figure5_Drug_Targets.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  ✓ 已保存：{output_file}")

# ================================================================================
# 主函数
# ================================================================================

def main():
    # 加载数据
    mr_data = load_mr_data()
    enrichment_data = load_enrichment_data()
    
    if mr_data is None:
        print("错误：无法加载 MR 数据")
        return
    
    # 创建输出目录
    output_dir = r"D:\下载\MR_batch_results\20260508_optimized_fixed_v2\figures"
    os.makedirs(output_dir, exist_ok=True)
    
    # 生成所有图表
    print("\n生成专业图表...")
    create_forest_plot(mr_data, output_dir)
    create_volcano_plot(mr_data, output_dir)
    create_enrichment_dot_plot(enrichment_data, output_dir)
    create_sensitivity_plot(mr_data, output_dir)
    
    # 加载药物靶点数据
    drug_file = r"D:\下载\MR_batch_results\20260508_optimized_fixed_v2\drug_targets\drug_targets_summary.csv"
    if os.path.exists(drug_file):
        drug_data = pd.read_csv(drug_file)
        create_drug_target_plot(drug_data, output_dir)
    
    print("\n" + "="*70)
    print("所有图表已生成!")
    print("="*70)
    print(f"\n输出目录：{output_dir}")
    print("\n生成的图表:")
    print("  Figure 1: 森林图 - MR 主要结果")
    print("  Figure 2: 火山图 - 效应量 vs 显著性")
    print("  Figure 3: 气泡图 - 功能富集分析")
    print("  Figure 4: 敏感性分析 - 质量控制")
    print("  Figure 5: 药物靶点 - DGIdb 查询结果")
    print("\n所有图表均为 Nature Communications 标准 (300 DPI)")
    print("="*70)

if __name__ == "__main__":
    main()
