"""
HGT模型预测结果可视化 - 石竹烯靶基因与铜死亡机制
生成出版级可视化图表
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
import warnings
warnings.filterwarnings('ignore')

# 设置中文字体支持
matplotlib.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial']
matplotlib.rcParams['axes.unicode_minus'] = False
matplotlib.rcParams['figure.dpi'] = 300
matplotlib.rcParams['savefig.dpi'] = 300
matplotlib.rcParams['font.size'] = 10

# 文件路径
BASE_DIR = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙"
OUTPUT_DIR = os.path.join(BASE_DIR, "hgt_analysis_results")
SCORES_FILE = os.path.join(BASE_DIR, "bridge_pathway_scores.csv")
CANDIDATES_FILE = os.path.join(BASE_DIR, "bridge_pathway_cirI_candidates.csv")

os.makedirs(OUTPUT_DIR, exist_ok=True)

# 铜死亡核心基因
CUPROPTOSIS_GENES = ['FDX1', 'LIAS', 'LIPT1', 'DLD', 'DLAT', 'PDHA1', 'PDHB', 
                     'MTF1', 'GLS', 'CDKN2A', 'SLC31A1', 'ATP7A', 'ATP7B', 
                     'DBT', 'DLST', 'PDHA2', 'GCSH']

# 加载数据
print("加载数据...")
df_main = pd.read_csv(SCORES_FILE)
df_candidates = pd.read_csv(CANDIDATES_FILE)

print(f"主数据: {len(df_main)} 条记录")
print(f"候选数据: {len(df_candidates)} 条记录")

# ======================== 图1: 分数分布直方图 ========================
print("\n生成图1: 分数分布直方图...")

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# 整体分布
axes[0].hist(df_main['score'], bins=100, color='#2196F3', edgecolor='white', alpha=0.8)
axes[0].axvline(df_main['score'].mean(), color='#FF5722', linestyle='--', linewidth=2, 
                label=f'Mean: {df_main["score"].mean():.4f}')
axes[0].axvline(df_main['score'].median(), color='#4CAF50', linestyle='--', linewidth=2,
                label=f'Median: {df_main["score"].median():.4f}')
axes[0].set_xlabel('Prediction Score', fontsize=12, fontweight='bold')
axes[0].set_ylabel('Frequency', fontsize=12, fontweight='bold')
axes[0].set_title('Overall Score Distribution', fontsize=14, fontweight='bold')
axes[0].legend(fontsize=11)
axes[0].grid(axis='y', alpha=0.3)

# Top基因分布
gene_stats = df_main.groupby('gene_symbol')['score'].agg(['mean', 'count']).reset_index()
gene_stats = gene_stats.sort_values('mean', ascending=False).head(30)

axes[1].barh(range(len(gene_stats)), gene_stats['mean'], color='#9C27B0', edgecolor='white')
axes[1].set_yticks(range(len(gene_stats)))
axes[1].set_yticklabels(gene_stats['gene_symbol'], fontsize=9)
axes[1].set_xlabel('Mean Prediction Score', fontsize=12, fontweight='bold')
axes[1].set_title('Top 30 Genes by Mean Score', fontsize=14, fontweight='bold')
axes[1].invert_yaxis()
axes[1].grid(axis='x', alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "fig1_score_distribution.pdf"), bbox_inches='tight', format='pdf')
plt.savefig(os.path.join(OUTPUT_DIR, "fig1_score_distribution.png"), bbox_inches='tight', dpi=300)
print("  ✓ 已保存 fig1_score_distribution")
plt.close()

# ======================== 图2: 铜死亡基因分析 ========================
print("\n生成图2: 铜死亡核心基因分析...")

# 找出存在的铜死亡基因
genes_in_data = set(df_main['gene_symbol'].unique())
cuproptosis_found = [g for g in CUPROPTOSIS_GENES if g in genes_in_data]

# 获取这些基因的分数统计
cuproptosis_data = df_main[df_main['gene_symbol'].isin(cuproptosis_found)]
cuproptosis_gene_stats = cuproptosis_data.groupby('gene_symbol')['score'].agg(['mean', 'median', 'std', 'max', 'count']).reset_index()
cuproptosis_gene_stats = cuproptosis_gene_stats.sort_values('mean', ascending=False)

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# 箱线图
genes_to_plot = cuproptosis_gene_stats['gene_symbol'].tolist()
data_to_plot = [cuproptosis_data[cuproptosis_data['gene_symbol'] == g]['score'].values for g in genes_to_plot]

bp = axes[0].boxplot(data_to_plot, labels=genes_to_plot, patch_artist=True, 
                      showfliers=False, widths=0.6)

# 设置颜色
colors = plt.cm.Set3(np.linspace(0, 1, len(genes_to_plot)))
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.8)

axes[0].set_ylabel('Prediction Score', fontsize=12, fontweight='bold')
axes[0].set_title('Cuproptosis Core Genes - Score Distribution', fontsize=14, fontweight='bold')
axes[0].tick_params(axis='x', rotation=45, labelsize=9)
axes[0].grid(axis='y', alpha=0.3)

# 柱状图 - 平均分
x = np.arange(len(cuproptosis_gene_stats))
bars = axes[1].bar(x, cuproptosis_gene_stats['mean'], 
                   yerr=cuproptosis_gene_stats['std'],
                   capsize=5, color='#FF9800', edgecolor='white', alpha=0.8)

axes[1].set_xticks(x)
axes[1].set_xticklabels(cuproptosis_gene_stats['gene_symbol'], rotation=45, ha='right', fontsize=9)
axes[1].set_ylabel('Mean Prediction Score', fontsize=12, fontweight='bold')
axes[1].set_title('Cuproptosis Genes - Mean Score ± SD', fontsize=14, fontweight='bold')
axes[1].grid(axis='y', alpha=0.3)

# 添加数值标签
for bar, mean_val in zip(bars, cuproptosis_gene_stats['mean']):
    axes[1].text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.01,
                f'{mean_val:.3f}', ha='center', va='bottom', fontsize=8, fontweight='bold')

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "fig2_cuproptosis_genes.pdf"), bbox_inches='tight', format='pdf')
plt.savefig(os.path.join(OUTPUT_DIR, "fig2_cuproptosis_genes.png"), bbox_inches='tight', dpi=300)
print("  ✓ 已保存 fig2_cuproptosis_genes")
plt.close()

# ======================== 图3: 机制相关通路热图 ========================
print("\n生成图3: 铜死亡机制相关通路热图...")

# 机制关键词
MECHANISM_KEYWORDS = ['TCA cycle', 'oxidative phosphorylation', 'metal ion', 
                      'ferroptosis', 'mitochondrial', 'oxidative stress',
                      'apoptosis', 'cell death']

# 为每个关键词找top基因
mechanism_top_genes = {}
for keyword in MECHANISM_KEYWORDS:
    matches = df_main[df_main['pathway_name'].str.contains(keyword, case=False, na=False)]
    if len(matches) > 0:
        top_genes = matches.groupby('gene_symbol')['score'].mean().sort_values(ascending=False).head(10)
        mechanism_top_genes[keyword] = top_genes

# 创建热图数据
all_genes = list(set().union(*[g.index.tolist() for g in mechanism_top_genes.values()]))
all_genes = sorted(all_genes)

heatmap_data = pd.DataFrame(index=all_genes, columns=MECHANISM_KEYWORDS)

for keyword in MECHANISM_KEYWORDS:
    if keyword in mechanism_top_genes:
        for gene in mechanism_top_genes[keyword].index:
            heatmap_data.loc[gene, keyword] = mechanism_top_genes[keyword][gene]

heatmap_data = heatmap_data.fillna(0).astype(float)

fig, ax = plt.subplots(figsize=(12, 8))
im = ax.imshow(heatmap_data.values, cmap='YlOrRd', aspect='auto')

ax.set_xticks(range(len(MECHANISM_KEYWORDS)))
ax.set_xticklabels(MECHANISM_KEYWORDS, rotation=45, ha='right', fontsize=10)
ax.set_yticks(range(len(all_genes)))
ax.set_yticklabels(all_genes, fontsize=9)

# 添加数值
for i in range(len(all_genes)):
    for j in range(len(MECHANISM_KEYWORDS)):
        val = heatmap_data.values[i, j]
        if val > 0:
            ax.text(j, i, f'{val:.3f}', ha='center', va='center', fontsize=6,
                   color='black' if val < 0.5 else 'white')

ax.set_title('Cuproptosis Mechanism-Related Pathways\nTop Gene-Keyword Associations', 
             fontsize=14, fontweight='bold')
plt.colorbar(im, ax=ax, label='Mean Prediction Score')
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "fig3_mechanism_heatmap.pdf"), bbox_inches='tight', format='pdf')
plt.savefig(os.path.join(OUTPUT_DIR, "fig3_mechanism_heatmap.png"), bbox_inches='tight', dpi=300)
print("  ✓ 已保存 fig3_mechanism_heatmap")
plt.close()

# ======================== 图4: CIRI候选基因分析 ========================
print("\n生成图4: CIRI候选基因与铜死亡基因交集...")

# 找出交集
candidate_genes = set(df_candidates['gene_symbol'].unique())
overlap_genes = candidate_genes.intersection(set(CUPROPTOSIS_GENES))

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Venn图风格的展示
categories = ['Cuproptosis\nGenes', 'CIRI\nCandidates', 'Overlap']
counts = [len(CUPROPTOSIS_GENES), len(candidate_genes), len(overlap_genes)]

colors = ['#FF5722', '#2196F3', '#4CAF50']
bars = axes[0].bar(categories, counts, color=colors, edgecolor='white', width=0.6)

for bar, count in zip(bars, counts):
    axes[0].text(bar.get_x() + bar.get_width()/2., bar.get_height() + 5,
                str(count), ha='center', va='bottom', fontsize=14, fontweight='bold')

axes[0].set_ylabel('Number of Genes', fontsize=12, fontweight='bold')
axes[0].set_title('Gene Set Overlap Analysis', fontsize=14, fontweight='bold')
axes[0].grid(axis='y', alpha=0.3)

# 交集基因的分数
if len(overlap_genes) > 0:
    overlap_data = df_candidates[df_candidates['gene_symbol'].isin(overlap_genes)]
    overlap_stats = overlap_data.groupby('gene_symbol')['score'].agg(['mean', 'max']).reset_index()
    overlap_stats = overlap_stats.sort_values('mean', ascending=False)
    
    x = np.arange(len(overlap_stats))
    bars = axes[1].bar(x, overlap_stats['mean'], color='#4CAF50', edgecolor='white', alpha=0.8)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(overlap_stats['gene_symbol'], rotation=45, ha='right', fontsize=11)
    axes[1].set_ylabel('Mean Score in CIRI Candidates', fontsize=12, fontweight='bold')
    axes[1].set_title('Overlap Genes - CIRI Candidate Scores', fontsize=14, fontweight='bold')
    axes[1].grid(axis='y', alpha=0.3)
    
    for bar, val in zip(bars, overlap_stats['mean']):
        axes[1].text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.01,
                    f'{val:.3f}', ha='center', va='bottom', fontsize=9, fontweight='bold')

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "fig4_ciri_overlap.pdf"), bbox_inches='tight', format='pdf')
plt.savefig(os.path.join(OUTPUT_DIR, "fig4_ciri_overlap.png"), bbox_inches='tight', dpi=300)
print("  ✓ 已保存 fig4_ciri_overlap")
plt.close()

# ======================== 图5: 关键基因通路排名 ========================
print("\n生成图5: 关键铜死亡基因Top通路...")

# 选择几个关键基因展示
key_genes = ['FDX1', 'DLAT', 'PDHA1', 'LIAS', 'CDKN2A', 'ATP7A']
key_genes_found = [g for g in key_genes if g in genes_in_data]

fig, axes = plt.subplots(2, 3, figsize=(16, 10))
axes = axes.flatten()

for idx, gene in enumerate(key_genes_found):
    gene_data = df_main[df_main['gene_symbol'] == gene].sort_values('score', ascending=False).head(10)
    
    axes[idx].barh(range(len(gene_data)), gene_data['score'].values, 
                   color='#E91E63', edgecolor='white', alpha=0.8)
    axes[idx].set_yticks(range(len(gene_data)))
    axes[idx].set_yticklabels([p[:40] + '...' if len(p) > 40 else p for p in gene_data['pathway_name'].values], 
                              fontsize=8)
    axes[idx].set_xlabel('Score', fontsize=10)
    axes[idx].set_title(f'{gene} - Top 10 Pathways', fontsize=12, fontweight='bold')
    axes[idx].invert_yaxis()
    axes[idx].grid(axis='x', alpha=0.3)

# 隐藏多余的子图
for idx in range(len(key_genes_found), len(axes)):
    axes[idx].axis('off')

plt.suptitle('Key Cuproptosis Genes - Top Pathway Predictions', fontsize=16, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "fig5_key_genes_pathways.pdf"), bbox_inches='tight', format='pdf')
plt.savefig(os.path.join(OUTPUT_DIR, "fig5_key_genes_pathways.png"), bbox_inches='tight', dpi=300)
print("  ✓ 已保存 fig5_key_genes_pathways")
plt.close()

print("\n" + "=" * 60)
print("所有可视化图表已生成完成！")
print(f"输出目录: {OUTPUT_DIR}")
print("=" * 60)
