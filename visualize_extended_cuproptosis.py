"""
扩展铜死亡基因集可视化 - 30个基因深度分析
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
import warnings
warnings.filterwarnings('ignore')

matplotlib.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial']
matplotlib.rcParams['axes.unicode_minus'] = False
matplotlib.rcParams['figure.dpi'] = 300
matplotlib.rcParams['savefig.dpi'] = 300
matplotlib.rcParams['font.size'] = 10

BASE_DIR = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙"
OUTPUT_DIR = os.path.join(BASE_DIR, "hgt_analysis_results_extended")
SCORES_FILE = os.path.join(BASE_DIR, "bridge_pathway_scores.csv")
CANDIDATES_FILE = os.path.join(BASE_DIR, "bridge_pathway_cirI_candidates.csv")

os.makedirs(OUTPUT_DIR, exist_ok=True)

EXTENDED_GENES = [
    'FDX1', 'LIAS', 'LIPT1', 'DLD', 'DLAT', 'PDHA1', 'PDHB', 
    'MTF1', 'GLS', 'CDKN2A', 'COX11', 'MFN2', 'TOMM20', 'NDUFB9', 
    'ATP6V1E1', 'NFE2L2', 'NLRP3', 'ATP7B', 'ATP7A', 'SLC31A1', 
    'LIPT2', 'DBT', 'GCSH', 'DLST', 'SURF1', 'NDUFB2', 'NDUFB6', 
    'NDUFA8', 'NDUFA1', 'NDUFC1'
]

GENE_CATEGORIES = {
    'FDX1': '核心调控', 'LIAS': '脂酰合成', 'LIPT1': '脂酰转移', 'LIPT2': '脂酰转移',
    'DLD': '脱氢酶', 'DLAT': '丙酮酸脱氢酶', 'PDHA1': '丙酮酸脱氢酶', 'PDHB': '丙酮酸脱氢酶',
    'DBT': '支链酮酸脱氢酶', 'DLST': 'α-酮戊二酸脱氢酶', 'GCSH': '甘氨酸裂解',
    'MTF1': '金属调控', 'SLC31A1': '铜转运', 'ATP7A': '铜转运', 'ATP7B': '铜转运',
    'COX11': '细胞色素c氧化酶', 'SURF1': '线粒体组装',
    'MFN2': '线粒体融合', 'TOMM20': '线粒体输入',
    'NDUFB9': '复合物I', 'NDUFB2': '复合物I', 'NDUFB6': '复合物I', 
    'NDUFA8': '复合物I', 'NDUFA1': '复合物I', 'NDUFC1': '复合物I',
    'ATP6V1E1': 'ATP酶',
    'NFE2L2': '抗氧化', 'NLRP3': '炎症小体', 'GLS': '谷氨酰胺代谢', 'CDKN2A': '细胞周期'
}

CATEGORY_COLORS = {
    '核心调控': '#E91E63', '脂酰合成': '#9C27B0', '脂酰转移': '#673AB7',
    '脱氢酶': '#3F51B5', '丙酮酸脱氢酶': '#2196F3', '金属调控': '#00BCD4',
    '铜转运': '#009688', '细胞色素c氧化酶': '#4CAF50', '线粒体组装': '#8BC34A',
    '线粒体融合': '#CDDC39', '线粒体输入': '#FFEB3B', '复合物I': '#FF9800',
    'ATP酶': '#FF5722', '抗氧化': '#795548', '炎症小体': '#607D8B',
    '谷氨酰胺代谢': '#E91E63', '细胞周期': '#F44336', '支链酮酸脱氢酶': '#3F51B5',
    'α-酮戊二酸脱氢酶': '#2196F3', '甘氨酸裂解': '#4CAF50'
}

print("加载数据...")
df_main = pd.read_csv(SCORES_FILE)
df_candidates = pd.read_csv(CANDIDATES_FILE)
df_detailed = pd.read_csv(os.path.join(OUTPUT_DIR, "gene_detailed_summary.csv"))
df_ranked = pd.read_csv(os.path.join(OUTPUT_DIR, "gene_composite_ranking.csv"))
df_mechanism = pd.read_csv(os.path.join(OUTPUT_DIR, "mechanism_keyword_search.csv"))
df_ciri = pd.read_csv(os.path.join(OUTPUT_DIR, "ciri_cuproptosis_overlap.csv"))

genes_in_main = set(df_main['gene_symbol'].unique())
genes_found = [g for g in EXTENDED_GENES if g in genes_in_main]

# ======================== 图1: 基因存在性与分类 ========================
print("\n生成图1: 基因存在性与功能分类...")

fig, ax = plt.subplots(figsize=(14, 8))

# 按类别统计
category_stats = df_detailed.groupby('category').agg({
    'gene': 'count',
    'mean_score': 'mean',
    'max_score': 'max'
}).reset_index()
category_stats.columns = ['category', 'count', 'mean_score', 'max_score']

x = np.arange(len(category_stats))
width = 0.35

bars1 = ax.bar(x - width/2, category_stats['count'], width, label='基因数量', 
               color='#2196F3', edgecolor='white', alpha=0.8)
bars2 = ax.bar(x + width/2, category_stats['mean_score'] * 100, width, label='平均分×100', 
               color='#FF9800', edgecolor='white', alpha=0.8)

ax.set_xticks(x)
ax.set_xticklabels(category_stats['category'], rotation=45, ha='right', fontsize=9)
ax.set_ylabel('Count / Score×100', fontsize=12, fontweight='bold')
ax.set_title('Extended Cuproptosis Genes - Category Distribution', fontsize=14, fontweight='bold')
ax.legend()
ax.grid(axis='y', alpha=0.3)

for bar in bars1:
    ax.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.3,
            str(int(bar.get_height())), ha='center', va='bottom', fontsize=9, fontweight='bold')

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "fig1_gene_categories.pdf"), bbox_inches='tight', format='pdf')
plt.savefig(os.path.join(OUTPUT_DIR, "fig1_gene_categories.png"), bbox_inches='tight', dpi=300)
print("  ✓ 已保存 fig1_gene_categories")
plt.close()

# ======================== 图2: 综合排名 ========================
print("\n生成图2: 综合排名...")

fig, ax = plt.subplots(figsize=(12, 10))

df_ranked_sorted = df_ranked.sort_values('composite_score', ascending=True)
colors = [CATEGORY_COLORS.get(row['category'], '#999999') for _, row in df_ranked_sorted.iterrows()]

bars = ax.barh(range(len(df_ranked_sorted)), df_ranked_sorted['composite_score'], 
               color=colors, edgecolor='white', alpha=0.85)

ax.set_yticks(range(len(df_ranked_sorted)))
labels = [f"{row['gene']} ({row['category']})" for _, row in df_ranked_sorted.iterrows()]
ax.set_yticklabels(labels, fontsize=10)

# CIRI标记
ciri_genes = set(df_ciri['gene'].unique()) if len(df_ciri) > 0 else set()
for i, (_, row) in enumerate(df_ranked_sorted.iterrows()):
    marker = " ★CIRI" if row['gene'] in ciri_genes else ""
    ax.text(row['composite_score'] + 0.01, i, 
            f"{row['composite_score']:.3f}{marker}", 
            va='center', fontsize=9, fontweight='bold')

ax.set_xlabel('Composite Score', fontsize=12, fontweight='bold')
ax.set_title('Extended Cuproptosis Genes - Composite Ranking\n(Mean 40% + Max 30% + CIRI 20% + Mechanism 10%)', 
             fontsize=14, fontweight='bold')
ax.grid(axis='x', alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "fig2_composite_ranking.pdf"), bbox_inches='tight', format='pdf')
plt.savefig(os.path.join(OUTPUT_DIR, "fig2_composite_ranking.png"), bbox_inches='tight', dpi=300)
print("  ✓ 已保存 fig2_composite_ranking")
plt.close()

# ======================== 图3: 机制关键词匹配 ========================
print("\n生成图3: 机制关键词匹配...")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# 扩展基因匹配数
df_mech_sorted = df_mechanism.sort_values('ext_gene_matches', ascending=True).tail(15)

axes[0].barh(range(len(df_mech_sorted)), df_mech_sorted['ext_gene_matches'], 
             color='#4CAF50', edgecolor='white', alpha=0.8)
axes[0].set_yticks(range(len(df_mech_sorted)))
axes[0].set_yticklabels(df_mech_sorted['keyword'], fontsize=10)
axes[0].set_xlabel('Extended Gene Matches', fontsize=12, fontweight='bold')
axes[0].set_title('Mechanism Keywords - Gene Match Count', fontsize=14, fontweight='bold')
axes[0].grid(axis='x', alpha=0.3)

# 最高分
axes[1].barh(range(len(df_mech_sorted)), df_mech_sorted['max_score'], 
             color='#FF5722', edgecolor='white', alpha=0.8)
axes[1].set_yticks(range(len(df_mech_sorted)))
axes[1].set_yticklabels(df_mech_sorted['keyword'], fontsize=10)
axes[1].set_xlabel('Max Score', fontsize=12, fontweight='bold')
axes[1].set_title('Mechanism Keywords - Max Score', fontsize=14, fontweight='bold')
axes[1].grid(axis='x', alpha=0.3)

for i, (bar, val) in enumerate(zip(axes[1].patches, df_mech_sorted['max_score'])):
    axes[1].text(val + 0.01, i, f'{val:.3f}', va='center', fontsize=9, fontweight='bold')

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "fig3_mechanism_keywords.pdf"), bbox_inches='tight', format='pdf')
plt.savefig(os.path.join(OUTPUT_DIR, "fig3_mechanism_keywords.png"), bbox_inches='tight', dpi=300)
print("  ✓ 已保存 fig3_mechanism_keywords")
plt.close()

# ======================== 图4: CIRI交集基因 ========================
print("\n生成图4: CIRI交集基因分析...")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# 交集展示
total_genes = len(EXTENDED_GENES)
ciri_overlap = len(df_ciri['gene'].unique()) if len(df_ciri) > 0 else 0
no_overlap = total_genes - ciri_overlap

categories = ['Extended\nCuproptosis\nGenes', 'CIRI\nOverlap', 'No\nOverlap']
counts = [total_genes, ciri_overlap, no_overlap]
colors = ['#2196F3', '#4CAF50', '#FF9800']

bars = axes[0].bar(categories, counts, color=colors, edgecolor='white', width=0.6)
for bar, count in zip(bars, counts):
    axes[0].text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.2,
                str(count), ha='center', va='bottom', fontsize=14, fontweight='bold')

axes[0].set_ylabel('Number of Genes', fontsize=12, fontweight='bold')
axes[0].set_title('CIRI Overlap Analysis', fontsize=14, fontweight='bold')
axes[0].grid(axis='y', alpha=0.3)

# 交集基因通路
if len(df_ciri) > 0:
    ciri_genes_list = df_ciri['gene'].unique()
    for gene in ciri_genes_list:
        gene_data = df_ciri[df_ciri['gene'] == gene].sort_values('score', ascending=False)
        axes[1].barh(range(len(gene_data)), gene_data['score'].values, 
                     edgecolor='white', alpha=0.8, label=gene)
        axes[1].set_yticks(range(len(gene_data)))
        axes[1].set_yticklabels([p[:35] + '...' if len(p) > 35 else p for p in gene_data['pathway'].values], 
                                fontsize=8)
        axes[1].invert_yaxis()
        axes[1].set_xlabel('Score', fontsize=10)
        axes[1].set_title('CIRI Overlap Genes - Pathway Scores', fontsize=12, fontweight='bold')
        axes[1].legend(fontsize=9)
        axes[1].grid(axis='x', alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "fig4_ciri_overlap.pdf"), bbox_inches='tight', format='pdf')
plt.savefig(os.path.join(OUTPUT_DIR, "fig4_ciri_overlap.png"), bbox_inches='tight', dpi=300)
print("  ✓ 已保存 fig4_ciri_overlap")
plt.close()

# ======================== 图5: 关键基因Top通路 ========================
print("\n生成图5: 关键基因Top通路...")

key_genes = ['FDX1', 'DLAT', 'CDKN2A', 'ATP7A', 'NDUFB9', 'COX11', 'MFN2', 'NFE2L2']
key_genes_found = [g for g in key_genes if g in genes_found]

fig, axes = plt.subplots(2, 4, figsize=(20, 10))
axes = axes.flatten()

for idx, gene in enumerate(key_genes_found):
    gene_data = df_main[df_main['gene_symbol'] == gene].sort_values('score', ascending=False).head(10)
    
    colors_gene = [CATEGORY_COLORS.get(GENE_CATEGORIES.get(gene, ''), '#2196F3')] * len(gene_data)
    axes[idx].barh(range(len(gene_data)), gene_data['score'].values, 
                   color=colors_gene, edgecolor='white', alpha=0.8)
    axes[idx].set_yticks(range(len(gene_data)))
    axes[idx].set_yticklabels([p[:30] + '...' if len(p) > 30 else p for p in gene_data['pathway_name'].values], 
                              fontsize=8)
    axes[idx].set_xlabel('Score', fontsize=10)
    axes[idx].set_title(f'{gene} - {GENE_CATEGORIES.get(gene, "")}', fontsize=12, fontweight='bold')
    axes[idx].invert_yaxis()
    axes[idx].grid(axis='x', alpha=0.3)

for idx in range(len(key_genes_found), len(axes)):
    axes[idx].axis('off')

plt.suptitle('Key Extended Cuproptosis Genes - Top 10 Pathway Predictions', 
             fontsize=16, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "fig5_key_genes_pathways.pdf"), bbox_inches='tight', format='pdf')
plt.savefig(os.path.join(OUTPUT_DIR, "fig5_key_genes_pathways.png"), bbox_inches='tight', dpi=300)
print("  ✓ 已保存 fig5_key_genes_pathways")
plt.close()

# ======================== 图6: 分数分布对比 ========================
print("\n生成图6: 不同类别基因分数分布对比...")

fig, ax = plt.subplots(figsize=(12, 6))

# 按类别分组
for category in df_detailed['category'].unique():
    cat_genes = df_detailed[df_detailed['category'] == category]['gene'].tolist()
    cat_data = df_main[df_main['gene_symbol'].isin(cat_genes)]['score']
    
    ax.hist(cat_data, bins=50, alpha=0.5, label=category, density=True)

ax.set_xlabel('Prediction Score', fontsize=12, fontweight='bold')
ax.set_ylabel('Density', fontsize=12, fontweight='bold')
ax.set_title('Score Distribution by Gene Category', fontsize=14, fontweight='bold')
ax.legend(fontsize=8, ncol=3)
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "fig6_score_distribution_by_category.pdf"), bbox_inches='tight', format='pdf')
plt.savefig(os.path.join(OUTPUT_DIR, "fig6_score_distribution_by_category.png"), bbox_inches='tight', dpi=300)
print("  ✓ 已保存 fig6_score_distribution_by_category")
plt.close()

print("\n" + "=" * 60)
print("所有可视化图表已生成！")
print(f"输出目录: {OUTPUT_DIR}")
print("=" * 60)
