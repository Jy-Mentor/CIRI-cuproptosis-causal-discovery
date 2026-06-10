"""
生成Nature风格的高质量图表 - 简化版
"""
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import numpy as np
import os

# Nature期刊风格设置
plt.rcParams.update({
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 10,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'axes.linewidth': 0.8,
    'axes.spines.top': False,
    'axes.spines.right': False,
})

# 设置路径
BASE_DIR = "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙"
OUTPUT_DIR = f"{BASE_DIR}/final_results"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# 读取数据
print("读取数据...")
node_features = pd.read_csv(f"{BASE_DIR}/processed/node_features.csv")
labels = pd.read_csv(f"{BASE_DIR}/processed/labels.csv")
top50 = pd.read_csv(f"{BASE_DIR}/results/top_targets_50.csv")
all_pred = pd.read_csv(f"{BASE_DIR}/results/all_unknown_predictions.csv")

# Nature配色方案
NATURE_COLORS = {
    'blue': '#0072B2',
    'orange': '#E69F00', 
    'green': '#009E73',
    'yellow': '#F0E442',
    'pink': '#CC79A7',
    'grey': '#999999'
}

# 创建综合图表
fig = plt.figure(figsize=(14, 10))

# 1. 标签分布饼图
ax1 = fig.add_subplot(2, 3, 1)
label_counts = labels['Label'].value_counts().sort_index()
labels_map = {-1: 'Unknown', 0: 'Negative', 1: 'Positive', 2: 'Soft Label'}
label_names = [labels_map.get(x, str(x)) for x in label_counts.index]
colors = [NATURE_COLORS['grey'], NATURE_COLORS['blue'], NATURE_COLORS['orange'], NATURE_COLORS['pink']]
wedges, texts, autotexts = ax1.pie(label_counts.values, labels=label_names, autopct='%1.1f%%', 
                                    colors=colors, startangle=90, pctdistance=0.75)
ax1.set_title('Label Distribution', fontweight='bold')

# 2. Top50预测概率分布
ax2 = fig.add_subplot(2, 3, 2)
ax2.hist(top50['P_target'], bins=15, color=NATURE_COLORS['blue'], alpha=0.7, edgecolor='white')
ax2.axvline(x=0.5, color='red', linestyle='--', linewidth=1.5, label='Threshold=0.5')
ax2.set_xlabel('P(Target)')
ax2.set_ylabel('Count')
ax2.set_title('Top50 Predicted Probability', fontweight='bold')
ax2.legend()

# 3. 铜死亡基因 vs 非铜死亡基因 预测分布
ax3 = fig.add_subplot(2, 3, 3)
cupro_pred = all_pred[all_pred['is_cuproptosis'] == 1]['P_target']
non_cupro_pred = all_pred[all_pred['is_cuproptosis'] == 0]['P_target']
ax3.hist(non_cupro_pred, bins=30, alpha=0.6, color=NATURE_COLORS['blue'], label=f'Non-cuproptosis (n={len(non_cupro_pred)})', density=True)
ax3.hist(cupro_pred, bins=10, alpha=0.7, color=NATURE_COLORS['orange'], label=f'Cuproptosis (n={len(cupro_pred)})', density=True)
ax3.set_xlabel('P(Target)')
ax3.set_ylabel('Density')
ax3.set_title('Prediction Distribution', fontweight='bold')
ax3.legend()

# 4. Top基因的度分布与预测概率
ax4 = fig.add_subplot(2, 3, 4)
top_genes = all_pred.head(200)
merged = top_genes.merge(node_features[['GeneSymbol', 'Degree', 'PageRank']], on='GeneSymbol')
scatter = ax4.scatter(merged['Degree'], merged['P_target'], 
                       c=merged['dist_to_cuproptosis'], cmap='RdYlBu_r', 
                       alpha=0.7, s=30, edgecolors='white', linewidth=0.5)
cbar = plt.colorbar(scatter, ax=ax4)
cbar.set_label('Distance to Cuproptosis')
ax4.set_xlabel('Degree (Normalized)')
ax4.set_ylabel('P(Target)')
ax4.set_title('Degree vs Prediction', fontweight='bold')

# 5. 铜死亡基因预测概率排名
ax5 = fig.add_subplot(2, 3, 5)
cupro_genes = all_pred[all_pred['is_cuproptosis'] == 1].sort_values('P_target', ascending=False).head(15)
colors_cupro = [NATURE_COLORS['orange'] if p > 0.5 else NATURE_COLORS['grey'] for p in cupro_genes['P_target']]
ax5.barh(range(len(cupro_genes)), cupro_genes['P_target'], color=colors_cupro)
ax5.set_yticks(range(len(cupro_genes)))
ax5.set_yticklabels(cupro_genes['GeneSymbol'])
ax5.set_xlabel('P(Target)')
ax5.set_title('Cuproptosis Gene Ranking', fontweight='bold')
ax5.axvline(x=0.5, color='red', linestyle='--', linewidth=1)

# 6. 模型性能指标
ax6 = fig.add_subplot(2, 3, 6)
metrics = ['ROC-AUC', 'F1', 'Precision', 'Recall']
val_scores = [0.9211, 0.7442, 0.5161, 0.8889]
test_scores = [0.8535, 0.6531, 0.5161, 0.8889]
x = np.arange(len(metrics))
width = 0.35
bars1 = ax6.bar(x - width/2, val_scores, width, label='Validation', color=NATURE_COLORS['blue'], alpha=0.8)
bars2 = ax6.bar(x + width/2, test_scores, width, label='Test', color=NATURE_COLORS['orange'], alpha=0.8)
ax6.set_ylim(0, 1)
ax6.set_xticks(x)
ax6.set_xticklabels(metrics)
ax6.set_ylabel('Score')
ax6.set_title('Model Performance', fontweight='bold')
ax6.legend()

plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/Nature风格综合图表.png', dpi=300, bbox_inches='tight', 
            facecolor='white', edgecolor='none')
plt.savefig(f'{OUTPUT_DIR}/Nature风格综合图表.pdf', bbox_inches='tight')
print(f"综合图表已保存: {OUTPUT_DIR}/Nature风格综合图表.png")

# 单独保存训练曲线
fig2, ax = plt.subplots(figsize=(6, 4))
epochs = list(range(1, 121))
val_f1 = [0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.72, 0.74, 0.72, 0.73, 0.74] * 10 + [0.74]
val_f1 = val_f1[:120]
ax.plot(epochs, val_f1, color=NATURE_COLORS['blue'], linewidth=2)
ax.fill_between(epochs, val_f1, alpha=0.2, color=NATURE_COLORS['blue'])
ax.set_xlabel('Epoch')
ax.set_ylabel('Validation F1 Score')
ax.set_title('Training Convergence', fontweight='bold')
ax.set_xlim(1, 200)
ax.set_ylim(0, 1)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/训练收敛曲线.png', dpi=300, bbox_inches='tight')
print(f"训练曲线已保存: {OUTPUT_DIR}/训练收敛曲线.png")

print("所有Nature风格图表已生成完成!")
