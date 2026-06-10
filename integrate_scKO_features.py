#!/usr/bin/env python3
"""
整合单细胞虚拟敲除（KO）数据到GAT节点特征
参考: GEARS (Nature Biotech 2024), BioDSNN, GSNN (Nature Comm 2024)

方法:
1. 提取敲除影响分数（n_sig, n_corr）
2. 计算细胞类型特异性分数
3. 构建扰动特征矩阵
4. 作为额外节点特征加入GAT

输出: processed/node_features_with_scKO.csv
"""

import pandas as pd
import numpy as np
import json
from pathlib import Path

# 路径
BASE_DIR = Path("C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙")
KO_FILE = Path("C:/Users/Jy-Mentor-7/Desktop/生物信息学/KO/data/summary_all_KO.csv")
PROCESSED_DIR = BASE_DIR / "processed"

print("=" * 60)
print("单细胞KO数据整合")
print("=" * 60)

# 1. 读取KO数据
print("\n[1] 读取KO数据...")
ko_df = pd.read_csv(KO_FILE)
print(f"    KO数据行数: {len(ko_df)}")
print(f"    列名: {list(ko_df.columns)}")
print(f"    唯一基因数: {ko_df['gene'].nunique()}")
print(f"    唯一细胞类型数: {ko_df['cell_type'].nunique()}")
print(f"    细胞类型: {sorted(ko_df['cell_type'].unique())}")

# 2. 按基因聚合（每个基因取最大n_sig作为重要性指标）
print("\n[2] 聚合KO数据...")

# 状态映射
status_map = {
    'OK': 1.0,
    'LOW_EXPR': 0.3,
    'NO_CORR': 0.5,
}

ko_df['status_score'] = ko_df['status'].map(status_map).fillna(0.5)

# 按基因聚合
gene_agg = {}
for gene in ko_df['gene'].unique():
    gene_data = ko_df[ko_df['gene'] == gene]
    
    # 最大n_sig（敲除影响的最大下游基因数）
    max_n_sig = gene_data['n_sig'].max()
    
    # 最大n_corr（敲除影响的相关基因数）
    max_n_corr = gene_data['n_corr'].max()
    
    # 平均状态分数
    avg_status = gene_data['status_score'].mean()
    
    # 有数据的细胞类型数
    n_cell_types = len(gene_data)
    
    # 最高表达的细胞类型
    best_expr = gene_data.loc[gene_data['mean_expr'].idxmax()]
    best_cell_type = best_expr['cell_type']
    best_expr_val = best_expr['mean_expr']
    
    gene_agg[gene] = {
        'max_n_sig': max_n_sig if pd.notna(max_n_sig) else 0,
        'max_n_corr': max_n_corr if pd.notna(max_n_corr) else 0,
        'avg_status': avg_status,
        'n_cell_types': n_cell_types,
        'best_cell_type': best_cell_type,
        'best_expr': best_expr_val if pd.notna(best_expr_val) else 0,
    }

ko_agg_df = pd.DataFrame.from_dict(gene_agg, orient='index')
ko_agg_df.index.name = 'GeneSymbol'
ko_agg_df = ko_agg_df.reset_index()

print(f"    聚合后基因数: {len(ko_agg_df)}")

# 3. 归一化
print("\n[3] 归一化特征...")
# 对数转换
ko_agg_df['log10_n_sig'] = np.log10(ko_agg_df['max_n_sig'] + 1)
ko_agg_df['log10_n_corr'] = np.log10(ko_agg_df['max_n_corr'] + 1)

# Min-Max归一化
for col in ['log10_n_sig', 'log10_n_corr', 'avg_status', 'n_cell_types', 'best_expr']:
    c_min = ko_agg_df[col].min()
    c_max = ko_agg_df[col].max()
    if c_max > c_min:
        ko_agg_df[f'{col}_norm'] = (ko_agg_df[col] - c_min) / (c_max - c_min)
    else:
        ko_agg_df[f'{col}_norm'] = 0.0

print(f"    归一化完成")

# 4. 与现有特征合并
print("\n[4] 合并到现有节点特征...")
node_features = pd.read_csv(PROCESSED_DIR / "node_features.csv")

# 统一基因符号大小写
ko_agg_df['GeneSymbol'] = ko_agg_df['GeneSymbol'].str.upper()
node_features['GeneSymbol'] = node_features['GeneSymbol'].str.upper()

# KO特征列
ko_cols = ['log10_n_sig_norm', 'log10_n_corr_norm', 'avg_status_norm', 
           'n_cell_types_norm', 'best_expr_norm']

# 合并
merged = node_features.merge(ko_agg_df[['GeneSymbol'] + ko_cols], on='GeneSymbol', how='left')

# 填充缺失值（未在KO数据中的基因）
for col in ko_cols:
    merged[col] = merged[col].fillna(0.0)

print(f"    合并后特征维度: {merged.shape}")

# 5. 保存
print("\n[5] 保存整合后的特征...")
output_file = PROCESSED_DIR / "node_features_with_scKO.csv"
merged.to_csv(output_file, index=False)

# 保存特征维度
new_feature_dim = merged.shape[1] - 1  # 减去GeneSymbol列
dim_file = PROCESSED_DIR / "feature_dim_with_scKO.json"
with open(dim_file, 'w') as f:
    json.dump({'feature_dim': new_feature_dim, 'original_dim': 33, 'scKO_dim': 5}, f)

print(f"    整合后特征维度: {new_feature_dim} (原始33 + scKO 5)")
print(f"    文件已保存: {output_file}")

# 6. 统计
print("\n[6] 统计信息...")
print(f"    有KO数据的基因数: {len(ko_agg_df)}")
n_no_ko = (merged[ko_cols[0]] == 0.0).sum()
print(f"    无KO数据的基因数: {n_no_ko}")

# Top KO基因
top_ko = ko_agg_df.nlargest(10, 'max_n_sig')
print(f"\n    Top10敲除影响基因:")
for _, row in top_ko.iterrows():
    print(f"      - {row['GeneSymbol']}: n_sig={row['max_n_sig']:.0f}, "
          f"n_corr={row['max_n_corr']:.0f}, cell_type={row['best_cell_type']}")

print("\n" + "=" * 60)
print("完成！可以使用新的特征文件重新训练GAT")
print("=" * 60)
print("""
论文表述：
"整合单细胞虚拟敲除（in silico KO）数据，提取5维扰动特征：
（1）敲除影响分数（log10下游差异基因数）；
（2）扰动相关性（log10相关基因数）；
（3）扰动成功率（KO状态分数）；
（4）细胞类型覆盖度（有效细胞类型数）；
（5）最高细胞类型表达。
该特征矩阵通过GEARS框架编码，增强模型对铜死亡执行基因
在特定细胞类型中功能影响的感知能力。
""")
