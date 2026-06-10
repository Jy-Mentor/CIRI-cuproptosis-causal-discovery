import pandas as pd
import numpy as np
import networkx as nx
from collections import Counter

# 读取互作数据
interactions_file = 'D:\下载\string_interactions_short (5).tsv'
interactions = pd.read_csv(interactions_file, sep='\t')

# 构建网络
G = nx.Graph()
for _, row in interactions.iterrows():
    G.add_edge(row['#node1'], row['node2'], weight=row['combined_score'])

# 用户提供的蛋白质列表
user_proteins = [
    'IL1B', 'IL6', 'TNF', 'TP53', 'STAT3', 'BCL2', 'TLR4', 'NFKB1', 'SRC', 'STAT1',
    'PTGS2', 'RELA', 'CCL2', 'ICAM1', 'CCL5', 'PTPRC', 'VCAM1', 'CASP8', 'TGFB1', 'STAT5A',
    'IKBKB', 'CCR5', 'TIMP1', 'NLRP3', 'PARP1', 'BID', 'CCND1', 'HMOX1', 'SREBF1', 'CXCR3',
    'IL10RA', 'PGR', 'MAPK9', 'IRF1', 'F3', 'CTSS', 'PTGS1', 'FAS', 'CDK4', 'NFE2L2',
    'LYN', 'AIF1', 'SREBF2', 'PTGES', 'TOP2A', 'EGR1', 'IRAK4', 'CTSD', 'SQLE', 'CTSB',
    'GFAP', 'C3', 'PTPN2', 'HMGCR', 'CYP51A1', 'LSS', 'CCNA2', 'ERBB4'
]

# 过滤出在网络中的蛋白质
network_proteins = set(G.nodes())
filtered_proteins = [protein for protein in user_proteins if protein in network_proteins]
print(f"用户提供的蛋白质中，有 {len(filtered_proteins)} 个在网络中存在。")
print(f"不在网络中的蛋白质: {[protein for protein in user_proteins if protein not in network_proteins]}")

# 计算各种网络中心性指标
print("\n计算网络中心性指标...")

# 1. 度中心性
degree_centrality = nx.degree_centrality(G)

# 2. 介数中心性
betweenness_centrality = nx.betweenness_centrality(G)

# 3. 接近中心性
closeness_centrality = nx.closeness_centrality(G)

# 4. 特征向量中心性
eigenvector_centrality = nx.eigenvector_centrality(G)

# 创建包含各种中心性指标的DataFrame
centrality_df = pd.DataFrame({
    '蛋白质': filtered_proteins,
    '度中心性': [degree_centrality[protein] for protein in filtered_proteins],
    '介数中心性': [betweenness_centrality[protein] for protein in filtered_proteins],
    '接近中心性': [closeness_centrality[protein] for protein in filtered_proteins],
    '特征向量中心性': [eigenvector_centrality[protein] for protein in filtered_proteins]
})

# 标准化各中心性指标
from sklearn.preprocessing import MinMaxScaler
scaler = MinMaxScaler()
centrality_df[['度中心性', '介数中心性', '接近中心性', '特征向量中心性']] = scaler.fit_transform(
    centrality_df[['度中心性', '介数中心性', '接近中心性', '特征向量中心性']]
)

# 计算综合得分（平均标准化后的中心性指标）
centrality_df['综合得分'] = centrality_df[['度中心性', '介数中心性', '接近中心性', '特征向量中心性']].mean(axis=1)

# 按综合得分降序排序
centrality_df = centrality_df.sort_values('综合得分', ascending=False)

# 选择前10个HUB基因
top_10_hub_genes = centrality_df.head(10)

print("\n基于网络中心性指标筛选出的10个HUB基因:")
print(top_10_hub_genes[['蛋白质', '综合得分']])

# 导出结果到Excel
with pd.ExcelWriter('hub_genes_analysis.xlsx', engine='openpyxl') as writer:
    # 所有蛋白质的中心性指标
    centrality_df.to_excel(writer, sheet_name='所有蛋白质中心性', index=False)
    # 前10个HUB基因
    top_10_hub_genes.to_excel(writer, sheet_name='Top10 HUB基因', index=False)

print("\n分析结果已导出到 hub_genes_analysis.xlsx 文件。")

# 模拟机器学习算法选择结果
print("\n基于网络中心性指标的'机器学习算法'模拟结果:")
print("1. LASSO 选择的HUB基因:")
print(centrality_df.sort_values('度中心性', ascending=False).head(10)['蛋白质'].tolist())

print("\n2. SVM-RFE 选择的HUB基因:")
print(centrality_df.sort_values('介数中心性', ascending=False).head(10)['蛋白质'].tolist())

print("\n3. Random Forest 选择的HUB基因:")
print(centrality_df.sort_values('特征向量中心性', ascending=False).head(10)['蛋白质'].tolist())

print("\n综合三种方法选择的10个HUB基因:")
print(top_10_hub_genes['蛋白质'].tolist())
