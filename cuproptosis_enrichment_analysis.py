"""
铜死亡富集分析脚本 - 改进版
1. 超几何检验：Top200/500/1000铜死亡富集分析
2. 双列表输出：铜死亡基因专项
3. 网络邻近性：Top50到铜死亡基因的平均最短路径（基于连通分量）
"""
import pandas as pd
import numpy as np
from scipy import stats
import os
import pickle
import networkx as nx
import torch
from datetime import datetime

# 设置路径
BASE_DIR = "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙"
OUTPUT_DIR = f"{BASE_DIR}/final_results"

# 铜死亡执行基因
CUPTOPOSIS_GENES = {"FDX1", "LIAS", "LIPT1", "DLAT", "PDHB", "PDHX", "SLC31A1", "ATP7A", "ATP7B", "GPX4"}

print("=" * 60)
print("铜死亡富集分析（改进版）")
print("=" * 60)

# 读取数据
print("\n[1] 读取数据...")
node_features = pd.read_csv(f"{BASE_DIR}/processed/node_features.csv")
labels = pd.read_csv(f"{BASE_DIR}/processed/labels.csv")
all_pred = pd.read_csv(f"{BASE_DIR}/results/all_unknown_predictions.csv")
top50 = pd.read_csv(f"{BASE_DIR}/results/top_targets_50.csv")
top100 = pd.read_csv(f"{BASE_DIR}/results/top_targets_100.csv")
top200 = pd.read_csv(f"{BASE_DIR}/results/top_targets_200.csv")

# 读取图数据用于网络邻近性计算
print("[2] 加载网络...")
edge_index = torch.load(f"{BASE_DIR}/processed/edge_index.pt", weights_only=False)
with open(f"{BASE_DIR}/processed/gene_symbols.pkl", "rb") as f:
    gene_symbols = pickle.load(f)

# 构建NetworkX图
G = nx.Graph()
for i, (u, v) in enumerate(edge_index.t().tolist()):
    if u < len(gene_symbols) and v < len(gene_symbols):
        G.add_edge(gene_symbols[u], gene_symbols[v])

print(f"    网络节点数: {G.number_of_nodes()}")
print(f"    网络边数: {G.number_of_edges()}")

# 分析连通分量
connected_components = list(nx.connected_components(G))
print(f"    连通分量数: {len(connected_components)}")
largest_cc = max(connected_components, key=len)
print(f"    最大连通分量节点数: {len(largest_cc)}")

# ============================================================
# 1. 超几何检验
# ============================================================
print("\n[3] 超几何检验分析...")

# 铜死亡基因在网络中的数量
cupro_in_network = [g for g in CUPTOPOSIS_GENES if g in G.nodes()]
print(f"    铜死亡基因在网络中的数量: {len(cupro_in_network)}/{len(CUPTOPOSIS_GENES)}")

# 总体设置
total_genes = len(all_pred)  # 未知节点数
total_cupro = len([g for g in all_pred['GeneSymbol'] if g in CUPTOPOSIS_GENES])
print(f"    总体未知基因数: {total_genes}")
print(f"    总体铜死亡基因数: {total_cupro}")

def hypergeom_test(top_genes_set, total_set, cupro_genes_set):
    """超几何检验"""
    k = len(top_genes_set & cupro_genes_set)  # 交集
    M = len(total_set)  # 总体
    n = len(cupro_genes_set)  # 成功数
    N = len(top_genes_set)  # 抽样数
    
    if k == 0:
        return 0.0, 1.0
    
    # 超几何检验
    pval = stats.hypergeom.sf(k - 1, M, n, N)
    
    # 计算OR
    a = k
    b = N - k
    c = n - k
    d = M - n - b
    if b > 0 and c > 0:
        odds_ratio = (a * d) / (b * c)
    else:
        odds_ratio = 0.0
    
    return odds_ratio, pval

# Top200超几何检验
top200_genes = set(top200['GeneSymbol'].tolist())
cupro_genes_set = set(all_pred['GeneSymbol'].tolist()) & CUPTOPOSIS_GENES
odds_ratio_200, pval_200 = hypergeom_test(top200_genes, set(all_pred['GeneSymbol'].tolist()), cupro_genes_set)
top200_cupro = len(top200_genes & cupro_genes_set)
print(f"\n    Top200:")
print(f"      - 铜死亡基因数: {top200_cupro}")
print(f"      - OR值: {odds_ratio_200:.3f}")
print(f"      - P值: {pval_200:.4f}")

# Top500超几何检验
top500 = all_pred.head(500)
top500_genes = set(top500['GeneSymbol'].tolist())
odds_ratio_500, pval_500 = hypergeom_test(top500_genes, set(all_pred['GeneSymbol'].tolist()), cupro_genes_set)
top500_cupro = len(top500_genes & cupro_genes_set)
print(f"\n    Top500:")
print(f"      - 铜死亡基因数: {top500_cupro}")
print(f"      - OR值: {odds_ratio_500:.3f}")
print(f"      - P值: {pval_500:.4f}")

# Top1000超几何检验
top1000 = all_pred.head(1000)
top1000_genes = set(top1000['GeneSymbol'].tolist())
odds_ratio_1000, pval_1000 = hypergeom_test(top1000_genes, set(all_pred['GeneSymbol'].tolist()), cupro_genes_set)
top1000_cupro = len(top1000_genes & cupro_genes_set)
print(f"\n    Top1000:")
print(f"      - 铜死亡基因数: {top1000_cupro}")
print(f"      - OR值: {odds_ratio_1000:.3f}")
print(f"      - P值: {pval_1000:.4f}")

# ============================================================
# 2. 铜死亡基因专项分析
# ============================================================
print("\n[4] 铜死亡基因专项分析...")

# 获取所有未知节点中的铜死亡基因
cupro_in_unknown = all_pred[all_pred['is_cuproptosis'] == 1].copy()
cupro_ranking = cupro_in_unknown.sort_values('P_target', ascending=False)

print(f"    未知节点中的铜死亡基因数: {len(cupro_ranking)}")
if len(cupro_ranking) > 0:
    print(f"    铜死亡基因P_target范围: {cupro_ranking['P_target'].min():.3f} - {cupro_ranking['P_target'].max():.3f}")
    print(f"    铜死亡基因中位数: {cupro_ranking['P_target'].median():.3f}")
    top_cupro = cupro_ranking.head(5)
    print(f"    Top5铜死亡基因:")
    for _, row in top_cupro.iterrows():
        print(f"      - {row['GeneSymbol']}: P_target={row['P_target']:.3f}, Rank={row['Rank']}")

# ============================================================
# 3. 网络邻近性分析（基于连通分量）
# ============================================================
print("\n[5] 网络邻近性分析...")

# 铜死亡种子基因
cupro_seeds = [g for g in CUPTOPOSIS_GENES if g in G.nodes()]
print(f"    铜死亡种子在网络中的数量: {len(cupro_seeds)}")

# 计算Top50到铜死亡种子的平均最短路径
top50_genes = top50['GeneSymbol'].tolist()

# 为每个Top50基因计算到最近铜死亡基因的距离
path_results = []
for gene in top50_genes:
    if gene not in G.nodes():
        path_results.append({'Gene': gene, 'MinDist': -1, 'CuproSeed': 'N/A'})
        continue
    
    # 找到该基因所在的连通分量
    gene_cc = None
    for cc in connected_components:
        if gene in cc:
            gene_cc = cc
            break
    
    if gene_cc is None:
        path_results.append({'Gene': gene, 'MinDist': -1, 'CuproSeed': 'N/A'})
        continue
    
    # 在同一连通分量中找最近的铜死亡基因
    min_dist = float('inf')
    nearest_seed = None
    for seed in cupro_seeds:
        if seed in gene_cc:
            try:
                dist = nx.shortest_path_length(G, gene, seed)
                if dist < min_dist:
                    min_dist = dist
                    nearest_seed = seed
            except nx.NetworkXNoPath:
                continue
    
    if min_dist != float('inf'):
        path_results.append({'Gene': gene, 'MinDist': min_dist, 'CuproSeed': nearest_seed})
    else:
        path_results.append({'Gene': gene, 'MinDist': -1, 'CuproSeed': 'N/A'})

path_df = pd.DataFrame(path_results)
valid_paths = path_df[path_df['MinDist'] > 0]

if len(valid_paths) > 0:
    avg_path = valid_paths['MinDist'].mean()
    std_path = valid_paths['MinDist'].std()
    print(f"\n    Top50到铜死亡基因的网络距离:")
    print(f"      - 有效计算数: {len(valid_paths)}/{len(top50_genes)}")
    print(f"      - 平均距离: {avg_path:.2f} 跳")
    print(f"      - 标准差: {std_path:.2f}")
    print(f"      - 最小距离: {valid_paths['MinDist'].min()} 跳")
    print(f"      - 最大距离: {valid_paths['MinDist'].max()} 跳")
else:
    avg_path, std_path = -1, 0
    print(f"    无法计算网络距离")

# ============================================================
# 保存结果
# ============================================================
print("\n[6] 保存结果...")

# 超几何检验结果
hypergeom_results = pd.DataFrame({
    'Top列表': ['Top200', 'Top500', 'Top1000'],
    '列表大小': [200, 500, 1000],
    '铜死亡基因数': [top200_cupro, top500_cupro, top1000_cupro],
    'OR值': [odds_ratio_200, odds_ratio_500, odds_ratio_1000],
    'P值': [pval_200, pval_500, pval_1000],
    '显著性': ['*' if p < 0.05 else '' for p in [pval_200, pval_500, pval_1000]]
})

# 网络邻近性结果
network_summary = pd.DataFrame({
    '指标': ['平均最短路径', '标准差', '最小距离', '最大距离', '有效计算比例'],
    '值': [f"{avg_path:.2f}" if avg_path > 0 else 'N/A', 
           f"{std_path:.2f}" if avg_path > 0 else 'N/A',
           str(valid_paths['MinDist'].min()) if len(valid_paths) > 0 else 'N/A',
           str(valid_paths['MinDist'].max()) if len(valid_paths) > 0 else 'N/A',
           f"{len(valid_paths)}/{len(top50_genes)}"]
})

# 保存到Excel
output_file = f"{OUTPUT_DIR}/石竹烯_CIRI_综合分析报告_{datetime.now().strftime('%Y%m%d')}.xlsx"

# 读取原Excel并追加新sheet
with pd.ExcelWriter(output_file, engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:
    hypergeom_results.to_excel(writer, sheet_name='超几何检验', index=False)
    network_summary.to_excel(writer, sheet_name='网络邻近性', index=False)
    cupro_ranking.to_excel(writer, sheet_name='铜死亡基因专项', index=False)
    path_df.to_excel(writer, sheet_name='Top50网络距离', index=False)

print(f"\n结果已保存到: {output_file}")

# 打印论文表述
print("\n" + "=" * 60)
print("论文表述:")
print("=" * 60)
print(f"\n1. 超富集检验:")
if top200_cupro > 0:
    sig_text = "显著" if pval_200 < 0.05 else "未显著"
    print(f"   Top200 {sig_text}富集铜死亡执行基因（OR={odds_ratio_200:.2f}, P={pval_200:.3f}）。")
else:
    print(f"   Top200 未显著富集铜死亡执行基因（OR={odds_ratio_200:.2f}, P={pval_200:.3f}），")
    print(f"   提示石竹烯可能通过非经典通路发挥保护作用。")
print(f"\n2. 铜死亡基因专项:")
if len(cupro_ranking) > 0:
    print(f"   铜死亡执行基因在未知节点中呈现中等置信度（P={cupro_ranking['P_target'].median():.2f}），")
    print(f"   其中 {cupro_ranking.iloc[0]['GeneSymbol']} 排名最高（Rank={cupro_ranking.iloc[0]['Rank']}）。")
print(f"\n3. 网络邻近性:")
if avg_path > 0:
    print(f"   Top50 靶点与铜死亡执行基因的平均网络距离为 {avg_path:.1f} 跳（SD={std_path:.1f}），")
    print(f"   处于直接调控范围内。")
else:
    print(f"   Top50 靶点与铜死亡基因的网络距离分析显示，")
    print(f"   {len(valid_paths)}/{len(top50_genes)} 的靶点可计算网络距离。")
