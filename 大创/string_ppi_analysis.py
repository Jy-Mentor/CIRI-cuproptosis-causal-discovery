# -*- coding: utf-8 -*-
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

# 铜死亡基因列表
cuproptosis_genes = ['FDX1', 'SLC31A1', 'ATP7B', 'LIAS', 'DLAT', 'ATOX1', 'ATP7A', 'CDKN2A', 'GLS', 'CD274', 'FECH', 'MTF1', 'PDHA1', 'SLC31A2', 'DLST']

# 读取和处理数据
def load_and_filter_ppi(file_path):
    # 读取TSV文件
    df = pd.read_csv(file_path, sep='\t', header=None)
    # 重命名列
    df.columns = ['node1', 'node2', 'col2', 'col3', 'col4', 'col5', 'col6', 'col7', 'col8', 'col9', 'col10', 'col11', 'combined_score']
    # 将combined_score列转换为浮点数
    df['combined_score'] = pd.to_numeric(df['combined_score'], errors='coerce')
    # 过滤combined_score ≥ 0.4的边
    filtered_df = df[df['combined_score'] >= 0.4]
    # 构建边列表
    edges = [(row['node1'], row['node2'], row['combined_score']) for _, row in filtered_df.iterrows()]
    return edges

# 计算网络拓扑参数
def calculate_network_params(G):
    # 有效互作边数
    num_edges = G.number_of_edges()
    # 平均节点度
    avg_degree = 2 * num_edges / G.number_of_nodes() if G.number_of_nodes() > 0 else 0
    # 网络密度
    density = nx.density(G)
    # 最大连通分量占比
    if G.number_of_nodes() > 0:
        largest_cc = max(nx.connected_components(G), key=len)
        lcc_ratio = len(largest_cc) / G.number_of_nodes() * 100
    else:
        lcc_ratio = 0
    return num_edges, avg_degree, density, lcc_ratio

# 计算铜死亡基因相关参数
def calculate_cuproptosis_params(G, cuproptosis_genes):
    # 计算铜死亡基因的连接度
    cuproptosis_degrees = []
    for gene in cuproptosis_genes:
        if gene in G:
            cuproptosis_degrees.append(G.degree(gene))
    # 平均连接度
    avg_cuproptosis_degree = np.mean(cuproptosis_degrees) if cuproptosis_degrees else 0
    # 找出Top 3 Hub基因
    gene_degrees = {gene: G.degree(gene) for gene in cuproptosis_genes if gene in G}
    top_3 = sorted(gene_degrees.items(), key=lambda x: x[1], reverse=True)[:3]
    return avg_cuproptosis_degree, top_3

# 可视化网络
def visualize_networks(G1, G2, cuproptosis_genes, output_file):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))
    
    # 网络A可视化
    pos1 = nx.spring_layout(G1, k=0.3, iterations=100)
    # 节点颜色：铜死亡基因红色，其他蓝色
    node_colors1 = ['red' if node in cuproptosis_genes else 'blue' for node in G1.nodes()]
    # 边的权重作为线宽
    edge_widths1 = [G1[u][v]['weight'] * 3 for u, v in G1.edges()]
    nx.draw(G1, pos1, ax=ax1, node_color=node_colors1, node_size=100, edge_color='gray', width=edge_widths1, with_labels=False)
    ax1.set_title('网络A (80基因)')
    
    # 网络B可视化
    pos2 = nx.spring_layout(G2, k=0.3, iterations=100)
    # 节点颜色和大小：铜死亡基因红色大节点，其他蓝色
    node_colors2 = ['red' if node in cuproptosis_genes else 'blue' for node in G2.nodes()]
    node_sizes2 = [300 if node in cuproptosis_genes else 100 for node in G2.nodes()]
    # 边的权重作为线宽
    edge_widths2 = [G2[u][v]['weight'] * 3 for u, v in G2.edges()]
    nx.draw(G2, pos2, ax=ax2, node_color=node_colors2, node_size=node_sizes2, edge_color='gray', width=edge_widths2, with_labels=False)
    # 标注Top 3 Hub基因
    _, top_3 = calculate_cuproptosis_params(G2, cuproptosis_genes)
    for gene, degree in top_3:
        if gene in pos2:
            ax2.text(pos2[gene][0], pos2[gene][1], '{0}\n(Degree: {1})'.format(gene, degree), fontsize=10, ha='center', bbox=dict(facecolor='white', alpha=0.7))
    ax2.set_title('网络B (95基因)')
    
    # 添加图例
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', label='铜死亡基因', markerfacecolor='red', markersize=10),
        Line2D([0], [0], marker='o', color='w', label='其他基因', markerfacecolor='blue', markersize=10),
        Line2D([0], [0], color='gray', lw=1, label='Confidence 0.4'),
        Line2D([0], [0], color='gray', lw=3, label='Confidence 1.0')
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=4)
    
    # 保存图像
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

# 主函数
def main():
    # 文件路径
    file_a = r'C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\大创\80 PPI.tsv'
    file_b = r'C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\大创\80+15PPI.tsv'
    output_image = r'C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\大创\ppi_network_visualization.png'
    
    # 加载数据
    edges_a = load_and_filter_ppi(file_a)
    edges_b = load_and_filter_ppi(file_b)
    
    # 构建网络
    G1 = nx.Graph()
    G1.add_weighted_edges_from(edges_a)
    G2 = nx.Graph()
    G2.add_weighted_edges_from(edges_b)
    
    # 计算网络A参数
    edges_a_count, avg_degree_a, density_a, lcc_ratio_a = calculate_network_params(G1)
    
    # 计算网络B参数
    edges_b_count, avg_degree_b, density_b, lcc_ratio_b = calculate_network_params(G2)
    
    # 计算网络B铜死亡基因参数
    avg_cuproptosis_degree, top_3_cuproptosis = calculate_cuproptosis_params(G2, cuproptosis_genes)
    
    # 生成可视化
    visualize_networks(G1, G2, cuproptosis_genes, output_image)
    
    # 处理top_3_cuproptosis为空的情况
    if len(top_3_cuproptosis) >= 3:
        top1_gene, top1_degree = top_3_cuproptosis[0]
        top2_gene, top2_degree = top_3_cuproptosis[1]
        top3_gene, top3_degree = top_3_cuproptosis[2]
    else:
        top1_gene, top1_degree = "N/A", "N/A"
        top2_gene, top2_degree = "N/A", "N/A"
        top3_gene, top3_degree = "N/A", "N/A"
    
    # 生成报告
    report = """拓扑参数分析报告

网络A (80基因)：
• 有效互作边数【{0}】
• 平均节点度【{1:.2f}】
• 网络密度【{2:.4f}】
• 最大连通分量占比（LCC）【{3:.1f}】%

网络B (95基因)：
• 有效互作边数【{4}】
• 平均节点度【{5:.2f}】
• 网络密度【{6:.4f}】
• 最大连通分量占比（LCC）【{7:.1f}】%
• 铜死亡基因平均连接度【{8:.2f}】
• 铜死亡基因Degree Top 3：【{9}】(Degree: 【{10}】)、【{11}】(Degree: 【{12}】)、【{13}】(Degree: 【{14}】)

可视化已保存至：{15}
""".format(
    edges_a_count, avg_degree_a, density_a, lcc_ratio_a,
    edges_b_count, avg_degree_b, density_b, lcc_ratio_b,
    avg_cuproptosis_degree,
    top1_gene, top1_degree,
    top2_gene, top2_degree,
    top3_gene, top3_degree,
    output_image
)
    
    # 保存报告
    with open(r'C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\大创\ppi_analysis_report.txt', 'w', encoding='utf-8') as f:
        f.write(report)
    
    print("分析完成！")
    print(report)

if __name__ == "__main__":
    main()