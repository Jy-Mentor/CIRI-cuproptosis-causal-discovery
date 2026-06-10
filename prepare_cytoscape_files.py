import pandas as pd
import os

# 设置路径
work_dir = "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙"
result_dir = os.path.join(work_dir, "String_Network_Systematic_Analysis")
output_dir = os.path.join(work_dir, "Cytoscape_Import_Files")

# 创建输出目录
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# 读取String网络边列表（预处理后，已删除孤立节点）
edges_file = os.path.join(result_dir, "05_cytoscape_edge_list.txt")
edges = pd.read_csv(edges_file, sep="\t")

print(f"读取边列表: {len(edges)} 条边")
print(f"列名: {edges.columns.tolist()}")

# 读取节点属性（包含中心性指标和RRA排名）
node_attr_file = os.path.join(result_dir, "05_cytoscape_node_attributes.txt")
node_attr = pd.read_csv(node_attr_file, sep="\t")

print(f"\n读取节点属性: {len(node_attr)} 个节点")
print(f"列名: {node_attr.columns.tolist()}")

# 读取RRA结果
rra_file = os.path.join(result_dir, "04_rra_corrected_all_nodes.txt")
rra_df = pd.read_csv(rra_file, sep="\t")

print(f"\n读取RRA结果: {len(rra_df)} 个节点")

# ==================== 1. 准备边列表文件（SIF格式） ====================
print("\n========================================")
print("1. 生成边列表文件 (SIF格式)")
print("========================================")

# SIF格式: node1 interaction_type node2
sif_edges = edges.copy()
sif_edges['interaction'] = 'pp'  # protein-protein interaction
sif_output = sif_edges[['source', 'interaction', 'target']]

sif_file = os.path.join(output_dir, "cytoscape_network.sif")
sif_output.to_csv(sif_file, sep='\t', header=False, index=False)
print(f"SIF格式边列表已保存: {sif_file}")
print(f"  共 {len(sif_output)} 条边")

# ==================== 2. 准备边属性文件（带String分数） ====================
print("\n========================================")
print("2. 生成边属性文件")
print("========================================")

# 读取原始String交互文件获取分数
original_edges_file = "C:/Users/Jy-Mentor-7/Downloads/string_interactions (4).tsv"
original_edges = pd.read_csv(original_edges_file, sep='\t', comment='#', header=None)

# 设置列名
if original_edges.shape[1] >= 13:
    original_edges.columns = ["node1", "node2", "node1_string_id", "node2_string_id", 
                              "neighborhood_on_chromosome", "gene_fusion", "phylogenetic_cooccurrence",
                              "homology", "coexpression", "experimentally_determined_interaction",
                              "database_annotated", "automated_textmining", "combined_score"]

# 合并边信息
edge_attr = edges.copy()
edge_attr = edge_attr.merge(
    original_edges[['node1', 'node2', 'combined_score', 'coexpression', 
                    'experimentally_determined_interaction', 'database_annotated']],
    left_on=['source', 'target'],
    right_on=['node1', 'node2'],
    how='left'
)

# 保存边属性
edge_attr_file = os.path.join(output_dir, "cytoscape_edge_attributes.txt")
edge_attr_output = edge_attr[['source', 'target', 'combined_score', 'coexpression', 
                               'experimentally_determined_interaction', 'database_annotated']]
edge_attr_output.to_csv(edge_attr_file, sep='\t', index=False)
print(f"边属性文件已保存: {edge_attr_file}")

# ==================== 3. 准备节点属性文件（完整版） ====================
print("\n========================================")
print("3. 生成节点属性文件（完整版）")
print("========================================")

# 合并节点属性和RRA结果
node_full = node_attr.copy()
node_full = node_full.merge(rra_df[['Name', 'Rank', 'Score']], 
                            left_on='Node', right_on='Name', how='left')

# 重命名RRA列为更清晰的名称
node_full.rename(columns={'Rank': 'RRA_Rank', 'Score': 'RRA_Score'}, inplace=True)
if 'Name' in node_full.columns:
    node_full.drop('Name', axis=1, inplace=True)

# 保存完整节点属性
node_full_file = os.path.join(output_dir, "cytoscape_node_attributes_full.txt")
node_full.to_csv(node_full_file, sep='\t', index=False)
print(f"完整节点属性文件已保存: {node_full_file}")
print(f"  共 {len(node_full)} 个节点")
print(f"  包含列: {node_full.columns.tolist()}")

# ==================== 4. 准备高置信度边子集（用于MCODE） ====================
print("\n========================================")
print("4. 生成高置信度边子集（combined_score >= 400）")
print("========================================")

# 筛选高置信度边（combined_score >= 0.4，即String的400分）
high_conf_edges = edge_attr[edge_attr['combined_score'] >= 0.4].copy()

high_conf_sif = high_conf_edges[['source', 'target']]
high_conf_sif.insert(1, 'interaction', 'pp')

high_conf_file = os.path.join(output_dir, "cytoscape_high_confidence_network.sif")
high_conf_sif.to_csv(high_conf_file, sep='\t', header=False, index=False)
print(f"高置信度网络已保存: {high_conf_file}")
print(f"  共 {len(high_conf_sif)} 条边 (combined_score >= 400)")

# ==================== 5. 准备Hub基因子网络 ====================
print("\n========================================")
print("5. 生成Top 20 Hub基因子网络")
print("========================================")

# 读取Top 20 Hub基因
top20_file = os.path.join(result_dir, "04_rra_top50_hub_nodes.txt")
top20 = pd.read_csv(top20_file, sep='\t', nrows=20)
hub_genes = top20['Name'].tolist()

# 提取Hub基因之间的边
hub_edges = edges[
    (edges['source'].isin(hub_genes)) & (edges['target'].isin(hub_genes))
].copy()

if len(hub_edges) > 0:
    hub_sif = hub_edges[['source', 'target']]
    hub_sif.insert(1, 'interaction', 'pp')
    
    hub_file = os.path.join(output_dir, "cytoscape_top20_hub_subnetwork.sif")
    hub_sif.to_csv(hub_file, sep='\t', header=False, index=False)
    print(f"Top 20 Hub子网络已保存: {hub_file}")
    print(f"  共 {len(hub_sif)} 条边")
else:
    print("Top 20 Hub基因之间没有直接连接边")

# 提取Hub基因的一阶邻居（Hub基因 + 直接连接的节点）
hub_neighbors = set(hub_genes)
for gene in hub_genes:
    neighbors = set(edges[edges['source'] == gene]['target'].tolist() + 
                    edges[edges['target'] == gene]['source'].tolist())
    hub_neighbors.update(neighbors)

# 提取Hub基因及其邻居的子网络
hub_neighbor_edges = edges[
    (edges['source'].isin(hub_neighbors)) & (edges['target'].isin(hub_neighbors))
].copy()

hub_neighbor_sif = hub_neighbor_edges[['source', 'target']]
hub_neighbor_sif.insert(1, 'interaction', 'pp')

hub_neighbor_file = os.path.join(output_dir, "cytoscape_hub_neighbors_subnetwork.sif")
hub_neighbor_sif.to_csv(hub_neighbor_file, sep='\t', header=False, index=False)
print(f"Hub基因邻居子网络已保存: {hub_neighbor_file}")
print(f"  包含 {len(hub_neighbors)} 个节点")
print(f"  共 {len(hub_neighbor_sif)} 条边")

# ==================== 6. 生成铜死亡基因子网络 ====================
print("\n========================================")
print("6. 生成铜死亡基因子网络")
print("========================================")

cup_genes = ["FDX1", "LIAS", "SLC31A1", "DLAT", "PDHB", "PDHX", 
             "GPX4", "CP", "ATP7A", "ATOX1", "HIF1A", "NFKB1"]

# 检查哪些铜死亡基因在网络中
cup_in_network = [g for g in cup_genes if g in node_full['Node'].values]
print(f"网络中的铜死亡基因: {len(cup_in_network)}/{len(cup_genes)}")
print(f"  {cup_in_network}")

# 提取铜死亡基因之间的边
cup_edges = edges[
    (edges['source'].isin(cup_in_network)) & (edges['target'].isin(cup_in_network))
].copy()

if len(cup_edges) > 0:
    cup_sif = cup_edges[['source', 'target']]
    cup_sif.insert(1, 'interaction', 'pp')
    
    cup_file = os.path.join(output_dir, "cytoscape_copper_death_subnetwork.sif")
    cup_sif.to_csv(cup_file, sep='\t', header=False, index=False)
    print(f"铜死亡基因子网络已保存: {cup_file}")
    print(f"  共 {len(cup_sif)} 条边")
else:
    print("铜死亡基因之间没有直接连接边")

# 提取铜死亡基因及其邻居
cup_neighbors = set(cup_in_network)
for gene in cup_in_network:
    neighbors = set(edges[edges['source'] == gene]['target'].tolist() + 
                    edges[edges['target'] == gene]['source'].tolist())
    cup_neighbors.update(neighbors)

cup_neighbor_edges = edges[
    (edges['source'].isin(cup_neighbors)) & (edges['target'].isin(cup_neighbors))
].copy()

cup_neighbor_sif = cup_neighbor_edges[['source', 'target']]
cup_neighbor_sif.insert(1, 'interaction', 'pp')

cup_neighbor_file = os.path.join(output_dir, "cytoscape_copper_neighbors_subnetwork.sif")
cup_neighbor_sif.to_csv(cup_neighbor_file, sep='\t', header=False, index=False)
print(f"铜死亡基因邻居子网络已保存: {cup_neighbor_file}")
print(f"  包含 {len(cup_neighbors)} 个节点")
print(f"  共 {len(cup_neighbor_sif)} 条边")

# ==================== 7. 生成Cytoscape导入说明 ====================
print("\n========================================")
print("7. 生成Cytoscape导入说明")
print("========================================")

guide_content = """# Cytoscape网络分析导入指南

## 文件清单

### 1. 完整网络文件
- **cytoscape_network.sif** - 完整网络（135节点，1944边）
- **cytoscape_edge_attributes.txt** - 边属性（包含String combined_score等）
- **cytoscape_node_attributes_full.txt** - 节点属性（包含中心性指标和RRA排名）

### 2. 子网络文件（用于特定分析）
- **cytoscape_high_confidence_network.sif** - 高置信度网络（combined_score >= 400）
- **cytoscape_top20_hub_subnetwork.sif** - Top 20 Hub基因子网络
- **cytoscape_hub_neighbors_subnetwork.sif** - Hub基因及其邻居子网络
- **cytoscape_copper_death_subnetwork.sif** - 铜死亡基因子网络
- **cytoscape_copper_neighbors_subnetwork.sif** - 铜死亡基因邻居子网络

## Cytoscape导入步骤

### 导入完整网络
1. File → Import → Network from File
2. 选择: cytoscape_network.sif
3. 列映射:
   - Source Interaction: Column 1
   - Interaction Type: Column 2 (pp)
   - Target Interaction: Column 3

### 导入边属性
1. File → Import → Table from File
2. 选择: cytoscape_edge_attributes.txt
3. 关键列: 选择 combined_score, coexpression 等

### 导入节点属性
1. File → Import → Table from File
2. 选择: cytoscape_node_attributes_full.txt
3. 关键列: DC, BC, CC, EC, RRA_Rank, RRA_Score

## MCODE聚类分析
1. Apps → MCODE
2. 参数设置:
   - Network Scoring: Degree Cutoff = 2, K-Core = 2
   - Cluster Finding: Node Score Cutoff = 0.2, Haircut = TRUE, Fluff = FALSE
3. 点击 Analyze Cluster

## cytoHubba中心性分析
1. Apps → cytoHubba
2. 选择网络 → 选择节点
3. 选择计算方法:
   - SC (Subgraph Centrality)
   - NC (Network Centrality)
   - LAC (Local Average Connectivity)
   - MCC (Maximum Clique Centrality)
4. 点击 Calculate

## 可视化建议
1. 使用 RRA_Rank 设置节点大小（排名越小节点越大）
2. 使用 combined_score 设置边的粗细
3. 使用 MCODE聚类结果设置节点颜色
4. 标注铜死亡基因（如FDX1, NFKB1等）
"""

guide_file = os.path.join(output_dir, "Cytoscape_Import_Guide.txt")
with open(guide_file, 'w') as f:
    f.write(guide_content)
print(f"导入说明已保存: {guide_file}")

# ==================== 总结 ====================
print("\n========================================")
print("         文件准备完成")
print("========================================")
print(f"\n输出目录: {output_dir}")
print("\n生成文件列表:")
for f in sorted(os.listdir(output_dir)):
    filepath = os.path.join(output_dir, f)
    size = os.path.getsize(filepath)
    print(f"  - {f} ({size} bytes)")
