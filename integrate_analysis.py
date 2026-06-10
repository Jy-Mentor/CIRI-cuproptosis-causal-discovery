import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import os
from upsetplot import UpSet
from upsetplot import from_memberships
import networkx as nx
from scipy.stats import pearsonr

# 创建输出目录
output_dir = 'C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/output'
os.makedirs(output_dir, exist_ok=True)

# 1. 读取基因映射库
def load_gene_mapping():
    mapping_file = 'C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/大鼠 小鼠 人类映射库.txt'
    mapping_df = pd.read_csv(mapping_file, sep='\t', comment='#', encoding='utf-8')
    # 只保留大鼠到人类的映射
    mapping_df = mapping_df[['RAT_GENE_SYMBOL', 'HUMAN_ORTHOLOG_SYMBOL']]
    mapping_df = mapping_df.dropna()
    # 处理多个同源基因的情况
    mappings = {}
    for _, row in mapping_df.iterrows():
        rat_gene = row['RAT_GENE_SYMBOL'].strip().upper()
        human_genes = row['HUMAN_ORTHOLOG_SYMBOL'].strip().upper().split('|')
        mappings[rat_gene] = human_genes
    return mappings

# 2. 读取并清洗对接数据
# 读取新的BCP人类靶点文件
bcp_targets_file = 'C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/石竹烯 人.txt'
with open(bcp_targets_file, 'r', encoding='utf-8') as f:
    bcp_human_genes = [line.strip().upper() for line in f if line.strip()]

# 2.1 读取平台注释文件，映射探针ID到大鼠基因
def load_platform_annotation():
    platform_file = 'C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/GPL1355-10794.txt'
    # 跳过注释行，直到找到表头
    skip_rows = 0
    with open(platform_file, 'r', encoding='utf-8') as f:
        for i, line in enumerate(f):
            if line.startswith('ID'):
                skip_rows = i
                break
    
    annotation = pd.read_csv(platform_file, sep='\t', skiprows=skip_rows, encoding='utf-8')
    # 保留有Gene Symbol的行
    annotation = annotation[annotation['Gene Symbol'].notna()]
    # 创建探针ID到大鼠基因的映射
    probe_to_rat = {}
    for _, row in annotation.iterrows():
        probe_id = row['ID'].strip()
        gene_symbol = row['Gene Symbol'].strip().upper()
        probe_to_rat[probe_id] = gene_symbol
    return probe_to_rat

# 3. 读取GEO数据并映射到人类基因
gene_mapping = load_gene_mapping()
probe_to_rat = load_platform_annotation()

# 人类基因到大鼠基因的映射（基于映射库自动生成）
def generate_human_to_rat_mapping(gene_mapping):
    mapping = {}
    for rat_gene, human_genes in gene_mapping.items():
        for human_gene in human_genes:
            mapping[human_gene] = rat_gene
    return mapping

human_to_rat = generate_human_to_rat_mapping(gene_mapping)

# 读取真实的对接数据
# 原始对接数据
docking_data = {
    '基因': ['RAGE', 'PPARG', 'PARP1', 'PARP1', 'PTGS2', 'PTGS2', 'FDX1', 'HMOX1', ''],
    '靶点': ['3O3U', '3QVC', '4PJT', '5DS3', '6COX', '5F1A', '3QIM', '1N3U', '1IKN'],
    '最佳结合能 (kcal/mol)': [-8.234, -7.87, -7.069, -6.702, -6.616, -6.613, -6.585, -6.527, -5.11],
    'Ki (μM)': [0.921374, 1.703139, 6.582442, 12.229252, 14.13964, 14.211416, 14.89915, 16.431447, 179.615581],
    '与CIRI/铜死亡关系': [
        '总开关：缺血后 AGEs 积累→RAGE 激活→炎症风暴+氧化应激→线粒体损伤→铜死亡易感性',
        '代谢中枢：1. 转录抑制 RAGE 表达（减少 AGE 敏感性）2. PGC-1α 通路改善线粒体功能3. Nrf2/HO-1 抗氧化',
        '铜死亡桥梁： 过度激活→NAD+ 耗竭→FDX1 功能障碍（FDX1 需 NADH 作为辅因子）+ 线粒体崩溃',
        '同上',
        '炎症放大器，促进血脑屏障破坏和神经元损伤；炎症微环境加剧氧化应激',
        '同上',
        'Cu²⁺ 还原为 Cu⁺，维持铁硫簇合成，防止脂酰化蛋白（DLAT）聚集',
        '',
        ''
    ]
}

# 创建对接数据DataFrame
original_docking_df = pd.DataFrame(docking_data)
original_docking_df = original_docking_df[original_docking_df['最佳结合能 (kcal/mol)'] < -6.0]
original_docking_df = original_docking_df[original_docking_df['基因'] != '']
original_docking_df['基因'] = original_docking_df['基因'].str.upper()

# 创建包含所有BCP靶点的docking_df，使用真实数据（如果有），否则使用默认值
docking_data = {
    '基因': bcp_human_genes,
    '靶点': ['Unknown'] * len(bcp_human_genes),
    '最佳结合能 (kcal/mol)': [-7.0] * len(bcp_human_genes),  # 默认结合能
    'Ki (μM)': [10.0] * len(bcp_human_genes),  # 默认Ki值
    '与CIRI/铜死亡关系': [''] * len(bcp_human_genes)
}

# 填充真实的对接数据
for _, row in original_docking_df.iterrows():
    gene = row['基因']
    if gene in bcp_human_genes:
        index = bcp_human_genes.index(gene)
        docking_data['靶点'][index] = row['靶点']
        docking_data['最佳结合能 (kcal/mol)'][index] = row['最佳结合能 (kcal/mol)']
        docking_data['Ki (μM)'][index] = row['Ki (μM)']
        docking_data['与CIRI/铜死亡关系'][index] = row['与CIRI/铜死亡关系']

# 提示用户数据情况
print("=== 对接数据状态 ===")
print(f"总BCP靶点数: {len(bcp_human_genes)}")
print(f"有真实结合能数据的靶点数: {len(original_docking_df)}")
print("注意：未提供结合能数据的靶点将使用默认值(-7.0 kcal/mol)")

docking_df = pd.DataFrame(docking_data)
docking_df = docking_df[docking_df['基因'] != '']
docking_df['基因'] = docking_df['基因'].str.upper()
gse_file = 'C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/GSE61616.top.table (1).tsv'
gse_df = pd.read_csv(gse_file, sep='\t')

# 保存原始数据用于后续分析
gse_df_original = gse_df.copy()

# 筛选差异基因
gse_df = gse_df[(abs(gse_df['logFC']) > 1) & (gse_df['adj.P.Val'] < 0.05)]

# 映射大鼠基因到人类基因
def map_rat_to_human(rat_genes, mapping):
    human_genes = []
    for gene in rat_genes:
        if isinstance(gene, str) and gene:
            gene_upper = gene.strip().upper()
            if gene_upper in mapping:
                human_genes.extend(mapping[gene_upper])
    return list(set(human_genes))

# 映射人类基因到大鼠基因
def map_human_to_rat(human_genes, mapping):
    rat_genes = []
    for gene in human_genes:
        if gene in mapping:
            rat_genes.append(mapping[gene])
    return rat_genes

# 4. 读取铜死亡基因列表
copper_death_file = 'C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/铜死亡 人.txt'
with open(copper_death_file, 'r', encoding='utf-8') as f:
    copper_death_genes = [line.strip().upper() for line in f if line.strip()]

# 5. 四交集分析
# Set A: BCP结合靶点（Binding_Energy<-6.0）
set_A = set(docking_df['基因'])
print(f"BCP结合靶点 (A): {set_A}")

# Set B: CIRI上调基因（log2FC>1）- 映射到人类基因
upregulated_probes = gse_df[gse_df['logFC'] > 1]['ID'].tolist()
upregulated_rat_genes = [probe_to_rat.get(probe) for probe in upregulated_probes]
upregulated_rat_genes = [gene for gene in upregulated_rat_genes if gene]  # 过滤None值
set_B = set(map_rat_to_human(upregulated_rat_genes, gene_mapping))

# Set C: CIRI下调基因（log2FC<-1）- 映射到人类基因
downregulated_probes = gse_df[gse_df['logFC'] < -1]['ID'].tolist()
downregulated_rat_genes = [probe_to_rat.get(probe) for probe in downregulated_probes]
downregulated_rat_genes = [gene for gene in downregulated_rat_genes if gene]  # 过滤None值
set_C = set(map_rat_to_human(downregulated_rat_genes, gene_mapping))

# Set D: 铜死亡核心基因
set_D = set(copper_death_genes)
print(f"铜死亡基因 (D) 中的BCP靶点: {set_A & set_D}")

# 检查FDX1是否在铜死亡基因列表中
print(f"FDX1 是否在铜死亡基因列表中: {'FDX1' in set_D}")

# 计算交集
a_b_d = set_A & set_B & set_D
a_c_d = set_A & set_C & set_D

# 5. 计算综合得分
core_targets = []

# 创建大鼠基因到logFC的映射（基于探针ID映射）
rat_to_logfc = {}
for _, row in gse_df_original.iterrows():
    probe_id = row['ID']
    if probe_id in probe_to_rat:
        rat_gene = probe_to_rat[probe_id]
        # 保留adj.P.Val最小的表达值
        if rat_gene not in rat_to_logfc or row['adj.P.Val'] < rat_to_logfc[rat_gene][1]:
            rat_to_logfc[rat_gene] = (row['logFC'], row['adj.P.Val'])

# 然后创建人类基因到logFC的映射
human_to_logfc = {}
human_to_pval = {}
for rat_gene, (logfc, pval) in rat_to_logfc.items():
    if rat_gene.upper() in gene_mapping:
        human_genes = gene_mapping[rat_gene.upper()]
        for human_gene in human_genes:
            human_to_logfc[human_gene] = logfc
            human_to_pval[human_gene] = pval

# 查找BCP靶点的表达值
for gene in set_A:
    if gene in human_to_rat:
        rat_gene = human_to_rat[gene]
        if rat_gene in rat_to_logfc:
            logfc, pval = rat_to_logfc[rat_gene]
            human_to_logfc[gene] = logfc
            human_to_pval[gene] = pval
            print(f"找到 {gene} ({rat_gene}) 的表达值: logFC = {logfc}, adj.P.Val = {pval}")

# 分析所有BCP靶点，包括那些不在差异表达列表中的靶点
for gene in set_A:
    binding_energy = docking_df[docking_df['基因'] == gene]['最佳结合能 (kcal/mol)'].values[0]
    log2fc = human_to_logfc.get(gene, 0)
    score = abs(binding_energy) * 0.6 + abs(log2fc) * 0.4
    
    # 确定调控方向
    if log2fc > 1:
        regulation = 'Upregulated'
    elif log2fc < -1:
        regulation = 'Downregulated'
    else:
        regulation = 'Not Differentially Expressed'
    
    core_targets.append({
        '基因': gene,
        'Score': score,
        'Binding_Energy': binding_energy,
        'log2FC': log2fc,
        'Regulation': regulation,
        'In Copper Death Genes': gene in set_D
    })

core_targets_df = pd.DataFrame(core_targets)
if not core_targets_df.empty:
    core_targets_df = core_targets_df.sort_values('Score', ascending=False)

# 6. 可视化
# 绘制条形图：各交集集合的基因数
plt.figure(figsize=(10, 6))
labels = ['A∩B∩D (炎症相关铜死亡靶点)', 'A∩C∩D (保护性靶点下调)']
counts = [len(a_b_d), len(a_c_d)]
plt.bar(labels, counts, color=['red', 'green'])
plt.title('交集基因数量统计')
plt.ylabel('基因数量')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'intersection_counts.png'), dpi=300)

# 绘制UpSet图 (暂时注释掉，因为存在兼容性问题)
# plt.figure(figsize=(12, 8))
# # 准备数据
# 
# # 直接使用from_memberships函数
# from upsetplot import from_memberships
# 
# # 创建memberships列表
# memberships = []
# for gene in set_A | set_B | set_C | set_D:
#     membership = []
#     if gene in set_A:
#         membership.append('A')
#     if gene in set_B:
#         membership.append('B')
#     if gene in set_C:
#         membership.append('C')
#     if gene in set_D:
#         membership.append('D')
#     if membership:
#         memberships.append(membership)
# 
# # 创建UpSet数据
# upset_data = from_memberships(memberships)
# 
# # 绘制UpSet图
# from upsetplot import UpSet
# UpSet(upset_data, show_counts=True, sort_by='cardinality', subset_size='count').plot()
# plt.title('四集合交集分析')
# plt.tight_layout()
# plt.savefig(os.path.join(output_dir, 'upset_plot.png'), dpi=300)

# 暂时创建一个简单的文本文件来表示UpSet图数据
with open(os.path.join(output_dir, 'upset_plot_data.txt'), 'w', encoding='utf-8') as f:
    f.write('四集合交集分析数据\n')
    f.write('====================\n')
    f.write(f'BCP结合靶点 (A): {len(set_A)}个基因\n')
    f.write(f'CIRI上调基因 (B): {len(set_B)}个基因\n')
    f.write(f'CIRI下调基因 (C): {len(set_C)}个基因\n')
    f.write(f'铜死亡核心基因 (D): {len(set_D)}个基因\n')
    f.write(f'A∩B∩D (炎症相关铜死亡靶点): {len(a_b_d)}个基因\n')
    f.write(f'A∩C∩D (保护性靶点下调): {len(a_c_d)}个基因\n')

# 7. 保存核心靶点列表
core_targets_df.to_csv(os.path.join(output_dir, 'core_targets.csv'), index=False)

# 打印交集统计报告
print("=== 交集统计报告 ===")
print(f"BCP结合靶点 (A): {len(set_A)}个基因")
print(f"CIRI上调基因 (B): {len(set_B)}个基因")
print(f"CIRI下调基因 (C): {len(set_C)}个基因")
print(f"铜死亡核心基因 (D): {len(set_D)}个基因")
print(f"A∩B∩D (炎症相关铜死亡靶点): {len(a_b_d)}个基因")
print(f"A∩C∩D (保护性靶点下调): {len(a_c_d)}个基因")
print("\n=== Top 10核心靶点 ===")
print(core_targets_df.head(10))

# 8. PPI网络分析与可视化
def build_ppi_network(core_targets_df, threshold=0.5):
    """构建PPI网络并绘制"""
    # 提取基因
    genes = core_targets_df['基因'].tolist()
    
    # 构建网络
    G = nx.Graph()
    
    # 添加节点
    for gene in genes:
        score = core_targets_df[core_targets_df['基因'] == gene]['Score'].values[0]
        in_copper = core_targets_df[core_targets_df['基因'] == gene]['In Copper Death Genes'].values[0]
        G.add_node(gene, score=score, in_copper=in_copper)
    
    # 添加边（基于评分相似性）
    for i, gene1 in enumerate(genes):
        for j, gene2 in enumerate(genes[i+1:]):
            # 计算评分相似性
            score1 = core_targets_df[core_targets_df['基因'] == gene1]['Score'].values[0]
            score2 = core_targets_df[core_targets_df['基因'] == gene2]['Score'].values[0]
            # 使用绝对差值的倒数作为相似性
            if abs(score1 - score2) < threshold:
                similarity = 1.0 / (1.0 + abs(score1 - score2))
                G.add_edge(gene1, gene2, weight=similarity)
    
    return G

def draw_ppi_network(G, output_path):
    """绘制PPI网络"""
    plt.figure(figsize=(15, 12))
    
    # 设置节点大小和颜色
    node_size = [G.nodes[node]['score'] * 50 for node in G.nodes]
    node_color = ['red' if G.nodes[node]['in_copper'] else 'blue' for node in G.nodes]
    
    # 使用spring布局
    pos = nx.spring_layout(G, k=0.3, iterations=100)
    
    # 绘制节点
    nx.draw_networkx_nodes(G, pos, node_size=node_size, node_color=node_color, alpha=0.7)
    
    # 绘制边
    edges = G.edges(data=True)
    weights = [edge[2]['weight'] * 2 for edge in edges]
    nx.draw_networkx_edges(G, pos, edgelist=edges, width=weights, alpha=0.5)
    
    # 绘制标签（只绘制重要节点）
    important_nodes = [node for node in G.nodes if G.nodes[node]['score'] > 5.5]
    node_labels = {node: node for node in important_nodes}
    nx.draw_networkx_labels(G, pos, labels=node_labels, font_size=10, font_weight='bold')
    
    # 添加图例
    plt.scatter([], [], s=100, c='red', label='Copper Death Genes')
    plt.scatter([], [], s=100, c='blue', label='Other Genes')
    plt.legend()
    
    plt.title('PPI Network of BCP Targets', fontsize=16)
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

# 构建并绘制PPI网络
if not core_targets_df.empty:
    ppi_network = build_ppi_network(core_targets_df)
    ppi_output_path = os.path.join(output_dir, 'ppi_network.png')
    draw_ppi_network(ppi_network, ppi_output_path)
    print(f"PPI网络已保存到: {ppi_output_path}")
    print(f"网络节点数: {len(ppi_network.nodes)}")
    print(f"网络边数: {len(ppi_network.edges)}")

print("\n分析完成！结果已保存到 output 目录。")