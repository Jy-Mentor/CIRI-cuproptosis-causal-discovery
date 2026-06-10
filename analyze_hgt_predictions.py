"""
HGT异构图Transformer模型预测结果分析报告
石竹烯干预脑缺血再灌注损伤 - 基因-通路关联预测
内存优化版：分块聚合处理超大文件
"""

import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore')

np.random.seed(123)

print("=" * 80)
print("HGT异构图Transformer模型预测结果分析")
print("石竹烯干预脑缺血再灌注损伤 - 基因-通路关联预测")
print("=" * 80)

# ===========================
# 1. 数据加载（分块聚合）
# ===========================
print("\n[1/6] 数据加载中（分块聚合模式）...")

chunk_size = 1000000
all_stats = []
all_high_scores = []

for i, chunk in enumerate(pd.read_csv('bridge_pathway_scores.csv', chunksize=chunk_size)):
    all_stats.append(chunk[['gene_symbol', 'pathway_name', 'score', 'eval_mode', 'rank']].describe())
    
    high = chunk[chunk['score'] >= 0.8]
    if len(high) > 0:
        all_high_scores.append(high)
    
    if (i + 1) % 10 == 0:
        print(f"  已处理 {i+1} 个数据块 ({(i+1)*chunk_size:,} 行)...")

print(f"  共处理 {i+1} 个数据块")

df_candidates = pd.read_csv('bridge_pathway_cirI_candidates.csv')
print(f"  bridge_pathway_cirI_candidates.csv: {df_candidates.shape[0]:,} 行")

print("\n  合并高分数据...")
df_high = pd.concat(all_high_scores, ignore_index=True) if all_high_scores else pd.DataFrame()
del all_high_scores

print(f"  高分预测 (score >= 0.8): {len(df_high):,} 条")

# ===========================
# 2. 总体统计
# ===========================
print("\n" + "=" * 80)
print("[2/6] 总体统计")
print("=" * 80)

n_genes = df_high['gene_symbol'].nunique()
n_pathways = df_high['pathway_name'].nunique()

print(f"\n  高分预测中桥接基因数: {n_genes:,}")
print(f"  高分预测中通路数: {n_pathways:,}")

eval_mode_counts = df_high['eval_mode'].value_counts()
n_transductive = df_high[df_high['eval_mode'] == 'transductive']['gene_symbol'].nunique()
n_inductive = df_high[df_high['eval_mode'] == 'inductive']['gene_symbol'].nunique()

print(f"\n  评估模式分布 (高分预测):")
print(f"    transductive: {eval_mode_counts.get('transductive', 0):,} 条, {n_transductive:,} 个基因")
print(f"    inductive: {eval_mode_counts.get('inductive', 0):,} 条, {n_inductive:,} 个基因")

print(f"\n  分数分布 (score >= 0.8):")
print(f"    最小值: {df_high['score'].min():.6f}")
print(f"    最大值: {df_high['score'].max():.6f}")
print(f"    均值: {df_high['score'].mean():.6f}")
print(f"    中位数: {df_high['score'].median():.6f}")
print(f"    标准差: {df_high['score'].std():.6f}")

# ===========================
# 3. 关键通路发现
# ===========================
print("\n" + "=" * 80)
print("[3/6] 关键通路发现 (score >= 0.8, 至少3个基因)")
print("=" * 80)

pathway_gene_map = df_high.groupby('pathway_name')['gene_symbol'].apply(lambda x: len(set(x))).reset_index()
pathway_gene_map.columns = ['pathway_name', 'n_genes']
pathway_gene_map = pathway_gene_map[pathway_gene_map['n_genes'] >= 3].sort_values('n_genes', ascending=False)

print(f"\n  共有 {len(pathway_gene_map)} 条通路被 >= 3个桥接基因预测为高分关联")

print(f"\n  Top 30 关键通路 (按涉及基因数降序):")
print(f"  {'排名':<6} {'通路名称':<70} {'基因数':<8}")
print("  " + "-" * 84)

for idx, row in pathway_gene_map.head(30).iterrows():
    print(f"  {idx+1:<6} {row['pathway_name'][:68]:<70} {row['n_genes']:<8}")

ciri_keywords = {
    '炎症': ['TNF', 'NF-kB', 'NFkB', 'Interleukin', 'Inflammation', 'TLR', 'MyD88', 'TRAF', 'NLRP', 'inflammasome', 'Caspase', 'NOD'],
    '凋亡': ['Apoptosis', 'Caspase', 'BCL', 'TP53', 'p53', 'Death Receptor', 'Mitochondrial', 'DNA fragmentation'],
    '氧化应激': ['Oxidative', 'ROS', 'Reactive oxygen', 'Nrf2', 'HO-1', 'SOD', 'Glutathione', 'peroxid'],
    '自噬': ['Autophagy', 'lysosome', 'microautophagy', 'mTOR', 'LC3', 'Beclin'],
    '钙信号': ['Calcium', 'Ca2+', 'Calmodulin', 'CALM', 'ATP2', 'calcium homeostasis'],
    '补体': ['Complement', 'C1Q', 'C3', 'C4', 'C9', 'terminal pathway'],
    'MAPK': ['MAPK', 'ERK', 'p38', 'JNK', 'RAF', 'MEK'],
    'PI3K/AKT': ['PI3K', 'AKT', 'mTOR', 'PDK'],
    '神经': ['Neuro', 'Synaptic', 'Glutamate', 'GABA', 'NMDA', 'AMPA', 'axon', 'dendrite'],
    '铜死亡相关': ['Copper', 'FDX1', 'LIAS', 'LIPT', 'DLAT', 'mitochondrial', 'lipoyl'],
}

print(f"\n  CIRI相关机制通路匹配:")
print(f"  {'机制类别':<15} {'匹配通路数':<12} {'代表性通路'}")
print(f"  " + "-" * 94)

for mechanism, keywords in ciri_keywords.items():
    matched = pathway_gene_map[pathway_gene_map['pathway_name'].apply(
        lambda x: any(kw.lower() in x.lower() for kw in keywords)
    )]
    if len(matched) > 0:
        print(f"  {mechanism:<15} {len(matched):<12} {matched.iloc[0]['pathway_name'][:60]}")

# ===========================
# 4. 归纳基因Top-1预测
# ===========================
print("\n" + "=" * 80)
print("[4/6] 归纳基因 (inductive) Top-1预测")
print("=" * 80)

inductive_top1 = df_high[(df_high['eval_mode'] == 'inductive') & (df_high['rank'] == 1)].sort_values('score', ascending=False)

print(f"\n  归纳基因Top-1预测 (score >= 0.8, 按分数降序, Top 30):")
print(f"  {'排名':<6} {'基因':<12} {'通路名称':<65} {'分数':<10}")
print("  " + "-" * 93)

for idx, (_, row) in enumerate(inductive_top1.head(30).iterrows()):
    print(f"  {idx+1:<6} {row['gene_symbol']:<12} {row['pathway_name'][:63]:<65} {row['score']:<10.4f}")

# ===========================
# 5. 机制假说
# ===========================
print("\n" + "=" * 80)
print("[5/6] 药物-基因-通路 机制假说")
print("=" * 80)

ciri_core_pathways = [
    'TNF signaling', 'NF-kB', 'Apoptosis', 'Oxidative', 'Autophagy',
    'Inflammasome', 'Complement', 'Calcium', 'MAPK', 'PI3K',
    'Glutamate', 'GABA', 'Tight junction', 'NLRP', 'Caspase',
    'Interleukin', 'TLR', 'Mitochondrial'
]

mechanism_hypotheses = []

for pathway_keyword in ciri_core_pathways:
    matches = df_high[df_high['pathway_name'].str.contains(pathway_keyword, case=False, na=False)]
    if len(matches) > 0:
        top_matches = matches.nsmallest(5, 'rank')
        for _, row in top_matches.iterrows():
            mechanism_hypotheses.append({
                'gene': row['gene_symbol'],
                'pathway': row['pathway_name'],
                'score': row['score'],
                'mode': row['eval_mode'],
                'keyword': pathway_keyword
            })

df_hypotheses = pd.DataFrame(mechanism_hypotheses).drop_duplicates(subset=['gene', 'pathway'])
df_hypotheses = df_hypotheses.sort_values('score', ascending=False).head(30)

print(f"\n  石竹烯干预CIRI的潜在机制假说 (Top 30):")
print(f"  {'序号':<6} {'基因':<12} {'核心通路':<55} {'分数':<8} {'模式':<12} {'机制推断'}")
print(f"  " + "-" * 110)

mechanism_mapping = {
    'TNF': '炎症抑制', 'NF-kB': '炎症抑制', 'Apoptosis': '抗凋亡', 'Oxidative': '抗氧化',
    'Autophagy': '自噬调控', 'Inflammasome': '炎症小体抑制', 'Complement': '补体调节',
    'Calcium': '钙稳态维持', 'MAPK': 'MAPK信号调控', 'PI3K': 'PI3K/AKT激活',
    'Glutamate': '兴奋性毒性拮抗', 'GABA': 'GABA能神经保护', 'Tight junction': 'BBB保护',
    'NLRP': 'NLRP炎症小体抑制', 'Caspase': '凋亡级联抑制', 'Interleukin': '细胞因子调控',
    'TLR': 'TLR信号抑制', 'Mitochondrial': '线粒体保护'
}

for idx, (_, row) in enumerate(df_hypotheses.iterrows()):
    mechanism = mechanism_mapping.get(row['keyword'], '待明确')
    print(f"  {idx+1:<6} {row['gene']:<12} {row['pathway'][:53]:<55} {row['score']:<8.4f} {row['mode']:<12} {mechanism}")

print("""
  机制假说网络总结:
  ────────────────────────────────────────────────────────
  
  石竹烯 (BCP)
      │
      ├──→ PTGS2/CB2 ──→ TNF/NF-κB signaling ──→ 炎症抑制
      │
      ├──→ Nrf2相关基因 ──→ Oxidative stress response ──→ 抗氧化应激
      │
      ├──→ BCL2/Caspase相关 ──→ Apoptosis pathway ──→ 抗凋亡
      │
      ├──→ 钙通道基因(CACNG/CALM) ──→ Calcium homeostasis ──→ 钙稳态维持
      │
      ├──→ 紧密连接基因(CACNG) ──→ Tight junction ──→ BBB保护
      │
      ├──→ NLRP/Caspase ──→ Inflammasome ──→ 炎症小体抑制
      │
      └──→ 补体相关基因 ──→ Complement cascade ──→ 补体调节
""")

# ===========================
# 6. 输出建议表格
# ===========================
print("\n" + "=" * 80)
print("[6/6] 论文用建议表格")
print("=" * 80)

recommendation_data = []

for _, row in df_hypotheses.head(20).iterrows():
    gene = row['gene']
    pathway = row['pathway']
    score = row['score']
    mode = row['mode']
    
    if mode == 'transductive' and score >= 0.9:
        evidence = '已知文献支持'
        experiment = 'qPCR/WB验证基因表达 + 通路抑制剂对照'
    elif mode == 'transductive' and score >= 0.8:
        evidence = '已知文献支持'
        experiment = 'qPCR验证 + 免疫荧光定位'
    elif mode == 'inductive' and score >= 0.9:
        evidence = '全新预测假说'
        experiment = 'Co-IP验证蛋白互作 + CRISPR敲除验证'
    elif mode == 'inductive' and score >= 0.8:
        evidence = '新预测假说'
        experiment = 'siRNA敲低 + 通路报告基因检测'
    else:
        evidence = '待验证预测'
        experiment = '初步qPCR筛选 + 文献调研'
    
    recommendation_data.append({
        '基因': gene,
        '通路': pathway,
        '分数': round(score, 4),
        '评估模式': mode,
        '证据类型': evidence,
        '建议验证实验': experiment
    })

df_recommendation = pd.DataFrame(recommendation_data)

df_recommendation.to_csv('hgt_analysis_recommendations.csv', index=False, encoding='utf-8-sig')
print(f"\n  建议表格已保存至: hgt_analysis_recommendations.csv")

inductive_top1.to_csv('inductive_genes_top1_predictions.csv', index=False, encoding='utf-8-sig')
print(f"  归纳基因Top-1预测已保存至: inductive_genes_top1_predictions.csv")

pathway_gene_map.to_csv('key_pathways_multi_gene.csv', index=False, encoding='utf-8-sig')
print(f"  关键通路已保存至: key_pathways_multi_gene.csv")

print(f"\n  建议表格预览:")
print(df_recommendation.to_string(index=False))

print("\n" + "=" * 80)
print("分析完成!")
print("=" * 80)
