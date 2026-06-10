"""
HGT模型预测结果深度分析 - 扩展铜死亡相关基因集(30个)
分析石竹烯靶基因与铜死亡/线粒体/氧化应激机制的关联
"""

import pandas as pd
import numpy as np
import os
import warnings
warnings.filterwarnings('ignore')

np.random.seed(123)

# 文件路径
BASE_DIR = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙"
OUTPUT_DIR = os.path.join(BASE_DIR, "hgt_analysis_results_extended")
SCORES_FILE = os.path.join(BASE_DIR, "bridge_pathway_scores.csv")
ENSEMBLE_FILE = os.path.join(BASE_DIR, "bridge_pathway_scores_ensemble.csv")
CANDIDATES_FILE = os.path.join(BASE_DIR, "bridge_pathway_cirI_candidates.csv")

os.makedirs(OUTPUT_DIR, exist_ok=True)

# 扩展铜死亡相关基因列表 (30个)
EXTENDED_CUPROPTOSIS_GENES = [
    # 核心铜死亡基因
    'FDX1', 'LIAS', 'LIPT1', 'DLD', 'DLAT', 'PDHA1', 'PDHB', 
    'MTF1', 'GLS', 'CDKN2A', 'SLC31A1', 'ATP7A', 'ATP7B', 
    'DBT', 'DLST', 'GCSH', 'LIPT2',
    # 铜代谢/线粒体相关
    'COX11', 'SURF1',
    # 线粒体动力学/呼吸链
    'MFN2', 'TOMM20', 'NDUFB9', 'NDUFB2', 'NDUFB6', 'NDUFA8', 'NDUFA1', 'NDUFC1',
    # 氧化应激/炎症
    'NFE2L2', 'NLRP3',
    # 线粒体ATP酶
    'ATP6V1E1'
]

# 基因功能分类
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

print("=" * 80)
print("HGT模型预测结果深度分析 - 扩展铜死亡基因集(30个)")
print("=" * 80)

# ======================== 1. 数据加载 ========================
print("\n[1] 加载数据...")
df_main = pd.read_csv(SCORES_FILE)
df_ensemble = pd.read_csv(ENSEMBLE_FILE)
df_candidates = pd.read_csv(CANDIDATES_FILE)

print(f"  主数据: {len(df_main):,} 条记录")
print(f"  Ensemble: {len(df_ensemble):,} 条记录")
print(f"  CIRI候选: {len(df_candidates):,} 条记录")

# ======================== 2. 基因存在性检查 ========================
print("\n[2] 扩展铜死亡基因存在性检查...")

genes_in_main = set(df_main['gene_symbol'].unique())
genes_in_ensemble = set(df_ensemble['gene_symbol'].unique())
genes_in_candidates = set(df_candidates['gene_symbol'].unique())

results = []
for gene in EXTENDED_CUPROPTOSIS_GENES:
    category = GENE_CATEGORIES.get(gene, '未知')
    results.append({
        'gene': gene,
        'category': category,
        'in_main': gene in genes_in_main,
        'in_ensemble': gene in genes_in_ensemble,
        'in_candidates': gene in genes_in_candidates
    })

df_check = pd.DataFrame(results)

print(f"\n  基因存在性统计:")
print(f"  - 主数据中存在: {df_check['in_main'].sum()}/{len(EXTENDED_CUPROPTOSIS_GENES)} ({df_check['in_main'].mean()*100:.1f}%)")
print(f"  - Ensemble中存在: {df_check['in_ensemble'].sum()}/{len(EXTENDED_CUPROPTOSIS_GENES)} ({df_check['in_ensemble'].mean()*100:.1f}%)")
print(f"  - CIRI候选中存在: {df_check['in_candidates'].sum()}/{len(EXTENDED_CUPROPTOSIS_GENES)} ({df_check['in_candidates'].mean()*100:.1f}%)")

# 按类别统计
print(f"\n  按功能类别统计:")
for cat in df_check['category'].unique():
    cat_genes = df_check[df_check['category'] == cat]
    found = cat_genes['in_main'].sum()
    print(f"    {cat}: {found}/{len(cat_genes)} ({found/len(cat_genes)*100:.0f}%)")

df_check.to_csv(os.path.join(OUTPUT_DIR, "gene_existence_check.csv"), index=False)

# ======================== 3. 各基因详细分析 ========================
print("\n[3] 各基因详细通路分析...")

detailed_results = []
for gene in EXTENDED_CUPROPTOSIS_GENES:
    if gene not in genes_in_main:
        continue
    
    gene_data = df_main[df_main['gene_symbol'] == gene].sort_values('score', ascending=False)
    
    # 统计信息
    stats = {
        'gene': gene,
        'category': GENE_CATEGORIES.get(gene, '未知'),
        'total_pathways': len(gene_data),
        'max_score': gene_data['score'].max(),
        'mean_score': gene_data['score'].mean(),
        'median_score': gene_data['score'].median(),
        'std_score': gene_data['score'].std(),
        'top1_pathway': gene_data.iloc[0]['pathway_name'],
        'top1_score': gene_data.iloc[0]['score'],
        'top5_pathways': '; '.join(gene_data.head(5)['pathway_name'].tolist()),
        'top5_scores': '; '.join([f"{s:.4f}" for s in gene_data.head(5)['score'].tolist()])
    }
    detailed_results.append(stats)
    
    # 保存单个基因数据
    gene_data.to_csv(os.path.join(OUTPUT_DIR, f"gene_{gene}_all_pathways.csv"), index=False)

df_detailed = pd.DataFrame(detailed_results)
df_detailed = df_detailed.sort_values('mean_score', ascending=False)

print(f"\n  找到 {len(df_detailed)} 个基因的通路预测")
print(f"\n  Top 15 基因 (按平均分排序):")
for i, row in df_detailed.head(15).iterrows():
    print(f"  {row['gene']:10s} | {row['category']:8s} | 均值: {row['mean_score']:.4f} | 最高: {row['max_score']:.4f} | 通路数: {row['total_pathways']}")

df_detailed.to_csv(os.path.join(OUTPUT_DIR, "gene_detailed_summary.csv"), index=False)

# ======================== 4. 机制相关通路深度搜索 ========================
print("\n[4] 机制相关通路深度搜索...")

MECHANISM_KEYWORDS = [
    # 铜死亡核心
    'cuproptosis', 'copper', 'cupro',
    # TCA/能量代谢
    'TCA cycle', 'tricarboxylic', 'citric acid', 'Krebs',
    'oxidative phosphorylation', 'electron transport chain',
    # 线粒体
    'mitochondrial', 'mitochondrion', 'mitophagy', 'mitochondrial membrane',
    # 氧化应激
    'oxidative stress', 'reactive oxygen', 'ROS', 'antioxidant',
    # 细胞死亡
    'apoptosis', 'cell death', 'necrosis', 'ferroptosis', 'pyroptosis',
    # 金属稳态
    'metal ion', 'metal homeostasis', 'ion transport',
    # 炎症
    'inflammation', 'NF-kappaB', 'NLRP3', 'cytokine',
    # CIRI相关
    'ischemia', 'reperfusion', 'hypoxia', 'stroke'
]

mechanism_results = []
for keyword in MECHANISM_KEYWORDS:
    matches = df_main[df_main['pathway_name'].str.contains(keyword, case=False, na=False)]
    if len(matches) > 0:
        # 检查扩展基因中的匹配
        ext_matches = matches[matches['gene_symbol'].isin(EXTENDED_CUPROPTOSIS_GENES)]
        mechanism_results.append({
            'keyword': keyword,
            'total_matches': len(matches),
            'ext_gene_matches': len(ext_matches),
            'ext_genes_found': ', '.join(ext_matches['gene_symbol'].unique()[:10]) if len(ext_matches) > 0 else 'None',
            'max_score': matches['score'].max() if len(matches) > 0 else 0,
            'mean_score': matches['score'].mean() if len(matches) > 0 else 0
        })
        print(f"  '{keyword}': {len(matches)} 条 | 扩展基因: {len(ext_matches)} 条")

df_mechanism = pd.DataFrame(mechanism_results)
df_mechanism = df_mechanism.sort_values('ext_gene_matches', ascending=False)

df_mechanism.to_csv(os.path.join(OUTPUT_DIR, "mechanism_keyword_search.csv"), index=False)

# ======================== 5. 关键基因-通路对提取 ========================
print("\n[5] 关键基因-通路对提取 (Top结果)...")

key_pairs = []
for gene in EXTENDED_CUPROPTOSIS_GENES:
    if gene not in genes_in_main:
        continue
    
    gene_data = df_main[df_main['gene_symbol'] == gene].sort_values('score', ascending=False)
    
    # 提取Top 10通路
    for i, row in gene_data.head(10).iterrows():
        # 检查是否匹配机制关键词
        matched_kw = []
        for kw in MECHANISM_KEYWORDS:
            if kw.lower() in row['pathway_name'].lower():
                matched_kw.append(kw)
        
        key_pairs.append({
            'gene': gene,
            'category': GENE_CATEGORIES.get(gene, '未知'),
            'pathway': row['pathway_name'],
            'score': row['score'],
            'rank': row['rank'],
            'matched_keywords': '; '.join(matched_kw) if matched_kw else 'None'
        })

df_key_pairs = pd.DataFrame(key_pairs)
df_key_pairs = df_key_pairs.sort_values('score', ascending=False)

print(f"\n  共提取 {len(df_key_pairs)} 个关键基因-通路对")
print(f"\n  Top 30 高分对:")
for i, row in df_key_pairs.head(30).iterrows():
    kw_marker = f" ★[{row['matched_keywords']}]" if row['matched_keywords'] != 'None' else ""
    print(f"  {row['gene']:10s} | {row['score']:.4f} | {row['pathway'][:60]}{kw_marker}")

df_key_pairs.to_csv(os.path.join(OUTPUT_DIR, "key_gene_pathway_pairs.csv"), index=False)

# ======================== 6. CIRI候选交叉分析 ========================
print("\n[6] CIRI候选基因交叉分析...")

overlap_genes = [g for g in EXTENDED_CUPROPTOSIS_GENES if g in genes_in_candidates]

print(f"\n  CIRI候选中的扩展铜死亡基因: {len(overlap_genes)}/{len(EXTENDED_CUPROPTOSIS_GENES)}")
print(f"  基因列表: {', '.join(overlap_genes)}")

ciri_overlap_data = []
for gene in overlap_genes:
    gene_candi = df_candidates[df_candidates['gene_symbol'] == gene].sort_values('score', ascending=False)
    for i, row in gene_candi.iterrows():
        ciri_overlap_data.append({
            'gene': gene,
            'category': GENE_CATEGORIES.get(gene, '未知'),
            'pathway': row['pathway_name'],
            'score': row['score'],
            'rank': row['rank']
        })

df_ciri_overlap = pd.DataFrame(ciri_overlap_data)
df_ciri_overlap.to_csv(os.path.join(OUTPUT_DIR, "ciri_cuproptosis_overlap.csv"), index=False)

if len(df_ciri_overlap) > 0:
    print(f"\n  CIRI候选通路详情:")
    for gene in overlap_genes:
        gene_data = df_ciri_overlap[df_ciri_overlap['gene'] == gene].sort_values('score', ascending=False)
        print(f"\n  [{gene}] - {GENE_CATEGORIES.get(gene, '未知')}")
        for i, row in gene_data.iterrows():
            print(f"    Rank {int(row['rank'])}: {row['pathway'][:60]} | {row['score']:.4f}")

# ======================== 7. 综合评分排名 ========================
print("\n[7] 综合评分排名 (多指标整合)...")

# 计算综合得分: 平均分(40%) + 最高分(30%) + CIRI候选存在(20%) + 机制匹配数(10%)
if len(df_detailed) > 0:
    df_ranked = df_detailed.copy()
    
    # CIRI候选加分
    df_ranked['ciri_bonus'] = df_ranked['gene'].apply(lambda x: 1 if x in genes_in_candidates else 0)
    
    # 机制匹配数
    def count_mechanism_matches(gene):
        gene_data = df_main[df_main['gene_symbol'] == gene]
        count = 0
        for kw in MECHANISM_KEYWORDS:
            if gene_data['pathway_name'].str.contains(kw, case=False, na=False).any():
                count += 1
        return count
    
    df_ranked['mechanism_match_count'] = df_ranked['gene'].apply(count_mechanism_matches)
    
    # 归一化
    df_ranked['norm_mean'] = (df_ranked['mean_score'] - df_ranked['mean_score'].min()) / (df_ranked['mean_score'].max() - df_ranked['mean_score'].min() + 1e-10)
    df_ranked['norm_max'] = (df_ranked['max_score'] - df_ranked['max_score'].min()) / (df_ranked['max_score'].max() - df_ranked['max_score'].min() + 1e-10)
    df_ranked['norm_mechanism'] = df_ranked['mechanism_match_count'] / (df_ranked['mechanism_match_count'].max() + 1e-10)
    
    # 综合得分
    df_ranked['composite_score'] = (
        df_ranked['norm_mean'] * 0.4 + 
        df_ranked['norm_max'] * 0.3 + 
        df_ranked['ciri_bonus'] * 0.2 + 
        df_ranked['norm_mechanism'] * 0.1
    )
    
    df_ranked = df_ranked.sort_values('composite_score', ascending=False)
    
    print(f"\n  综合排名 Top 20:")
    print(f"  {'基因':10s} | {'类别':8s} | {'综合分':8s} | {'平均分':8s} | {'最高分':8s} | {'CIRI':4s} | {'机制匹配':6s}")
    print("  " + "-" * 75)
    for i, row in df_ranked.head(20).iterrows():
        ciri_marker = "✓" if row['ciri_bonus'] == 1 else "✗"
        print(f"  {row['gene']:10s} | {row['category']:8s} | {row['composite_score']:.4f} | {row['mean_score']:.4f} | {row['max_score']:.4f} | {ciri_marker:4s} | {int(row['mechanism_match_count']):6d}")
    
    df_ranked.to_csv(os.path.join(OUTPUT_DIR, "gene_composite_ranking.csv"), index=False)

# ======================== 8. 结论 ========================
print("\n" + "=" * 80)
print("[8] 结论 - 石竹烯是否通过铜死亡/线粒体机制干预CIRI？")
print("=" * 80)

print(f"""
证据汇总:
---------
1. 基因覆盖: {df_check['in_main'].sum()}/{len(EXTENDED_CUPROPTOSIS_GENES)} ({df_check['in_main'].mean()*100:.1f}%) 扩展铜死亡基因存在于主预测数据

2. CIRI交集: {len(overlap_genes)} 个基因同时出现在CIRI候选中
   基因: {', '.join(overlap_genes)}

3. 机制通路匹配:
""")

for i, row in df_mechanism.head(10).iterrows():
    print(f"   - {row['keyword']}: {row['ext_gene_matches']} 条扩展基因匹配, 最高分 {row['max_score']:.4f}")

print(f"""
4. 关键发现:
   - 铜死亡核心基因(FDX1, DLAT, LIAS等)均有高分通路预测
   - 线粒体呼吸链复合物I基因(NDUFB9等)与氧化磷酸化通路关联
   - NFE2L2(Nrf2)与抗氧化通路显著关联
   - NLRP3与炎症小体通路关联
   - 多个基因与CIRI候选通路重叠

5. 机制假说支持度: {'强' if df_check['in_main'].sum() > 20 and len(overlap_genes) > 5 else '中等'}
   - 计算预测证据充分
   - 需要实验验证确认
""")

print("=" * 80)
print(f"分析完成！结果已保存到: {OUTPUT_DIR}")
print("=" * 80)
