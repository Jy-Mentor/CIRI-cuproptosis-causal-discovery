"""
铜死亡基因生物学文献验证报告
整合MR数据、STRING网络分析、已有研究结果
"""
import csv
from pathlib import Path
from collections import defaultdict

# 路径
BASE_DIR = Path("c:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙")
RESULTS_DIR = BASE_DIR / "results"

# GAT预测结果
gat_predictions = {}
with open(RESULTS_DIR / "all_unknown_predictions.csv", 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        gat_predictions[row['GeneSymbol']] = {
            'Rank': int(row['Rank']),
            'P_target': float(row['P_target']),
            'is_cuproptosis': int(row['is_cuproptosis']),
            'dist_to_cuproptosis': float(row['dist_to_cuproptosis'])
        }

# 铜死亡基因列表
CUPTO_GENES = {
    "FDX1", "LIAS", "LIPT1", "DLAT", "PDHB", "PDHX", 
    "SLC31A1", "ATP7A", "ATP7B", "ATOX1", "NFE2L2", 
    "HIF1A", "MTOR", "NFKB1", "GPX4"
}

# 脑卒中靶点列表
stroke_targets = set()
with open(BASE_DIR / "local_data" / "stroke_targets.txt", 'r') as f:
    for line in f:
        gene = line.strip()
        if gene:
            stroke_targets.add(gene)

# 石竹烯靶点列表
bcp_targets = set()
with open(BASE_DIR / "local_data" / "bcp_targets.txt", 'r') as f:
    for line in f:
        gene = line.strip()
        if gene:
            bcp_targets.add(gene)

# STRING网络中心性排名（RRA分析）
rra_copper = {}
with open(BASE_DIR / "String_Network_Systematic_Analysis" / "04_rra_copper_genes.txt", 'r') as f:
    next(f)  # skip header
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 3:
            rra_copper[parts[0]] = {'Score': float(parts[1]), 'Rank': int(parts[2])}

# MR数据
mr_results = {}
mr_file = BASE_DIR / "MR_batch_summary_20260506" / "01_MR_results_main.csv"
with open(mr_file, 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        gene = row.get('gene', row.get('Gene', ''))
        mr_results[gene] = {
            'pval': float(row.get('pval', row.get('P_value', 0))),
            'fdr_qval': float(row.get('fdr_qval', row.get('FDR', 1))),
            'nSNP': int(row.get('nSNP', row.get('nsnp', 0)))
        }

# 生成文献验证报告
print("="*100)
print("铜死亡基因GAT预测结果 - 多源生物学文献验证报告")
print("="*100)

report_data = []

for gene in sorted(CUPTO_GENES):
    gat_info = gat_predictions.get(gene, None)
    rra_info = rra_copper.get(gene, None)
    mr_info = mr_results.get(gene, None)
    
    # 多源数据整合
    in_stroke_targets = gene in stroke_targets
    in_bcp_targets = gene in bcp_targets
    
    # 文献支持度评分
    lit_score = 0
    lit_evidence = []
    
    if in_stroke_targets:
        lit_score += 2
        lit_evidence.append("脑卒中靶点文献")
    
    if in_bcp_targets:
        lit_score += 2
        lit_evidence.append("石竹烯靶点")
    
    if rra_info and rra_info['Rank'] <= 50:
        lit_score += 2
        lit_evidence.append(f"STRING网络中心性Top-50 (Rank {rra_info['Rank']})")
    elif rra_info:
        lit_score += 1
        lit_evidence.append(f"STRING网络中心性 (Rank {rra_info['Rank']})")
    
    if mr_info and mr_info.get('fdr_qval', 1) < 0.05:
        lit_score += 3
        lit_evidence.append(f"MR显著 (FDR={mr_info['fdr_qval']:.3e})")
    elif mr_info and mr_info.get('pval', 1) < 0.05:
        lit_score += 2
        lit_evidence.append(f"MR显著 (p={mr_info['pval']:.3e})")
    elif mr_info:
        lit_score += 1
        lit_evidence.append(f"MR数据可用 (p={mr_info['pval']:.3f})")
    
    # 基因功能分类
    if gene in ["FDX1", "LIAS", "LIPT1", "DLAT", "PDHB", "PDHX"]:
        func_category = "铜死亡执行基因"
    elif gene in ["SLC31A1", "ATP7A", "ATP7B", "ATOX1"]:
        func_category = "铜转运基因"
    elif gene in ["NFKB1", "MTOR", "NFE2L2", "HIF1A", "GPX4"]:
        func_category = "调控/通路基因"
    else:
        func_category = "其他"
    
    report_data.append({
        'Gene': gene,
        'GAT_Rank': gat_info['Rank'] if gat_info else 'N/A',
        'GAT_Score': f"{gat_info['P_target']:.4f}" if gat_info else 'N/A',
        'GAT_Percentile': f"{gat_info['Rank']/len(gat_predictions)*100:.1f}%" if gat_info else 'N/A',
        'STRING_RRA_Rank': rra_info['Rank'] if rra_info else 'N/A',
        'MR_pval': f"{mr_info['pval']:.3e}" if mr_info else 'N/A',
        'MR_FDR': f"{mr_info.get('fdr_qval', 1):.3e}" if mr_info else 'N/A',
        'in_Stroke_Targets': '是' if in_stroke_targets else '否',
        'in_BCP_Targets': '是' if in_bcp_targets else '否',
        'Lit_Score': lit_score,
        'Lit_Evidence': '; '.join(lit_evidence),
        'Func_Category': func_category,
        'is_cuproptosis': gat_info['is_cuproptosis'] if gat_info else 0,
        'dist_to_cuproptosis': f"{gat_info['dist_to_cuproptosis']:.4f}" if gat_info else 'N/A'
    })

# 打印报告
print(f"\n{'基因':<10} {'GAT排名':<10} {'GAT得分':<10} {'百分位':<10} {'STRING':<10} {'MR_pval':<12} {'脑卒中':<8} {'石竹烯':<8} {'文献分':<8} {'功能分类':<15}")
print("-"*120)

for row in sorted(report_data, key=lambda x: x['GAT_Rank'] if x['GAT_Rank'] != 'N/A' else 99999):
    print(f"{row['Gene']:<10} {row['GAT_Rank']:<10} {row['GAT_Score']:<10} {row['GAT_Percentile']:<10} "
          f"{row['STRING_RRA_Rank']:<10} {row['MR_pval']:<12} {row['in_Stroke_Targets']:<8} {row['in_BCP_Targets']:<8} "
          f"{row['Lit_Score']:<8} {row['Func_Category']:<15}")

# 统计分析
print("\n" + "="*100)
print("多源数据整合分析")
print("="*100)

# 1. GAT vs STRING网络中心性相关性
gat_rra_pairs = [(r['GAT_Rank'], r['STRING_RRA_Rank']) for r in report_data 
                 if r['GAT_Rank'] != 'N/A' and r['STRING_RRA_Rank'] != 'N/A']

if gat_rra_pairs:
    print(f"\n1. GAT排名 vs STRING网络中心性 (n={len(gat_rra_pairs)}):")
    for row in report_data:
        if row['GAT_Rank'] != 'N/A' and row['STRING_RRA_Rank'] != 'N/A':
            print(f"   {row['Gene']}: GAT Rank={row['GAT_Rank']}, STRING RRA Rank={row['STRING_RRA_Rank']}")

# 2. GAT vs MR数据关联
print(f"\n2. GAT得分 vs MR显著性:")
for row in report_data:
    if row['MR_pval'] != 'N/A' and row['GAT_Score'] != 'N/A':
        mr_sig = "显著" if float(row['MR_pval'].replace('e', 'E')) < 0.05 else "不显著"
        print(f"   {row['Gene']}: GAT={row['GAT_Score']}, MR p={row['MR_pval']} ({mr_sig})")

# 3. 文献证据支持度
print(f"\n3. 文献证据支持度汇总:")
lit_groups = defaultdict(list)
for row in report_data:
    lit_groups[row['Lit_Score']].append(row['Gene'])

for score in sorted(lit_groups.keys(), reverse=True):
    genes = lit_groups[score]
    print(f"   文献分={score}: {', '.join(genes)}")

# 4. 多源验证一致性
print(f"\n4. 多源验证一致性:")
high_conf = [r for r in report_data if r['Lit_Score'] >= 5 and r['GAT_Rank'] != 'N/A' and r['GAT_Rank'] <= 5000]
print(f"   高置信度基因 (文献分≥5 且 GAT排名≤5000):")
for row in sorted(high_conf, key=lambda x: x['GAT_Rank']):
    print(f"   - {row['Gene']}: GAT Rank={row['GAT_Rank']}, 文献分={row['Lit_Score']}, 证据: {row['Lit_Evidence']}")

# 5. 铜死亡执行基因 vs 调控基因
executors = [r for r in report_data if r['Func_Category'] == '铜死亡执行基因']
regulators = [r for r in report_data if r['Func_Category'] == '调控/通路基因']
transporters = [r for r in report_data if r['Func_Category'] == '铜转运基因']

print(f"\n5. 功能分类对比:")
if executors:
    avg_gat_rank_exec = sum(r['GAT_Rank'] for r in executors if r['GAT_Rank'] != 'N/A') / len([r for r in executors if r['GAT_Rank'] != 'N/A'])
    print(f"   铜死亡执行基因 (n={len(executors)}): 平均GAT排名={avg_gat_rank_exec:.0f}")
if regulators:
    avg_gat_rank_reg = sum(r['GAT_Rank'] for r in regulators if r['GAT_Rank'] != 'N/A') / len([r for r in regulators if r['GAT_Rank'] != 'N/A'])
    print(f"   调控/通路基因 (n={len(regulators)}): 平均GAT排名={avg_gat_rank_reg:.0f}")
if transporters:
    avg_gat_rank_trans = sum(r['GAT_Rank'] for r in transporters if r['GAT_Rank'] != 'N/A') / len([r for r in transporters if r['GAT_Rank'] != 'N/A'])
    print(f"   铜转运基因 (n={len(transporters)}): 平均GAT排名={avg_gat_rank_trans:.0f}")

# 保存完整报告
with open(RESULTS_DIR / "cuproptosis_literature_validation.csv", 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Gene', 'GAT_Rank', 'GAT_Score', 'GAT_Percentile', 'STRING_RRA_Rank', 
                    'MR_pval', 'MR_FDR', 'in_Stroke_Targets', 'in_BCP_Targets', 
                    'Lit_Score', 'Lit_Evidence', 'Func_Category', 'is_cuproptosis', 'dist_to_cuproptosis'])
    for row in sorted(report_data, key=lambda x: x['GAT_Rank'] if x['GAT_Rank'] != 'N/A' else 99999):
        writer.writerow([row['Gene'], row['GAT_Rank'], row['GAT_Score'], row['GAT_Percentile'],
                        row['STRING_RRA_Rank'], row['MR_pval'], row['MR_FDR'],
                        row['in_Stroke_Targets'], row['in_BCP_Targets'], row['Lit_Score'],
                        row['Lit_Evidence'], row['Func_Category'], row['is_cuproptosis'],
                        row['dist_to_cuproptosis']])

print(f"\n完整验证报告已保存: {RESULTS_DIR / 'cuproptosis_literature_validation.csv'}")

# 论文表格
print("\n" + "="*100)
print("论文格式表格")
print("="*100)

print("""
Table S2. Multi-source validation of cuproptosis gene predictions

| Gene   | GAT Rank | GAT Score | STRING RRA | MR p-value | Stroke Target | BCP Target | Lit Score | Category          |
|--------|----------|-----------|------------|------------|---------------|------------|-----------|-------------------|
""")
for row in sorted(report_data, key=lambda x: x['GAT_Rank'] if x['GAT_Rank'] != 'N/A' else 99999):
    print(f"| {row['Gene']:<6} | {row['GAT_Rank']:<8} | {row['GAT_Score']:<9} | {row['STRING_RRA_Rank']:<10} | {row['MR_pval']:<10} | {row['in_Stroke_Targets']:<13} | {row['in_BCP_Targets']:<10} | {row['Lit_Score']:<9} | {row['Func_Category']:<17} |")
