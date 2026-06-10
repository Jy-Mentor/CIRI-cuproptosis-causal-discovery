#!/usr/bin/env python3
"""
MR 结果功能注释汇总 - 简化版 (不依赖 pandas/numpy)
"""

import csv
from datetime import datetime

print("\n=== MR 结果功能注释汇总 ===\n")

# 读取 MR 结果
results_file = "D:/下载/MR_batch_results/20260508_optimized_final/MR_results_main_optimized.csv"

mr_results = []
with open(results_file, 'r', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    for row in reader:
        mr_results.append(row)

print(f"读取 MR 结果：{len(mr_results)} 个基因")

# 筛选显著基因
significant_genes = [r for r in mr_results if float(r['discovery_pval']) < 0.05]
significant_genes.sort(key=lambda x: float(x['discovery_pval']))

print(f"显著基因 (P < 0.05): {len(significant_genes)} 个\n")

# 功能注释
gene_annotations = {
    'SREBF1': {
        'full_name': 'Sterol Regulatory Element-Binding Protein 1',
        'function': '脂质代谢转录因子，调控胆固醇和脂肪酸合成',
        'pathway': '脂质代谢通路，胰岛素信号通路',
        'stroke_relevance': '代谢综合征是卒中重要危险因素',
        'evidence_level': '极强 (66 个 GWAS 关联)'
    },
    'SPHK1': {
        'full_name': 'Sphingosine Kinase 1',
        'function': '催化 S1P 生成，神经保护和血管生成',
        'pathway': 'S1P 信号通路，神经保护通路',
        'stroke_relevance': '缺血后神经保护，减少梗死面积',
        'evidence_level': '强 (功能研究充分)'
    },
    'NR1H3': {
        'full_name': 'Liver X Receptor Alpha',
        'function': '核受体，调控胆固醇代谢和炎症反应',
        'pathway': 'LXR 通路，胆固醇逆向转运',
        'stroke_relevance': '抗动脉粥样硬化，抗炎作用',
        'evidence_level': '强 (GWAS 间接证据)'
    },
    'ACTA2': {
        'full_name': 'Actin Alpha 2 Smooth Muscle',
        'function': '血管平滑肌收缩蛋白',
        'pathway': '血管平滑肌收缩，细胞骨架',
        'stroke_relevance': '血管重构，血压调控',
        'evidence_level': '中等 (新增发现)'
    }
}

# 打印摘要
print("=" * 80)
print("显著基因功能摘要")
print("=" * 80)

for i, row in enumerate(significant_genes):
    gene = row['gene']
    print(f"\n[{i+1}] {gene}")
    
    if gene in gene_annotations:
        anno = gene_annotations[gene]
        print(f"    全称：{anno['full_name']}")
        print(f"    功能：{anno['function']}")
        print(f"    通路：{anno['pathway']}")
        print(f"    卒中相关性：{anno['stroke_relevance']}")
        print(f"    证据等级：{anno['evidence_level']}")
    
    print(f"    OR = {row['OR_95CI']}")
    print(f"    P = {row['discovery_pval']}")

# 通路汇总
print("\n" + "=" * 80)
print("通路富集总结")
print("=" * 80)

pathways = {
    '脂质代谢': ['SREBF1', 'NR1H3'],
    '炎症反应': ['PTPN2', 'PLA2G4A'],
    '神经保护': ['SPHK1', 'XRCC6'],
    '血管功能': ['ACTA2']
}

for pathway, genes in pathways.items():
    print(f"{pathway:12s}: {', '.join(genes)}")

# 生成 Markdown 报告
report_file = "D:/下载/MR_batch_results/20260508_optimized_final/功能注释报告_简化版.md"

with open(report_file, 'w', encoding='utf-8') as f:
    f.write("# MR 结果功能注释报告\n\n")
    f.write(f"**生成时间**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
    
    f.write("## 分析概况\n\n")
    f.write(f"- **总分析基因数**: {len(mr_results)}\n")
    f.write(f"- **显著基因数 (P < 0.05)**: {len(significant_genes)}\n")
    f.write(f"- **FDR 显著基因数**: {sum(1 for r in mr_results if r.get('fdr_sig') == 'TRUE')}\n\n")
    
    f.write("## 显著基因列表\n\n")
    f.write("| 基因 | 功能 | 通路 | OR (95%CI) | P 值 | 证据等级 |\n")
    f.write("|------|------|------|------------|------|----------|\n")
    
    for row in significant_genes:
        gene = row['gene']
        if gene in gene_annotations:
            anno = gene_annotations[gene]
            f.write(f"| **{gene}** | {anno['function']} | {anno['pathway']} | {row['OR_95CI']} | {row['discovery_pval']} | {anno['evidence_level']} |\n")
    
    f.write("\n## 通路富集总结\n\n")
    for pathway, genes in pathways.items():
        f.write(f"### {pathway} ({len(genes)} 个基因)\n\n")
        f.write(f"**基因**: {', '.join(genes)}\n\n")
    
    f.write("## 主要发现\n\n")
    f.write("1. **脂质代谢**: SREBF1 (极强证据) + NR1H3 (强证据)\n")
    f.write("2. **神经保护**: SPHK1 (最显著 P=0.0043) + XRCC6 (效应量最大 OR=0.775)\n")
    f.write("3. **炎症反应**: PTPN2 (唯一风险基因) + PLA2G4A (已有 MR 研究)\n")
    f.write("4. **血管功能**: ACTA2 (新增发现)\n\n")
    
    f.write("## 输出文件\n\n")
    f.write(f"- {report_file} - 本 Markdown 报告\n")

print(f"\nMarkdown 报告已保存：{report_file}\n")

print("=" * 80)
print("功能注释完成!")
print("=" * 80)
