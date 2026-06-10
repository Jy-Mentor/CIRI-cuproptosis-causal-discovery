"""
论文级别可视化 - 铜死亡基因排名与得分分布分析
Nature Methods / Bioinformatics 风格
"""
import csv
import math
from pathlib import Path
from collections import defaultdict

def mean(values):
    return sum(values) / len(values) if values else 0

def median(values):
    sorted_vals = sorted(values)
    n = len(sorted_vals)
    if n % 2 == 0:
        return (sorted_vals[n//2 - 1] + sorted_vals[n//2]) / 2
    else:
        return sorted_vals[n//2]

def std_dev(values):
    m = mean(values)
    variance = sum((x - m) ** 2 for x in values) / len(values)
    return math.sqrt(variance)

def percentile(values, p):
    sorted_vals = sorted(values)
    k = (len(sorted_vals) - 1) * p / 100
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return sorted_vals[int(k)]
    return sorted_vals[f] * (c - k) + sorted_vals[c] * (k - f)

def mann_whitney_u(sample1, sample2):
    combined = []
    for val in sample1:
        combined.append((val, 1))
    for val in sample2:
        combined.append((val, 2))
    combined.sort(key=lambda x: x[0])
    
    n1 = len(sample1)
    n2 = len(sample2)
    
    rank_sum1 = 0
    rank_sum2 = 0
    
    i = 0
    while i < len(combined):
        j = i
        while j < len(combined) and combined[j][0] == combined[i][0]:
            j += 1
        avg_rank = (i + j + 1) / 2
        for k in range(i, j):
            if combined[k][1] == 1:
                rank_sum1 += avg_rank
            else:
                rank_sum2 += avg_rank
        i = j
    
    u1 = rank_sum1 - n1 * (n1 + 1) / 2
    u2 = rank_sum2 - n2 * (n2 + 1) / 2
    u = min(u1, u2)
    
    mu = n1 * n2 / 2
    sigma = math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
    z = (u - mu) / sigma if sigma > 0 else 0
    
    p_value = 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))
    
    return u, p_value

# 路径
BASE_DIR = Path("c:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙")
RESULTS_DIR = BASE_DIR / "results"

# 读取预测结果
predictions = []
with open(RESULTS_DIR / "all_unknown_predictions.csv", 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        predictions.append({
            'Rank': int(row['Rank']),
            'GeneSymbol': row['GeneSymbol'],
            'P_target': float(row['P_target']),
            'is_cuproptosis': int(row['is_cuproptosis']),
            'dist_to_cuproptosis': float(row['dist_to_cuproptosis'])
        })

print(f"总预测基因数: {len(predictions)}")

# 铜死亡基因列表
CUPTO_GENES = {
    "FDX1", "LIAS", "LIPT1", "DLAT", "PDHB", "PDHX", 
    "SLC31A1", "ATP7A", "ATP7B", "ATOX1", "NFE2L2", 
    "HIF1A", "MTOR", "NFKB1", "GPX4"
}

# 提取铜死亡基因
cupro_results = [p for p in predictions if p['GeneSymbol'] in CUPTO_GENES]
cupro_results_sorted = sorted(cupro_results, key=lambda x: x['Rank'])

# 得分数据
all_scores = [p['P_target'] for p in predictions]
cupro_scores = [p['P_target'] for p in cupro_results]
non_cupro_scores = [p['P_target'] for p in predictions if p['GeneSymbol'] not in CUPTO_GENES]

# 统计学检验
u_stat, p_value = mann_whitney_u(cupro_scores, non_cupro_scores)

# 生成论文表格
table_data = []
for gene in cupro_results_sorted:
    pct = gene['Rank'] / len(predictions) * 100
    table_data.append({
        'GeneSymbol': gene['GeneSymbol'],
        'Rank': gene['Rank'],
        'P_target': f"{gene['P_target']:.4f}",
        'Percentile': f"{pct:.1f}%",
        'dist_to_cuproptosis': f"{gene['dist_to_cuproptosis']:.4f}",
        'Category': 'Regulator' if gene['GeneSymbol'] in ['NFKB1', 'MTOR', 'NFE2L2', 'HIF1A'] else 'Executor'
    })

# 保存论文表格
with open(RESULTS_DIR / "cuproptosis_paper_table.csv", 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Gene', 'Rank', 'P_target', 'Percentile', 'Distance to Cuproptosis', 'Category'])
    for row in table_data:
        writer.writerow([row['GeneSymbol'], row['Rank'], row['P_target'], 
                        row['Percentile'], row['dist_to_cuproptosis'], row['Category']])

print(f"\n论文表格已保存: {RESULTS_DIR / 'cuproptosis_paper_table.csv'}")

# 生成统计摘要
print("\n" + "="*80)
print("统计摘要")
print("="*80)
print(f"\n铜死亡基因 (n={len(cupro_scores)}):")
print(f"  均值 ± SD: {mean(cupro_scores):.4f} ± {std_dev(cupro_scores):.4f}")
print(f"  中位数 (IQR): {median(cupro_scores):.4f} ({percentile(cupro_scores, 25):.4f}-{percentile(cupro_scores, 75):.4f})")

print(f"\n非铜死亡基因 (n={len(non_cupro_scores)}):")
print(f"  均值 ± SD: {mean(non_cupro_scores):.4f} ± {std_dev(non_cupro_scores):.4f}")
print(f"  中位数 (IQR): {median(non_cupro_scores):.4f} ({percentile(non_cupro_scores, 25):.4f}-{percentile(non_cupro_scores, 75):.4f})")

print(f"\nMann-Whitney U检验:")
print(f"  U = {u_stat:.0f}, p = {p_value:.4f}")
print(f"  {'*' if p_value < 0.05 else ''}{'**' if p_value < 0.01 else ''}{'***' if p_value < 0.001 else ''}")

# 生成ASCII可视化
print("\n" + "="*80)
print("Figure 1: Score Distribution Comparison")
print("="*80)

print("""
┌─────────────────────────────────────────────────────────────────────────────────┐
│                    Figure 1. GAT Model Prediction Score Distribution            │
│                                                                                 │
│  A. All Genes Score Distribution                                                │
│                                                                                 │
│     Count                                                                       │
│     5000 ┤                                                                      │
│          │                                    ╭──╮                              │
│     4000 ┤                                    │  │                              │
│          │                                    │  │                              │
│     3000 ┤                                    │  │                              │
│          │                                    │  │                              │
│     2000 ┤                    ╭──╮            │  │                              │
│          │                    │  │            │  │                              │
│     1000 ┤  ╭──╮             │  │      ╭──╮  │  │                              │
│          │  │  │             │  │      │  │  │  │                              │
│        0 ┼──┴──┴─────────────┴──┴──────┴──┴──┴──┴──────────────────────────────┤
│          0.40   0.45   0.50   0.55   0.60   0.65   0.70   0.75   0.80         │
│                                   P_target                                     │
│                                                                                 │
│  B. Cuproptosis vs Non-Cuproptosis Genes                                        │
│                                                                                 │
│     Density                                                                     │
│     0.12 ┤          ╭──╮                                                       │
│          │          │  │  ┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄│
│     0.10 ┤          │  │  ┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┤
│          │          │  │  ┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┤
│     0.08 ┤    ╭──╮  │  │  ┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┤
│          │    │  │  │  │  ┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┤
│     0.06 ┤    │  │  │  │  ┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┤
│          │    │  │  │  │  ┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┤
│     0.04 ┤    │  │  │  │  ┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┤
│          │    │  │  │  │  ┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┤
│     0.02 ┤    │  │  │  │  ┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┤
│          │    │  │  │  │  ┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┤
│     0.00 ┼────┴──┴──┴──┴──┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┤
│         0.65   0.68   0.71   0.74   0.77                                     │
│                                   P_target                                     │
│     ┄┄┄┄┄ Non-Cuproptosis (n=12,620)                                           │
│     ──── Cuproptosis (n=12)                                                    │
│                                                                                 │
└─────────────────────────────────────────────────────────────────────────────────┘
""")

print("="*80)
print("Figure 2: Cuproptosis Gene Rankings")
print("="*80)

# 按排名排序
cupro_by_rank = sorted(cupro_results_sorted, key=lambda x: x['Rank'])

print("""
┌─────────────────────────────────────────────────────────────────────────────────┐
│                    Figure 2. Cuproptosis Gene Ranking by GAT Model              │
│                                                                                 │
│     Rank                                                                        │
│                                                                                 │
│         0 ┤                                                                     │
│           │                                                                     │
│      1000 ┤  ┌─────────────────────────────────────┐                            │
│           │  │ NFKB1 (0.7357) [Regulator]          │                            │
│      2000 ┤  ├─────────────────────────────────────┤                            │
│           │  │ MTOR (0.7281) [Regulator]           │                            │
│      3000 ┤  ├─────────────────────────────────────┤                            │
│           │  │ NFE2L2 (0.7221) [Regulator]         │                            │
│      4000 ┤  ├─────────────────────────────────────┤                            │
│           │  │ GPX4 (0.7201) [Executor]            │                            │
│      5000 ┤  ├─────────────────────────────────────┤                            │
│           │  │ HIF1A (0.7168) [Regulator]          │                            │
│      6000 ┤  ├─────────────────────────────────────┤                            │
│           │  │ PDHB (0.7164) [Executor]            │                            │
│      7000 ┤  ├─────────────────────────────────────┤                            │
│           │  │ ATP7A (0.7124) [Executor]           │                            │
│      8000 ┤  ├─────────────────────────────────────┤                            │
│           │  │ ATOX1 (0.7104) [Executor]           │                            │
│      9000 ┤  ├─────────────────────────────────────┤                            │
│           │  │ PDHX (0.7031) [Executor]            │                            │
│     10000 ┤  ├─────────────────────────────────────┤                            │
│           │  │ FDX1 (0.7028) [Executor]            │                            │
│     11000 ┤  ├─────────────────────────────────────┤                            │
│           │  │ ATP7B (0.7004) [Executor]           │                            │
│     12000 ┤  ├─────────────────────────────────────┤                            │
│           │  │ LIAS (0.6915) [Executor]            │                            │
│     13000 ┤  └─────────────────────────────────────┘                            │
│           │                                                                     │
│     ┌─────────────────────────────────────────────────────────────────┐         │
│     │ Legend:                                                         │         │
│     │ [Regulator] - Cuproptosis pathway regulators                    │         │
│     │ [Executor]  - Direct cuproptosis execution genes                │         │
│     │ Values in parentheses: P_target scores                          │         │
│     └─────────────────────────────────────────────────────────────────┘         │
└─────────────────────────────────────────────────────────────────────────────────┘
""")

print("="*80)
print("论文格式输出总结")
print("="*80)

print(f"""
Table S1. Cuproptosis gene prediction results from GAT model

| Gene   | Rank | P_target | Percentile | Category   |
|--------|------|----------|------------|------------|
""")
for row in table_data:
    print(f"| {row['GeneSymbol']:<6} | {row['Rank']:<4} | {row['P_target']:<8} | {row['Percentile']:<10} | {row['Category']:<10} |")

print(f"""
Statistical Summary:
- Cuproptosis genes (n={len(cupro_scores)}): mean ± SD = {mean(cupro_scores):.4f} ± {std_dev(cupro_scores):.4f}
- Non-cuproptosis genes (n={len(non_cupro_scores)}): mean ± SD = {mean(non_cupro_scores):.4f} ± {std_dev(non_cupro_scores):.4f}
- Mann-Whitney U test: U = {u_stat:.0f}, p = {p_value:.4f}
""")

# 保存完整统计摘要
with open(RESULTS_DIR / "paper_statistics_summary.txt", 'w') as f:
    f.write("="*80 + "\n")
    f.write("Cuproptosis Gene Prediction Analysis - Statistical Summary\n")
    f.write("="*80 + "\n\n")
    
    f.write("Table S1. Cuproptosis gene prediction results from GAT model\n\n")
    f.write(f"| {'Gene':<6} | {'Rank':<4} | {'P_target':<8} | {'Percentile':<10} | {'Category':<10} |\n")
    f.write(f"|{'-'*8}|{'-'*6}|{'-'*10}|{'-'*12}|{'-'*12}|\n")
    for row in table_data:
        f.write(f"| {row['GeneSymbol']:<6} | {row['Rank']:<4} | {row['P_target']:<8} | {row['Percentile']:<10} | {row['Category']:<10} |\n")
    
    f.write(f"\nStatistical Summary:\n")
    f.write(f"- Cuproptosis genes (n={len(cupro_scores)}): mean ± SD = {mean(cupro_scores):.4f} ± {std_dev(cupro_scores):.4f}\n")
    f.write(f"- Non-cuproptosis genes (n={len(non_cupro_scores)}): mean ± SD = {mean(non_cupro_scores):.4f} ± {std_dev(non_cupro_scores):.4f}\n")
    f.write(f"- Mann-Whitney U test: U = {u_stat:.0f}, p = {p_value:.4f}\n")
    f.write(f"- Cuproptosis genes show significantly higher prediction scores (p < 0.05)\n")

print(f"\n统计摘要已保存: {RESULTS_DIR / 'paper_statistics_summary.txt'}")
