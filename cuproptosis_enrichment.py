"""
铜死亡基因在共识未知基因集中的超几何检验富集分析
"""
import pandas as pd
from scipy.stats import hypergeom

cuproptosis_genes = [
    "FDX1", "LIAS", "LIPT1", "LIPT2", "DLD", "DLAT", "DLST", "DBT", "GCSH",
    "PDHA1", "PDHB", "SLC31A1", "ATP7A", "ATP7B", "SLC11A2", "ATOX1", "CCS",
    "COX17", "SCO1", "SCO2", "COX11", "MT1A", "MT2A", "MTF1", "GLS", "CDKN2A",
    "SOD1", "SOD3", "CP", "TYR", "LOX", "LOXL1", "LOXL2", "LOXL3", "LOXL4",
    "ALB", "HAMP", "NLRP3", "NFE2L2",
]

df_consensus = pd.read_csv("consensus_unknown_genes.csv")
all_genes = df_consensus["gene_symbol"].tolist()
all_genes_set = set(all_genes)

n_total = 4414
n_cuproptosis_total = len(cuproptosis_genes)
n_consensus = len(df_consensus)

cupr_in_consensus = [g for g in cuproptosis_genes if g in all_genes_set]
n_cupr_in = len(cupr_in_consensus)

p_hyper = hypergeom.sf(n_cupr_in - 1, n_total, n_cuproptosis_total, n_consensus)

df_consensus["rank"] = range(1, len(df_consensus) + 1)

print("=" * 60)
print("  铜死亡基因富集分析 — 超几何检验")
print("=" * 60)
print(f"  背景总基因数 (完全未知):       {n_total}")
print(f"  共识基因数 (全部未知基因):      {n_consensus}")
print(f"  铜死亡基因总数:                 {n_cuproptosis_total}")
print(f"  铜死亡基因出现在共识集中:       {n_cupr_in}")
print(f"  期望值 (E):                     {n_cuproptosis_total * n_consensus / n_total:.2f}")
print(f"  P-value (超几何检验):           {p_hyper:.6e}")
print()

if n_cupr_in > 0:
    print(f"  命中铜死亡基因 ({n_cupr_in}个):")
    df_cupr = df_consensus[df_consensus["gene_symbol"].isin(cupr_in_consensus)].sort_values(
        "consensus_score", ascending=False
    )
    for _, row in df_cupr.iterrows():
        print(
            f"    排名 {int(row['rank']):>4d}  {row['gene_symbol']:10s}  "
            f"consensus={row['consensus_score']:.2f}  "
            f"GAT={row['gat_rank_pct']:.2f}  LGBM={row['lgbm_rank_pct']:.2f}"
        )
else:
    print("  铜死亡基因未出现在共识基因集中")

print()
print("-" * 60)
print("  Top-K 截断富集分析")
print("-" * 60)
for k in [50, 100, 200, 500, 1000, 2000]:
    df_topk = df_consensus.head(k)
    topk_set = set(df_topk["gene_symbol"].tolist())
    cupr_topk = [g for g in cuproptosis_genes if g in topk_set]
    n_cupr_k = len(cupr_topk)
    p_k = hypergeom.sf(n_cupr_k - 1, n_total, n_cuproptosis_total, k)
    exp_k = n_cuproptosis_total * k / n_total
    sig = ""
    if p_k < 0.001:
        sig = " ***"
    elif p_k < 0.01:
        sig = " **"
    elif p_k < 0.05:
        sig = " *"
    print(f"  Top-{k:>4d}: {n_cupr_k}个铜死亡基因 (期望 {exp_k:.1f}), P={p_k:.4e}{sig}")
    if n_cupr_k > 0 and k <= 200:
        for g in cupr_topk:
            rank_val = int(df_topk.loc[df_topk["gene_symbol"] == g, "rank"].values[0])
            print(f"            #{rank_val} {g}")

not_found = [g for g in cuproptosis_genes if g not in all_genes_set]
if not_found:
    print()
    print(f"  不在共识未知基因集中: {len(not_found)} 个")
    print(f"    {not_found}")

print()
print("=" * 60)