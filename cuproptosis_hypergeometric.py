# -*- coding: utf-8 -*-
"""铜死亡基因超几何富集检验"""

import os
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
GAT_DIR = r"C:\Users\Jy-Mentor-7\Desktop\GAT"

# 用户提供的铜死亡相关基因列表
CUPROPTOSIS_GENES = {
    "FDX1", "LIAS", "LIPT1", "LIPT2", "DLD", "DLAT", "DLST", "DBT",
    "GCSH", "PDHA1", "PDHB", "SLC31A1", "ATP7A", "ATP7B", "SLC11A2",
    "ATOX1", "CCS", "COX17", "SCO1", "SCO2", "COX11", "MT1A", "MT2A",
    "MTF1", "GLS", "CDKN2A", "SOD1", "SOD3", "CP", "TYR", "LOX",
    "LOXL1", "LOXL2", "LOXL3", "LOXL4", "ALB", "HAMP", "NLRP3", "NFE2L2"
}

def load_genes(filepath):
    with open(filepath) as f:
        return set(line.strip() for line in f if line.strip())

def hypergeom_test(N, K, n, k):
    """
    超几何检验 (单尾, 右尾: 检验是否显著多于期望)
    N: 总体大小
    K: 总体中成功数
    n: 抽样数
    k: 抽样中成功数
    返回 p-value, fold_enrichment, expected
    """
    if n == 0 or K == 0:
        return 1.0, 0.0, 0.0
    expected = n * K / N
    fold = k / expected if expected > 0 else float('inf')
    rv = hypergeom(N, K, n)
    p_value = rv.sf(k - 1)  # P(X >= k)
    return p_value, fold, expected


def main():
    print("=" * 70)
    print("  铜死亡基因超几何富集检验")
    print("=" * 70)

    # ---- 加载数据 ----
    all_subgraph_genes = load_genes(os.path.join(GAT_DIR, "subgraph_genes.txt"))
    drug_targets = load_genes(os.path.join(GAT_DIR, "drug_targets.txt"))
    disease_genes = load_genes(os.path.join(GAT_DIR, "disease_genes.txt"))
    known_genes = drug_targets | disease_genes

    # 读取 GAT 预测结果
    bridge_df = pd.read_csv(os.path.join(SCRIPT_DIR, "all_bridge_genes.csv"))
    predicted_genes = set(bridge_df["gene_symbol"].values)

    # ---- 铜死亡基因在图中的分布 ----
    cupro_in_graph = CUPROPTOSIS_GENES & all_subgraph_genes
    cupro_not_in_graph = CUPROPTOSIS_GENES - all_subgraph_genes

    print(f"\n[铜死亡基因总览]")
    print(f"  列表总数: {len(CUPROPTOSIS_GENES)}")
    print(f"  在图中: {len(cupro_in_graph)}")
    print(f"  不在图中 (无特征): {len(cupro_not_in_graph)}")
    if cupro_not_in_graph:
        print(f"  → {', '.join(sorted(cupro_not_in_graph))}")

    # 分类
    cupro_known = cupro_in_graph & known_genes
    cupro_unknown = cupro_in_graph - known_genes
    cupro_drug_target = cupro_in_graph & drug_targets
    cupro_disease = cupro_in_graph & disease_genes

    print(f"\n[铜死亡基因分类]")
    print(f"  已知靶点 (drug_target): {len(cupro_drug_target)} → {', '.join(sorted(cupro_drug_target)) if cupro_drug_target else '无'}")
    print(f"  已知疾病基因 (disease): {len(cupro_disease)} → {', '.join(sorted(cupro_disease)) if cupro_disease else '无'}")
    print(f"  已知 (drug ∪ disease): {len(cupro_known)}")
    print(f"  未知基因 (候选池): {len(cupro_unknown)}")
    print(f"  → {', '.join(sorted(cupro_unknown)) if cupro_unknown else '无'}")

    # ---- 背景计算 ----
    total_in_graph = len(all_subgraph_genes)
    total_known = len(known_genes & all_subgraph_genes)
    total_unknown = total_in_graph - total_known
    print(f"\n[背景统计]")
    print(f"  全图基因: {total_in_graph}")
    print(f"  已知基因 (drug ∪ disease): {total_known}")
    print(f"  未知基因 (候选池): {total_unknown}")
    print(f"  铜死亡基因在图中: {len(cupro_in_graph)}")

    # ========================================================================
    # 检验1: 铜死亡基因在未知基因池中是否显著富集
    # ========================================================================
    print(f"\n{'='*70}")
    print(f"[检验1] 全部未知基因中，铜死亡基因是否显著富集？")
    print(f"{'='*70}")
    N1 = total_in_graph          # 总体: 图中所有基因
    K1 = len(cupro_in_graph)      # 总体中铜死亡基因数
    n1 = total_unknown            # 抽样: 未知基因数
    k1 = len(cupro_unknown)       # 未知基因中铜死亡基因数
    p1, fold1, exp1 = hypergeom_test(N1, K1, n1, k1)

    print(f"  背景 N = {N1} (全图基因)")
    print(f"  铜死亡总数 K = {K1}")
    print(f"  未知基因数 n = {n1}")
    print(f"  未知中铜死亡数 k = {k1}")
    print(f"  期望值 E[k] = {exp1:.2f}")
    print(f"  富集倍数 (Fold Enrichment) = {fold1:.3f}")
    print(f"  p-value (超几何, 右尾) = {p1:.6f}")
    print(f"  显著性: {'*** 极显著' if p1 < 0.001 else '** 非常显著' if p1 < 0.01 else '* 显著' if p1 < 0.05 else 'ns 不显著'}")

    # ========================================================================
    # 检验2: GAT预测的桥梁靶点中，铜死亡基因是否显著富集
    # ========================================================================
    cupro_predicted = cupro_in_graph & predicted_genes
    cupro_in_top20 = cupro_in_graph & set(bridge_df.head(20)["gene_symbol"])
    cupro_in_top50 = cupro_in_graph & set(bridge_df.head(50)["gene_symbol"])
    cupro_in_top100 = cupro_in_graph & set(bridge_df.head(100)["gene_symbol"])

    print(f"\n{'='*70}")
    print(f"[检验2] GAT预测的桥梁靶点中，铜死亡基因是否显著富集？")
    print(f"{'='*70}")

    tests = []
    # 检验2a: Top-20
    N2 = total_unknown                # 候选池作为背景
    K2 = len(cupro_unknown)           # 候选池中铜死亡数
    n_top20 = 20
    k_top20 = len(cupro_in_top20)
    p_top20, fold_top20, exp_top20 = hypergeom_test(N2, K2, n_top20, k_top20)
    tests.append(("Top-20", p_top20, fold_top20))
    print(f"\n  --- Top-20 ---")
    print(f"  背景 N = {N2} (未知基因)")
    print(f"  铜死亡在未知中 K = {K2}")
    print(f"  Top-20 中铜死亡数 k = {k_top20}")
    print(f"  期望值 E[k] = {exp_top20:.2f}")
    print(f"  富集倍数 = {fold_top20:.3f}")
    print(f"  p-value = {p_top20:.6f}")

    # 检验2b: Top-50
    n_top50 = 50
    k_top50 = len(cupro_in_top50)
    p_top50, fold_top50, exp_top50 = hypergeom_test(N2, K2, n_top50, k_top50)
    tests.append(("Top-50", p_top50, fold_top50))
    print(f"\n  --- Top-50 ---")
    print(f"  Top-50 中铜死亡数 k = {k_top50}")
    print(f"  期望值 E[k] = {exp_top50:.2f}")
    print(f"  富集倍数 = {fold_top50:.3f}")
    print(f"  p-value = {p_top50:.6f}")

    # 检验2c: Top-100
    n_top100 = 100
    k_top100 = len(cupro_in_top100)
    p_top100, fold_top100, exp_top100 = hypergeom_test(N2, K2, n_top100, k_top100)
    tests.append(("Top-100", p_top100, fold_top100))
    print(f"\n  --- Top-100 ---")
    print(f"  Top-100 中铜死亡数 k = {k_top100}")
    print(f"  期望值 E[k] = {exp_top100:.2f}")
    print(f"  富集倍数 = {fold_top100:.3f}")
    print(f"  p-value = {p_top100:.6f}")

    # 检验2d: 全部预测 (4414)
    n_all = len(predicted_genes)
    k_all = len(cupro_predicted)
    p_all, fold_all, exp_all = hypergeom_test(N2, K2, n_all, k_all)
    tests.append(("全部预测", p_all, fold_all))
    print(f"\n  --- 全部预测基因 (n={n_all}) ---")
    print(f"  预测基因中铜死亡数 k = {k_all}")
    print(f"  期望值 E[k] = {exp_all:.2f}")
    print(f"  富集倍数 = {fold_all:.3f}")
    print(f"  p-value = {p_all:.6f}")

    # Multiple-testing correction (BH)
    pvals = [t[1] for t in tests]
    _, pvals_corrected, _, _ = multipletests(pvals, method='fdr_bh')
    print(f"\n  --- 多重检验校正 (Benjamini-Hochberg FDR) ---")
    for (name, _, fold), p_corr in zip(tests, pvals_corrected):
        sig = "***" if p_corr < 0.001 else "**" if p_corr < 0.01 else "*" if p_corr < 0.05 else "ns"
        print(f"  {name:12s}: FDR={p_corr:.6f} {sig}  Fold={fold:.3f}")

    # ========================================================================
    # 详细列表: 铜死亡基因在 GAT 预测中的排名
    # ========================================================================
    print(f"\n{'='*70}")
    print(f"[铜死亡基因在 GAT 桥梁预测中的排名]")
    print(f"{'='*70}")

    bridge_df_reset = bridge_df.reset_index(drop=True)
    bridge_df_reset["rank"] = bridge_df_reset.index + 1
    cupro_in_pred = bridge_df_reset[bridge_df_reset["gene_symbol"].isin(cupro_in_graph)]
    if len(cupro_in_pred) > 0:
        for _, row in cupro_in_pred.iterrows():
            print(f"  Rank {int(row['rank']):4d} | {row['gene_symbol']:12s} | "
                  f"DT={row['drug_target_score']:.4f} | "
                  f"TD={row['target_disease_score']:.4f} | "
                  f"Combined={row['combined_score']:.4f}")
    else:
        print("  无铜死亡基因在预测列表中")

    # 不在预测列表中的铜死亡基因
    cupro_not_pred = cupro_unknown - cupro_predicted
    if cupro_not_pred:
        print(f"\n  未出现在预测列表的铜死亡基因: {', '.join(sorted(cupro_not_pred))}")

    # ========================================================================
    # 汇总
    # ========================================================================
    print(f"\n{'='*70}")
    print(f"[汇总]")
    print(f"{'='*70}")
    print(f"  铜死亡基因总数: {len(CUPROPTOSIS_GENES)}")
    print(f"  在图中: {len(cupro_in_graph)} | 不在图: {len(cupro_not_in_graph)}")
    print(f"  已知 (drug/disease): {len(cupro_known)}")
    print(f"  未知 (候选): {len(cupro_unknown)}")
    print(f"  预测为桥梁 (all): {len(cupro_predicted)}")
    print(f"  预测为桥梁 (Top-20): {len(cupro_in_top20)}")
    print(f"  预测为桥梁 (Top-50): {len(cupro_in_top50)}")
    print(f"  预测为桥梁 (Top-100): {len(cupro_in_top100)}")
    print(f"\n  关键解释:")
    print(f"  - 检验1 p={p1:.4f}: 铜死亡基因是否倾向于'未知'而非'已知'")
    print(f"  - 检验2 检验不同排名阈值下铜死亡基因的预测富集程度")
    print(f"  - Fold > 1: 富集 (观察值 > 期望值) | Fold < 1: 贫化")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
