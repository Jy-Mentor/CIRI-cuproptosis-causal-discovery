"""Step 1-4: RF vs GAT 交叉验证分析"""
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RF_DIR = r"C:\Users\Jy-Mentor-7\Desktop\随机森林"

rf_pred = pd.read_csv(os.path.join(RF_DIR, "bridge_rf_predictions.csv"))
gat_scores = pd.read_csv(os.path.join(SCRIPT_DIR, "all_bridge_genes.csv"))
gene_table = pd.read_csv(os.path.join(RF_DIR, "gene_features_table.csv"))

print("=" * 70)
print("Step 1: RF Top 300 标签组成")
print("=" * 70)
rf_top300 = rf_pred.head(300)
tag_counts = rf_top300[['is_drug_target', 'is_disease_gene', 'is_bridge_gene']].sum()
print(f"  RF Top 300:")
print(f"    is_drug_target   = {int(tag_counts['is_drug_target'])}")
print(f"    is_disease_gene  = {int(tag_counts['is_disease_gene'])}")
print(f"    is_bridge_gene   = {int(tag_counts['is_bridge_gene'])}")

dt_only = ((rf_top300['is_drug_target'] == 1) & (rf_top300['is_disease_gene'] == 0) & (rf_top300['is_bridge_gene'] == 0)).sum()
dg_only = ((rf_top300['is_drug_target'] == 0) & (rf_top300['is_disease_gene'] == 1) & (rf_top300['is_bridge_gene'] == 0)).sum()
both_known = ((rf_top300['is_drug_target'] == 1) & (rf_top300['is_disease_gene'] == 1)).sum()
novel = ((rf_top300['is_drug_target'] == 0) & (rf_top300['is_disease_gene'] == 0)).sum()
print(f"    桥基因 (DT+DG)    = {int(both_known)}")
print(f"    仅 DT              = {int(dt_only)}")
print(f"    仅 DG              = {int(dg_only)}")
print(f"    全新无标签          = {int(novel)}")

print()
print("=" * 70)
print("Step 2: 限制未知基因候选池，重新计算 RF-GAT 交集")
print("=" * 70)

unknown_mask = (rf_pred['is_drug_target'] == 0) & (rf_pred['is_disease_gene'] == 0)
rf_unknown = rf_pred[unknown_mask].head(300)
print(f"  RF 未知基因候选池: {unknown_mask.sum():,} genes total")
print(f"  GAT 总基因数: {len(gat_scores):,}")

gat_top300 = gat_scores.head(300)
rf_unknown_genes = set(rf_unknown['gene_symbol'])
gat_top300_genes = set(gat_top300['gene_symbol'])
overlap = rf_unknown_genes & gat_top300_genes
print(f"  RF unknown Top 300: {len(rf_unknown_genes)} genes")
print(f"  GAT Top 300:        {len(gat_top300_genes)} genes")
print(f"  交集数量:           {len(overlap)}")
if len(overlap) <= 15:
    print(f"  交集基因: {sorted(overlap)}")
else:
    print(f"  前15个交集基因: {sorted(overlap)[:15]} ...")

print()
print("=" * 70)
print("Step 3: 全局未知基因 Spearman 相关性")
print("=" * 70)

merged = pd.merge(rf_pred[unknown_mask][['gene_symbol', 'bridge_probability']],
                  gat_scores[['gene_symbol', 'combined_score', 'drug_target_score', 'target_disease_score']],
                  on='gene_symbol', how='inner')
print(f"  可匹配未知基因数: {len(merged):,}")

corr, pval = spearmanr(merged['bridge_probability'], merged['combined_score'])
print(f"  RF vs GAT combined:    Spearman r = {corr:.4f}  (p = {pval:.4e})")

corr_dt, pval_dt = spearmanr(merged['bridge_probability'], merged['drug_target_score'])
print(f"  RF vs GAT drug_target: Spearman r = {corr_dt:.4f}  (p = {pval_dt:.4e})")

corr_td, pval_td = spearmanr(merged['bridge_probability'], merged['target_disease_score'])
print(f"  RF vs GAT target->dis: Spearman r = {corr_td:.4f}  (p = {pval_td:.4e})")

print()
print("=" * 70)
print("Step 4: GAT 对已知桥基因的评分验证")
print("=" * 70)

bridge_genes = set(gene_table[gene_table['is_bridge_gene'] == 1]['gene_symbol'])
print(f"  已知桥基因总数: {len(bridge_genes)}")

bridge_in_gat = gat_scores[gat_scores['gene_symbol'].isin(bridge_genes)].copy()
print(f"  GAT 中可匹配: {len(bridge_in_gat)}")

if len(bridge_in_gat) > 0:
    bridge_median = bridge_in_gat['combined_score'].median()
    global_median = gat_scores['combined_score'].median()
    bridge_mean = bridge_in_gat['combined_score'].mean()
    global_mean = gat_scores['combined_score'].mean()
    bridge_rank_median = bridge_in_gat.index.to_series().median() + 1
    print(f"  桥基因 combined_score 中位数: {bridge_median:.6f}")
    print(f"  全局 combined_score 中位数:   {global_median:.6f}")
    print(f"  桥基因 combined_score 均值:   {bridge_mean:.6f}")
    print(f"  全局 combined_score 均值:     {global_mean:.6f}")
    print(f"  桥基因排名中位数 (1-indexed): {bridge_rank_median:.0f} / {len(gat_scores):,}")
    print()
    
    top_pct = (bridge_in_gat.index < len(gat_scores) * 0.1).sum()
    print(f"  进入 GAT Top 10% 的桥基因:  {top_pct} / {len(bridge_in_gat)}")
    
    if len(bridge_in_gat) <= 20:
        for _, r in bridge_in_gat.sort_values('combined_score', ascending=False).iterrows():
            rank = bridge_in_gat.index.get_loc(r.name) + 1 if r.name in bridge_in_gat.index else '?'
            print(f"    {r['gene_symbol']:<12s} combined={r['combined_score']:.5f}  rank≈{r.name+1}")

print()
print("=" * 70)
print("RF 对已知桥基因的评分 (训练集内)")
print("=" * 70)

bridge_in_rf = rf_pred[rf_pred['gene_symbol'].isin(bridge_genes)]
print(f"  已知桥基因数: {len(bridge_in_rf)}")
print(f"  桥基因概率中位数: {bridge_in_rf['bridge_probability'].median():.4f}")
print(f"  桥基因概率均值:   {bridge_in_rf['bridge_probability'].mean():.4f}")
print(f"  排名中位数 (1-indexed): {bridge_in_rf.index.to_series().median()+1:.0f}")
top10_in_top50 = (bridge_in_rf.index < 50).sum()
print(f"  Top 50 中桥基因数: {top10_in_top50}")

print()
print("=" * 70)
print("综合解读")
print("=" * 70)
if corr < 0.1:
    print("  Spearman r ≈ 0: RF 与 GAT 在未知基因上几乎独立运作，无共识。")
    print("  两模型利用的信息维度正交（网络嵌入 vs 手工特征+序列），")
    print("  最终候选需两模型交叉验证或取交集作为高置信集。")
elif corr < 0.3:
    print(f"  Spearman r = {corr:.3f}: RF 与 GAT 存在弱正相关。")
    print("  交集小可能仅是排序差异的尾部效应。")
else:
    print(f"  Spearman r = {corr:.3f}: RF 与 GAT 有明显共识。")
print()