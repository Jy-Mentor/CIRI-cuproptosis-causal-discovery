"""
共识桥接基因分析 + 通路富集分析
HeteroGAT + LightGBM 双模型共识筛选完全未知桥接基因
"""
import pandas as pd
import numpy as np
from scipy.stats import percentileofscore
import os

# ============ 文件路径配置 ============
GAT_FILE = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\all_bridge_genes.csv"
LGBM_FILE = r"C:\Users\Jy-Mentor-7\Desktop\随机森林\bridge_predictions_v8_final.csv"
OUTPUT_DIR = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙"
CONSENSUS_OUTPUT = os.path.join(OUTPUT_DIR, "consensus_unknown_genes.csv")
TOP50_OUTPUT = os.path.join(OUTPUT_DIR, "consensus_top50.csv")
KEGG_OUTPUT = os.path.join(OUTPUT_DIR, "enrichment_kegg.csv")
GOBP_OUTPUT = os.path.join(OUTPUT_DIR, "enrichment_gobp.csv")

# ============ 1. 数据加载 ============
print("=" * 60)
print("加载数据...")

df_gat = pd.read_csv(GAT_FILE)
df_lgbm = pd.read_csv(LGBM_FILE)

print(f"  GAT 模型基因数: {len(df_gat)}")
print(f"  LightGBM 模型基因数: {len(df_lgbm)}")

# GAT 文件列: gene_symbol, drug_target_score, target_disease_score, combined_score
# LightGBM 文件列: gene_symbol, is_drug_target, is_disease_gene, is_bridge_gene,
#                  prob_drug_target, prob_disease_gene, bridge_probability

# 将 LightGBM 的标签列合并到 GAT 数据中（GAT 文件缺少标签列）
df_labels = df_lgbm[["gene_symbol", "is_drug_target", "is_disease_gene", "is_bridge_gene"]]
df_gat = df_gat.merge(df_labels, on="gene_symbol", how="inner")

print(f"  合并后 GAT 模型基因数: {len(df_gat)}")

# ============ 2. 筛选完全未知基因 ============
unknown_mask_gat = (df_gat["is_drug_target"] == 0) & (df_gat["is_disease_gene"] == 0)
unknown_mask_lgbm = (df_lgbm["is_drug_target"] == 0) & (df_lgbm["is_disease_gene"] == 0)

df_gat_unknown = df_gat[unknown_mask_gat].copy()
df_lgbm_unknown = df_lgbm[unknown_mask_lgbm].copy()

print(f"\n  完全未知基因 (GAT 中): {len(df_gat_unknown)}")
print(f"  完全未知基因 (LightGBM 中): {len(df_lgbm_unknown)}")

# 取两个模型共有的未知基因
common_genes = set(df_gat_unknown["gene_symbol"]) & set(df_lgbm_unknown["gene_symbol"])
print(f"  两模型共有未知基因: {len(common_genes)}")

df_gat_unknown = df_gat_unknown[df_gat_unknown["gene_symbol"].isin(common_genes)].copy()
df_lgbm_unknown = df_lgbm_unknown[df_lgbm_unknown["gene_symbol"].isin(common_genes)].copy()

# ============ 3. 计算排名百分位数 ============
# 在未知基因子集中计算排名百分位数
# 值越大 → 排名越靠前 → 百分位越高
gat_scores = df_gat_unknown["combined_score"].values
lgbm_scores = df_lgbm_unknown["bridge_probability"].values

# 使用 scipy.percentileofscore 计算百分位排名
df_gat_unknown["gat_rank_pct"] = [
    percentileofscore(gat_scores, v, kind="rank") for v in gat_scores
]
df_lgbm_unknown["lgbm_rank_pct"] = [
    percentileofscore(lgbm_scores, v, kind="rank") for v in lgbm_scores
]

# 合并
consensus = df_gat_unknown[["gene_symbol", "gat_rank_pct"]].merge(
    df_lgbm_unknown[["gene_symbol", "lgbm_rank_pct"]],
    on="gene_symbol",
    how="inner",
)

# 共识分数 = 两个百分位数的平均值
consensus["consensus_score"] = (
    consensus["gat_rank_pct"] + consensus["lgbm_rank_pct"]
) / 2.0

# 按共识分数降序排列
consensus = consensus.sort_values("consensus_score", ascending=False).reset_index(drop=True)

print(f"\n  共识基因数: {len(consensus)}")
print(f"  Top 10 共识基因:")
print(consensus.head(10).to_string(index=False))

# ============ 4. 输出文件 ============
consensus.to_csv(CONSENSUS_OUTPUT, index=False)
print(f"\n  共识分析结果已保存: {CONSENSUS_OUTPUT}")

top50 = consensus.head(50)
top50.to_csv(TOP50_OUTPUT, index=False)
print(f"  Top 50 共识基因已保存: {TOP50_OUTPUT}")

top50_genes = top50["gene_symbol"].tolist()
print(f"\n  **共识 Top 50 基因**:")
for i, g in enumerate(top50_genes, 1):
    print(f"    {i:2d}. {g}")

# ============ 5. 通路富集分析 ============
print("\n" + "=" * 60)
print("通路富集分析 (gseapy → Enrichr API)...")

BG_GENE_COUNT = 15389  # 背景基因总数

try:
    import gseapy as gp

    gene_list = top50_genes

    # KEGG 2021 Human
    print("\n  [1/2] KEGG 富集分析...")
    try:
        kegg_results = gp.enrichr(
            gene_list=gene_list,
            gene_sets="KEGG_2021_Human",
            organism="human",
            outdir=None,
            no_plot=True,
            cutoff=0.05,
        )
        kegg_df = kegg_results.results.copy()
        # 计算 overlap 基因数
        kegg_df["overlap_genes"] = kegg_df["Overlap"].apply(
            lambda x: x.split("/")[0] if isinstance(x, str) else ""
        )
        kegg_df = kegg_df.rename(columns={
            "Term": "term",
            "Overlap": "overlap",
            "P-value": "p_value",
            "Adjusted P-value": "adj_p_value",
            "Genes": "genes",
        })
        kegg_out = kegg_df[["term", "overlap", "p_value", "adj_p_value", "genes"]]
        kegg_out.to_csv(KEGG_OUTPUT, index=False)
        print(f"    KEGG 显著通路数: {len(kegg_out)}")
        print(f"    已保存: {KEGG_OUTPUT}")
        print(kegg_out.head(10).to_string(index=False))
    except Exception as e:
        print(f"    KEGG 富集失败: {e}")
        kegg_out = None

    # GO Biological Process 2023
    print("\n  [2/2] GO Biological Process 富集分析...")
    try:
        gobp_results = gp.enrichr(
            gene_list=gene_list,
            gene_sets="GO_Biological_Process_2023",
            organism="human",
            outdir=None,
            no_plot=True,
            cutoff=0.05,
        )
        gobp_df = gobp_results.results.copy()
        gobp_df["overlap_genes"] = gobp_df["Overlap"].apply(
            lambda x: x.split("/")[0] if isinstance(x, str) else ""
        )
        gobp_df = gobp_df.rename(columns={
            "Term": "term",
            "Overlap": "overlap",
            "P-value": "p_value",
            "Adjusted P-value": "adj_p_value",
            "Genes": "genes",
        })
        gobp_out = gobp_df[["term", "overlap", "p_value", "adj_p_value", "genes"]]
        gobp_out.to_csv(GOBP_OUTPUT, index=False)
        print(f"    GO BP 显著通路数: {len(gobp_out)}")
        print(f"    已保存: {GOBP_OUTPUT}")
        print(gobp_out.head(10).to_string(index=False))
    except Exception as e:
        print(f"    GO BP 富集失败: {e}")
        gobp_out = None

except ImportError:
    print("  gseapy 未安装，尝试安装中...")
    import subprocess
    subprocess.check_call(["pip", "install", "gseapy", "-q"])
    print("  请重新运行脚本以使用 gseapy。")

print("\n" + "=" * 60)
print("分析完成！")
print(f"  共识文件: {CONSENSUS_OUTPUT}")
print(f"  Top 50:   {TOP50_OUTPUT}")
if kegg_out is not None:
    print(f"  KEGG:     {KEGG_OUTPUT}")
if gobp_out is not None:
    print(f"  GO BP:    {GOBP_OUTPUT}")
