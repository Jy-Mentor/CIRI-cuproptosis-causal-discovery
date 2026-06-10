# ============================================================
# 论文图表数据准备脚本
# 从管道结果中提取6张图所需数据
# ============================================================

import pandas as pd
import numpy as np
import os
import sys

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
RESULT_DIR = os.path.join(BASE_DIR, "铜死亡和CIRI", "1", "results")
OUTPUT_DIR = os.path.join(BASE_DIR, "figure_data")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================
# 图1: 铜死亡GSVA评分箱线图 (Stroke vs Control)
# 数据源: cuproptosis_gsva_scores.csv
# ============================================================
print("[1/6] 准备图1数据: GSVA评分箱线图...")

gsva_scores_path = os.path.join(RESULT_DIR, "cuproptosis_gsva", "cuproptosis_gsva_scores.csv")
gsva_stats_path = os.path.join(RESULT_DIR, "cuproptosis_gsva", "cuproptosis_gsva_stats.csv")

if os.path.exists(gsva_scores_path):
    df_gsva = pd.read_csv(gsva_scores_path)
    # 只保留Stroke和Control组 (GSE61616)
    df_fig1 = df_gsva[df_gsva["Group"].isin(["Stroke", "Control"])].copy()
    df_fig1 = df_fig1.rename(columns={"Cuproptosis_Score": "Cuproptosis_Score"})
    df_fig1.to_csv(os.path.join(OUTPUT_DIR, "data_fig1_gsva_score.csv"), index=False)
    print(f"  图1数据: {df_fig1.shape[0]} 样本 (Stroke={sum(df_fig1['Group']=='Stroke')}, Control={sum(df_fig1['Group']=='Control')})")
else:
    print(f"  ⚠ GSVA评分文件不存在: {gsva_scores_path}")

# ============================================================
# 图2: GSEA富集分析结果条形图
# 数据源: cuproptosis_gsea_summary.csv
# ============================================================
print("[2/6] 准备图2数据: GSEA富集条形图...")

gsea_summary_path = os.path.join(RESULT_DIR, "cuproptosis_gsea", "cuproptosis_gsea_summary.csv")

if os.path.exists(gsea_summary_path):
    df_gsea = pd.read_csv(gsea_summary_path)
    # 计算-log10(FDR)用于可视化
    df_gsea["neg_log10_FDR"] = -np.log10(df_gsea["FDR"].clip(lower=1e-300))
    df_gsea["neg_log10_P"] = -np.log10(df_gsea["P_value"].clip(lower=1e-300))
    # 添加方向标签
    df_gsea["Direction"] = df_gsea["NES"].apply(lambda x: "Activated" if x > 0 else "Suppressed")
    df_gsea.to_csv(os.path.join(OUTPUT_DIR, "data_fig2_gsea.csv"), index=False)
    print(f"  图2数据: {df_gsea.shape[0]} 基因集")
else:
    print(f"  ⚠ GSEA汇总文件不存在: {gsea_summary_path}")

# ============================================================
# 图3: 单细胞铜死亡评分小提琴图
# 数据源: cuproptosis_by_celltype.csv
# ============================================================
print("[3/6] 准备图3数据: 单细胞小提琴图...")

sc_path = os.path.join(RESULT_DIR, "cuproptosis_singlecell", "cuproptosis_by_celltype.csv")

if os.path.exists(sc_path):
    df_sc = pd.read_csv(sc_path)
    # 排除Unknown类型
    df_sc = df_sc[df_sc["cell_type"] != "Unknown"].copy()
    # 计算-log10(p)用于标注显著性
    # 由于没有原始p值，用delta和样本量估算
    df_sc.to_csv(os.path.join(OUTPUT_DIR, "data_fig3_single_cell.csv"), index=False)
    print(f"  图3数据: {df_sc.shape[0]} 细胞类型")
else:
    print(f"  ⚠ 单细胞文件不存在: {sc_path}")

# ============================================================
# 图4: 铜死亡基因-炎症因子相关性气泡图
# 数据源: cuproptosis_cytokine_correlation.csv
# ============================================================
print("[4/6] 准备图4数据: 基因-炎症因子相关性气泡图...")

cyto_path = os.path.join(RESULT_DIR, "cuproptosis_immunology", "cuproptosis_cytokine_correlation.csv")

if os.path.exists(cyto_path):
    df_cyto = pd.read_csv(cyto_path)
    # 取Top 30相关性对
    df_cyto["abs_corr"] = df_cyto["Correlation"].abs()
    df_cyto = df_cyto.sort_values("abs_corr", ascending=False).head(30)
    df_cyto = df_cyto.drop(columns=["abs_corr"])
    df_cyto.to_csv(os.path.join(OUTPUT_DIR, "data_fig4_cytokine_corr.csv"), index=False)
    print(f"  图4数据: Top {df_cyto.shape[0]} 相关性对")
else:
    print(f"  ⚠ 炎症因子相关性文件不存在: {cyto_path}")

# ============================================================
# 图5: WGCNA模块关联 + 免疫浸润
# 数据源: cuproptosis_module_enrichment_extended.csv + immune_infiltration_estimates.csv
# ============================================================
print("[5/6] 准备图5数据: WGCNA模块 + 免疫浸润...")

wgcna_path = os.path.join(RESULT_DIR, "cuproptosis_wgcna", "cuproptosis_module_enrichment_extended.csv")
immune_path = os.path.join(RESULT_DIR, "cuproptosis_immunology", "immune_infiltration_estimates.csv")
immune_corr_path = os.path.join(RESULT_DIR, "cuproptosis_immunology", "cuproptosis_immune_correlation.csv")

if os.path.exists(wgcna_path):
    df_wgcna = pd.read_csv(wgcna_path)
    df_wgcna.to_csv(os.path.join(OUTPUT_DIR, "data_fig5_wgcna_modules.csv"), index=False)
    print(f"  图5 WGCNA数据: {df_wgcna.shape[0]} 模块")
else:
    print(f"  ⚠ WGCNA文件不存在: {wgcna_path}")

if os.path.exists(immune_corr_path):
    df_immune = pd.read_csv(immune_corr_path)
    df_immune.to_csv(os.path.join(OUTPUT_DIR, "data_fig5_immune.csv"), index=False)
    print(f"  图5 免疫数据: {df_immune.shape[0]} 细胞类型")
else:
    print(f"  ⚠ 免疫相关性文件不存在: {immune_corr_path}")

# ============================================================
# 图6: PPI网络拓扑 + 多维度融合评分
# 数据源: ppi_hub_genes.csv + ppi_neighbor_per_gene_stats.csv
# ============================================================
print("[6/6] 准备图6数据: PPI拓扑 + 融合评分...")

ppi_hub_path = os.path.join(RESULT_DIR, "cuproptosis_ppi", "ppi_hub_genes.csv")
ppi_neighbor_path = os.path.join(RESULT_DIR, "cuproptosis_ppi", "ppi_neighbor_per_gene_stats.csv")

# 合并PPI数据
ppi_data = {}

if os.path.exists(ppi_hub_path):
    df_hub = pd.read_csv(ppi_hub_path)
    ppi_data["hub_genes"] = df_hub
    print(f"  图6 Hub基因: {df_hub.shape[0]} 个")

if os.path.exists(ppi_neighbor_path):
    df_neighbor = pd.read_csv(ppi_neighbor_path)
    ppi_data["neighbor_stats"] = df_neighbor
    print(f"  图6 邻居统计: {df_neighbor.shape[0]} 个铜死亡基因")

# 创建融合数据: 铜死亡基因的PPI参数 + 融合评分
if "hub_genes" in ppi_data and "neighbor_stats" in ppi_data:
    df_hub = ppi_data["hub_genes"]
    df_neighbor = ppi_data["neighbor_stats"]
    
    # 为每个铜死亡基因计算融合评分
    fusion_rows = []
    for _, row in df_neighbor.iterrows():
        gene = row["Cuproptosis_Gene"]
        n_neighbors = row["N_Neighbors"]
        n_sig = row["N_Significant"]
        n_up = row["N_Up"]
        n_down = row["N_Down"]
        
        # 在hub基因中查找
        hub_match = df_hub[df_hub["Gene"] == gene]
        if len(hub_match) > 0:
            degree = hub_match.iloc[0]["Degree"]
            hub_score = hub_match.iloc[0]["HubScore"]
        else:
            degree = 0
            hub_score = 0
        
        # 计算综合评分
        sig_ratio = n_sig / max(n_neighbors, 1)
        fusion_score = 0.4 * sig_ratio + 0.3 * min(degree / 1000, 1.0) + 0.3 * hub_score
        
        fusion_rows.append({
            "Gene": gene,
            "Degree": degree,
            "N_Neighbors": n_neighbors,
            "N_Significant": n_sig,
            "N_Up": n_up,
            "N_Down": n_down,
            "Sig_Ratio": round(sig_ratio, 4),
            "HubScore": round(hub_score, 4),
            "Fusion_Score": round(fusion_score, 4)
        })
    
    df_fusion = pd.DataFrame(fusion_rows)
    df_fusion = df_fusion.sort_values("Fusion_Score", ascending=False)
    df_fusion.to_csv(os.path.join(OUTPUT_DIR, "data_fig6_ppi_fusion.csv"), index=False)
    print(f"  图6 融合数据: {df_fusion.shape[0]} 个铜死亡基因")
else:
    print(f"  ⚠ PPI数据不完整")

print(f"\n所有数据已保存至: {OUTPUT_DIR}")
print("完成!")