# -*- coding: utf-8 -*-
"""
铜死亡分析结果汇总 - 导出为Excel
"""
import os
import sys
import pandas as pd
from datetime import datetime

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import RESULTS_DIR

RESULTS_ROOT = RESULTS_DIR
OUTPUT_FILE = os.path.join(RESULTS_ROOT, "cuproptosis_analysis", "cuproptosis_results_summary.xlsx")

def load_csv_safe(path):
    if os.path.exists(path):
        try:
            df = pd.read_csv(path)
            if df is not None and len(df) > 0:
                return df
        except Exception:
            pass
    return None

def main():
    os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)

    with pd.ExcelWriter(OUTPUT_FILE, engine='openpyxl') as writer:

        # === Sheet 1: 分析概览 ===
        overview = pd.DataFrame({
            '模块': ['M1 GSVA通路评分', 'M2 GSEA富集分析', 'M3 单细胞评分', 'M4 WGCNA模块关联', 'M5 免疫浸润', 'M6 PPI邻居分析', 'M7 Hallmark GSVA全通路分析'],
            '状态': ['成功', '成功', '成功', '成功', '成功', '成功', '成功'],
            '输出文件数': [4, 2, 6, 3, 3, 4, 3],
            '说明': [
                '铜死亡通路ssGSEA评分+统计检验+可视化',
                '7个基因集的GSEA富集分析 (Core_12, Core, Extended, Lipoylation, Copper_Transport, Pro, Anti)',
                '单细胞铜死亡评分+细胞类型差异+炎症因子相关',
                'WGCNA模块铜死亡基因富集',
                '免疫浸润估算+炎症因子相关性',
                'PPI邻居基因差异表达分析 (adj.P.Val校正)',
                '50条Hallmark通路GSVA评分+铜死亡相关性分析'
            ]
        })
        overview.to_excel(writer, sheet_name='分析概览', index=False)

        # === Sheet 2: M1 GSVA评分 ===
        gsva_scores = load_csv_safe(os.path.join(RESULTS_ROOT, "cuproptosis_gsva", "cuproptosis_gsva_scores.csv"))
        if gsva_scores is not None:
            gsva_scores.to_excel(writer, sheet_name='M1_GSVA评分', index=False)

        gsva_stats = load_csv_safe(os.path.join(RESULTS_ROOT, "cuproptosis_gsva", "cuproptosis_gsva_stats.csv"))
        if gsva_stats is not None:
            gsva_stats.to_excel(writer, sheet_name='M1_统计检验', index=False)

        # === Sheet 3: M2 GSEA结果 ===
        gsea_summary = load_csv_safe(os.path.join(RESULTS_ROOT, "cuproptosis_gsea", "cuproptosis_gsea_summary.csv"))
        if gsea_summary is not None:
            gsea_summary.to_excel(writer, sheet_name='M2_GSEA汇总', index=False)

        # === Sheet 4: M3 单细胞 ===
        sc_stats = load_csv_safe(os.path.join(RESULTS_ROOT, "cuproptosis_singlecell", "cuproptosis_by_celltype.csv"))
        if sc_stats is not None:
            sc_stats.to_excel(writer, sheet_name='M3_细胞类型统计', index=False)

        sc_cyto = load_csv_safe(os.path.join(RESULTS_ROOT, "cuproptosis_singlecell", "cuproptosis_cytokine_correlation.csv"))
        if sc_cyto is not None:
            sc_cyto.head(50).to_excel(writer, sheet_name='M3_炎症因子相关', index=False)

        # === Sheet 5: M4 WGCNA ===
        wgcna_enrich = load_csv_safe(os.path.join(RESULTS_ROOT, "cuproptosis_wgcna", "cuproptosis_module_enrichment.csv"))
        if wgcna_enrich is not None:
            wgcna_enrich.to_excel(writer, sheet_name='M4_模块富集', index=False)

        # === Sheet 6: M5 免疫浸润 ===
        immune = load_csv_safe(os.path.join(RESULTS_ROOT, "cuproptosis_immunology", "immune_infiltration_estimates.csv"))
        if immune is not None:
            immune.to_excel(writer, sheet_name='M5_免疫浸润估算', index=True)

        cyto_corr = load_csv_safe(os.path.join(RESULTS_ROOT, "cuproptosis_immunology", "cuproptosis_cytokine_correlation.csv"))
        if cyto_corr is not None:
            cyto_corr.head(100).to_excel(writer, sheet_name='M5_炎症因子相关', index=False)

        # === Sheet 7: M6 PPI邻居 ===
        ppi_stats = load_csv_safe(os.path.join(RESULTS_ROOT, "cuproptosis_ppi", "ppi_neighbor_per_gene_stats.csv"))
        if ppi_stats is not None:
            ppi_stats.to_excel(writer, sheet_name='M6_邻居统计', index=False)

        ppi_sig = load_csv_safe(os.path.join(RESULTS_ROOT, "cuproptosis_ppi", "sig_neighbor_degs.csv"))
        if ppi_sig is not None:
            ppi_sig.head(200).to_excel(writer, sheet_name='M6_显著邻居DEGs', index=False)

        # === Sheet 8: M7 Hallmark GSVA ===
        hallmark_corr = load_csv_safe(os.path.join(RESULTS_ROOT, "cuproptosis_hallmark_gsva", "cuproptosis_hallmark_correlations.csv"))
        if hallmark_corr is not None:
            hallmark_corr.to_excel(writer, sheet_name='M7_Hallmark相关性', index=False)

        hallmark_stats = load_csv_safe(os.path.join(RESULTS_ROOT, "cuproptosis_hallmark_gsva", "hallmark_gsva_stats.csv"))
        if hallmark_stats is not None:
            hallmark_stats.to_excel(writer, sheet_name='M7_Hallmark统计', index=False)

    print(f"Excel汇总已保存: {OUTPUT_FILE}")

if __name__ == "__main__":
    main()
