# -*- coding: utf-8 -*-
"""
============================================================================
主控脚本: 铜死亡模块分析管线
============================================================================
目的: 一键运行所有铜死亡相关分析模块
输入: 现有管线输出数据
输出: 所有模块的分析结果和可视化

模块列表:
1. GSVA通路活性评分 (cuproptosis_gsva.py)
2. GSEA富集分析 (cuproptosis_gsea.py)
3. 单细胞铜死亡评分 (cuproptosis_singlecell.py)
4. WGCNA模块关联 (cuproptosis_wgcna.py)
5. 免疫浸润相关性 (cuproptosis_immunology.py)
6. PPI邻居差异表达 (cuproptosis_ppi_neighbors.py)
7. Hallmark GSVA全通路分析 (cuproptosis_hallmark_gsva.py)

参考文献:
- Tsvetkov P, et al. Science 2022 (PMID:35298263) - 铜死亡核心基因集
- Hanzelmann et al. BMC Bioinformatics 2013 - GSVA方法
- Liberzon et al. Cell Systems 2015 - MSigDB Hallmark基因集
- Liu H, et al. Am J Cancer Res 2022 (PMID:36119826) - 12核心铜死亡基因
============================================================================
"""

import os
import sys
import time
import logging
import warnings
from datetime import datetime
from typing import Dict, List, Optional

os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import BASE_DIR, RESULTS_DIR
from scripts.utils import setup_logger, ensure_dir

warnings.filterwarnings('ignore')

STAGE_DIR = os.path.join(RESULTS_DIR, "cuproptosis_analysis")
ensure_dir(STAGE_DIR)

logger = setup_logger("cuproptosis_main", os.path.join(STAGE_DIR, "cuproptosis_analysis.log"))


def run_module(module_name: str, module_func, *args, **kwargs):
    """运行单个模块并记录日志"""
    logger.info("\n" + "=" * 70)
    logger.info(f"  开始运行: {module_name}")
    logger.info("=" * 70)
    
    start_time = time.time()
    
    try:
        result = module_func(*args, **kwargs)
        elapsed = time.time() - start_time
        logger.info(f"  {module_name} 完成! 耗时: {elapsed:.1f}秒")
        return result
    except Exception as e:
        elapsed = time.time() - start_time
        logger.error(f"  {module_name} 失败! 耗时: {elapsed:.1f}秒")
        logger.error(f"  错误: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return None


def main():
    """主函数 - 运行所有铜死亡分析模块"""
    start_time = time.time()
    
    logger.info("\n" + "=" * 70)
    logger.info("  铜死亡模块分析管线启动")
    logger.info(f"  时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info("=" * 70)
    
    results = {}
    
    # 模块1: GSVA通路活性评分
    logger.info("\n" + "-" * 70)
    logger.info("  模块1: GSVA通路活性评分")
    logger.info("-" * 70)
    try:
        from scripts.cuproptosis_gsva import main as gsva_main
        results['gsva'] = run_module("GSVA通路活性评分", gsva_main)
    except Exception as e:
        logger.error(f"模块1导入失败: {e}")
        results['gsva'] = None
    
    # 模块2: GSEA富集分析
    logger.info("\n" + "-" * 70)
    logger.info("  模块2: GSEA富集分析")
    logger.info("-" * 70)
    try:
        from scripts.cuproptosis_gsea import main as gsea_main
        results['gsea'] = run_module("GSEA富集分析", gsea_main)
    except Exception as e:
        logger.error(f"模块2导入失败: {e}")
        results['gsea'] = None
    
    # 模块3: 单细胞铜死亡评分 (可选，无数据时跳过)
    logger.info("\n" + "-" * 70)
    logger.info("  模块3: 单细胞铜死亡评分")
    logger.info("-" * 70)
    sc_data_path = os.path.join(RESULTS_DIR, "stage2_single_cell", "sc_adata.h5ad")
    if os.path.exists(sc_data_path):
        try:
            from scripts.cuproptosis_singlecell import main as sc_main
            results['singlecell'] = run_module("单细胞铜死亡评分", sc_main)
        except Exception as e:
            logger.error(f"模块3运行失败: {e}")
            results['singlecell'] = None
    else:
        logger.warning("单细胞数据不存在，跳过模块3")
        results['singlecell'] = None
    
    # 模块4: WGCNA模块关联
    logger.info("\n" + "-" * 70)
    logger.info("  模块4: WGCNA模块关联")
    logger.info("-" * 70)
    try:
        from scripts.cuproptosis_wgcna import main as wgcna_main
        results['wgcna'] = run_module("WGCNA模块关联", wgcna_main)
    except Exception as e:
        logger.error(f"模块4导入失败: {e}")
        results['wgcna'] = None
    
    # 模块5: 免疫浸润相关性
    logger.info("\n" + "-" * 70)
    logger.info("  模块5: 免疫浸润相关性")
    logger.info("-" * 70)
    try:
        from scripts.cuproptosis_immunology import main as imm_main
        results['immunology'] = run_module("免疫浸润相关性", imm_main)
    except Exception as e:
        logger.error(f"模块5导入失败: {e}")
        results['immunology'] = None
    
    # 模块6: PPI邻居差异表达
    logger.info("\n" + "-" * 70)
    logger.info("  模块6: PPI邻居差异表达")
    logger.info("-" * 70)
    try:
        from scripts.cuproptosis_ppi_neighbors import main as ppi_main
        results['ppi'] = run_module("PPI邻居差异表达", ppi_main)
    except Exception as e:
        logger.error(f"模块6导入失败: {e}")
        results['ppi'] = None
    
    # 模块7: Hallmark GSVA全通路分析
    logger.info("\n" + "-" * 70)
    logger.info("  模块7: Hallmark GSVA全通路分析")
    logger.info("-" * 70)
    try:
        from scripts.cuproptosis_hallmark_gsva import main as hallmark_main
        results['hallmark_gsva'] = run_module("Hallmark GSVA全通路分析", hallmark_main)
    except Exception as e:
        logger.error(f"模块7导入失败: {e}")
        results['hallmark_gsva'] = None
    
    # 汇总
    total_time = time.time() - start_time
    
    logger.info("\n" + "=" * 70)
    logger.info("  铜死亡模块分析管线完成!")
    logger.info(f"  总耗时: {total_time:.1f}秒 ({total_time/60:.1f}分钟)")
    logger.info("  模块运行状态:")
    for module, result in results.items():
        status = "成功" if result is not None else "失败"
        logger.info(f"    {module:15s}: {status}")
    logger.info("=" * 70)
    
    # v12 FIX: 跨模块一致性检验 (检查清单2.10要求)
    _run_cross_module_validation(results)
    
    return results


def _run_cross_module_validation(results: Dict):
    """
    v12 FIX: 跨模块一致性检验 (检查清单2.10要求)
    
    验证铜死亡基因在不同模块中的表现一致性
    """
    logger.info("\n" + "=" * 70)
    logger.info("跨模块一致性检验 (检查清单2.10)")
    logger.info("=" * 70)
    
    import pandas as pd
    from scipy.stats import spearmanr
    
    # 1. GSVA vs GSEA一致性
    gsva_file = os.path.join(RESULTS_DIR, "cuproptosis_gsva", "cuproptosis_gsva_stats.csv")
    gsea_file = os.path.join(RESULTS_DIR, "cuproptosis_gsea", "cuproptosis_gsea_summary.csv")
    
    if os.path.exists(gsva_file) and os.path.exists(gsea_file):
        try:
            gsva_stats = pd.read_csv(gsva_file)
            gsea_stats = pd.read_csv(gsea_file)
            
            logger.info(f"  GSVA结果: {gsva_stats.shape[0]} 组比较")
            logger.info(f"  GSEA结果: {gsea_stats.shape[0]} 基因集")
            
            if 'logFC' in gsva_stats.columns and 'NES' in gsea_stats.columns:
                gsea_stats['NES_num'] = pd.to_numeric(gsea_stats['NES'], errors='coerce')
                gsva_neg_logfc = (gsva_stats['logFC'] < 0).sum()
                gsea_neg_nes = (gsea_stats['NES_num'] < 0).sum()
                logger.info(f"  GSVA负向logFC: {gsva_neg_logfc}/{len(gsva_stats)}")
                logger.info(f"  GSEA负向NES: {gsea_neg_nes}/{len(gsea_stats)}")
        except Exception as e:
            logger.warning(f"  GSVA-GSEA一致性检验失败: {e}")
    
    # 2. WGCNA富集 vs PPI中心性重叠检验
    wgcna_file = os.path.join(RESULTS_DIR, "cuproptosis_wgcna", "cuproptosis_module_enrichment.csv")
    ppi_file = os.path.join(RESULTS_DIR, "cuproptosis_ppi", "sig_neighbor_degs.csv")
    
    if os.path.exists(wgcna_file) and os.path.exists(ppi_file):
        try:
            wgcna_enrich = pd.read_csv(wgcna_file)
            ppi_degs = pd.read_csv(ppi_file)
            
            wgcna_genes = set()
            if 'Cuproptosis_Genes' in wgcna_enrich.columns:
                for genes_str in wgcna_enrich['Cuproptosis_Genes'].dropna():
                    wgcna_genes.update(str(genes_str).split(','))
            
            ppi_genes = set()
            if 'Gene' in ppi_degs.columns:
                ppi_genes = set(ppi_degs['Gene'].dropna().str.upper())
            
            overlap = wgcna_genes & ppi_genes
            logger.info(f"  WGCNA富集基因: {len(wgcna_genes)}")
            logger.info(f"  PPI显著基因: {len(ppi_genes)}")
            logger.info(f"  重叠基因: {len(overlap)}")
            
            if len(overlap) > 0:
                logger.info(f"  重叠基因: {', '.join(sorted(overlap)[:10])}")
        except Exception as e:
            logger.warning(f"  WGCNA-PPI重叠检验失败: {e}")
    
    # 3. 免疫浸润与铜死亡评分相关性
    immuno_file = os.path.join(RESULTS_DIR, "cuproptosis_immunology", "immune_cytokine_correlation.csv")
    
    if os.path.exists(immuno_file):
        try:
            immuno_corr = pd.read_csv(immuno_file)
            fdr_col = [c for c in immuno_corr.columns if 'FDR' in c.upper()]
            
            if fdr_col:
                sig_corrs = immuno_corr[immuno_corr[fdr_col[0]] < 0.05]
                
                if len(sig_corrs) > 0:
                    logger.info(f"  免疫-铜死亡显著相关对: {len(sig_corrs)}")
                    r_col = [c for c in immuno_corr.columns if 'R' in c.upper() or 'CORR' in c.upper()]
                    if r_col:
                        for _, row in sig_corrs.head(5).iterrows():
                            corr_type = row.get('Cell_Type', row.get('Gene', 'Unknown'))
                            r_val = row[r_col[0]]
                            logger.info(f"    {corr_type}: R={r_val:.3f}")
        except Exception as e:
            logger.warning(f"  免疫相关性检验失败: {e}")
    
    logger.info("=" * 70)


if __name__ == "__main__":
    main()
