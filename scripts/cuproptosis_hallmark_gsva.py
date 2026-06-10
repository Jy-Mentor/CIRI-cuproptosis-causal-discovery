# -*- coding: utf-8 -*-
"""模块7: Hallmark通路GSVA分析 (V4修复版)"""

import os
import sys
import logging
import warnings
import numpy as np
import pandas as pd
from scipy import stats
from typing import Dict, List, Optional

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import BASE_DIR, RESULTS_DIR, FIG_FORMAT, FIG_DPI
from scripts.utils import setup_logger, ensure_dir

warnings.filterwarnings('ignore')

STAGE_DIR = os.path.join(RESULTS_DIR, "cuproptosis_hallmark_gsva")
ensure_dir(STAGE_DIR)
logger = setup_logger("cuproptosis_hallmark_v4", os.path.join(STAGE_DIR, "cuproptosis_hallmark_v4.log"))

CUPROPTOSIS_GENES = ["FDX1","LIAS","LIPT1","DLAT","PDHA1","PDHB","MTF1","GLS",
                     "CDKN2A","SLC31A1","ATP7A","ATP7B","DLD","DBT","DLST","PDHA2","GCSH"]


# ============================================================
# FIX P4: Hallmark铜死亡评分使用GSVA而非简单均值
# ============================================================

def compute_cuproptosis_gsva_score(expr: pd.DataFrame,
                                    cupro_genes: List[str],
                                    logger=None) -> pd.Series:
    """
    使用GSVA计算铜死亡通路评分 (替代简单均值)
    优势: 考虑基因的相对排序, 非参数化, 对异常值更稳健
    依据: Hanzelmann S, et al. BMC Bioinformatics 2013
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    matched = [g for g in cupro_genes if g in expr.index]
    if len(matched) < 5:
        logger.info(f"  铜死亡基因不足5个 ({len(matched)}), 回退到均值法")
        return expr.loc[matched].mean(axis=0) if matched else pd.Series(dtype=float)

    try:
        import gseapy as gp
        gsva_result = gp.gsva(
            data=expr.loc[matched],
            gene_sets={'Cuproptosis': matched},
            mx_diff=True, tau=1.0,
            min_size=3, max_size=500,
            outdir=None, no_plot=True, seed=42, verbose=False
        )
        if isinstance(gsva_result, pd.DataFrame):
            scores = gsva_result.loc['Cuproptosis'] if 'Cuproptosis' in gsva_result.index else gsva_result.iloc[0]
            logger.info(f"  铜死亡评分 (GSVA): mean={scores.mean():.3f}")
            return scores
        elif hasattr(gsva_result, 'res2d') and gsva_result.res2d is not None:
            scores = gsva_result.res2d.set_index('Term')['ES']
            if 'Cuproptosis' in scores.index:
                logger.info(f"  铜死亡评分 (GSVA): mean={scores.mean():.3f}")
                return scores
    except Exception as e:
        logger.warning(f"  GSVA铜死亡评分失败 ({e}), 回退到均值法")

    return expr.loc[matched].mean(axis=0)


def _get_builtin_hallmark_sets() -> Dict[str, List[str]]:
    """内置Hallmark基因集 (简化版, 优先从gseapy加载完整版)"""
    return {
        'HALLMARK_TNFA_SIGNALING_VIA_NFKB': ['TNF','IL1B','IL6','CXCL8','CXCL1','CXCL2','CCL2','PTGS2','ICAM1','VCAM1','SELE','NFKB1','RELA','NFKBIA'],
        'HALLMARK_HYPOXIA': ['HIF1A','VEGFA','EGLN1','EGLN3','CA9','SLC2A1','PGK1','ENO1','LDHA','BNIP3','NOS2'],
        'HALLMARK_APOPTOSIS': ['BAX','BAK1','CASP3','CASP8','CASP9','FAS','FASLG','BCL2','BCL2L1','BAD','BID','BIM','PUMA','NOXA'],
        'HALLMARK_INFLAMMATORY_RESPONSE': ['IL1B','IL6','TNF','CCL2','CCL5','CXCL10','CXCL8','PTGS2','ICAM1','VCAM1','SELE','NFKB1'],
        'HALLMARK_OXIDATIVE_PHOSPHORYLATION': ['NDUFA1','NDUFA2','NDUFA3','NDUFA4','NDUFA5','NDUFA6','NDUFA7','NDUFA8','NDUFA9','NDUFA10','SDHA','SDHB','SDHC','SDHD','UQCRC1','UQCRC2','COX4I1','COX5A','COX5B','COX6A1','COX6B1','COX7A1','ATP5F1A','ATP5F1B','ATP5F1C','ATP5F1D'],
        'HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY': ['SOD1','SOD2','SOD3','CAT','GPX1','GPX2','GPX3','GPX4','PRDX1','PRDX2','PRDX3','PRDX4','PRDX5','PRDX6'],
        'HALLMARK_P53_PATHWAY': ['TP53','CDKN1A','MDM2','GADD45A','BAX','FAS','CASP3','CASP8','CASP9','PERP','SESN1','SESN2'],
        'HALLMARK_PI3K_AKT_MTOR_SIGNALING': ['AKT1','AKT2','AKT3','MTOR','RPTOR','RICTOR','PIK3CA','PIK3CB','PIK3CD','PIK3R1','PIK3R2','PTEN','TSC1','TSC2','RHEB','RPS6KB1','EIF4EBP1'],
        'HALLMARK_E2F_TARGETS': ['E2F1','E2F2','E2F3','E2F4','E2F5','E2F6','E2F7','E2F8','CCNA2','CCNB1','CCNB2','CCNE1','CCNE2','CDK1','CDK2','MCM2','MCM3','MCM4','MCM5','MCM6','MCM7','PCNA','RFC1','RFC2','RFC3','RFC4','RFC5'],
        'HALLMARK_MYC_TARGETS_V1': ['MYC','MYCN','MYCL','NPM1','NCL','FBL','NOP56','NOP58','RRP1','RRP9','UTP6','UTP15','BOP1','WDR12','WDR43','WDR75'],
    }


def _load_hallmark_sets() -> Dict[str, List[str]]:
    """加载Hallmark基因集, 优先gseapy完整版"""
    try:
        import gseapy as gp
        hallmark = gp.get_library(name='MSigDB_Hallmark_2020')
        if hallmark and len(hallmark) >= 10:
            logger.info(f"  从gseapy加载Hallmark: {len(hallmark)} 条通路")
            return hallmark
    except Exception:
        pass
    logger.info("  使用内置Hallmark基因集 (10条)")
    return _get_builtin_hallmark_sets()


def run_hallmark_gsva(expr: pd.DataFrame,
                      hallmark_sets: Dict[str, List[str]] = None,
                      output_dir: str = None,
                      logger=None) -> pd.DataFrame:
    """V4 Hallmark GSVA主函数"""
    if logger is None:
        logger = logging.getLogger(__name__)
    if output_dir is None:
        output_dir = STAGE_DIR
    if hallmark_sets is None:
        hallmark_sets = _load_hallmark_sets()

    import gseapy as gp

    # 1. 铜死亡评分 (GSVA)
    logger.info("[Hallmark V4] 铜死亡评分 (GSVA)...")
    cupro_scores = compute_cuproptosis_gsva_score(expr, CUPROPTOSIS_GENES, logger)

    # 2. Hallmark GSVA
    logger.info(f"[Hallmark V4] 运行GSVA: {len(hallmark_sets)} 条通路...")
    gsva_result = gp.gsva(
        data=expr,
        gene_sets=hallmark_sets,
        mx_diff=True, tau=1.0,
        min_size=5, max_size=500,
        outdir=None, no_plot=True, seed=42, verbose=False
    )

    if hasattr(gsva_result, 'res2d') and gsva_result.res2d is not None:
        hallmark_scores = gsva_result.res2d.set_index('Term')['ES']
    else:
        hallmark_scores = pd.Series(dtype=float)

    # 3. 相关性分析
    correlations = []
    for pathway in hallmark_scores.index:
        r, p = stats.spearmanr(cupro_scores, hallmark_scores[pathway])
        correlations.append({
            'Pathway': pathway,
            'Spearman_R': r,
            'P_value': p,
            'Significant': p < 0.05,
        })

    corr_df = pd.DataFrame(correlations).sort_values('P_value')
    corr_df.to_csv(os.path.join(output_dir, "cuproptosis_hallmark_correlations_v4.csv"), index=False)

    logger.info(f"[Hallmark V4] 完成: {len(corr_df)} 条通路, {(corr_df['Significant']).sum()} 个显著")
    return corr_df


def main():
    logger.info("=" * 60)
    logger.info("模块7 (V4): Hallmark通路GSVA分析")
    logger.info("修复: P4 铜死亡评分使用GSVA (非简单均值)")
    logger.info("=" * 60)

    expr_file = os.path.join(RESULTS_DIR, "cuproptosis_gsva", "expr_matrix_gene_symbol.csv")
    if not os.path.exists(expr_file):
        logger.error(f"表达矩阵不存在: {expr_file}")
        return 1

    expr = pd.read_csv(expr_file, index_col=0)
    logger.info(f"[1/2] 加载表达矩阵: {expr.shape}")

    logger.info("[2/2] V4 Hallmark GSVA分析...")
    corr_df = run_hallmark_gsva(expr, output_dir=STAGE_DIR, logger=logger)

    logger.info(f"\n{'='*60}\n模块7 (V4) 完成!\n{'='*60}")
    return 0


if __name__ == "__main__":
    sys.exit(main())