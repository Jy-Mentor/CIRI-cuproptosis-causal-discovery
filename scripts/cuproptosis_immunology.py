# -*- coding: utf-8 -*-
"""模块5: 免疫浸润与炎症因子分析 (V4修复版)"""

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

STAGE_DIR = os.path.join(RESULTS_DIR, "cuproptosis_immunology")
ensure_dir(STAGE_DIR)
logger = setup_logger("cuproptosis_immuno_v4", os.path.join(STAGE_DIR, "cuproptosis_immunology_v4.log"))

CUPROPTOSIS_GENES = ["FDX1","LIAS","LIPT1","DLAT","PDHA1","PDHB","MTF1","GLS",
                     "CDKN2A","SLC31A1","ATP7A","ATP7B","DLD","DBT","DLST","PDHA2","GCSH"]


# ============================================================
# FIX P2: 免疫浸润 — 多方法ssGSEA, 不强制归一化
# ============================================================

def estimate_immune_infiltration_v4(expr: pd.DataFrame,
                                     cell_markers: Optional[Dict] = None,
                                     species: Optional[str] = None,
                                     methods: List[str] = None) -> Dict[str, pd.DataFrame]:
    """
    多方法免疫浸润估算 (V4新增)
    methods: 优先级列表, e.g. ['xcell', 'ssgsea', 'mean']
    FIX P2: 保留原始ssGSEA分数, 不强制归一化为和=1
    Returns: {method_name: scores_df}
    """
    if methods is None:
        methods = ['xcell', 'ssgsea', 'mean']

    if species is None:
        species = _detect_species(expr.index)

    if cell_markers is None:
        cell_markers = _get_markers_by_species(species)

    results = {}

    for method in methods:
        try:
            if method == 'xcell':
                scores = _run_xcell(expr)
            elif method == 'ssgsea':
                scores = _run_ssgsea_infiltration(expr, cell_markers)
            elif method == 'mean':
                scores = _run_mean_infiltration(expr, cell_markers)
            else:
                continue

            if scores is not None and not scores.empty:
                results[method] = scores
                logger.info(f"[免疫浸润] {method}: {scores.shape[1]} 种细胞类型")
                break
        except Exception as e:
            logger.warning(f"[免疫浸润] {method} failed: {e}")
            continue

    if not results:
        logger.error("[免疫浸润] 所有方法均失败")

    return results


def _run_xcell(expr: pd.DataFrame) -> Optional[pd.DataFrame]:
    """xCell enrichment-based deconvolution (ssGSEA-based)"""
    try:
        import gseapy as gp
        xcell_genesets = _get_xcell_genesets()
        if not xcell_genesets:
            raise ValueError("xCell gene sets not available")

        res = gp.ssgsea(data=expr, gene_sets=xcell_genesets,
                       outdir=None, no_plot=True, seed=42)
        scores = res.res2d.set_index('Name').T if hasattr(res, 'res2d') else pd.DataFrame()
        return scores.fillna(0)
    except Exception as e:
        raise RuntimeError(f"xCell failed: {e}")


def _run_ssgsea_infiltration(expr: pd.DataFrame, cell_markers: Dict) -> pd.DataFrame:
    """ssGSEA-based免疫浸润 (V4首选)"""
    import gseapy as gp

    valid_sets = {}
    for cell_type, markers in cell_markers.items():
        matched = [m for m in markers if m in expr.index]
        if len(matched) >= 5:
            valid_sets[cell_type] = matched

    if not valid_sets:
        raise ValueError("No valid marker sets")

    res = gp.ssgsea(data=expr, gene_sets=valid_sets,
                   outdir=None, no_plot=True, seed=42)

    if hasattr(res, 'res2d') and res.res2d is not None:
        scores = res.res2d.set_index('Name').T
        scores.columns = scores.columns.str.strip()
        return scores.fillna(0)
    raise RuntimeError("ssGSEA returned invalid format")


def _run_mean_infiltration(expr: pd.DataFrame, cell_markers: Dict) -> pd.DataFrame:
    """简单均值法 (V4回退)"""
    scores = {}
    for cell_type, markers in cell_markers.items():
        matched = [m for m in markers if m in expr.index]
        if matched:
            scores[cell_type] = expr.loc[matched].mean(axis=0)
    return pd.DataFrame(scores).fillna(0)


def _detect_species(expr_index) -> str:
    """根据基因名格式检测物种"""
    import re
    genes = list(expr_index)[:100]
    mouse_n = sum(1 for g in genes if isinstance(g, str) and re.match(r'^[A-Z][a-z0-9_]+$', g))
    human_n = sum(1 for g in genes if isinstance(g, str) and re.match(r'^[A-Z][A-Z0-9_]+$', g))
    total = mouse_n + human_n
    if total == 0:
        return 'human'
    return 'mouse' if (mouse_n / total) > 0.5 else 'human'


def _get_markers_by_species(species: str) -> Dict:
    """获取物种特异性免疫细胞marker"""
    if species == 'mouse':
        return {
            'Microglia': ['Cx3cr1', 'Tmem119', 'P2ry12', 'Csf1r', 'Aif1', 'Trem2'],
            'Macrophage': ['Cd68', 'Cd163', 'Cd14', 'Fcgr1', 'Mrc1'],
            'Neutrophil': ['Fcgr3', 'Cxcr2', 'S100a8', 'S100a9'],
            'T_cell': ['Cd3d', 'Cd3e', 'Cd4', 'Cd8a'],
            'B_cell': ['Cd19', 'Cd79a', 'Ms4a1'],
            'NK_cell': ['Ncam1', 'Nkg7', 'Klrd1'],
            'Astrocyte': ['Gfap', 'Aqp4', 'Slc1a3'],
            'Endothelial': ['Cldn5', 'Pecam1', 'Flt1'],
        }
    return {
        'Microglia': ['CX3CR1', 'TMEM119', 'P2RY12', 'CSF1R', 'AIF1', 'TREM2'],
        'Macrophage': ['CD68', 'CD163', 'CD14', 'FCGR1A', 'MRC1'],
        'Neutrophil': ['FCGR3B', 'CXCR2', 'S100A8', 'S100A9'],
        'T_cell': ['CD3D', 'CD3E', 'CD4', 'CD8A'],
        'B_cell': ['CD19', 'CD79A', 'MS4A1'],
        'NK_cell': ['NCAM1', 'NKG7', 'KLRD1'],
        'Astrocyte': ['GFAP', 'AQP4', 'SLC1A3'],
        'Endothelial': ['CLDN5', 'PECAM1', 'FLT1'],
    }


def _get_xcell_genesets() -> Dict:
    """xCell基因集 (简化版, 实际使用时应加载完整基因集)"""
    return {}


# ============================================================
# 兼容旧接口
# ============================================================

def estimate_immune_infiltration(expr, cell_markers=None, logger=None):
    """V4兼容接口"""
    if logger is None:
        logger = logging.getLogger(__name__)

    try:
        results = estimate_immune_infiltration_v4(expr, cell_markers,
                                                   methods=['xcell', 'ssgsea', 'mean'])
        for method, scores in results.items():
            logger.info(f"  {method}: {scores.shape[1]} 种细胞类型")
            return scores
    except Exception as e:
        logger.warning(f"  V4多方法估算失败 ({e}), 回退到ssGSEA")

    return _run_ssgsea_infiltration(expr, cell_markers or _get_markers_by_species('human'), logger)


def main():
    logger.info("=" * 60)
    logger.info("模块5 (V4): 免疫浸润与炎症因子分析")
    logger.info("修复: P2 多方法ssGSEA, 不强制归一化")
    logger.info("=" * 60)

    expr_file = os.path.join(RESULTS_DIR, "cuproptosis_gsva", "expr_matrix_gene_symbol.csv")
    if not os.path.exists(expr_file):
        logger.error(f"表达矩阵不存在: {expr_file}")
        return 1

    expr = pd.read_csv(expr_file, index_col=0)
    logger.info(f"[1/2] 加载表达矩阵: {expr.shape}")

    logger.info("[2/2] V4免疫浸润估算...")
    results = estimate_immune_infiltration_v4(expr, methods=['xcell', 'ssgsea', 'mean'])

    for method, scores in results.items():
        scores.to_csv(os.path.join(STAGE_DIR, f"immune_infiltration_{method}_v4.csv"))
        logger.info(f"  {method}: {scores.shape}")

    logger.info(f"\n{'='*60}\n模块5 (V4) 完成!\n{'='*60}")
    return 0


if __name__ == "__main__":
    sys.exit(main())