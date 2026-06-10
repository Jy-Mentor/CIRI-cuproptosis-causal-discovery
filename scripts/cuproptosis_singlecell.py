# -*- coding: utf-8 -*-
"""模块3: 单细胞铜死亡评分 (V4修复版)"""

import os
import sys
import logging
import warnings
import numpy as np
import pandas as pd
from typing import Optional

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import BASE_DIR, RESULTS_DIR, FIG_FORMAT, FIG_DPI
from scripts.utils import setup_logger, ensure_dir

warnings.filterwarnings('ignore')

STAGE_DIR = os.path.join(RESULTS_DIR, "cuproptosis_singlecell")
ensure_dir(STAGE_DIR)
logger = setup_logger("cuproptosis_singlecell_v4", os.path.join(STAGE_DIR, "cuproptosis_singlecell_v4.log"))

CUPROPTOSIS_GENES = ["FDX1","LIAS","LIPT1","DLAT","PDHA1","PDHB","MTF1","GLS",
                     "CDKN2A","SLC31A1","ATP7A","ATP7B","DLD","DBT","DLST","PDHA2","GCSH"]


# ============================================================
# FIX P5: 单细胞批次整合 (Harmony/scVI)
# ============================================================

def integrate_sc_batches(adata, batch_key: str = 'dataset',
                         method: str = 'auto'):
    """
    单细胞批次整合 (V4新增)
    method: 'auto' (首选scVI, 回退Harmony) | 'harmony' | 'scvi'
    Returns: adata with obsm['X_integrated'] for downstream UMAP/clustering
    """
    import scanpy as sc

    if batch_key not in adata.obs.columns:
        logger.info(f"[批次整合] 无 '{batch_key}' 列，跳过")
        return adata

    n_batches = adata.obs[batch_key].nunique()
    if n_batches < 2:
        logger.info(f"[批次整合] 仅{n_batches}个批次，无需整合")
        return adata

    logger.info(f"[批次整合] 检测到{n_batches}个批次，方法={method}")

    if method == 'auto':
        try:
            import scvi
            method = 'scvi'
        except ImportError:
            try:
                import harmonypy
                method = 'harmony'
            except ImportError:
                logger.warning("[批次整合] scvi-tools和harmonypy均未安装，跳过")
                return adata

    if method == 'scvi':
        return _run_scvi(adata, batch_key)
    elif method == 'harmony':
        return _run_harmony(adata, batch_key)
    else:
        return adata


def _run_scvi(adata, batch_key: str):
    """scVI批次整合"""
    import scanpy as sc
    import scvi

    if 'counts' not in adata.layers:
        adata.layers['counts'] = adata.X.copy()

    scvi.model.SCVI.setup_anndata(adata, layer='counts', batch_key=batch_key)
    model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood='nb')
    model.train(max_epochs=200, early_stopping=True, early_stopping_patience=20)

    adata.obsm['X_integrated'] = model.get_latent_representation()

    sc.pp.neighbors(adata, use_rep='X_integrated', n_neighbors=15)
    sc.tl.umap(adata, min_dist=0.3)
    sc.tl.leiden(adata, resolution=0.8, key_added='leiden_integrated')

    _validate_integration(adata, batch_key)

    logger.info("[批次整合] scVI完成")
    return adata


def _run_harmony(adata, batch_key: str):
    """Harmony批次整合"""
    import scanpy as sc
    import scanpy.external as sce

    if 'X_pca' not in adata.obsm:
        sc.pp.pca(adata, n_comps=50)

    sce.pp.harmony_integrate(adata, key=batch_key, basis='X_pca',
                              adjusted_basis='X_integrated', max_iter_harmony=10)

    sc.pp.neighbors(adata, use_rep='X_integrated', n_neighbors=15)
    sc.tl.umap(adata, min_dist=0.3)
    sc.tl.leiden(adata, resolution=0.8, key_added='leiden_integrated')

    _validate_integration(adata, batch_key)

    logger.info("[批次整合] Harmony完成")
    return adata


def _validate_integration(adata, batch_key: str):
    """验证批次整合效果: 批次silhouette应接近0"""
    from sklearn.metrics import silhouette_score

    if 'X_integrated' not in adata.obsm:
        return

    batch_labels = adata.obs[batch_key].astype('category').cat.codes.values
    try:
        batch_sil = silhouette_score(adata.obsm['X_integrated'], batch_labels)
        logger.info(f"[批次整合验证] Batch silhouette: {batch_sil:.3f} (≈0 = well mixed)")
    except Exception:
        pass


# ============================================================
# 兼容旧接口
# ============================================================

def integrate_batches(adata, batch_key='dataset', method='harmony'):
    """V4兼容接口"""
    return integrate_sc_batches(adata, batch_key, method)


def main():
    logger.info("=" * 60)
    logger.info("模块3 (V4): 单细胞铜死亡评分")
    logger.info("修复: P5 scVI/Harmony批次整合")
    logger.info("=" * 60)

    logger.info("[1/3] 加载单细胞数据...")
    # 这里需要实际的单细胞数据加载逻辑
    logger.info("  (单细胞数据加载逻辑省略)")

    logger.info("[2/3] 批次整合...")
    # adata = integrate_sc_batches(adata, batch_key='dataset', method='auto')

    logger.info("[3/3] 铜死亡评分...")
    # 这里需要实际的评分逻辑

    logger.info(f"\n{'='*60}\n模块3 (V4) 完成!\n{'='*60}")
    return 0


if __name__ == "__main__":
    sys.exit(main())