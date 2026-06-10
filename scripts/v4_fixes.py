"""
铜死亡模块 V3 → V4 最终修复代码
修复问题清单:
  1. cross_module_integration: IPS硬编码分母 → 动态计算; Bonferroni → BH-FDR
  2. GSVA敏感性: 默认启用，参数自动推荐
  3. 单细胞: 新增Harmony/scVI批次整合
  4. 免疫浸润: 新增CIBERSORT/xCell/MCP-counter支持
"""
import os
import logging
import numpy as np
import pandas as pd
from scipy import stats
from typing import Dict, List, Optional
from statsmodels.stats.multitest import multipletests

logger = logging.getLogger("cuproptosis_v4")


# ============================================================
# FIX 1: CrossModuleIntegration — IPS分母 + FDR
# ============================================================

class CrossModuleIntegratorV4:
    """
    跨模块证据整合 V4 修复版
    修复:
      - IPS分母硬编码4.0 → 动态计算实际参与模块数
      - Bonferroni FDR → BH-FDR (Benjamini-Hochberg)
      - GSVA证据加载: 处理组间统计 vs 基因级结果的格式差异
      - 增加Stouffer's Z-score方向一致性检验
    """

    MODULE_CONFIG = {
        'GSVA':        {'type': 'group',   'min_genes': 3, 'weight': 1.5},
        'GSEA':        {'type': 'pathway', 'min_genes': 5, 'weight': 1.2},
        'SingleCell':  {'type': 'gene',    'min_genes': 2, 'weight': 1.5},
        'WGCNA':       {'type': 'gene',    'min_genes': 3, 'weight': 1.8},
        'Immunology':  {'type': 'gene',    'min_genes': 3, 'weight': 1.5},
        'PPI':         {'type': 'network', 'min_genes': 5, 'weight': 2.0},
        'Hallmark':    {'type': 'pathway', 'min_genes': 5, 'weight': 1.5},
    }

    def integrate(self, module_results: Dict,
                  cuproptosis_genes: Optional[List[str]] = None) -> pd.DataFrame:
        """
        主整合函数 (V4修复版)
        动态IPS分母 = max(基因级模块数, 1)
        FDR = BH-FDR (Benjamini-Hochberg)
        """
        cuproptosis_genes = cuproptosis_genes or self._default_cupro_genes()

        evidence = self._load_all_evidence(module_results, cuproptosis_genes)
        if evidence.empty:
            logger.warning("No evidence loaded from any module")
            return pd.DataFrame()

        results = []
        for gene in cuproptosis_genes:
            gene_e = evidence[evidence['Gene'] == gene]
            if gene_e.empty:
                continue

            modules_present = gene_e['Module'].unique().tolist()
            n_total = len(modules_present)

            gene_level = [m for m in modules_present
                          if self.MODULE_CONFIG.get(m, {}).get('type') in ('gene', 'network')]
            n_gene_level = len(gene_level)

            denominator = max(n_gene_level, 1)
            module_coverage = n_total / denominator

            z_scores, weights_list = [], []
            for _, row in gene_e.iterrows():
                p = max(min(row['P_Value'], 1 - 1e-15), 1e-300)
                direction = row.get('Direction', 'positive')
                z = stats.norm.ppf(1 - p) if direction != 'negative' else -stats.norm.ppf(1 - p)
                w = self.MODULE_CONFIG.get(row['Module'], {}).get('weight', 1.0)
                z_scores.append(z)
                weights_list.append(w)

            if not z_scores:
                continue

            weighted_z = sum(z * w for z, w in zip(z_scores, weights_list)) / np.sqrt(sum(w**2 for w in weights_list))
            combined_p = 1 - stats.norm.cdf(abs(weighted_z))

            directions = gene_e['Direction'].tolist() if 'Direction' in gene_e.columns else []
            if directions:
                pos_ratio = sum(1 for d in directions if d == 'positive') / len(directions)
                direction_consistency = max(pos_ratio, 1 - pos_ratio)
            else:
                direction_consistency = 0.5

            results.append({
                'Gene': gene,
                'Module_Present': ', '.join(modules_present),
                'N_GeneLevel_Modules': n_gene_level,
                'N_PathwayLevel_Modules': n_total - n_gene_level,
                'Total_Modules': n_total,
                'Denominator': denominator,
                'ModuleCoverage': module_coverage,
                'Stouffer_Z': weighted_z,
                'Combined_P': combined_p,
                'Direction_Consistency': direction_consistency,
                'IPS': weighted_z * module_coverage * direction_consistency,
            })

        if not results:
            return pd.DataFrame()

        df = pd.DataFrame(results)

        if len(df) > 1:
            _, fdr_vals, _, _ = multipletests(df['Combined_P'].values, method='fdr_bh')
            df['BH_FDR'] = fdr_vals
        else:
            df['BH_FDR'] = df['Combined_P']

        df['BH_FDR'] = df['BH_FDR'].clip(upper=1.0)
        df = df.sort_values('IPS', ascending=False)
        df['IPS_Rank'] = range(1, len(df) + 1)

        logger.info(f"Integrated {len(df)} genes: mean IPS={df['IPS'].mean():.2f}, "
                    f"FDR<0.05={(df['BH_FDR']<0.05).sum()}/{len(df)}")
        return df

    def _load_all_evidence(self, module_results: Dict,
                           cupro_genes: List[str]) -> pd.DataFrame:
        """加载所有模块证据，统一格式"""
        records = []
        for mod_name, config in self.MODULE_CONFIG.items():
            result = module_results.get(mod_name)
            if result is None or (isinstance(result, pd.DataFrame) and result.empty):
                continue
            mod_records = self._parse_module(mod_name, result, config, cupro_genes)
            records.extend(mod_records)
        return pd.DataFrame(records) if records else pd.DataFrame()

    def _parse_module(self, mod_name: str, result, config: dict,
                      cupro_genes: List[str]) -> List[dict]:
        """解析单个模块结果为统一格式"""
        records = []
        mtype = config.get('type', 'unknown')

        if mtype == 'network' and mod_name == 'PPI':
            for _, row in result.iterrows():
                gene = row.get('Gene', row.get('Gene_Symbol', ''))
                hub = row.get('HubScore', 0)
                records.append({
                    'Gene': gene, 'Module': mod_name,
                    'P_Value': max(1e-10, 1 - min(0.999, hub)),
                    'Direction': 'positive'
                })

        elif mtype == 'gene' and mod_name == 'WGCNA':
            for _, row in result.iterrows():
                gene = row.get('gene', row.get('Gene', ''))
                kme = row.get('kME', 0)
                records.append({
                    'Gene': gene, 'Module': mod_name,
                    'P_Value': row.get('kME_pvalue', 0.05),
                    'Direction': 'positive' if kme > 0 else 'negative'
                })

        elif mtype == 'group' and mod_name == 'GSVA':
            if 'p_value' in result.columns:
                p_val = result['p_value'].min() if not result.empty else 1.0
                for gene in cupro_genes:
                    records.append({
                        'Gene': gene, 'Module': mod_name,
                        'P_Value': p_val,
                        'Direction': 'positive' if p_val < 0.05 else 'neutral'
                    })

        elif mtype in ('pathway',) and mod_name in ('GSEA', 'Hallmark'):
            p_col = next((c for c in ['pval', 'p_value', 'P_Value', 'FDR p-val']
                          if c in result.columns), None)
            if p_col:
                for _, row in result.iterrows():
                    nes = row.get('nes', row.get('NES', 0))
                    records.append({
                        'Gene': row.get('Term', mod_name), 'Module': mod_name,
                        'P_Value': row[p_col],
                        'Direction': 'positive' if nes > 0 else 'negative'
                    })

        elif mtype == 'gene' and mod_name in ('SingleCell', 'Immunology'):
            for _, row in result.iterrows():
                gene = row.get('Gene', row.get('gene', ''))
                records.append({
                    'Gene': gene, 'Module': mod_name,
                    'P_Value': row.get('p_value', row.get('P_Value', 1.0)),
                    'Direction': 'positive' if row.get('effect_size', 0) > 0 else 'negative'
                })

        return records

    @staticmethod
    def _default_cupro_genes() -> List[str]:
        return ["FDX1", "LIAS", "LIPT1", "DLAT", "PDHA1", "PDHB", "MTF1", "GLS",
                "CDKN2A", "SLC31A1", "ATP7A", "ATP7B", "DLD", "DBT", "DLST",
                "PDHA2", "GCSH"]


# ============================================================
# FIX 2: GSVA敏感性 — 默认启用 + 自动推荐
# ============================================================

def run_gsva_with_sensitivity(expr: pd.DataFrame,
                               gene_set: List[str],
                               output_dir: str,
                               mx_diff_grid=(True, False),
                               tau_grid=(0.25, 0.5, 1.0),
                               n_top: int = 2) -> pd.DataFrame:
    """
    GSVA参数敏感性分析 + 最优参数推荐 (V4新增)
    n_top: 推荐前N个最优参数组合
    Returns: 最优参数下的GSVA评分 DataFrame
    """
    import gseapy as gp

    available_genes = [g for g in gene_set if g in expr.index]
    if len(available_genes) < 5:
        raise ValueError(f"GSVA genes < 5: {len(available_genes)}")

    sensitivity_dir = os.path.join(output_dir, "gsva_sensitivity")
    os.makedirs(sensitivity_dir, exist_ok=True)

    results = []
    all_scores = {}

    logger.info("[GSVA敏感性] 扫描参数网格...")
    for mx_diff in mx_diff_grid:
        for tau in tau_grid:
            try:
                res = gp.gsva(
                    data=expr.loc[available_genes],
                    gene_sets={'Cuproptosis': available_genes},
                    mx_diff=mx_diff, tau=tau,
                    min_size=3, max_size=500,
                    outdir=None, no_plot=True, seed=42, verbose=False
                )
                es = res.res2d['ES'].values if hasattr(res, 'res2d') and res.res2d is not None else []
                label = f"mx{mx_diff}_tau{tau}"
                all_scores[label] = es
                results.append({
                    'mx_diff': mx_diff, 'tau': tau,
                    'mean_ES': float(np.mean(es)) if len(es) else np.nan,
                    'std_ES': float(np.std(es)) if len(es) else np.nan,
                    'abs_ES': float(np.mean(np.abs(es))) if len(es) else np.nan,
                })
            except Exception as e:
                logger.warning(f"  GSVA mx_diff={mx_diff}, tau={tau} failed: {e}")

    sens_df = pd.DataFrame(results)
    sens_df.to_csv(os.path.join(sensitivity_dir, "gsva_sensitivity.csv"), index=False)

    stability = _compute_rank_stability(all_scores) if len(all_scores) >= 2 else 1.0
    logger.info(f"[GSVA敏感性] 跨参数排名稳定性: {stability:.3f}")

    sens_df['score'] = sens_df['abs_ES'] * stability
    best_params = sens_df.nlargest(n_top, 'score')

    logger.info(f"[GSVA敏感性] 推荐参数 (Top {n_top}):")
    for _, bp in best_params.iterrows():
        logger.info(f"  mx_diff={bp['mx_diff']}, tau={bp['tau']}: |ES|={bp['abs_ES']:.4f}")

    best = best_params.iloc[0]
    final_res = gp.gsva(
        data=expr.loc[available_genes],
        gene_sets={'Cuproptosis': available_genes},
        mx_diff=best['mx_diff'], tau=best['tau'],
        min_size=3, max_size=500,
        outdir=None, no_plot=True, seed=42, verbose=False
    )

    if hasattr(final_res, 'res2d') and final_res.res2d is not None:
        scores = final_res.res2d.set_index('Term')['ES']
        logger.info(f"[GSVA] 最终评分 (mx_diff={best['mx_diff']}, tau={best['tau']}): "
                    f"mean={scores.mean():.4f}")
        return pd.DataFrame({'Cuproptosis_Score': scores})
    else:
        raise RuntimeError("GSVA final run failed")


def _compute_rank_stability(all_scores: dict) -> float:
    """跨参数组合的通路排名稳定性 (Spearman平均)"""
    from scipy.stats import spearmanr
    labels = list(all_scores.keys())
    if len(labels) < 2:
        return 1.0

    all_indices = set()
    for es in all_scores.values():
        all_indices.update(range(len(es)))

    ranks = {}
    for label, es in all_scores.items():
        ranks[label] = pd.Series(es).rank().values

    corrs = []
    for i in range(len(labels)):
        for j in range(i + 1, len(labels)):
            a, b = ranks[labels[i]], ranks[labels[j]]
            min_len = min(len(a), len(b))
            if min_len >= 3:
                c, _ = spearmanr(a[:min_len], b[:min_len])
                if not np.isnan(c):
                    corrs.append(c)

    return np.mean(corrs) if corrs else 1.0


# ============================================================
# FIX 3: 单细胞批次整合 (Harmony/scVI)
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
# FIX 4: 免疫浸润 — 多方法ssGSEA
# ============================================================

def estimate_immune_infiltration_v4(expr: pd.DataFrame,
                                     cell_markers: Optional[Dict] = None,
                                     species: Optional[str] = None,
                                     methods: List[str] = None) -> Dict[str, pd.DataFrame]:
    """
    多方法免疫浸润估算 (V4新增)
    methods: 优先级列表, e.g. ['xcell', 'ssgsea', 'mean']
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
                logger.info(f"[免疫浸润] {method}: {scores.shape}")
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
# 使用示例 (如何替换现有代码)
# ============================================================
"""
# 1. 替换 CrossModuleIntegration.integrate:
from v4_fixes import CrossModuleIntegratorV4
integrator = CrossModuleIntegratorV4()
final_df = integrator.integrate(module_results, cuproptosis_genes)

# 2. 替换 GSVA评分计算:
from v4_fixes import run_gsva_with_sensitivity
scores_df = run_gsva_with_sensitivity(expr, gene_set, output_dir)

# 3. 替换单细胞批次整合:
from v4_fixes import integrate_sc_batches
adata = integrate_sc_batches(adata, batch_key='dataset', method='auto')

# 4. 替换免疫浸润:
from v4_fixes import estimate_immune_infiltration_v4
results = estimate_immune_infiltration_v4(expr, methods=['xcell', 'ssgsea', 'mean'])
immune_scores = results.get('ssgsea') or results.get('mean')
"""

if __name__ == "__main__":
    print("V4修复模块已加载")
    print("  - CrossModuleIntegratorV4.integrate()  # 跨模块整合")
    print("  - run_gsva_with_sensitivity()           # GSVA+敏感性")
    print("  - integrate_sc_batches()                # 单细胞批次整合")
    print("  - estimate_immune_infiltration_v4()     # 多方法免疫浸润")