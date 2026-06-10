# -*- coding: utf-8 -*-
"""模块8: 跨模块证据融合 (V4修复版)"""

import sys
import os
import logging
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STAGE_DIR = os.path.join(BASE_DIR, "results", "cross_module_integration")
os.makedirs(STAGE_DIR, exist_ok=True)

logger = logging.getLogger(__name__)
FIG_FORMAT = 'png'
FIG_DPI = 300


def setup_logging():
    for h in logger.handlers[:]:
        logger.removeHandler(h)
    logger.setLevel(logging.INFO)
    fh = logging.FileHandler(os.path.join(STAGE_DIR, 'cross_module_integration_v4.log'), mode='w', encoding='utf-8')
    fh.setLevel(logging.INFO)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    fmt = logging.Formatter('%(asctime)s | %(levelname)-8s | %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    fh.setFormatter(fmt)
    ch.setFormatter(fmt)
    logger.addHandler(fh)
    logger.addHandler(ch)


# ============================================================
# FIX P6: Stouffer's Z - 方向由signs参数传入
# ============================================================

def stouffer_z(p_values, weights=None, signs=None):
    """
    Stouffer's Z-score method (V4修复版)
    FIX P6: 方向不再由p值推导 (np.sign(0.5 - valid_p) 总是+1),
    而是由调用方通过 signs 参数传入真实的效应方向 (+1/-1)
    Z = sum(w_i * z_i) / sqrt(sum(w_i^2))
    其中 z_i = Phi^{-1}(1 - p_i/2) * sign_i
    """
    valid_p, valid_w, valid_s = [], [], []

    for i, p in enumerate(p_values):
        if p is None or np.isnan(p) or p <= 0 or p >= 1:
            continue
        valid_p.append(np.clip(p, 1e-300, 1.0))
        valid_w.append(weights[i] if weights is not None else 1.0)
        s = signs[i] if signs is not None and i < len(signs) and signs[i] is not None else 1.0
        valid_s.append(float(s))

    if len(valid_p) == 0:
        return None, None, 0

    valid_p = np.array(valid_p)
    valid_w = np.array(valid_w)
    valid_s = np.array(valid_s)

    z_scores = norm.ppf(1 - valid_p / 2) * valid_s

    w_sum = np.sum(valid_w)
    if w_sum > 0:
        valid_w = valid_w / w_sum

    combined_z = np.sum(valid_w * z_scores) / np.sqrt(np.sum(valid_w**2))
    combined_p = 2 * (1 - norm.cdf(abs(combined_z)))
    combined_p = np.clip(combined_p, 1e-300, 1.0)

    return combined_z, combined_p, len(valid_p)


# ============================================================
# V4 基因级整合器 (IPS动态分母 + BH-FDR)
# ============================================================

class CrossModuleIntegratorV4:
    """V4跨模块整合: 动态IPS分母 + BH-FDR + 方向一致性"""

    MODULE_CONFIG = {
        'GSVA':        {'type': 'group',   'min_genes': 3, 'weight': 1.5},
        'GSEA':        {'type': 'pathway', 'min_genes': 5, 'weight': 1.2},
        'SingleCell':  {'type': 'gene',    'min_genes': 2, 'weight': 1.5},
        'WGCNA':       {'type': 'gene',    'min_genes': 3, 'weight': 1.8},
        'Immunology':  {'type': 'gene',    'min_genes': 3, 'weight': 1.5},
        'PPI':         {'type': 'network', 'min_genes': 5, 'weight': 2.0},
        'Hallmark':    {'type': 'pathway', 'min_genes': 5, 'weight': 1.5},
    }

    def integrate(self, module_results: dict, cuproptosis_genes: list = None) -> pd.DataFrame:
        cuproptosis_genes = cuproptosis_genes or [
            "FDX1","LIAS","LIPT1","DLAT","PDHA1","PDHB","MTF1","GLS",
            "CDKN2A","SLC31A1","ATP7A","ATP7B","DLD","DBT","DLST","PDHA2","GCSH"
        ]

        evidence = self._load_all_evidence(module_results, cuproptosis_genes)
        if evidence.empty:
            logger.warning("无证据加载")
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

            # V4: 动态分母 = max(实际基因级模块数, 1)
            denominator = max(n_gene_level, 1)
            module_coverage = n_total / denominator

            pvals, weights_list, sign_list = [], [], []
            for _, row in gene_e.iterrows():
                p = max(min(row['P_Value'], 1 - 1e-15), 1e-300)
                direction = row.get('Direction', 'positive')
                sign = 1.0 if direction != 'negative' else -1.0
                w = self.MODULE_CONFIG.get(row['Module'], {}).get('weight', 1.0)
                pvals.append(p)
                weights_list.append(w)
                sign_list.append(sign)

            if not pvals:
                continue

            # V4 FIX P6: signs参数传入方向
            z, p_combined, n_src = stouffer_z(pvals, weights_list, sign_list)
            if p_combined is None or p_combined >= 1.0:
                continue

            # 方向一致性
            directions = gene_e['Direction'].tolist() if 'Direction' in gene_e.columns else []
            if directions:
                pos_ratio = sum(1 for d in directions if d == 'positive') / len(directions)
                dir_cons = max(pos_ratio, 1 - pos_ratio)
            else:
                dir_cons = 0.5

            results.append({
                'Gene': gene,
                'Module_Present': ', '.join(modules_present),
                'N_GeneLevel_Modules': n_gene_level,
                'N_PathwayLevel_Modules': n_total - n_gene_level,
                'Total_Modules': n_total,
                'Denominator': denominator,
                'ModuleCoverage': module_coverage,
                'Stouffer_Z': z,
                'Combined_P': p_combined,
                'Direction_Consistency': dir_cons,
                'IPS': z * module_coverage * dir_cons,
            })

        if not results:
            return pd.DataFrame()

        df = pd.DataFrame(results)

        # V4: BH-FDR
        if len(df) > 1:
            _, fdr_vals, _, _ = multipletests(df['Combined_P'].values, method='fdr_bh')
            df['BH_FDR'] = fdr_vals
        else:
            df['BH_FDR'] = df['Combined_P']
        df['BH_FDR'] = df['BH_FDR'].clip(upper=1.0)

        df = df.sort_values('IPS', ascending=False)
        df['IPS_Rank'] = range(1, len(df) + 1)

        n_sig = (df['BH_FDR'] < 0.05).sum()
        logger.info(f"整合完成: {len(df)} 基因, FDR<0.05: {n_sig}/{len(df)}, 平均IPS={df['IPS'].mean():.2f}")
        return df

    def _load_all_evidence(self, module_results, cupro_genes):
        records = []
        for mod_name, config in self.MODULE_CONFIG.items():
            result = module_results.get(mod_name)
            if result is None or (isinstance(result, pd.DataFrame) and result.empty):
                continue
            records.extend(self._parse_module(mod_name, result, config, cupro_genes))
        return pd.DataFrame(records) if records else pd.DataFrame()

    def _parse_module(self, mod_name, result, config, cupro_genes):
        records = []
        mtype = config.get('type', 'unknown')

        if mtype == 'network' and mod_name == 'PPI':
            for _, row in result.iterrows():
                gene = row.get('Gene', row.get('Gene_Symbol', ''))
                hub = row.get('HubScore', 0)
                records.append({'Gene': gene, 'Module': mod_name,
                               'P_Value': max(1e-10, 1 - min(0.999, hub)),
                               'Direction': 'positive'})

        elif mtype == 'gene' and mod_name == 'WGCNA':
            for _, row in result.iterrows():
                gene = row.get('gene', row.get('Gene', ''))
                kme = row.get('kME', 0)
                records.append({'Gene': gene, 'Module': mod_name,
                               'P_Value': row.get('kME_pvalue', 0.05),
                               'Direction': 'positive' if kme > 0 else 'negative'})

        elif mtype == 'group' and mod_name == 'GSVA':
            if 'p_value' in result.columns:
                p_val = result['p_value'].min() if not result.empty else 1.0
                for gene in cupro_genes:
                    records.append({'Gene': gene, 'Module': mod_name,
                                   'P_Value': p_val,
                                   'Direction': 'positive' if p_val < 0.05 else 'neutral'})

        elif mtype == 'pathway' and mod_name in ('GSEA', 'Hallmark'):
            p_col = next((c for c in ['pval','p_value','P_Value','FDR p-val'] if c in result.columns), None)
            if p_col:
                for _, row in result.iterrows():
                    nes = row.get('nes', row.get('NES', 0))
                    records.append({'Gene': row.get('Term', mod_name), 'Module': mod_name,
                                   'P_Value': row[p_col],
                                   'Direction': 'positive' if nes > 0 else 'negative'})

        elif mtype == 'gene' and mod_name in ('SingleCell', 'Immunology'):
            for _, row in result.iterrows():
                gene = row.get('Gene', row.get('gene', ''))
                records.append({'Gene': gene, 'Module': mod_name,
                               'P_Value': row.get('p_value', row.get('P_Value', 1.0)),
                               'Direction': 'positive' if row.get('effect_size', 0) > 0 else 'negative'})
        return records


# ============================================================
# 通路级整合 (沿用Stouffer)
# ============================================================

def integrate_pathway_evidence(evidence: dict) -> pd.DataFrame:
    pathways = evidence.get('pathways', {})
    if not pathways:
        return pd.DataFrame()

    results = []
    for pathway, mod_evidence in pathways.items():
        pvals, weights, sign_list, modules_present = [], [], [], []
        effect_sizes = []

        for mod_name, val in mod_evidence.items():
            if '_pval' in mod_name and isinstance(val, (int, float)) and 0 < val < 1:
                pvals.append(val)
                modules_present.append(mod_name.replace('_pval', ''))
                w = 1.5 if 'Hallmark' in mod_name else 1.2 if 'GSEA' in mod_name else 1.0
                weights.append(w)
                if '_NES' in mod_name.replace('_pval', '_NES'):
                    sign_list.append(1.0)
                elif mod_evidence.get(mod_name.replace('_pval', '_R'), 0) > 0:
                    sign_list.append(1.0)
                elif mod_evidence.get(mod_name.replace('_pval', '_R'), 0) < 0:
                    sign_list.append(-1.0)
                else:
                    sign_list.append(1.0)

                for sfx in ['_NES', '_R', '_ratio']:
                    k = mod_name.replace('_pval', sfx)
                    if k in mod_evidence:
                        effect_sizes.append(abs(mod_evidence[k]))

        if len(pvals) == 0:
            continue

        z, p_combined, n_mod = stouffer_z(pvals, weights, sign_list)
        results.append({
            'Pathway': pathway, 'ModuleCount': len(modules_present),
            'Modules': ','.join(modules_present), 'Stouffer_Z': z,
            'Combined_P': p_combined,
            'Avg_Effect': np.mean(effect_sizes) if effect_sizes else 0,
            'N_Source': len(pvals),
        })

    df = pd.DataFrame(results)
    if not df.empty:
        df = df.sort_values('Combined_P')
        df['Rank'] = range(1, len(df) + 1)
    return df


def load_module_evidence() -> dict:
    """加载各模块证据 (兼容旧格式)"""
    evidence = {'genes': {}, 'pathways': {}}
    RESULTS = os.path.join(BASE_DIR, "results")

    f = os.path.join(RESULTS, "cuproptosis_gsva", "cuproptosis_gsva_stats.csv")
    if os.path.exists(f):
        try:
            df = pd.read_csv(f)
            pval = float(df.iloc[0]['p_value']) if 'p_value' in df.columns else 1.0
            cohens_d = float(df.iloc[0]['cohens_d']) if 'cohens_d' in df.columns else 0
            evidence['pathways']['Cuproptosis_GSVA'] = {'GSVA_pval': pval, 'GSVA_effect': cohens_d}
        except Exception as e:
            logger.warning(f"GSVA加载失败: {e}")

    f = os.path.join(RESULTS, "cuproptosis_gsea", "cuproptosis_gsea_summary.csv")
    if os.path.exists(f):
        try:
            df = pd.read_csv(f)
            for _, row in df.iterrows():
                sn = row['Gene_Set']
                pval = row.get('P_value', 1)
                nes = row.get('NES', 0)
                evidence['pathways'][sn] = evidence['pathways'].get(sn, {})
                evidence['pathways'][sn]['GSEA_pval'] = float(pval) if pd.notna(pval) else 1.0
                evidence['pathways'][sn]['GSEA_NES'] = float(nes) if pd.notna(nes) else 0
        except Exception as e:
            logger.warning(f"GSEA加载失败: {e}")

    f = os.path.join(RESULTS, "cuproptosis_wgcna", "cuproptosis_module_enrichment.csv")
    if os.path.exists(f):
        try:
            df = pd.read_csv(f)
            for _, row in df.iterrows():
                pval = row.get('P_value', 1)
                module = row.get('Module', 'unknown')
                evidence['pathways'][module] = evidence['pathways'].get(module, {})
                evidence['pathways'][module]['WGCNA_pval'] = float(pval) if pd.notna(pval) else 1.0
        except Exception as e:
            logger.warning(f"WGCNA加载失败: {e}")

    f = os.path.join(RESULTS, "cuproptosis_ppi", "ppi_hub_genes.csv")
    if os.path.exists(f):
        try:
            df = pd.read_csv(f)
            for _, row in df.iterrows():
                gene = str(row['Gene']).upper()
                pval = row.get('adj_P_Val', 1)
                score = row.get('HubScore', 0)
                evidence['genes'][gene] = evidence['genes'].get(gene, {})
                evidence['genes'][gene]['PPI_pval'] = float(pval) if pd.notna(pval) else 1.0
                evidence['genes'][gene]['PPI_HubScore'] = float(score) if pd.notna(score) else 0
        except Exception as e:
            logger.warning(f"PPI加载失败: {e}")

    f = os.path.join(RESULTS, "cuproptosis_hallmark_gsva", "cuproptosis_hallmark_correlations.csv")
    if os.path.exists(f):
        try:
            df = pd.read_csv(f)
            for _, row in df.iterrows():
                pathway = str(row['Pathway'])
                pval = row.get('P_value', 1)
                r = row.get('Spearman_R', 0)
                evidence['pathways'][pathway] = evidence['pathways'].get(pathway, {})
                evidence['pathways'][pathway]['Hallmark_pval'] = float(pval) if pd.notna(pval) else 1.0
                evidence['pathways'][pathway]['Hallmark_R'] = float(r) if pd.notna(r) else 0
        except Exception as e:
            logger.warning(f"Hallmark加载失败: {e}")

    return evidence


def _load_module_results() -> dict:
    """加载模块结果为DataFrame格式"""
    RESULTS = os.path.join(BASE_DIR, "results")
    module_results = {}

    for mod, path_key in [
        ('PPI', ("cuproptosis_ppi", "ppi_hub_genes.csv")),
        ('WGCNA', ("cuproptosis_wgcna", "cuproptosis_module_enrichment.csv")),
        ('GSVA', ("cuproptosis_gsva", "cuproptosis_gsva_stats.csv")),
        ('Hallmark', ("cuproptosis_hallmark_gsva", "cuproptosis_hallmark_correlations.csv")),
        ('GSEA', ("cuproptosis_gsea", "cuproptosis_gsea_summary.csv")),
    ]:
        f = os.path.join(RESULTS, path_key[0], path_key[1])
        if os.path.exists(f):
            try:
                module_results[mod] = pd.read_csv(f)
                logger.info(f"  {mod}: {len(module_results[mod])} 行")
            except Exception as e:
                logger.warning(f"  {mod}加载失败: {e}")

    return module_results


def plot_integration_heatmap(gene_df, pathway_df, output_dir):
    if gene_df.empty and pathway_df.empty:
        return
    fig, axes = plt.subplots(1, 2, figsize=(16, 8))

    if not gene_df.empty:
        top = gene_df.head(15).sort_values('IPS')
        axes[0].barh(range(len(top)), top['IPS'], color='#2E86AB', alpha=0.8)
        axes[0].set_yticks(range(len(top)))
        axes[0].set_yticklabels(top['Gene'], fontsize=9)
        axes[0].set_xlabel('Integrated Priority Score (IPS)', fontsize=11)
        axes[0].set_title('Top Cross-Module Genes (V4)', fontsize=12, fontweight='bold')
        for i, (_, row) in enumerate(top.iterrows()):
            axes[0].text(row['IPS'] + 0.01, i, f"M:{row['ModuleCount']}, FDR:{row['BH_FDR']:.2e}",
                        va='center', fontsize=7)

    if not pathway_df.empty:
        top = pathway_df.head(15).sort_values('Stouffer_Z', ascending=False)
        axes[1].barh(range(len(top)), top['Stouffer_Z'], color='#E8575A', alpha=0.8)
        axes[1].set_yticks(range(len(top)))
        axes[1].set_yticklabels([p[:25] for p in top['Pathway']], fontsize=9)
        axes[1].set_xlabel("Stouffer's Z-score", fontsize=11)
        axes[1].set_title('Top Cross-Module Pathways', fontsize=12, fontweight='bold')

    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, f"cross_module_integration_v4.{FIG_FORMAT}"),
                dpi=FIG_DPI, bbox_inches='tight')
    plt.close()


def main():
    setup_logging()
    logger.info("=" * 60)
    logger.info("模块8 (V4): 跨模块证据融合")
    logger.info("修复: P6 Stouffer方向逻辑 + 动态IPS分母 + BH-FDR")
    logger.info("=" * 60)

    logger.info("[1/4] 加载各模块结果...")
    module_results = _load_module_results()
    logger.info(f"  加载了 {len(module_results)} 个模块")

    logger.info("[2/4] V4基因级整合 (动态IPS + BH-FDR + 方向一致性)...")
    integrator = CrossModuleIntegratorV4()
    cupro_genes = ["FDX1","LIAS","LIPT1","DLAT","PDHA1","PDHB","MTF1","GLS",
                   "CDKN2A","SLC31A1","ATP7A","ATP7B","DLD","DBT","DLST","PDHA2","GCSH"]
    gene_df = integrator.integrate(module_results, cupro_genes)
    if not gene_df.empty:
        gene_df.to_csv(os.path.join(STAGE_DIR, "integrated_gene_ranking_v4.csv"), index=False)
        logger.info(f"  Top 10 基因:")
        for _, row in gene_df.head(10).iterrows():
            logger.info(f"    R{row['IPS_Rank']:2d}: {row['Gene']:10s} IPS={row['IPS']:.3f} Z={row['Stouffer_Z']:.2f} FDR={row['BH_FDR']:.2e} DirCons={row['Direction_Consistency']:.2f}")

    logger.info("[3/4] 通路级整合...")
    evidence = load_module_evidence()
    pathway_df = integrate_pathway_evidence(evidence)
    if not pathway_df.empty:
        pathway_df.to_csv(os.path.join(STAGE_DIR, "integrated_pathway_ranking_v4.csv"), index=False)
        logger.info(f"  Top 10 通路:")
        for _, row in pathway_df.head(10).iterrows():
            logger.info(f"    R{row['Rank']:2d}: {row['Pathway'][:30]:30s} Z={row['Stouffer_Z']:.2f} p={row['Combined_P']:.2e}")

    logger.info("[4/4] 可视化...")
    plot_integration_heatmap(gene_df, pathway_df, STAGE_DIR)

    logger.info(f"\n{'='*60}\n模块8 (V4) 完成!\n{'='*60}")
    return 0


if __name__ == "__main__":
    sys.exit(main())