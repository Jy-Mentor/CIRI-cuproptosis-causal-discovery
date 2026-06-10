# -*- coding: utf-8 -*-
"""
GEO数据统一处理模块
====================
提供GEO数据集解析、探针注释映射、差异表达分析等统一接口

统一封装了 deg_temporal_integration.py 和 gat_network_pharmacology_pipeline.py
中重复的GEO数据解析和DEG分析逻辑，消除DRY违规。

功能:
  - GEO Series Matrix文件解析
  - GPL1355平台探针注释解析
  - 探针→基因折叠
  - 芯片数据差异分析 (Welch t-test + BH校正)
  - RNA-seq数据差异分析 (PyDESeq2)

版本: v1.0 | 日期: 2026-05-28
"""

import os
import gzip
import io
import logging
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

logger = logging.getLogger(__name__)


class GEODataProcessor:
    """GEO数据统一处理器"""

    # ============================================================
    # 文件查找
    # ============================================================

    @staticmethod
    def find_file(dir_path, patterns):
        for f in os.listdir(dir_path):
            for pat in patterns:
                if pat.lower() in f.lower():
                    return os.path.join(dir_path, f)
        return None

    # ============================================================
    # GEO Series Matrix 解析
    # ============================================================

    @staticmethod
    def parse_series_matrix(filepath, return_meta=False):
        open_func = gzip.open if filepath.endswith('.gz') else open
        with open_func(filepath, 'rt', encoding='latin-1') as f:
            content = f.read()

        meta = {}
        if return_meta:
            for line in content.splitlines():
                if line.startswith('!Sample_title'):
                    parts = line.split('\t')
                    meta['sample_titles'] = [p.strip('"') for p in parts[1:]]
                elif line.startswith('!Sample_geo_accession'):
                    parts = line.split('\t')
                    meta['sample_geo'] = [p.strip('"') for p in parts[1:]]

        data_start = content.find('!series_matrix_table_begin')
        data_end = content.find('!series_matrix_table_end')
        if data_start == -1 or data_end == -1:
            raise ValueError(f'无法在 {filepath} 中找到 series_matrix_table 标记')

        table_text = content[data_start:data_end]
        table_text = table_text.replace('!series_matrix_table_begin', '').strip()

        df = pd.read_csv(io.StringIO(table_text), sep='\t', quoting=1,
                         dtype=str, low_memory=False)
        if 'ID_REF' in df.columns:
            df = df.set_index('ID_REF')
        else:
            df = df.set_index(df.columns[0])

        df = df.apply(pd.to_numeric, errors='coerce')

        if return_meta:
            return df, meta
        return df

    # ============================================================
    # GPL1355 平台注释解析
    # ============================================================

    @staticmethod
    def parse_gpl1355_annotation(filepath):
        mapping = {}
        with open(filepath, 'r', encoding='latin-1') as f:
            for line in f:
                if line.startswith('#') or line.strip() == '':
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 11:
                    continue
                probe_id = parts[0].strip()
                gene_symbol = parts[10].strip()
                if gene_symbol and gene_symbol != '---':
                    mapping[probe_id] = gene_symbol.split('///')[0].strip()
        return mapping

    # ============================================================
    # 探针→基因折叠
    # ============================================================

    @staticmethod
    def collapse_probes_to_genes(expr_df, probe_to_gene):
        mapped_probes = set(probe_to_gene.keys()) & set(expr_df.index)
        expr_mapped = expr_df.loc[list(mapped_probes)]
        probe_to_gene_sub = {p: probe_to_gene[p] for p in mapped_probes}

        gene_rows = []
        for gene in set(probe_to_gene_sub.values()):
            probes = [p for p in mapped_probes if probe_to_gene_sub[p] == gene]
            if len(probes) == 1:
                gene_rows.append((gene, expr_mapped.loc[probes[0]]))
            else:
                sub = expr_mapped.loc[probes]
                mean_expr = sub.mean(axis=1)
                best_probe = mean_expr.idxmax()
                gene_rows.append((gene, expr_mapped.loc[best_probe]))

        result = pd.DataFrame([row[1] for row in gene_rows],
                              index=[r[0] for r in gene_rows])
        return result

    # ============================================================
    # 芯片数据差异分析 (Welch t-test + BH校正)
    # ============================================================

    @staticmethod
    def deg_microarray_t_test(expr_df, case_samples, control_samples):
        results = []
        case = expr_df[case_samples].values
        control = expr_df[control_samples].values

        for i, gene in enumerate(expr_df.index):
            c = control[i, :].astype(float)
            t = case[i, :].astype(float)
            if np.all(np.isnan(c)) or np.all(np.isnan(t)):
                continue
            c = c[~np.isnan(c)]
            t = t[~np.isnan(t)]
            if len(c) < 2 or len(t) < 2:
                continue
            log2fc = np.mean(t) - np.mean(c)
            stat, pval = stats.ttest_ind(t, c, equal_var=False)
            results.append({
                'gene_symbol': gene,
                'log2FoldChange': log2fc,
                'stat': stat,
                'pvalue': pval
            })

        res_df = pd.DataFrame(results)
        if res_df.empty:
            return res_df

        reject, padj, _, _ = multipletests(res_df['pvalue'], method='fdr_bh')
        res_df['padj'] = padj
        res_df = res_df.sort_values('pvalue')
        return res_df

    # ============================================================
    # RNA-seq数据差异分析 (PyDESeq2)
    # ============================================================

    @staticmethod
    def deg_rnaseq_pydeseq2(counts_df, metadata, case_label, control_label):
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats

        samples = metadata.index.tolist()
        common_genes = counts_df.index.tolist()

        counts_sub = counts_df.loc[common_genes, samples].astype(int)
        counts_sub = counts_sub[~(counts_sub == 0).all(axis=1)]

        dds = DeseqDataSet(
            counts=counts_sub.T,
            metadata=metadata,
            design='~condition',
        )
        dds.deseq2()

        stat_res = DeseqStats(dds, contrast=['condition', case_label, control_label])
        stat_res.summary()

        result = stat_res.results_df.copy()
        result = result.reset_index()
        first_col = result.columns[0]
        if first_col != 'gene_symbol':
            result = result.rename(columns={first_col: 'gene_symbol'})
        result = result[['gene_symbol', 'log2FoldChange', 'pvalue', 'padj']].dropna()
        return result