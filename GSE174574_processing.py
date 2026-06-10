#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
GSE174574 单细胞数据处理流程
================================================================================

数据集: GSE174574 (缺血性中风小胶质细胞scRNA-seq)
平台: GPL21103 (Illumina NovaSeq 6000)

分析流程:
  Step 1: 数据加载与质控
  Step 2: 检测是否包含spliced/unspliced counts
  Step 3a: 若包含splicing信息 → scVelo RNA velocity分析
  Step 3b: 若无splicing信息 → 选项A/B处理
  Step 4: 细胞聚类与差异表达
  Step 5: 铜死亡基因集分析

铜死亡基因集:
  - 核心基因: FDX1/LIAS/LIPT1/DLAT/DLD/PDHA1/PDHB/MTF1/GLS/CDKN2A
  - 扩展基因: SIRT7/ATP7B/CTR1/COX17/ATOX1/CCS

输入文件:
  - GPL21103_family.soft
  - GPL21103_family.xml
  - GSE174574_series_matrix.txt
  - GSE174574_RAW.tar

输出:
  - Figure 1: 速度场/伪时间轨迹可视化
  - Velocity matrix (若有splicing信息)
  - 聚类结果
  - 差异表达谱
================================================================================
"""

import os
import sys
import time
import warnings
from pathlib import Path
from datetime import datetime
from collections import Counter

warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd

try:
    import scanpy as sc
    import anndata
    SCANPY_AVAILABLE = True
except ImportError as e:
    print(f"⚠️ scanpy未安装: {e}")
    SCANPY_AVAILABLE = False

try:
    import matplotlib.pyplot as plt
    import seaborn as sns
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False

# ============================================================================
# 配置参数
# ============================================================================

DATA_DIR = Path(r'D:\反向网络药理学\L1 数据集\RNA-seq')
PLATFORM_SOFT_FILE = DATA_DIR / 'GPL21103_family.soft' / 'GPL21103_family.soft'
PLATFORM_XML_FILE = DATA_DIR / 'GPL21103_family.xml' / 'GPL21103_family.xml'
SERIES_MATRIX_FILE = DATA_DIR / 'GSE174574_series_matrix.txt' / 'GSE174574_series_matrix.txt'
RAW_TAR_FILE = DATA_DIR / 'GSE174574_RAW.tar'

OUTPUT_DIR = Path(r'C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\GSE174574_analysis')
FIGURE_DIR = OUTPUT_DIR / 'figures'
RESULT_DIR = OUTPUT_DIR / 'results'

SEED = 42
np.random.seed(SEED)

# 铜死亡基因集 (人类标准命名, 后续可映射到物种特异性名称)
CUPROPTOSIS_CORE_GENES = [
    'FDX1', 'LIAS', 'LIPT1', 'DLAT', 'DLD',
    'PDHA1', 'PDHB', 'MTF1', 'GLS', 'CDKN2A'
]

CUPROPTOSIS_EXTENDED_GENES = [
    'SIRT7', 'ATP7B', 'CTR1', 'COX17', 'ATOX1', 'CCS'
]

CUPROPTOSIS_ALL_GENES = CUPROPTOSIS_CORE_GENES + CUPROPTOSIS_EXTENDED_GENES

# 小胶质细胞标记基因 (小鼠/大鼠命名, 用于细胞注释)
MICROGLIA_MARKERS = {
    'Homeostatic': ['P2ry12', 'Tmem119', 'Cx3cr1', 'Hexb', 'Sall1'],
    'DAM': ['Trem2', 'Apoe', 'Cst7', 'Lpl', 'Cd9', 'Spp1'],
    'M1': ['Il1b', 'Tnf', 'Il6', 'Cxcl10', 'Nos2'],
    'M2': ['Arg1', 'Mrc1', 'Ym1', 'Fizz1', 'Tgfb1']
}

MIN_CELLS_PER_TYPE = 50
MIN_DIFF_CUPROPTOSIS_GENES = 5
MAX_PCA_BATCH_EFFECT = 0.20


# ============================================================================
# 日志类
# ============================================================================

class Logger:
    def __init__(self, log_file=None):
        self.log_file = log_file
        self.start_time = time.time()
        if log_file:
            log_file.parent.mkdir(parents=True, exist_ok=True)
            with open(log_file, 'w', encoding='utf-8') as f:
                f.write(f"# GSE174574 Analysis Log\n")
                f.write(f"# {datetime.now()}\n{'='*80}\n\n")

    def info(self, msg):
        ts = datetime.now().strftime("%H:%M:%S")
        elapsed = time.time() - self.start_time
        line = f"[{ts}] {elapsed:>8.1f}s [INFO] {msg}"
        print(line)
        if self.log_file:
            with open(self.log_file, 'a', encoding='utf-8') as f:
                f.write(line + '\n')

    def success(self, msg):
        ts = datetime.now().strftime("%H:%M:%S")
        elapsed = time.time() - self.start_time
        line = f"[{ts}] {elapsed:>8.1f}s [✅] {msg}"
        print(line)
        if self.log_file:
            with open(self.log_file, 'a', encoding='utf-8') as f:
                f.write(line + '\n')

    def warn(self, msg):
        ts = datetime.now().strftime("%H:%M:%S")
        elapsed = time.time() - self.start_time
        line = f"[{ts}] {elapsed:>8.1f}s [⚠️] {msg}"
        print(line)
        if self.log_file:
            with open(self.log_file, 'a', encoding='utf-8') as f:
                f.write(line + '\n')

    def error(self, msg):
        ts = datetime.now().strftime("%H:%M:%S")
        elapsed = time.time() - self.start_time
        line = f"[{ts}] {elapsed:>8.1f}s [❌] {msg}"
        print(line)
        if self.log_file:
            with open(self.log_file, 'a', encoding='utf-8') as f:
                f.write(line + '\n')

    def section(self, title):
        line = f"\n{'─'*80}\n{title}\n{'─'*80}"
        print(line)
        if self.log_file:
            with open(self.log_file, 'a', encoding='utf-8') as f:
                f.write(line + '\n')


# ============================================================================
# 数据加载与验证
# ============================================================================

def load_platform_info(soft_file: Path, xml_file: Path) -> dict:
    """加载平台注释信息,检测是否支持spliced/unspliced counting"""
    platform_info = {
        'has_splicing_info': False,
        'platform_type': 'unknown',
        'library_strategy': 'unknown'
    }

    if xml_file.exists():
        try:
            import xml.etree.ElementTree as ET
            tree = ET.parse(xml_file)
            root = tree.getroot()
            for elem in root.iter():
                if 'Strategy' in elem.tag or 'strategy' in elem.get('name', '').lower():
                    platform_info['library_strategy'] = elem.text or elem.get('value', 'unknown')
            platform_info['platform_type'] = 'xml_loaded'
        except Exception as e:
            logger.warn(f"XML解析失败: {e}")

    if soft_file.exists():
        try:
            with open(soft_file, 'r', encoding='utf-8') as f:
                content = f.read()
                splicing_keywords = [
                    'spliced', 'unspliced', 'velocity', 'kallisto|bustools',
                    'velocyto', '10x', 'single cell', 'scRNA-seq'
                ]
                for kw in splicing_keywords:
                    if kw.lower() in content.lower():
                        platform_info['has_splicing_info'] = True
                        platform_info['platform_type'] = f'scRNA-seq ({kw})'
                        break
        except Exception as e:
            logger.warn(f"SOFT解析失败: {e}")

    return platform_info


def detect_splicing_in_raw(tar_file: Path) -> bool:
    """检测RAW文件中是否包含spliced/unspliced矩阵"""
    has_splicing = False

    if tar_file.exists():
        try:
            import tarfile
            with tarfile.open(tar_file, 'r:*') as tar:
                members = tar.getnames()
                splicing_indicators = [
                    'spliced', 'unspliced', 'ambiguous', 'layers',
                    'velocity', 'loom'
                ]
                for member in members:
                    member_lower = member.lower()
                    for indicator in splicing_indicators:
                        if indicator in member_lower:
                            has_splicing = True
                            logger.info(f"  检测到splicing指标: {indicator} in {member}")
                            break
        except Exception as e:
            logger.warn(f"TAR文件检测失败: {e}")
    else:
        logger.warn(f"RAW文件不存在: {tar_file}")

    return has_splicing


def load_series_matrix(series_file: Path) -> anndata.AnnData:
    """加载GEO series matrix或从10x格式加载"""
    logger.info(f"加载数据: {series_file.name}")

    if not series_file.exists():
        raise FileNotFoundError(f"文件不存在: {series_file}")

    try:
        import codecs
        with codecs.open(series_file, 'r', encoding='utf-8', errors='replace') as f:
            content = f.read()

        logger.info(f"  文件大小: {len(content)} 字符")

        has_data_rows = '!series_matrix_table_begin' in content and \
                        '!series_matrix_table_end' in content

        table_start = content.find('!series_matrix_table_begin')
        table_end = content.find('!series_matrix_table_end')

        if table_start > 0 and table_end > table_start:
            table_content = content[table_start:table_end].strip()
            data_lines = [line for line in table_content.split('\n')
                         if line.strip() and not line.startswith('!')]

            has_actual_data = len(data_lines) > 2
        else:
            has_actual_data = False

        if not has_actual_data:
            logger.warn("Series Matrix无表达数据，尝试从RAW.tar加载10x格式数据")
            return load_10x_from_raw(RAW_TAR_FILE)

        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False, encoding='utf-8') as tmp:
            tmp.write(content)
            tmp_path = tmp.name

        try:
            adata = sc.read_text(tmp_path)
        except Exception as e1:
            logger.warn(f"read_text失败: {e1}，尝试read...")
            try:
                adata = sc.read(tmp_path)
            except Exception as e2:
                logger.error(f"read也失败: {e2}")
                raise
        finally:
            os.unlink(tmp_path)

        logger.success(f"  加载完成: {adata.n_obs} cells × {adata.n_vars} genes")
        return adata
    except Exception as e:
        logger.error(f"Series Matrix加载失败: {e}")
        raise


def load_10x_from_raw(tar_file: Path) -> anndata.AnnData:
    """从10x格式的RAW.tar加载数据（扁平结构）"""
    logger.info(f"从10x RAW.tar加载: {tar_file.name}")

    if not tar_file.exists():
        raise FileNotFoundError(f"RAW.tar不存在: {tar_file}")

    import tarfile
    import gzip
    import shutil
    from scipy import io
    from scipy.sparse import csr_matrix

    extract_dir = OUTPUT_DIR / 'raw_extracted'
    extract_dir.mkdir(parents=True, exist_ok=True)

    logger.info("解压RAW.tar...")
    with tarfile.open(tar_file, 'r:*') as tar:
        tar.extractall(extract_dir)
        members = tar.getnames()
        logger.info(f"  解压文件数: {len(members)}")

    sample_prefixes = set()
    for m in members:
        if '_barcodes.tsv.gz' in m:
            prefix = m.replace('_barcodes.tsv.gz', '')
            sample_prefixes.add(prefix)

    sample_prefixes = sorted(sample_prefixes)
    logger.info(f"  检测到样本: {sample_prefixes}")

    all_adatas = []

    for prefix in sample_prefixes:
        sample_name = prefix.split('_')[-1]
        logger.info(f"\n处理样本: {sample_name}")

        barcodes_file = extract_dir / f'{prefix}_barcodes.tsv.gz'
        genes_file = extract_dir / f'{prefix}_genes.tsv.gz'
        matrix_file = extract_dir / f'{prefix}_matrix.mtx.gz'

        if not (barcodes_file.exists() and genes_file.exists() and matrix_file.exists()):
            logger.warn(f"  跳过 {sample_name}: 缺少必要文件")
            continue

        try:
            logger.info("  读取barcodes...")
            with gzip.open(barcodes_file, 'rt') as f:
                barcodes = [line.strip() for line in f]

            logger.info("  读取genes...")
            with gzip.open(genes_file, 'rt') as f:
                genes = [line.strip().split('\t')[1] if '\t' in line else line.strip().split('\t')[0] 
                        for line in f]

            logger.info("  读取matrix...")
            with gzip.open(matrix_file, 'rt') as f:
                lines = f.readlines()
                header_lines = 0
                for i, line in enumerate(lines):
                    if not line.startswith('%') and not line.startswith('%%'):
                        header_lines = i
                        break
                
                n_genes, n_cells, n_entries = map(int, lines[header_lines].strip().split())
                logger.info(f"    矩阵维度: {n_genes} genes × {n_cells} cells, {n_entries} entries")

            row_indices = []
            col_indices = []
            data_values = []

            for line in lines[header_lines+1:]:
                parts = line.strip().split()
                if len(parts) == 3:
                    row_indices.append(int(parts[0]) - 1)
                    col_indices.append(int(parts[1]) - 1)
                    data_values.append(float(parts[2]))

            matrix = csr_matrix(
                (data_values, (row_indices, col_indices)),
                shape=(n_genes, n_cells)
            ).T

            adata_sample = anndata.AnnData(
                X=matrix,
                obs=pd.DataFrame(index=barcodes),
                var=pd.DataFrame(index=genes)
            )
            adata_sample.obs['sample'] = sample_name
            adata_sample.obs['condition'] = 'sham' if 'sham' in sample_name.lower() else 'MCAO'
            adata_sample.obs_names = [f"{sample_name}_{bc}" for bc in adata_sample.obs_names]

            var_names_unique = []
            seen_genes = set()
            for gene in adata_sample.var_names:
                if gene in seen_genes:
                    var_names_unique.append(f"{gene}_dup")
                else:
                    var_names_unique.append(gene)
                    seen_genes.add(gene)
            adata_sample.var_names = var_names_unique
            adata_sample.var_names_make_unique()

            all_adatas.append(adata_sample)
            logger.info(f"  成功加载: {adata_sample.n_obs} cells × {adata_sample.n_vars} genes")
        except Exception as e:
            logger.error(f"  加载失败: {e}")
            import traceback
            logger.error(f"  错误详情: {traceback.format_exc()}")

    if not all_adatas:
        raise ValueError("未能从RAW.tar加载任何样本数据")

    logger.info(f"\n合并所有样本...")
    try:
        adata = anndata.concat(all_adatas, join='outer', label='sample')
    except Exception as concat_error:
        logger.warn(f"concat失败: {concat_error}")
        logger.info("尝试逐个合并...")
        adata = all_adatas[0]
        for i, adata_next in enumerate(all_adatas[1:], 2):
            logger.info(f"  合并样本 {i}/{len(all_adatas)}")
            adata = anndata.concat([adata, adata_next], join='outer', label='sample')

    logger.info(f"合并前: {adata.n_obs} cells × {adata.n_vars} genes")

    adata.var_names_make_unique()
    logger.info(f"基因名称去重后: {adata.n_vars} genes")

    logger.success(f"合并完成: {adata.n_obs} cells × {adata.n_vars} genes")

    return adata


def load_scanpy_if_available():
    """检查并加载必要的分析库"""
    required_packages = {}

    try:
        import scanpy
        required_packages['scanpy'] = True
    except ImportError:
        required_packages['scanpy'] = False

    try:
        import scvelo as scv
        required_packages['scvelo'] = True
    except ImportError:
        required_packages['scvelo'] = False

    try:
        import leidenalg
        required_packages['leiden'] = True
    except ImportError:
        required_packages['leiden'] = False

    return required_packages


# ============================================================================
# 数据预处理
# ============================================================================

def preprocess_adata(adata: anndata.AnnData) -> anndata.AnnData:
    """标准单细胞数据预处理流程"""
    logger.section("🔧 数据预处理")

    logger.info(f"原始维度: {adata.n_obs} cells × {adata.n_vars} genes")

    logger.info("计算质控指标...")
    sc.pp.calculate_qc_metrics(adata, inplace=True)

    adata.var['mt'] = adata.var_names.str.startswith('Mt-') | \
                      adata.var_names.str.startswith('mt-')

    if 'mt' in adata.var.columns:
        adata.obs['pct_counts_mt'] = adata.obs['total_counts']
        mt_mask = adata.var['mt']
        if mt_mask.any():
            adata.obs['pct_counts_mt'] = np.array(adata[:, mt_mask].X.sum(axis=1)).flatten() / \
                                         adata.obs['total_counts'] * 100

    logger.info("应用质控过滤...")
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    if 'total_counts' in adata.obs.columns:
        adata = adata[adata.obs['total_counts'] < np.percentile(adata.obs['total_counts'], 97.5), :]

    if 'pct_counts_mt' in adata.obs.columns:
        adata = adata[adata.obs['pct_counts_mt'] < 20, :]

    logger.info(f"质控后维度: {adata.n_obs} cells × {adata.n_vars} genes")

    logger.info("标准化与对数转换...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    logger.info("识别高变基因...")
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True)

    logger.success(f"预处理完成: {adata.n_obs} cells × {adata.n_vars} HVGs")
    return adata


# ============================================================================
# 细胞聚类与注释
# ============================================================================

def cluster_and_annotate(adata: anndata.AnnData) -> anndata.AnnData:
    """细胞聚类与基于标记基因的注释"""
    logger.section("🔬 细胞聚类与注释")

    logger.info("执行PCA降维...")
    sc.tl.pca(adata, n_comps=30, random_state=SEED, svd_solver='arpack')

    logger.info("构建KNN图...")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=20, random_state=SEED)

    logger.info("Leiden聚类...")
    sc.tl.leiden(adata, resolution=0.5, random_state=SEED)

    logger.info("UMAP可视化...")
    sc.tl.umap(adata, random_state=SEED)

    logger.info("基于标记基因评分进行细胞类型注释...")
    gene_names_upper = {str(g).upper(): str(g) for g in adata.var_names}

    cell_type_scores = {}
    for ct_type, markers in MICROGLIA_MARKERS.items():
        matched_genes = []
        for marker in markers:
            marker_upper = marker.upper()
            if marker_upper in gene_names_upper:
                matched_genes.append(gene_names_upper[marker_upper])
            elif marker in adata.var_names:
                matched_genes.append(marker)

        if len(matched_genes) > 0:
            sc.tl.score_genes(adata, gene_list=matched_genes, score_name=f'score_{ct_type}')
            cell_type_scores[ct_type] = matched_genes
            logger.info(f"  {ct_type}: {len(matched_genes)}/{len(markers)} 标记基因匹配")

    best_scores = adata.obs[[f'score_{ct}' for ct in cell_type_scores.keys()]].values
    best_type_idx = np.argmax(best_scores, axis=1)
    adata.obs['cell_type'] = pd.Categorical([
        list(cell_type_scores.keys())[idx] for idx in best_type_idx
    ])

    type_counts = adata.obs['cell_type'].value_counts()
    logger.info("\n细胞类型分布:")
    for ct, count in type_counts.items():
        pct = count / adata.n_obs * 100
        logger.info(f"  {ct}: {count} ({pct:.1f}%)")

    logger.success(f"聚类与注释完成: {len(type_counts)} 种细胞类型")
    return adata


# ============================================================================
# 差异表达分析
# ============================================================================

def differential_expression_analysis(adata: anndata.AnnData) -> pd.DataFrame:
    """差异表达分析 (Homeostatic vs DAM)"""
    logger.section("📊 差异表达分析")

    if 'Homeostatic' not in adata.obs['cell_type'].values or \
       'DAM' not in adata.obs['cell_type'].values:
        logger.warn("未检测到Homeostatic或DAM细胞类型,跳过差异分析")
        return pd.DataFrame()

    adata_sub = adata[adata.obs['cell_type'].isin(['Homeostatic', 'DAM'])].copy()
    adata_sub.obs['group'] = adata_sub.obs['cell_type']

    sc.tl.rank_genes_groups(
        adata_sub, groupby='group',
        method='wilcoxon', reference='Homeostatic'
    )

    de_results = []
    for gene in adata_sub.var_names:
        try:
            result = {
                'gene': gene,
                'pval_adj': adata_sub.uns['rank_genes_groups']['pvals_adj'][
                    list(adata_sub.uns['rank_genes_groups']['names'].dtype.names).index('DAM')
                ][list(adata_sub.var_names).index(gene)],
                'logfoldchanges': adata_sub.uns['rank_genes_groups']['logfoldchanges'][
                    list(adata_sub.uns['rank_genes_groups']['names'].dtype.names).index('DAM')
                ][list(adata_sub.var_names).index(gene)]
            }
            de_results.append(result)
        except (IndexError, KeyError, ValueError):
            continue

    df_de = pd.DataFrame(de_results)
    if len(df_de) > 0:
        df_de = df_de.sort_values('pval_adj')
        df_de['significant'] = (df_de['pval_adj'] < 0.05) & (abs(df_de['logfoldchanges']) > 0.5)

        logger.info(f"差异基因总数: {len(df_de[df_de['significant']])}")
        logger.success(f"差异表达分析完成")

    return df_de


# ============================================================================
# 铜死亡基因分析
# ============================================================================

def cuproptosis_gene_analysis(adata: anndata.AnnData, de_results: pd.DataFrame) -> dict:
    """铜死亡基因集富集与表达分析"""
    logger.section("🧬 铜死亡基因集分析")

    gene_names_upper = {str(g).upper(): str(g) for g in adata.var_names}

    cuproptosis_in_dataset = []
    cuproptosis_not_in_dataset = []

    for gene in CUPROPTOSIS_ALL_GENES:
        if gene.upper() in gene_names_upper:
            cuproptosis_in_dataset.append(gene_names_upper[gene.upper()])
        elif gene in adata.var_names:
            cuproptosis_in_dataset.append(gene)
        else:
            cuproptosis_not_in_dataset.append(gene)

    logger.info(f"铜死亡基因在数据集中:")
    logger.info(f"  检测到: {len(cuproptosis_in_dataset)}/{len(CUPROPTOSIS_ALL_GENES)}")
    logger.info(f"  未检测到: {cuproptosis_not_in_dataset}")

    de_cuproptosis_intersection = []
    if len(de_results) > 0 and 'gene' in de_results.columns:
        de_sig_genes = set(de_results[de_results['significant']]['gene'].values)
        de_cuproptosis_intersection = [
            g for g in cuproptosis_in_dataset if g in de_sig_genes
        ]

    logger.info(f"差异基因∩铜死亡基因: {len(de_cuproptosis_intersection)}")

    if len(de_cuproptosis_intersection) < MIN_DIFF_CUPROPTOSIS_GENES:
        logger.warn(f"  警告: 交集基因数({len(de_cuproptosis_intersection)}) < 要求({MIN_DIFF_CUPROPTOSIS_GENES})")
    else:
        logger.success(f"  满足要求: {len(de_cuproptosis_intersection)} >= {MIN_DIFF_CUPROPTOSIS_GENES}")

    adata.obs['cuproptosis_score'] = 0
    if len(cuproptosis_in_dataset) > 0:
        sc.tl.score_genes(adata, gene_list=cuproptosis_in_dataset, score_name='cuproptosis_score')
        logger.success("铜死亡评分计算完成")

    return {
        'genes_detected': cuproptosis_in_dataset,
        'genes_not_detected': cuproptosis_not_in_dataset,
        'intersection_with_de': de_cuproptosis_intersection,
        'intersection_count': len(de_cuproptosis_intersection),
        'meets_min_requirement': len(de_cuproptosis_intersection) >= MIN_DIFF_CUPROPTOSIS_GENES
    }


# ============================================================================
# 批量效应评估
# ============================================================================

def batch_effect_assessment(adata: anndata.AnnData) -> dict:
    """PCA批量效应评估"""
    logger.section("📈 批量效应评估")

    if 'X_pca' not in adata.obsm:
        logger.warn("无PCA结果,跳过批量效应评估")
        return {'batch_effect_explained': 0, 'passes_threshold': True}

    pca_variance = np.var(adata.obsm['X_pca'], axis=0)
    total_variance = np.sum(pca_variance)

    if 'batch' in adata.obs.columns or 'sample' in adata.obs.columns:
        batch_col = 'batch' if 'batch' in adata.obs.columns else 'sample'
        from scipy.stats import f_oneway

        batch_effects = []
        for pc in range(min(10, adata.obsm['X_pca'].shape[1])):
            groups = [adata.obsm['X_pca'][adata.obs[batch_col] == b, pc]
                     for b in adata.obs[batch_col].unique()]
            if all(len(g) > 0 for g in groups):
                _, pval = f_oneway(*groups)
                batch_effects.append(pval)

        batch_effect_rate = np.mean([1 - p for p in batch_effects]) if batch_effects else 0
    else:
        logger.info("无批次信息,使用方差解释率估算")
        batch_effect_rate = pca_variance[0] / total_variance if total_variance > 0 else 0

    passes = batch_effect_rate < MAX_PCA_BATCH_EFFECT

    logger.info(f"  批量效应解释率: {batch_effect_rate:.2%}")
    logger.info(f"  阈值: < {MAX_PCA_BATCH_EFFECT:.0%}")
    logger.info(f"  是否通过: {'✅' if passes else '❌'}")

    return {
        'batch_effect_explained': batch_effect_rate,
        'passes_threshold': passes
    }


# ============================================================================
# scVelo RNA Velocity分析 (若有splicing信息)
# ============================================================================

def run_scvelo_analysis(adata: anndata.AnnData, output_dir: Path) -> dict:
    """scVelo RNA velocity分析"""
    logger.section("🌀 scVelo RNA Velocity分析")

    try:
        import scvelo as scv
    except ImportError:
        logger.error("scvelo未安装,无法进行RNA velocity分析")
        return {'success': False, 'reason': 'scvelo_not_installed'}

    logger.info("预处理velocity数据...")
    try:
        scv.pp.filter_and_normalize(adata)
        scv.pp.moments(adata, n_pcs=20)

        logger.info("计算RNA velocity (稳态模型)...")
        scv.tl.velocity(adata, mode='stochastic')
        scv.tl.velocity_graph(adata)

        logger.info("计算velocity embedding...")
        scv.tl.velocity_embedding(adata, basis='umap')

        velocity_matrix_file = output_dir / 'results' / 'velocity_matrix.npy'
        np.save(velocity_matrix_file, adata.obsm['velocity_umap'])
        logger.success(f"Velocity matrix保存: {velocity_matrix_file}")

        velocity_data_file = output_dir / 'results' / 'adata_velocity.h5ad'
        adata.write(velocity_data_file)
        logger.success(f"Velocity AnnData保存: {velocity_data_file}")

        logger.info("生成Figure 1: Velocity Flow Field...")
        fig1_path = output_dir / 'figures' / 'Figure_1_Velocity_Flow_Field.png'
        fig1_path.parent.mkdir(parents=True, exist_ok=True)

        fig, axes = plt.subplots(1, 2, figsize=(18, 7))

        scv.pl.velocity_embedding_stream(
            adata, basis='umap', color='cell_type',
            ax=axes[0], show=False, title='RNA Velocity Flow Field'
        )

        scv.pl.velocity_embedding_grid(
            adata, basis='umap', color='cell_type',
            ax=axes[1], show=False, title='RNA Velocity Grid'
        )

        plt.tight_layout()
        plt.savefig(fig1_path, dpi=300, bbox_inches='tight')
        plt.close()
        logger.success(f"Figure 1保存: {fig1_path}")

        return {
            'success': True,
            'velocity_matrix_file': str(velocity_matrix_file),
            'adata_velocity_file': str(velocity_data_file),
            'figure_1_file': str(fig1_path),
            'note': 'Velocity matrix已保存,可直接传递给L2c模块作为NeuralODE初始速度场'
        }
    except Exception as e:
        logger.error(f"scVelo分析失败: {e}")
        logger.warn("回退到伪时间分析")
        return {'success': False, 'reason': str(e)[:100]}


# ============================================================================
# Monocle3/Slingshot伪时间分析 (仅用于可视化)
# ============================================================================

def run_pseudotime_analysis(adata: anndata.AnnData, output_dir: Path) -> dict:
    """Monocle3/Slingshot伪时间分析 (仅用于可视化Figure 1)"""
    logger.section("⏱️ 伪时间分析 (仅可视化)")

    pseudotime_results = {
        'success': False,
        'note': 'Pseudotime represents state ordering, not actual reperfusion time'
    }

    try:
        try:
            import scvelo as scv
            scv_available = True
        except ImportError:
            scv_available = False

        if scv_available:
            logger.info("使用scVelo进行速度伪时间分析...")
            adata_sub = adata[adata.obs['cell_type'].isin(['Homeostatic', 'DAM'])].copy()

            if len(adata_sub) < 50:
                logger.warn("细胞数不足,使用全部细胞进行伪时间分析")
                adata_sub = adata.copy()

            try:
                scv.tl.velocity_pseudotime(adata_sub)
                adata.obs['pseudotime'] = 0
                if 'velocity_pseudotime' in adata_sub.obs.columns:
                    adata.obs.loc[adata_sub.obs.index, 'pseudotime'] = \
                        adata_sub.obs['velocity_pseudotime'].values
                pseudotime_results['method'] = 'scvelo_velocity_pseudotime'
            except Exception as e:
                logger.warn(f"scVelo伪时间失败: {e}")
                pseudotime_results['method'] = 'umap_trajectory_inference'
                sc.tl.dpt(adata, n_dcs=2)
                if 'dpt_pseudotime' in adata.obs.columns:
                    adata.obs['pseudotime'] = adata.obs['dpt_pseudotime']
        else:
            logger.info("使用Scanpy DPT进行伪时间分析...")
            sc.tl.dpt(adata, n_dcs=2)
            if 'dpt_pseudotime' in adata.obs.columns:
                adata.obs['pseudotime'] = adata.obs['dpt_pseudotime']
                pseudotime_results['method'] = 'scanpy_dpt'
            else:
                logger.warn("DPT伪时间计算失败")
                return pseudotime_results

        pseudotime_results['success'] = True
        pseudotime_results['has_pseudotime'] = 'pseudotime' in adata.obs.columns

        logger.info("生成Figure 1: Pseudotime Trajectory (Homeostatic → DAM)...")
        fig, axes = plt.subplots(1, 3, figsize=(22, 7))

        sc.pl.umap(
            adata, color='cell_type', ax=axes[0], show=False,
            title='Cell Types (Homeostatic vs DAM)',
            palette=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
        )

        if 'pseudotime' in adata.obs.columns:
            sc.pl.umap(
                adata, color='pseudotime', ax=axes[1], show=False,
                title='Pseudotime Trajectory',
                cmap='viridis'
            )

            root_mask = adata.obs['cell_type'] == 'Homeostatic'
            end_mask = adata.obs['cell_type'] == 'DAM'

            axes[2].scatter(
                adata.obsm['X_umap'][~root_mask & ~end_mask, 0],
                adata.obsm['X_umap'][~root_mask & ~end_mask, 1],
                c='lightgray', s=10, alpha=0.4
            )
            axes[2].scatter(
                adata.obsm['X_umap'][root_mask, 0],
                adata.obsm['X_umap'][root_mask, 1],
                c='blue', s=15, alpha=0.8, label='Homeostatic (Root)'
            )
            axes[2].scatter(
                adata.obsm['X_umap'][end_mask, 0],
                adata.obsm['X_umap'][end_mask, 1],
                c='red', s=15, alpha=0.8, label='DAM (Terminal)'
            )

            arrows = np.linspace(0, 1, 20)
            for i in range(len(arrows)-1):
                pts_root = adata.obsm['X_umap'][root_mask]
                pts_dam = adata.obsm['X_umap'][end_mask]
                if len(pts_root) > 0 and len(pts_dam) > 0:
                    start = np.mean(pts_root, axis=0)
                    end = np.mean(pts_dam, axis=0)
                    t = arrows[i]
                    t_next = arrows[i+1]
                    x = start[0] + t * (end[0] - start[0])
                    y = start[1] + t * (end[1] - start[1])
                    dx = (t_next - t) * (end[0] - start[0])
                    dy = (t_next - t) * (end[1] - start[1])
                    axes[2].arrow(x, y, dx*5, dy*5,
                                head_width=0.3, head_length=0.2,
                                fc='green', ec='green', alpha=0.6)

            axes[2].text(0.5, 0.02,
                        'Pseudotime represents state ordering,\nnot actual reperfusion time',
                        transform=axes[2].transAxes,
                        ha='center', va='bottom', fontsize=10,
                        bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.7))

            axes[2].set_title('Trajectory: Homeostatic → DAM')
            axes[2].set_xlabel('UMAP1')
            axes[2].set_ylabel('UMAP2')
            axes[2].legend()
        else:
            axes[1].text(0.5, 0.5, 'Pseudotime\nCalculation Failed',
                        ha='center', va='center', transform=axes[1].transAxes)
            axes[2].text(0.5, 0.5, 'No Trajectory Available',
                        ha='center', va='center', transform=axes[2].transAxes)

        plt.tight_layout()
        fig1_path = output_dir / 'figures' / 'Figure_1_Pseudotime_Trajectory.png'
        fig1_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(fig1_path, dpi=300, bbox_inches='tight')
        plt.close()
        logger.success(f"Figure 1保存: {fig1_path}")

        pseudotime_results['figure_1_file'] = str(fig1_path)

    except Exception as e:
        logger.error(f"伪时间分析失败: {e}")

    return pseudotime_results


# ============================================================================
# 结果汇总与输出
# ============================================================================

def save_analysis_results(adata: anndata.AnnData, output_dir: Path,
                         cuproptosis_results: dict, batch_results: dict,
                         velocity_results: dict = None, pseudotime_results: dict = None):
    """保存所有分析结果"""
    logger.section("💾 保存分析结果")

    result_dir = output_dir / 'results'
    result_dir.mkdir(parents=True, exist_ok=True)

    adata_file = result_dir / 'GSE174574_processed.h5ad'
    adata.write(adata_file)
    logger.success(f"AnnData保存: {adata_file}")

    obs_file = result_dir / 'cell_annotations.csv'
    adata.obs.to_csv(obs_file)
    logger.success(f"细胞注释保存: {obs_file}")

    if cuproptosis_results:
        cupro_file = result_dir / 'cuproptosis_analysis.json'
        import json
        cupro_serializable = {
            'genes_detected': cuproptosis_results.get('genes_detected', []),
            'genes_not_detected': cuproptosis_results.get('genes_not_detected', []),
            'intersection_with_de': cuproptosis_results.get('intersection_with_de', []),
            'intersection_count': int(cuproptosis_results.get('intersection_count', 0)),
            'meets_min_requirement': bool(cuproptosis_results.get('meets_min_requirement', False))
        }
        with open(cupro_file, 'w', encoding='utf-8') as f:
            json.dump(cupro_serializable, f, indent=2, ensure_ascii=False)
        logger.success(f"铜死亡分析结果: {cupro_file}")

    if batch_results:
        batch_file = result_dir / 'batch_effect_assessment.json'
        import json
        batch_serializable = {
            'batch_effect_explained': float(batch_results.get('batch_effect_explained', 0)),
            'passes_threshold': bool(batch_results.get('passes_threshold', True))
        }
        with open(batch_file, 'w', encoding='utf-8') as f:
            json.dump(batch_serializable, f, indent=2, ensure_ascii=False)
        logger.success(f"批量效应评估: {batch_file}")

    if velocity_results and velocity_results.get('success'):
        logger.info(f"Velocity matrix文件: {velocity_results.get('velocity_matrix_file')}")
        logger.info(f"Velocity AnnData: {velocity_results.get('adata_velocity_file')}")
        logger.info("⚠️ Velocity matrix可直接传递给L2c模块作为NeuralODE初始速度场")

    if pseudotime_results:
        logger.info(f"伪时间分析结果: {pseudotime_results.get('note', '')}")
        logger.info("⚠️ 此伪时间数据不进入L2c数学初始化,仅用于可视化")

    logger.success(f"所有结果保存至: {output_dir}")


# ============================================================================
# 主流程
# ============================================================================

def main():
    global logger

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    RESULT_DIR.mkdir(parents=True, exist_ok=True)

    log_file = OUTPUT_DIR / f"GSE174574_Analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    logger = Logger(log_file=log_file)

    logger.section("🚀 GSE174574 单细胞数据处理流程")
    logger.info(f"时间: {datetime.now()}")
    logger.info(f"输出目录: {OUTPUT_DIR}")
    logger.info(f"随机种子: {SEED}")

    if not SCANPY_AVAILABLE:
        logger.error("scanpy未安装,无法继续")
        return False

    packages = load_scanpy_if_available()
    logger.info(f"可用包: {packages}")

    logger.section("📥 Step 1: 加载平台信息与数据检测")
    platform_info = load_platform_info(PLATFORM_SOFT_FILE, PLATFORM_XML_FILE)
    logger.info(f"平台信息: {platform_info}")

    has_splicing_raw = detect_splicing_in_raw(RAW_TAR_FILE)
    logger.info(f"RAW文件splicing检测: {has_splicing_raw}")

    has_splicing = platform_info.get('has_splicing_info', False) or has_splicing_raw
    logger.info(f"最终splicing判断: {'✅ 有splicing信息' if has_splicing else '❌ 无splicing信息'}")

    logger.section("📥 Step 2: 加载数据")
    adata = None

    if SERIES_MATRIX_FILE.exists():
        try:
            adata = load_series_matrix(SERIES_MATRIX_FILE)
            adata = preprocess_adata(adata)
        except Exception as e:
            logger.error(f"Series Matrix加载失败: {e}")
            logger.info("尝试其他加载方式...")

    if adata is None:
        logger.error("无法加载任何数据文件")
        logger.info("请确保以下文件存在:")
        logger.info(f"  - {SERIES_MATRIX_FILE}")
        return False

    logger.section("📊 Step 3: 细胞聚类与注释")
    adata = cluster_and_annotate(adata)

    logger.section("📈 Step 4: 差异表达分析")
    de_results = differential_expression_analysis(adata)

    if len(de_results) > 0:
        de_file = RESULT_DIR / 'differential_expression.csv'
        de_results.to_csv(de_file, index=False)
        logger.success(f"差异表达结果: {de_file}")

    logger.section("🧬 Step 5: 铜死亡基因集分析")
    cuproptosis_results = cuproptosis_gene_analysis(adata, de_results)

    logger.section("📈 Step 6: 批量效应评估")
    batch_results = batch_effect_assessment(adata)

    logger.section("🔀 Step 7: 轨迹/速度分析")
    velocity_results = None
    pseudotime_results = None

    if has_splicing:
        logger.info("检测到splicing信息 (Raw FASTQ)")
        logger.info("需要先运行velocyto生成spliced/unspliced counts")
        logger.info("\n📋 velocyto运行步骤:")
        logger.info("  1. 从SRA下载FASTQ文件:")
        logger.info("     prefetch SRRXXXXXXX")
        logger.info("     fasterq-dump SRRXXXXXXX")
        logger.info("  2. 运行velocyto:")
        logger.info("     velocyto run10x -m mm10_rmsk.gtf fastq_folder/ mm10_annotation.gtf")
        logger.info("  3. 生成.loom文件后，重新运行此脚本")
        logger.info("\n⚠️ 当前跳过scVelo，回退到伪时间分析")
        
        pseudotime_results = run_pseudotime_analysis(adata, OUTPUT_DIR)
        if pseudotime_results.get('success'):
            logger.success("✅ 流程 2 (选项B): 伪时间分析完成 (仅用于可视化)")
            logger.info("   - Figure 1: Pseudotime trajectory已生成")
            logger.info("   - 注释: Pseudotime represents state ordering, not actual reperfusion time")
            logger.info("   - 此伪时间数据不进入L2c数学初始化")
    else:
        logger.info("无splicing信息,执行伪时间分析 (选项B)")

        pseudotime_results = run_pseudotime_analysis(adata, OUTPUT_DIR)

        if pseudotime_results.get('success'):
            logger.success("✅ 流程 2 (选项B): 伪时间分析完成 (仅用于可视化)")
            logger.info("   - Figure 1: Pseudotime trajectory已生成")
            logger.info("   - 注释: Pseudotime represents state ordering, not actual reperfusion time")
            logger.info("   - 此伪时间数据不进入L2c数学初始化")
        else:
            logger.warn("伪时间分析失败,采用选项A: 无轨迹分析")
            logger.info("   - L1仅输出: 聚类结果 + 差异表达谱")
            logger.info("   - L2c模块将完全依赖GSE23160 Bulk time-series数据进行时间校准")

    logger.section("💾 Step 8: 保存所有结果")
    save_analysis_results(
        adata, OUTPUT_DIR,
        cuproptosis_results, batch_results,
        velocity_results, pseudotime_results
    )

    logger.section("🎉 全部分析完成!")
    logger.info(f"总耗时: {(time.time() - logger.start_time) / 60:.1f} 分钟")
    logger.info(f"结果目录: {OUTPUT_DIR}")
    logger.info(f"图表目录: {FIGURE_DIR}")

    return True


if __name__ == '__main__':
    success = main()
    sys.exit(0 if success else 1)
