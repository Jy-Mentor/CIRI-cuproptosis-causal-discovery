# ============================================================
# 阶段2: GSE210986 单细胞RNA-seq真实分析 (Scanpy)
# 数据: 小鼠MCAO模型, Sham vs MCAO
# 来源: GSE210986 (GEO)
# 输出: 聚类结果、细胞类型注释、铜死亡评分、差异表达
# ============================================================

import os
import sys
import logging
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

# 路径配置 - 使用统一导入接口
from pipeline_scripts import (
    RESULTS_DIR, DATA_DIR,
    CUPROPTOSIS_GENES, CUPROPTOSIS_RELATED,
    SC_QC_MIN_GENES, SC_QC_MAX_GENES, SC_QC_MIN_COUNTS, SC_QC_MAX_COUNTS,
    SC_QC_MAX_MITO, SC_N_HVGS, SC_N_PCS, SC_KNN_NEIGHBORS, SC_RESOLUTION,
)

logging.basicConfig(level=logging.INFO, format='%(asctime)s | %(levelname)-8s | %(message)s')
logger = logging.getLogger(__name__)

sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=150, facecolor='white')

STAGE_DIR = os.path.join(RESULTS_DIR, "stage2_single_cell")
os.makedirs(STAGE_DIR, exist_ok=True)

SC_DATA_DIR = os.path.join(DATA_DIR, "GSE210986")

# ---- 小鼠脑细胞marker基因 (文献来源) ----
BRAIN_MARKERS = {
    'Neuron': ['Rbfox3', 'Tubb3', 'Map2', 'Syp', 'Snap25', 'Syt1'],
    'Microglia': ['Cx3cr1', 'Tmem119', 'P2ry12', 'Csf1r', 'Aif1', 'Itgam'],
    'Astrocyte': ['Gfap', 'Aqp4', 'Slc1a3', 'Aldh1l1', 'S100b'],
    'Oligodendrocyte': ['Mbp', 'Mog', 'Plp1', 'Cnp', 'Olig1', 'Olig2'],
    'OPC': ['Pdgfra', 'Cspg4', 'Sox10', 'Nkx2-2'],
    'Endothelial': ['Cldn5', 'Pecam1', 'Flt1', 'Tek', 'Cdh5'],
    'Pericyte': ['Pdgfrb', 'Rgs5', 'Cspg4', 'Anpep'],
    'Ependymal': ['Foxj1', 'Tmem212', 'Cfap54'],
}

# ---- 铜死亡基因 (小鼠同源) ----
CUPROPTOSIS_MOUSE = [
    "Fdx1", "Lias", "Lipt1", "Dlat", "Pdha1", "Pdhb",
    "Mtf1", "Gls", "Cdkn2a", "Slc31a1", "Atp7a", "Atp7b",
    "Dld", "Dbt", "Dlst", "Pdha2", "Gcsh"
]

CUPROPTOSIS_RELATED_MOUSE = [
    "Sod1", "Sod2", "Cat", "Gpx1", "Gsr", "Hmox1", "Nfe2l2", "Keap1",
    "Mt1", "Mt2", "Cox17", "Ccs", "Atox1", "Cox11", "Sco1", "Sco2",
    "Cp", "Steap1", "Steap2", "Steap3", "Steap4",
    "Commd1", "Xiap", "Ccdc22"
]


def load_count_matrix(tsv_path):
    """加载gzip压缩的TSV计数矩阵"""
    logger.info(f"  加载: {os.path.basename(tsv_path)}")
    df = pd.read_csv(tsv_path, sep='\t', compression='gzip', index_col=0)
    logger.info(f"    基因数: {df.shape[0]}, 细胞数: {df.shape[1]}")
    return df


def create_anndata(count_df, condition):
    """从计数矩阵创建AnnData对象"""
    adata = sc.AnnData(
        X=sparse.csr_matrix(count_df.values.T),
        obs=pd.DataFrame(index=count_df.columns),
        var=pd.DataFrame(index=count_df.index)
    )
    adata.obs['condition'] = condition
    adata.obs_names = [f"{condition}_{cn}" for cn in count_df.columns]
    adata.var_names_make_unique()
    return adata


def annotate_cell_types(adata):
    """基于marker基因注释细胞类型"""
    logger.info("  细胞类型注释 (基于marker基因)...")
    
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon', n_genes=50)
    
    cell_type_map = {}
    for cluster in adata.obs['leiden'].cat.categories:
        cluster_genes = set(adata.uns['rank_genes_groups']['names'][cluster][:50])
        
        best_type = 'Unknown'
        best_score = 0
        for ct, markers in BRAIN_MARKERS.items():
            markers_in_cluster = [m for m in markers if m in cluster_genes]
            score = len(markers_in_cluster) / len(markers)
            if score > best_score:
                best_score = score
                best_type = ct
        
        cell_type_map[cluster] = best_type if best_score > 0.1 else 'Unknown'
        logger.info(f"    簇 {cluster}: {best_type} (score={best_score:.2f})")
    
    adata.obs['cell_type'] = adata.obs['leiden'].map(cell_type_map)
    return cell_type_map


def compute_cuproptosis_score(adata):
    """计算每个细胞的铜死亡评分"""
    logger.info("  计算铜死亡评分...")
    
    cupro_genes_present = [g for g in CUPROPTOSIS_MOUSE if g in adata.var_names]
    cupro_related_present = [g for g in CUPROPTOSIS_RELATED_MOUSE if g in adata.var_names]
    
    logger.info(f"    核心铜死亡基因: {len(cupro_genes_present)}/{len(CUPROPTOSIS_MOUSE)}")
    logger.info(f"    铜死亡相关基因: {len(cupro_related_present)}/{len(CUPROPTOSIS_RELATED_MOUSE)}")
    
    if len(cupro_genes_present) > 0:
        sc.tl.score_genes(adata, gene_list=cupro_genes_present, 
                          score_name='cuproptosis_score', ctrl_size=min(50, len(cupro_genes_present)))
    
    if len(cupro_related_present) > 0:
        sc.tl.score_genes(adata, gene_list=cupro_related_present,
                          score_name='cuproptosis_related_score', ctrl_size=min(50, len(cupro_related_present)))
    
    return cupro_genes_present, cupro_related_present


def analyze_cuproptosis_by_celltype(adata):
    """按细胞类型分析铜死亡评分差异"""
    logger.info("  按细胞类型分析铜死亡评分...")
    
    results = []
    for ct in adata.obs['cell_type'].unique():
        ct_mask = adata.obs['cell_type'] == ct
        sham_mask = adata.obs['condition'] == 'Sham'
        mcao_mask = adata.obs['condition'] == 'MCAO'
        
        sham_score = adata[ct_mask & sham_mask].obs.get('cuproptosis_score', pd.Series(dtype=float))
        mcao_score = adata[ct_mask & mcao_mask].obs.get('cuproptosis_score', pd.Series(dtype=float))
        
        results.append({
            'cell_type': ct,
            'n_sham': len(sham_score),
            'n_mcao': len(mcao_score),
            'sham_mean': sham_score.mean() if len(sham_score) > 0 else np.nan,
            'mcao_mean': mcao_score.mean() if len(mcao_score) > 0 else np.nan,
            'delta': (mcao_score.mean() - sham_score.mean()) if (len(sham_score) > 0 and len(mcao_score) > 0) else np.nan,
        })
    
    df = pd.DataFrame(results).sort_values('delta', ascending=False)
    for _, row in df.iterrows():
        logger.info(f"    {row['cell_type']:20s}: Sham={row['sham_mean']:.3f}, "
                   f"MCAO={row['mcao_mean']:.3f}, Δ={row['delta']:+.3f}")
    
    return df


def find_effect_cell_types(cupro_df, delta_threshold=0.02, top_n=5):
    """识别铜死亡评分显著变化的效应细胞类型（排除Unknown，取Top-N）"""
    known_df = cupro_df[cupro_df['cell_type'] != 'Unknown'].copy()
    known_df = known_df[
        (known_df['n_sham'] >= 10) & 
        (known_df['n_mcao'] >= 10)
    ]
    known_df['abs_delta'] = known_df['delta'].abs()
    known_df = known_df.sort_values('abs_delta', ascending=False)
    
    effect_types = known_df.head(top_n)['cell_type'].tolist()
    
    logger.info(f"  效应细胞类型 (Top-{top_n} |Δ|, 排除Unknown):")
    for _, row in known_df.head(top_n).iterrows():
        logger.info(f"    {row['cell_type']:20s}: Δ={row['delta']:+.4f}")
    
    return effect_types


def extract_marker_genes(adata, effect_types):
    """提取效应细胞类型的marker基因"""
    logger.info("  提取效应细胞类型marker基因...")
    
    all_markers = []
    for ct in effect_types:
        ct_cells = adata[adata.obs['cell_type'] == ct]
        if ct_cells.n_obs < 10:
            continue
        
        sc.tl.rank_genes_groups(ct_cells, 'condition', method='wilcoxon',
                                groups=['MCAO'], reference='Sham', n_genes=100)
        
        deg_df = sc.get.rank_genes_groups_df(ct_cells, group='MCAO')
        deg_df = deg_df[deg_df['pvals_adj'] < 0.05]
        deg_df['cell_type'] = ct
        
        top_genes = deg_df.head(30)['names'].tolist()
        all_markers.extend(top_genes)
        logger.info(f"    {ct}: {len(deg_df)} DEGs, top: {top_genes[:5]}")
    
    unique_markers = list(set(all_markers))
    logger.info(f"  效应细胞marker基因总数: {len(unique_markers)}")
    return unique_markers


def analyze_ligand_receptor(adata, effect_types):
    """分析效应细胞类型的配体-受体对（基于原始计数）"""
    logger.info("  分析配体-受体互作...")
    
    lr_pairs = [
        ('Cx3cl1', 'Cx3cr1'), ('Ccl2', 'Ccr2'), ('Ccl5', 'Ccr5'),
        ('Il1b', 'Il1r1'), ('Tnf', 'Tnfrsf1a'), ('Il6', 'Il6st'),
        ('Tgfb1', 'Tgfbr1'), ('Vegfa', 'Flt1'), ('Bdnf', 'Ntrk2'),
        ('Sema3a', 'Nrp1'), ('Ephb2', 'Ephb1'), ('Wnt3a', 'Fzd1'),
        ('Shh', 'Ptch1'), ('Notch1', 'Dll1'), ('Csf1', 'Csf1r'),
        ('Apoe', 'Lrp1'), ('Spp1', 'Cd44'), ('C1qa', 'Cd93'),
    ]
    
    ligands = []
    receptors = []
    for lig, rec in lr_pairs:
        if lig in adata.var_names:
            pct = (adata[:, lig].X.toarray().flatten() > 0).mean()
            if pct > 0.05:
                ligands.append(lig)
        if rec in adata.var_names:
            pct = (adata[:, rec].X.toarray().flatten() > 0).mean()
            if pct > 0.05:
                receptors.append(rec)
    
    logger.info(f"  表达的配体: {len(ligands)}: {ligands}")
    logger.info(f"  表达的受体: {len(receptors)}: {receptors}")
    return ligands, receptors


def main():
    logger.info("=" * 60)
    logger.info("阶段2: GSE210986 单细胞RNA-seq分析 (真实数据)")
    logger.info("=" * 60)
    
    # ---- 1. 加载数据 ----
    logger.info("[1/7] 加载计数矩阵...")
    
    sham_counts = load_count_matrix(
        os.path.join(SC_DATA_DIR, "GSM6443690_sham.counts.tsv.gz"))
    mcao_counts = load_count_matrix(
        os.path.join(SC_DATA_DIR, "GSM6443691_MCAO.counts.tsv.gz"))
    
    # ---- 2. 创建AnnData ----
    logger.info("[2/7] 创建AnnData对象...")
    
    adata_sham = create_anndata(sham_counts, 'Sham')
    adata_mcao = create_anndata(mcao_counts, 'MCAO')
    
    adata = adata_sham.concatenate(adata_mcao, batch_key='batch', 
                                    index_unique=None)
    logger.info(f"  合并后: {adata.n_obs} 细胞 × {adata.n_vars} 基因")
    
    # ---- 3. QC过滤 ----
    logger.info("[3/7] 质量控制...")
    
    adata.var['mt'] = adata.var_names.str.startswith('mt-')
    adata.var['ribo'] = adata.var_names.str.startswith(('Rps', 'Rpl'))
    
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True)
    
    logger.info(f"  QC前: {adata.n_obs} 细胞")
    
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_cells(adata, max_genes=5000)
    sc.pp.filter_genes(adata, min_cells=3)
    
    adata = adata[adata.obs.pct_counts_mt < 20, :].copy()
    
    logger.info(f"  QC后: {adata.n_obs} 细胞 × {adata.n_vars} 基因")
    logger.info(f"    Sham: {(adata.obs['condition']=='Sham').sum()}, "
                f"MCAO: {(adata.obs['condition']=='MCAO').sum()}")
    
    # ---- 4. 归一化与降维 ----
    logger.info("[4/7] 归一化与降维...")
    
    # 保存原始counts供scTenifoldKnk使用
    adata.raw = adata
    
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key='batch')
    logger.info(f"  高变基因: {adata.var.highly_variable.sum()}")
    
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=50, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.8)
    
    logger.info(f"  Leiden聚类数: {adata.obs['leiden'].nunique()}")
    
    # ---- 5. 细胞类型注释 ----
    logger.info("[5/7] 细胞类型注释...")
    cell_type_map = annotate_cell_types(adata)
    
    # ---- 6. 铜死亡评分 ----
    logger.info("[6/7] 铜死亡评分...")
    cupro_genes, cupro_related = compute_cuproptosis_score(adata)
    cupro_df = analyze_cuproptosis_by_celltype(adata)
    effect_types = find_effect_cell_types(cupro_df)
    
    # ---- 7. 提取贡献基因 ----
    logger.info("[7/7] 提取种子池贡献基因...")
    marker_genes = extract_marker_genes(adata, effect_types)
    ligands, receptors = analyze_ligand_receptor(adata, effect_types)
    
    # ---- 8. 导出count矩阵供scTenifoldKnk虚拟敲除使用 ----
    logger.info("[8/8] 导出count矩阵供scTenifoldKnk (mtx格式)...")
    
    # 保存原始count矩阵（QC过滤后，归一化前）
    raw_counts = adata.raw[:, adata.var_names].X if adata.raw is not None else adata.X
    
    # 过滤低表达基因（至少在50个细胞中表达）
    if sparse.issparse(raw_counts):
        # CSR矩阵: indices是列索引, bincount得到每列(基因)的非零计数
        n_cells_expressing = np.bincount(raw_counts.indices, minlength=raw_counts.shape[1])
    else:
        n_cells_expressing = np.asarray((raw_counts > 0).sum(axis=0)).flatten()
    
    keep_genes = n_cells_expressing >= 50
    raw_counts_filtered = raw_counts[:, keep_genes].T.tocsc()  # 转置为(基因×细胞)并转为CSC
    kept_genes = adata.var_names[keep_genes].tolist()
    
    logger.info(f"  scTenifoldKnk count矩阵: {len(kept_genes)} 基因 x {raw_counts_filtered.shape[1]} 细胞")
    
    # 保存为Matrix Market格式 (稀疏矩阵，比CSV小100倍)
    mtx_path = os.path.join(STAGE_DIR, "sc_count_matrix.mtx")
    genes_path = os.path.join(STAGE_DIR, "sc_genes.tsv")
    barcodes_path = os.path.join(STAGE_DIR, "sc_barcodes.tsv")
    
    from scipy.io import mmwrite
    # mmwrite不支持中文路径，先切到目标目录再写
    old_cwd = os.getcwd()
    os.chdir(STAGE_DIR)
    mmwrite("sc_count_matrix.mtx", raw_counts_filtered)
    os.chdir(old_cwd)
    
    # 保存基因名和细胞barcode
    pd.DataFrame(kept_genes, columns=['gene']).to_csv(genes_path, sep='\t', index=False, header=False)
    pd.DataFrame(adata.obs_names.tolist(), columns=['cell']).to_csv(barcodes_path, sep='\t', index=False, header=False)
    
    logger.info(f"  mtx已保存: sc_count_matrix.mtx ({len(kept_genes)}基因 x {raw_counts_filtered.shape[1]}细胞)")
    logger.info(f"  基因列表: sc_genes.tsv")
    logger.info(f"  细胞barcode: sc_barcodes.tsv")
    
    # 保存细胞注释
    cell_annot = adata.obs[['condition', 'cell_type', 'leiden']].copy()
    cell_annot.to_csv(os.path.join(STAGE_DIR, "sc_cell_annotations.csv"))
    logger.info(f"  细胞注释已保存: {len(cell_annot)} 细胞")
    
    # ---- 保存结果 ----
    logger.info("保存结果...")
    
    adata.write(os.path.join(STAGE_DIR, "sc_adata.h5ad"), compression='gzip')
    
    cupro_df.to_csv(os.path.join(STAGE_DIR, "cuproptosis_by_celltype.csv"), index=False)
    
    with open(os.path.join(STAGE_DIR, "sc_marker_genes.txt"), 'w', encoding='utf-8') as f:
        f.write("# 效应细胞类型marker基因\n")
        for g in sorted(marker_genes):
            f.write(f"{g}\n")
    
    with open(os.path.join(STAGE_DIR, "cellchat_ligands.txt"), 'w', encoding='utf-8') as f:
        f.write("# 表达的配体\n")
        for g in sorted(ligands):
            f.write(f"{g}\n")
    
    with open(os.path.join(STAGE_DIR, "cellchat_receptors.txt"), 'w', encoding='utf-8') as f:
        f.write("# 表达的受体\n")
        for g in sorted(receptors):
            f.write(f"{g}\n")
    
    with open(os.path.join(STAGE_DIR, "effect_cell_types.txt"), 'w', encoding='utf-8') as f:
        f.write("# 效应细胞类型\n")
        for ct in effect_types:
            f.write(f"{ct}\n")
    
    # ---- 输出摘要 ----
    logger.info("\n" + "=" * 60)
    logger.info("阶段2完成! 摘要:")
    logger.info(f"  细胞总数: {adata.n_obs}")
    logger.info(f"  细胞类型: {list(cell_type_map.values())}")
    logger.info(f"  效应细胞类型: {effect_types}")
    logger.info(f"  Marker基因: {len(marker_genes)}")
    logger.info(f"  配体: {len(ligands)}, 受体: {len(receptors)}")
    logger.info(f"  scTenifoldKnk输入: {len(kept_genes)} 基因 x {raw_counts_filtered.shape[1]} 细胞")
    logger.info(f"  铜死亡评分变化最大的细胞类型: {cupro_df.iloc[0]['cell_type']} "
                f"(Δ={cupro_df.iloc[0]['delta']:+.3f})")
    logger.info("=" * 60)
    
    return adata, cupro_df, effect_types, marker_genes, ligands, receptors


if __name__ == "__main__":
    main()