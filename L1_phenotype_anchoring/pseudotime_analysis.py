# ==================== L1 拟时序分析：CIRI-铜死亡轨迹推断 ====================
# 参考：Monocle 3 方法学（微信公众号推文）
# 输入数据：GSE174574（24h MCAO vs Sham，scRNA-seq 10X格式）
# 方案：Scanpy预处理 → PAGA轨迹图 → DPT扩散拟时序 → scVelo RNA速率 → 铜死亡基因沿轨迹表达
# 输出：拟时序轨迹图、铜死亡基因伪时间表达热图、RNA速率图

import scanpy as sc
import anndata as ann
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
from scipy import stats
from scipy.sparse import issparse
import os
import gzip
import shutil
import warnings
warnings.filterwarnings('ignore')

sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=100, frameon=False, figsize=(6, 6), facecolor='white')

SEED = 42
np.random.seed(SEED)

# ==================== 铜死亡基因集 ====================
CUPROPTOSIS_CORE = ['FDX1', 'LIAS', 'LIPT1', 'DLD', 'DLAT', 'PDHA1', 'PDHB', 'MTF1', 'GLS', 'CDKN2A']
CUPROPTOSIS_EXTENDED = ['SIRT7', 'ATP7B', 'SLC31A1', 'COX17', 'ATOX1', 'CCS']
CUPROPTOSIS_ALL = CUPROPTOSIS_CORE + CUPROPTOSIS_EXTENDED

# Mouse gene names (capitalize first letter, rest lowercase)
CUPROPTOSIS_MOUSE = [g.capitalize() for g in CUPROPTOSIS_ALL]

DATA_DIR = r'D:\反向网络药理学\L1 数据集\RNA-seq'
RAW_DIR = os.path.join(DATA_DIR, 'GSE174574_extracted')
TENX_DIR = os.path.join(DATA_DIR, 'GSE174574_10X_organized')

OUTPUT_DIR = '../results/L1_pseudotime'
FIGURE_DIR = '../figures/L1_pseudotime'
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(FIGURE_DIR, exist_ok=True)

MARKER_GENES = {
    'Microglia': ['Ptprc', 'Aif1', 'Cx3cr1', 'Tmem119', 'P2ry12'],
    'Neuron': ['Snap25', 'Syt1', 'Nefl', 'Nefm', 'Rbfox3'],
    'Astrocyte': ['Gfap', 'Aqp4', 'Slc1a3', 'Aldh1l1'],
    'Endothelial': ['Pecam1', 'Vwf', 'Cldn5', 'Cdh5', 'Flt1'],
    'Oligodendrocyte': ['Mbp', 'Plp1', 'Mog', 'Mag', 'Opalin'],
    'OPC': ['Pdgfra', 'Vcan', 'Cspg4', 'Olig1', 'Olig2'],
    'Pericyte': ['Rgs5', 'Pdgfrb', 'Acta2', 'Mcam', 'Cspg4'],
}

SAMPLE_INFO = [
    ('GSM5319987_sham1', 'Sham', 'Sham_1'),
    ('GSM5319988_sham2', 'Sham', 'Sham_2'),
    ('GSM5319989_sham3', 'Sham', 'Sham_3'),
    ('GSM5319990_MCAO1', 'MCAO', 'MCAO_1'),
    ('GSM5319991_MCAO2', 'MCAO', 'MCAO_2'),
    ('GSM5319992_MCAO3', 'MCAO', 'MCAO_3'),
]


def organize_and_load_samples():
    print("=== 步骤0-1: 加载scRNA-seq数据 ===")
    adatas = []

    for file_prefix, condition, sample_name in SAMPLE_INFO:
        barcodes_src = os.path.join(RAW_DIR, f'{file_prefix}_barcodes.tsv.gz')
        genes_src = os.path.join(RAW_DIR, f'{file_prefix}_genes.tsv.gz')
        matrix_src = os.path.join(RAW_DIR, f'{file_prefix}_matrix.mtx.gz')

        if not os.path.exists(matrix_src):
            print(f"  警告: {sample_name} 数据不完整，跳过")
            continue

        try:
            with gzip.open(barcodes_src, 'rt') as f:
                barcodes = [line.strip() for line in f.readlines()]
            with gzip.open(genes_src, 'rt') as f:
                gene_info = [line.strip().split('\t') for line in f.readlines()]
            gene_ids = [g[0] for g in gene_info]
            gene_names = [g[1] if len(g) > 1 else g[0] for g in gene_info]

            adata = sc.read_mtx(matrix_src)
            adata = adata.T
            adata.obs_names = barcodes
            adata.var_names = gene_ids
            adata.var['gene_name'] = gene_names
            adata.var_names_make_unique()

            adata.obs['sample'] = sample_name
            adata.obs['condition'] = condition
            adata.obs_names = [f"{sample_name}_{bc}" for bc in adata.obs_names]

            adatas.append(adata)
            print(f"  {sample_name} ({condition}): {adata.n_obs} cells, {adata.n_vars} genes")
        except Exception as e:
            import traceback
            print(f"  错误加载 {sample_name}: {e}")
            traceback.print_exc()

    if len(adatas) == 0:
        raise FileNotFoundError("没有成功加载任何样本数据")

    adata_merged = ann.concat(adatas, axis=0, join='outer')
    print(f"\n合并后: {adata_merged.n_obs} cells x {adata_merged.n_vars} genes")
    return adata_merged


def preprocess_data(adata):
    print("\n=== 步骤2: 数据预处理 ===")
    adata.var['mt'] = adata.var_names.str.startswith('mt-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    adata = adata[adata.obs['n_genes_by_counts'] < 4000, :].copy()
    adata = adata[adata.obs['pct_counts_mt'] < 20, :].copy()

    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    adata = adata[adata.obs['total_counts'] > 500, :].copy()

    print(f"QC后: {adata.n_obs} cells x {adata.n_vars} genes")

    adata.layers['counts'] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata.copy()

    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat_v3')

    if 'gene_name' in adata.var.columns:
        cupro_ensg = []
        name_lower = adata.var['gene_name'].str.lower()
        for g in CUPROPTOSIS_MOUSE:
            matches = adata.var_names[name_lower == g.lower()]
            if len(matches) > 0:
                cupro_ensg.extend(matches.tolist())
        adata.var.loc[adata.var_names.isin(cupro_ensg), 'highly_variable'] = True
        print(f"  保留 {len(cupro_ensg)} 个铜死亡基因在高变基因集中")

    adata = adata[:, adata.var['highly_variable']].copy()

    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack', n_comps=50, random_state=SEED)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, random_state=SEED)
    sc.tl.umap(adata, random_state=SEED)
    sc.tl.leiden(adata, resolution=0.5, random_state=SEED)

    return adata


def annotate_cell_types(adata):
    print("\n=== 步骤3: 细胞类型注释 ===")

    score_adata = adata.copy()
    for ct, genes in MARKER_GENES.items():
        valid_genes = [g for g in genes if g in score_adata.var_names]
        if len(valid_genes) >= 2:
            try:
                sc.tl.score_genes(score_adata, gene_list=valid_genes, score_name=f'{ct}_score')
                adata.obs[f'{ct}_score'] = score_adata.obs[f'{ct}_score'].values
            except Exception as e:
                print(f"  {ct} 评分失败: {e}")

    adata.obs['cell_type'] = 'Unannotated'
    for ct in MARKER_GENES.keys():
        score_col = f'{ct}_score'
        if score_col in adata.obs.columns:
            threshold = adata.obs[score_col].quantile(0.75)
            high_score = adata.obs[score_col] > threshold
            adata.obs.loc[high_score & (adata.obs['cell_type'] == 'Unannotated'), 'cell_type'] = ct

    n_per_type = adata.obs['cell_type'].value_counts()
    if n_per_type.get('Unannotated', 0) == adata.n_obs:
        print("  Marker基因注释未成功，使用Leiden聚类结果")
        adata.obs['cell_type'] = 'Cluster_' + adata.obs['leiden'].astype(str)

    print("细胞类型分布:")
    print(adata.obs['cell_type'].value_counts())
    return adata


def trajectory_inference(adata):
    print("\n=== 步骤4: 轨迹推断（DPT拟时序 + PAGA）===")

    n_cell_types = adata.obs['cell_type'].nunique()
    if n_cell_types > 1:
        try:
            sc.tl.paga(adata, groups='cell_type')
            sc.pl.paga(adata, plot=False)
            sc.tl.umap(adata, init_pos='paga', random_state=SEED)
        except Exception as e:
            print(f"  PAGA失败（{e}），使用标准UMAP")

    sham_indices = np.flatnonzero(adata.obs['condition'] == 'Sham')
    if len(sham_indices) > 0:
        adata.uns['iroot'] = sham_indices[0]
    else:
        adata.uns['iroot'] = 0

    sc.tl.diffmap(adata, random_state=SEED)
    sc.tl.dpt(adata, n_dcs=15)

    print(f"DPT拟时序范围: {adata.obs['dpt_pseudotime'].min():.3f} - {adata.obs['dpt_pseudotime'].max():.3f}")
    print(f"Sham平均拟时序: {adata[adata.obs['condition'] == 'Sham'].obs['dpt_pseudotime'].mean():.3f}")
    print(f"MCAO平均拟时序: {adata[adata.obs['condition'] == 'MCAO'].obs['dpt_pseudotime'].mean():.3f}")

    return adata


def rna_velocity_analysis(adata_raw):
    print("\n=== 步骤5: RNA速率分析（scVelo）===")
    try:
        import scvelo as scv
    except ImportError:
        print("scVelo未安装，跳过RNA速率分析")
        return adata_raw

    if 'spliced' not in adata_raw.layers and 'unspliced' not in adata_raw.layers:
        print("数据缺少 spliced/unspliced 层，跳过RNA速率分析（需要 loom 格式的原始数据）")
        return adata_raw

    try:
        adata_raw_for_velo = adata_raw.copy()
        scv.pp.filter_and_normalize(adata_raw_for_velo, min_shared_counts=20, n_top_genes=2000)
        scv.pp.moments(adata_raw_for_velo, n_pcs=30, n_neighbors=30)
        scv.tl.recover_dynamics(adata_raw_for_velo, n_jobs=4)
        scv.tl.velocity(adata_raw_for_velo, mode='dynamical')
        scv.tl.velocity_graph(adata_raw_for_velo)
        scv.tl.latent_time(adata_raw_for_velo)

        adata_raw.obsm['velocity_umap'] = adata_raw_for_velo.obsm.get('velocity_umap', None)
        adata_raw.obs['velocity_latent_time'] = adata_raw_for_velo.obs.get('latent_time', np.nan)
        adata_raw.obs['velocity_pseudotime'] = adata_raw_for_velo.obs.get('velocity_pseudotime', np.nan)
        print("RNA速率分析完成")
    except Exception as e:
        print(f"RNA速率分析失败: {e}")

    return adata_raw


def cuproptosis_pseudotime_analysis(adata):
    print("\n=== 步骤6: 铜死亡基因沿拟时序表达分析 ===")

    gene_name_map = None
    if 'gene_name' in adata.var.columns:
        gene_name_map = dict(zip(adata.var['gene_name'].str.lower(), adata.var_names))
        cupro_genes_in_data = []
        for g in CUPROPTOSIS_MOUSE:
            g_lower = g.lower()
            if g_lower in gene_name_map:
                cupro_genes_in_data.append(gene_name_map[g_lower])
            elif g in adata.var_names:
                cupro_genes_in_data.append(g)

        display_names = []
        for ensg in cupro_genes_in_data:
            for gene_name, eid in gene_name_map.items():
                if eid == ensg:
                    display_names.append(gene_name.capitalize())
                    break
            else:
                display_names.append(ensg)
    else:
        cupro_genes_in_data = [g for g in CUPROPTOSIS_MOUSE if g in adata.var_names]
        display_names = cupro_genes_in_data

    print(f"检测到的铜死亡基因: {len(cupro_genes_in_data)}/{len(CUPROPTOSIS_MOUSE)}")
    print(f"基因列表: {display_names}")

    if len(cupro_genes_in_data) == 0:
        print("警告: 未检测到铜死亡基因，使用前10个高变基因进行分析")
        cupro_genes_in_data = adata.var_names[:10].tolist()
        display_names = cupro_genes_in_data

    X = adata[:, cupro_genes_in_data].X
    if issparse(X):
        X = X.toarray()

    pseudotime = adata.obs['dpt_pseudotime'].values
    condition = adata.obs['condition'].values
    cell_type = adata.obs['cell_type'].values

    n_bins = 20
    bins = np.linspace(pseudotime.min(), pseudotime.max(), n_bins + 1)
    bin_centers = (bins[:-1] + bins[1:]) / 2

    trend_data = {}
    for i, gene in enumerate(cupro_genes_in_data):
        gene_expr = X[:, i]
        sham_trend = []
        mcao_trend = []

        for j in range(n_bins):
            mask_sham = (pseudotime >= bins[j]) & (pseudotime < bins[j + 1]) & (condition == 'Sham')
            mask_mcao = (pseudotime >= bins[j]) & (pseudotime < bins[j + 1]) & (condition == 'MCAO')

            sham_trend.append(np.mean(gene_expr[mask_sham]) if mask_sham.sum() > 0 else np.nan)
            mcao_trend.append(np.mean(gene_expr[mask_mcao]) if mask_mcao.sum() > 0 else np.nan)

        display_name = display_names[i] if i < len(display_names) else gene
        trend_data[display_name] = {'bin_centers': bin_centers, 'Sham': sham_trend, 'MCAO': mcao_trend}

    results_df = []
    for ct in np.unique(cell_type):
        if ct == 'Unannotated':
            continue
        ct_mask = cell_type == ct
        if ct_mask.sum() < 10:
            continue
        for i, gene in enumerate(cupro_genes_in_data):
            gene_expr = X[ct_mask, i]
            pt_ct = pseudotime[ct_mask]
            cond_ct = condition[ct_mask]

            sham_expr = gene_expr[cond_ct == 'Sham']
            mcao_expr = gene_expr[cond_ct == 'MCAO']

            sham_mean = np.mean(sham_expr) if len(sham_expr) > 0 else np.nan
            mcao_mean = np.mean(mcao_expr) if len(mcao_expr) > 0 else np.nan

            pt_corr = np.corrcoef(gene_expr, pt_ct)[0, 1] if len(gene_expr) > 3 else np.nan

            display_name = display_names[i] if i < len(display_names) else gene
            results_df.append({
                'cell_type': ct,
                'gene': display_name,
                'sham_mean': sham_mean,
                'mcao_mean': mcao_mean,
                'log2FC': np.log2(mcao_mean + 0.01) - np.log2(sham_mean + 0.01) if sham_mean > 0 else np.nan,
                'pseudotime_corr': pt_corr,
                'n_cells': ct_mask.sum()
            })

    results_df = pd.DataFrame(results_df)
    results_df.to_csv(os.path.join(OUTPUT_DIR, 'cuproptosis_pseudotime_summary.csv'), index=False)
    print(f"铜死亡拟时序分析结果已保存: {len(results_df)} 条记录")

    return display_names, trend_data, results_df


def visualize_results(adata, cupro_genes, trend_data, results_df):
    print("\n=== 步骤7: 可视化 ===")

    # 图1: UMAP按条件和细胞类型着色
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    sc.pl.umap(adata, color='condition', ax=axes[0], show=False, title='Condition (Sham vs MCAO)',
               palette={'Sham': '#4472C4', 'MCAO': '#ED7D31'})
    sc.pl.umap(adata, color='cell_type', ax=axes[1], show=False, title='Cell Type',
               palette='tab10')
    sc.pl.umap(adata, color='dpt_pseudotime', ax=axes[2], show=False, title='DPT Pseudotime',
               cmap='viridis')
    plt.tight_layout()
    fig.savefig(os.path.join(FIGURE_DIR, 'pseudotime_umap_overview.pdf'), dpi=300, bbox_inches='tight')
    plt.close()
    print("  UMAP概览图已保存")

    # 图2: PAGA轨迹图
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    sc.pl.paga(adata, color='cell_type', ax=ax, show=False, title='PAGA Trajectory Graph',
               node_size_scale=1.5, edge_width_scale=0.5)
    fig.savefig(os.path.join(FIGURE_DIR, 'paga_trajectory_graph.pdf'), dpi=300, bbox_inches='tight')
    plt.close()
    print("  PAGA轨迹图已保存")

    # 图3: 拟时序UMAP
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    sc.pl.umap(adata, color='dpt_pseudotime', ax=axes[0], show=False, title='DPT Pseudotime on UMAP',
               cmap='plasma')
    sc.pl.umap(adata, color='condition', ax=axes[1], show=False, title='Condition',
               palette={'Sham': '#4472C4', 'MCAO': '#ED7D31'})
    plt.tight_layout()
    fig.savefig(os.path.join(FIGURE_DIR, 'pseudotime_on_umap.pdf'), dpi=300, bbox_inches='tight')
    plt.close()
    print("  拟时序UMAP图已保存")

    # 图4: 铜死亡基因沿拟时序的表达趋势
    n_genes = len(cupro_genes)
    n_cols = min(3, n_genes)
    n_rows = int(np.ceil(n_genes / n_cols))

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 4 * n_rows))
    if n_genes == 1:
        axes = np.array([axes])
    axes = axes.flatten()

    for i, gene in enumerate(cupro_genes):
        ax = axes[i]
        td = trend_data[gene]
        ax.plot(td['bin_centers'], td['Sham'], 'o-', color='#4472C4', label='Sham', linewidth=2, markersize=6)
        ax.plot(td['bin_centers'], td['MCAO'], 's-', color='#ED7D31', label='MCAO', linewidth=2, markersize=6)
        ax.set_title(gene, fontsize=12, fontweight='bold')
        ax.set_xlabel('Pseudotime')
        ax.set_ylabel('Mean Expression')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    plt.tight_layout()
    fig.savefig(os.path.join(FIGURE_DIR, 'cuproptosis_pseudotime_trends.pdf'), dpi=300, bbox_inches='tight')
    plt.close()
    print("  铜死亡基因拟时序趋势图已保存")

    # 图5: 铜死亡基因拟时序热图
    significant_genes = results_df[
        (results_df['log2FC'].abs() > 0.25) & (results_df['pseudotime_corr'].abs() > 0.1)
    ]['gene'].unique()

    if len(significant_genes) > 2:
        pivot = results_df.pivot_table(
            values='pseudotime_corr', index='gene', columns='cell_type', aggfunc='mean'
        )
        pivot = pivot.loc[pivot.index.isin(significant_genes)]

        fig, ax = plt.subplots(figsize=(max(8, len(pivot.columns) * 1.2), max(6, len(pivot) * 0.4)))
        sns.heatmap(pivot, annot=True, fmt='.2f', cmap='RdBu_r', center=0,
                    linewidths=0.5, ax=ax, cbar_kws={'label': 'Pseudotime Correlation'})
        ax.set_title('Cuproptosis Gene - Pseudotime Correlation by Cell Type', fontsize=13, fontweight='bold')
        plt.tight_layout()
        fig.savefig(os.path.join(FIGURE_DIR, 'cuproptosis_pseudotime_heatmap.pdf'), dpi=300, bbox_inches='tight')
        plt.close()
        print("  铜死亡基因拟时序热图已保存")

    # 图6: 铜死亡基因UMAP表达
    top_genes = cupro_genes[:min(6, len(cupro_genes))]
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    axes = axes.flatten()
    for i, gene in enumerate(top_genes):
        sc.pl.umap(adata, color=gene, ax=axes[i], show=False, title=gene,
                   cmap='Reds', vmax='p99')
    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)
    plt.tight_layout()
    fig.savefig(os.path.join(FIGURE_DIR, 'cuproptosis_genes_umap.pdf'), dpi=300, bbox_inches='tight')
    plt.close()
    print("  铜死亡基因UMAP表达图已保存")

    # 图7: 拟时序密度分布
    fig, ax = plt.subplots(figsize=(10, 5))
    for cond, color in [('Sham', '#4472C4'), ('MCAO', '#ED7D31')]:
        subset = adata[adata.obs['condition'] == cond].obs['dpt_pseudotime']
        ax.hist(subset, bins=30, alpha=0.5, label=cond, color=color, density=True)
    ax.set_xlabel('DPT Pseudotime')
    ax.set_ylabel('Density')
    ax.set_title('Pseudotime Distribution: Sham vs MCAO')
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.savefig(os.path.join(FIGURE_DIR, 'pseudotime_density.pdf'), dpi=300, bbox_inches='tight')
    plt.close()
    print("  拟时序密度分布图已保存")

    # 图8: RNA速率图（如果可用）
    if 'velocity_umap' in adata.obsm and adata.obsm['velocity_umap'] is not None:
        try:
            import scvelo as scv
            fig, axes = plt.subplots(1, 2, figsize=(14, 5))
            scv.pl.velocity_embedding_stream(adata, basis='umap', color='cell_type',
                                              ax=axes[0], show=False, title='RNA Velocity Stream',
                                              palette='tab10', legend_loc='right margin')
            scv.pl.velocity_embedding_stream(adata, basis='umap', color='dpt_pseudotime',
                                              ax=axes[1], show=False, title='RNA Velocity + Pseudotime',
                                              cmap='viridis', legend_loc='right margin')
            plt.tight_layout()
            fig.savefig(os.path.join(FIGURE_DIR, 'rna_velocity.pdf'), dpi=300, bbox_inches='tight')
            plt.close()
            print("  RNA速率图已保存")
        except Exception as e:
            print(f"  RNA速率图生成失败: {e}")

    print(f"\n所有可视化已保存至: {FIGURE_DIR}")


def save_results(adata, results_df):
    print("\n=== 步骤8: 保存结果 ===")

    adata.write(os.path.join(OUTPUT_DIR, 'pseudotime_analysis.h5ad'), compression='gzip')
    print(f"  AnnData已保存: pseudotime_analysis.h5ad")

    pseudotime_stats = adata.obs.groupby(['condition', 'cell_type'])['dpt_pseudotime'].agg(
        ['mean', 'std', 'min', 'max', 'count']
    ).reset_index()
    pseudotime_stats.to_csv(os.path.join(OUTPUT_DIR, 'pseudotime_stats_by_celltype.csv'), index=False)
    print(f"  拟时序统计已保存")

    top_changes = results_df.sort_values('log2FC', key=abs, ascending=False).head(20)
    top_changes.to_csv(os.path.join(OUTPUT_DIR, 'cuproptosis_top_changes.csv'), index=False)
    print(f"  Top20铜死亡变化基因已保存")


def main():
    print("=" * 70)
    print("L1 拟时序分析：CIRI-铜死亡轨迹推断（参考Monocle 3方法学）")
    print("=" * 70)

    # 0-1. 整理并加载数据
    adata_raw = organize_and_load_samples()

    # 2. 预处理
    adata = preprocess_data(adata_raw)

    # 3. 细胞类型注释
    adata = annotate_cell_types(adata)

    # 4. 轨迹推断（PAGA + DPT）
    adata = trajectory_inference(adata)

    # 5. RNA速率分析
    adata = rna_velocity_analysis(adata)

    # 6. 铜死亡基因拟时序分析
    cupro_genes, trend_data, results_df = cuproptosis_pseudotime_analysis(adata)

    # 7. 可视化
    visualize_results(adata, cupro_genes, trend_data, results_df)

    # 8. 保存结果
    save_results(adata, results_df)

    print("\n" + "=" * 70)
    print("自检标准")
    print("=" * 70)

    n_cells = adata.n_obs
    n_cupro_genes = len(cupro_genes)
    pt_range = adata.obs['dpt_pseudotime'].max() - adata.obs['dpt_pseudotime'].min()
    sig_genes = len(results_df[results_df['log2FC'].abs() > 0.25])

    print(f"✓ 细胞总数: {n_cells} (要求 ≥ 100)")
    print(f"✓ 检测铜死亡基因数: {n_cupro_genes}/{len(CUPROPTOSIS_MOUSE)} (要求 ≥ 5)")
    print(f"✓ 拟时序跨度: {pt_range:.2f} (要求 > 0)")
    print(f"✓ 显著变化铜死亡基因: {sig_genes} (要求 ≥ 3)")

    checks_pass = n_cells >= 100 and n_cupro_genes >= 5 and pt_range > 0 and sig_genes >= 3

    if checks_pass:
        print("\n✓ 所有自检通过！拟时序分析完成。")
    else:
        print("\n⚠ 部分自检未通过，请检查数据质量。")

    return checks_pass


if __name__ == '__main__':
    success = main()
    exit(0 if success else 1)