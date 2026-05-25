# ==================== L1 QualTCA: scVelo RNA Velocity 分析 ====================
# 用途：为 Monocle 3 拟时序提供 RNA velocity 验证
# 输入：GSE174574 预处理后的 h5ad 文件
# 输出：scVelo velocity 图 + 每个细胞的瞬时转录方向

import scanpy as sc
import scvelo as scv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import warnings
warnings.filterwarnings('ignore')

SEED = 42
np.random.seed(SEED)

BASE_DIR = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\ciri-cuproptosis-causal-discovery"
OUTPUT_DIR = os.path.join(BASE_DIR, "results", "L1_QualTCA")
FIGURE_DIR = os.path.join(BASE_DIR, "figures", "L1_QualTCA")
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(FIGURE_DIR, exist_ok=True)

scv.settings.verbosity = 2
scv.settings.set_figure_params('scvelo', dpi=150)

DATA_10X_DIR = r"D:\反向网络药理学\L1 数据集\RNA-seq\GSE174574_10X_organized"
EXTRACT_DIR = r"D:\反向网络药理学\L1 数据集\RNA-seq\GSE174574_extracted"

# ==================== 步骤1: 加载/构建 scRNA 数据 ====================
print("=" * 60)
print("L1 QualTCA: scVelo RNA Velocity 分析")
print("=" * 60)

# 尝试从已有的 h5ad 文件加载
existing_h5ad = [
    os.path.join(BASE_DIR, "results", "L1_phenotype_anchoring", "sc_adata_cuproptosis.h5ad"),
    os.path.join(BASE_DIR, "results", "L1_phenotype_anchoring", "pseudotime_adata.h5ad"),
]
loaded = False
for h5ad_path in existing_h5ad:
    if os.path.exists(h5ad_path):
        print(f"\n从现有 h5ad 加载: {h5ad_path}")
        adata = sc.read_h5ad(h5ad_path)
        loaded = True
        break

if not loaded:
    print("\n未找到预处理 h5ad，从原始 10X 数据构建...")

    # 加载各个样本
    adatas = []
    sample_dirs = []
    if os.path.exists(DATA_10X_DIR):
        sample_dirs = [os.path.join(DATA_10X_DIR, d) for d in os.listdir(DATA_10X_DIR)
                       if os.path.isdir(os.path.join(DATA_10X_DIR, d))]
    elif os.path.exists(EXTRACT_DIR):
        # 从提取目录按前缀分组
        prefixes = set()
        for f in os.listdir(EXTRACT_DIR):
            if f.endswith("_matrix.mtx.gz"):
                prefix = f.replace("_matrix.mtx.gz", "")
                prefixes.add(prefix)
        sample_dirs = list(prefixes)

    if len(sample_dirs) == 0:
        print("错误：未找到 10X 数据目录")
        print("请先运行 scrna_analysis.py 生成预处理文件")
        import sys
        sys.exit(1)

    for sample_dir in sample_dirs:
        sample_name = os.path.basename(sample_dir) if os.path.isdir(sample_dir) else sample_dir
        print(f"\n  加载样本: {sample_name}")

        try:
            if os.path.isdir(sample_dir):
                adata_sample = sc.read_10x_mtx(
                    sample_dir,
                    var_names='gene_symbols',
                    cache=False
                )
            else:
                # 从提取目录加载
                import gzip
                from scipy.io import mmread, mmwrite
                from scipy.sparse import csr_matrix

                matrix_path = os.path.join(EXTRACT_DIR, f"{sample_dir}_matrix.mtx.gz")
                genes_path = os.path.join(EXTRACT_DIR, f"{sample_dir}_genes.tsv.gz")
                barcodes_path = os.path.join(EXTRACT_DIR, f"{sample_dir}_barcodes.tsv.gz")

                if not all(os.path.exists(p) for p in [matrix_path, genes_path, barcodes_path]):
                    continue

                with gzip.open(genes_path, 'rt') as f:
                    gene_names = [line.strip().split('\t')[1] if '\t' in line else line.strip()
                                  for line in f]
                with gzip.open(barcodes_path, 'rt') as f:
                    barcodes = [line.strip() for line in f]

                matrix = mmread(matrix_path).tocsr().T
                adata_sample = sc.AnnData(matrix)
                adata_sample.var_names = gene_names
                adata_sample.obs_names = [f"{bc}_{sample_dir}" for bc in barcodes]

            adata_sample.obs['sample'] = sample_name
            condition = 'MCAO' if 'MCAO' in sample_name.upper() else 'Sham'
            adata_sample.obs['condition'] = condition
            adatas.append(adata_sample)
            print(f"    {adata_sample.n_obs} cells, {adata_sample.n_vars} genes")

        except Exception as e:
            print(f"    加载失败: {e}")
            continue

    if len(adatas) == 0:
        print("错误：未能加载任何样本")
        import sys
        sys.exit(1)

    # 合并样本
    for adata in adatas:
        adata.obs_names_make_unique()
        adata.var_names_make_unique()
    adata = adatas[0].concatenate(adatas[1:], batch_key='batch', index_unique=None)
    print(f"\n合并后: {adata.n_obs} cells x {adata.n_vars} genes")

    # QC 过滤
    adata.var['mt'] = adata.var_names.str.startswith('mt-')
    adata.var['ribo'] = adata.var_names.str.match('^(Rps|Rpl|Mrps|Mrpl)')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], inplace=True)

    adata = adata[adata.obs.n_genes_by_counts >= 200, :]
    adata = adata[adata.obs.n_genes_by_counts < 5000, :]
    adata = adata[adata.obs.pct_counts_mt < 20, :]
    adata = adata[adata.obs.total_counts > 500, :]
    sc.pp.filter_genes(adata, min_cells=3)
    print(f"QC 后: {adata.n_obs} cells x {adata.n_vars} genes")

    # 标准化 + 高变基因
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]

    # PCA + 邻居 + UMAP
    sc.pp.scale(adata, max_value=10)
    sc.pp.pca(adata, n_comps=50, random_state=SEED)
    sc.pp.neighbors(adata, n_pcs=30, random_state=SEED)
    sc.tl.umap(adata, random_state=SEED)
    sc.tl.leiden(adata, resolution=0.5, random_state=SEED)

    # 细胞类型注释
    cell_type_markers = {
        'Microglia': ['Ptprc', 'Aif1', 'Cx3cr1', 'Tmem119', 'P2ry12', 'C1qa', 'C1qb'],
        'Neuron': ['Snap25', 'Syt1', 'Nefl', 'Rbfox3', 'Syn1', 'Nefm'],
        'Astrocyte': ['Gfap', 'Aqp4', 'Slc1a3', 'Aldh1l1', 'Slc1a2'],
        'Endothelial': ['Pecam1', 'Vwf', 'Cldn5', 'Cdh5', 'Flt1'],
        'Oligodendrocyte': ['Mbp', 'Plp1', 'Mog', 'Mag', 'Opalin'],
        'OPC': ['Pdgfra', 'Vcan', 'Cspg4', 'Olig1', 'Olig2'],
    }

    for ct, genes in cell_type_markers.items():
        present_genes = [g for g in genes if g in adata.var_names]
        if len(present_genes) >= 2:
            sc.tl.score_genes(adata, gene_list=present_genes, score_name=f'score_{ct}')

    score_cols = [f'score_{ct}' for ct in cell_type_markers if f'score_{ct}' in adata.obs.columns]
    if score_cols:
        adata.obs['cell_type'] = adata.obs[score_cols].idxmax(axis=1).str.replace('score_', '')

    print(f"\n细胞类型分布: {adata.obs['cell_type'].value_counts().to_dict()}")

# ==================== 步骤2: scVelo RNA Velocity ====================
print("\n========== scVelo RNA Velocity 分析 ==========")

scv_h5ad_path = os.path.join(OUTPUT_DIR, "scvelo_adata.h5ad")

# 重新加载原始 counts 用于 velocity
# 需要 spliced/unspliced 矩阵
# 注：10X 数据需要 loom 文件或 BAM 文件才能获取 spliced/unspliced counts
# 由于当前仅有 gene expression matrix，使用 scVelo 的 dynamical modeling
# 基于基因表达矩阵的分化潜力估算

loom_dir = r"D:\反向网络药理学\L1 数据集\RNA-seq"

# 方法 A: 尝试从 loom 文件加载
loom_files = []
if os.path.exists(loom_dir):
    for f in os.listdir(loom_dir):
        if f.endswith('.loom'):
            loom_files.append(os.path.join(loom_dir, f))

if len(loom_files) > 0:
    print(f"找到 {len(loom_files)} 个 loom 文件")
    ldata = scv.read(loom_files[0], cache=True)
    for lf in loom_files[1:]:
        ldata = ldata.concatenate(scv.read(lf, cache=True))

    # 合并 spliced/unspliced 到 adata
    adata = scv.utils.merge(adata, ldata)
    print(f"合并 loom 后: {adata.n_obs} cells, layers: {list(adata.layers.keys())}")

else:
    print("\n未找到 loom 文件，使用替代方法估算 RNA velocity...")
    print("将基于 CytoTRACE 基因特征估算分化方向（替代 velocity）")

    # 替代方案：使用扩散图 + CytoTRACE 启发式基因
    # 基于 pre-mRNA 特征基因代替 spliced/unspliced

    # 使用 CytoTRACE 特征基因进行"pseudo-velocity"估计
    cyto_genes = [
        'Ccnd2', 'Cenpf', 'Mki67', 'Mybl2', 'Tyms', 'Top2a', 'Hmgb2',
        'Cks1b', 'Cks2', 'Ccnb1', 'Aurkb', 'Bub3', 'Cdc20', 'Cenpa',
        'Smc4', 'Ube2c', 'Nusap1', 'Cdkn3', 'Prc1', 'Rrm2'
    ]
    cyto_present = [g for g in cyto_genes if g in adata.var_names]
    print(f"  CytoTRACE 特征基因可用: {len(cyto_present)}/{len(cyto_genes)}")

    if len(cyto_present) >= 5:
        sc.tl.score_genes(adata, gene_list=cyto_present, score_name='cyto_score')

    # 基于扩散图估算方向
    sc.tl.diffmap(adata, n_comps=15)

try:
    # 如果有 spliced/unspliced，运行完整的 scVelo 流程
    if 'spliced' in adata.layers or 'Ms' in adata.layers:
        print("\n运行 scVelo 标准流程...")

        scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
        scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

        scv.tl.velocity(adata, mode='stochastic')
        scv.tl.velocity_graph(adata)

        # Velocity UMAP
        scv.pl.velocity_embedding_stream(
            adata, basis='umap', color='cell_type',
            save=f'{FIGURE_DIR}/scvelo_stream_celltype.png',
            show=False, dpi=150, title='RNA Velocity (scVelo)'
        )

        scv.pl.velocity_embedding_grid(
            adata, basis='umap', color='cell_type',
            save=f'{FIGURE_DIR}/scvelo_grid_celltype.png',
            show=False, dpi=150
        )

        scv.pl.velocity_embedding_stream(
            adata, basis='umap', color='condition',
            save=f'{FIGURE_DIR}/scvelo_stream_condition.png',
            show=False, dpi=150
        )

        # 速度一致性
        scv.tl.velocity_confidence(adata)
        scv.pl.scatter(
            adata, basis='umap', color='velocity_confidence',
            cmap='coolwarm', perc=[2, 98],
            save=f'{FIGURE_DIR}/scvelo_confidence.png',
            show=False, dpi=150
        )

        # 提取 velocity 向量
        velocity_vectors = adata.obsm['velocity_umap']
        velocity_df = pd.DataFrame(
            velocity_vectors,
            index=adata.obs_names,
            columns=['v_umap1', 'v_umap2']
        )
        velocity_df['velocity_magnitude'] = np.sqrt(
            velocity_df['v_umap1']**2 + velocity_df['v_umap2']**2
        )
        velocity_df.to_csv(os.path.join(OUTPUT_DIR, "scvelo_vectors.csv"))
        print(f"  Velocity 向量已保存: {len(velocity_df)} cells")

        # 按细胞类型统计 velocity 方向
        velocity_df['cell_type'] = adata.obs['cell_type'].values
        for ct in velocity_df['cell_type'].unique():
            ct_data = velocity_df[velocity_df['cell_type'] == ct]
            print(f"  {ct}: mean velocity magnitude = {ct_data['velocity_magnitude'].mean():.4f}")

    else:
        # 无 spliced/unspliced → 使用替代方案
        print("\n运行替代 Velocity 分析（Diffusion Pseudotime + CytoTRACE）...")

        # DPT 扩散拟时序
        sc.tl.dpt(adata, n_branchings=1)

        # 基于 DPT 方向信息的 pseudo-velocity 图
        fig, axes = plt.subplots(1, 2, figsize=(16, 7))

        sc.pl.umap(adata, color='dpt_pseudotime', ax=axes[0],
                   show=False, title='Diffusion Pseudotime', cmap='viridis')

        sc.pl.umap(adata, color='cell_type', ax=axes[1],
                   show=False, title='Cell Types', legend_loc='right margin')

        plt.tight_layout()
        plt.savefig(os.path.join(FIGURE_DIR, 'pseudo_velocity_dpt.pdf'),
                    dpi=300, bbox_inches='tight')
        plt.close()

        # 保存 DPT 结果
        dpt_df = pd.DataFrame({
            'cell_barcode': adata.obs_names,
            'dpt_pseudotime': adata.obs['dpt_pseudotime'].values,
            'cell_type': adata.obs['cell_type'].values,
            'condition': adata.obs['condition'].values if 'condition' in adata.obs else '',
        })
        dpt_df.to_csv(os.path.join(OUTPUT_DIR, "dpt_pseudotime.csv"), index=False)
        print(f"  DPT 拟时序已保存: {len(dpt_df)} cells")

        # 按细胞类型的拟时序分布
        print("\n  细胞类型拟时序分布:")
        for ct in sorted(adata.obs['cell_type'].unique()):
            ct_pt = adata[adata.obs['cell_type'] == ct].obs['dpt_pseudotime']
            print(f"    {ct}: mean={ct_pt.mean():.3f}, median={ct_pt.median():.3f}, "
                  f"range=[{ct_pt.min():.3f}, {ct_pt.max():.3f}]")

except Exception as e:
    print(f"  Velocity 分析出错: {e}")
    import traceback
    traceback.print_exc()

# ==================== 步骤3: E/M/L 分期映射 ====================
print("\n========== E/M/L 分期映射 ==========")

for ct in ['Microglia', 'Neuron', 'Astrocyte']:
    if ct in adata.obs['cell_type'].values:
        ct_mask = adata.obs['cell_type'] == ct
        ct_adata = adata[ct_mask].copy()

        if 'dpt_pseudotime' in ct_adata.obs:
            pt_vals = ct_adata.obs['dpt_pseudotime'].values
            pt_vals = pt_vals[np.isfinite(pt_vals)]

            if len(pt_vals) > 30:
                # 三等分
                breaks = np.percentile(pt_vals, [0, 33.33, 66.67, 100])
                labels = ['E', 'M', 'L']
                stage = pd.cut(pt_vals, bins=breaks, labels=labels, include_lowest=True)

                stage_counts = pd.Series(stage).value_counts()
                print(f"  {ct}: E={stage_counts.get('E', 0)}, M={stage_counts.get('M', 0)}, "
                      f"L={stage_counts.get('L', 0)}")

                # 保存分期结果
                stage_df = pd.DataFrame({
                    'cell_barcode': ct_adata.obs_names[ct_adata.obs['dpt_pseudotime'].notna()],
                    'cell_type': ct,
                    'pseudotime': pt_vals,
                    'stage': stage.values
                })
                stage_out = os.path.join(OUTPUT_DIR, f"pseudotime_stages_{ct}.csv")
                stage_df.to_csv(stage_out, index=False)

# 汇总所有细胞类型
stage_files = [f for f in os.listdir(OUTPUT_DIR) if f.startswith('pseudotime_stages_')]
all_stages = []
for sf in stage_files:
    df = pd.read_csv(os.path.join(OUTPUT_DIR, sf))
    all_stages.append(df)

if all_stages:
    combined_stages = pd.concat(all_stages, ignore_index=True)
    combined_stages.to_csv(os.path.join(OUTPUT_DIR, "pseudotime_stages_all.csv"), index=False)
    print(f"\n  综合分期结果: {len(combined_stages)} cells")
    print(f"  E: {sum(combined_stages['stage']=='E')}, "
          f"M: {sum(combined_stages['stage']=='M')}, "
          f"L: {sum(combined_stages['stage']=='L')}")

# ==================== 步骤4: 铜死亡基因沿拟时序表达 ====================
print("\n========== 铜死亡基因沿拟时序表达趋势 ==========")

CUPROPTOSIS_GENES = [
    'Fdx1', 'Lias', 'Lipt1', 'Lipt2', 'Dld', 'Dlat', 'Dlst',
    'Pdha1', 'Pdhb', 'Dbt', 'Mtf1', 'Nfe2l2', 'Nlrp3', 'Gls',
    'Cdkn2a', 'Slc31a1', 'Atp7a', 'Atp7b', 'Gcsh', 'Atox1',
    'Ccs', 'Cox17', 'Cox11', 'Sco1', 'Sco2', 'Mt1', 'Mt2',
    'Cp', 'Commd1', 'Sod1', 'Sod3', 'Slc31a2', 'Slc11a2', 'Steap3'
]

cupro_present = [g for g in CUPROPTOSIS_GENES if g in adata.var_names]
print(f"  铜死亡基因可用: {len(cupro_present)}/{len(CUPROPTOSIS_GENES)}")

if 'dpt_pseudotime' in adata.obs and len(cupro_present) >= 5:
    # 分箱分析
    n_bins = 10
    adata.obs['pseudotime_bin'] = pd.cut(
        adata.obs['dpt_pseudotime'], bins=n_bins, labels=range(1, n_bins + 1)
    )

    # 每个 bin 的铜死亡基因平均表达
    bin_means = {}
    for gene in cupro_present:
        if gene in adata.var_names:
            bin_expr = adata[:, gene].X.toarray().flatten() if hasattr(
                adata[:, gene].X, 'toarray') else adata[:, gene].X.flatten()
            bin_means[gene] = pd.Series(bin_expr).groupby(
                adata.obs['pseudotime_bin'].values
            ).mean()

    bin_df = pd.DataFrame(bin_means)
    bin_df.to_csv(os.path.join(OUTPUT_DIR, "cuproptosis_genes_pseudotime_bins.csv"))
    print(f"  铜死亡基因拟时序分箱表达: {bin_df.shape[1]} genes × {n_bins} bins")

    # 热图
    if bin_df.shape[1] >= 5:
        fig, ax = plt.subplots(figsize=(14, 10))
        from matplotlib.colors import TwoSlopeNorm
        norm = TwoSlopeNorm(vmin=-2, vcenter=0, vmax=2)
        im = ax.imshow(bin_df.T.values, aspect='auto', cmap='RdBu_r', norm=norm)
        ax.set_xticks(range(n_bins))
        ax.set_xticklabels(range(1, n_bins + 1))
        ax.set_xlabel('Pseudotime Bin (E → L)')
        ax.set_yticks(range(len(bin_df.columns)))
        ax.set_yticklabels(bin_df.columns, fontsize=8)
        ax.set_title('Cuproptosis Genes along Pseudotime')
        plt.colorbar(im, ax=ax, label='Mean Expression (z-score)')
        plt.tight_layout()
        plt.savefig(os.path.join(FIGURE_DIR, 'cuproptosis_pseudotime_heatmap.pdf'),
                    dpi=300, bbox_inches='tight')
        plt.close()
        print("  铜死亡拟时序热图已保存")

# ==================== 保存 ====================
print("\n========== 保存结果 ==========")
adata.write(os.path.join(OUTPUT_DIR, "scvelo_adata.h5ad"), compression='gzip')

# 保存到 R 可读的 CSV
if 'dpt_pseudotime' in adata.obs:
    sc_summary = pd.DataFrame({
        'cell_barcode': adata.obs_names,
        'cell_type': adata.obs['cell_type'].values,
        'condition': adata.obs['condition'].values if 'condition' in adata.obs else 'Unknown',
        'pseudotime': adata.obs['dpt_pseudotime'].values,
    })
    sc_summary.to_csv(os.path.join(OUTPUT_DIR, "scRNA_pseudotime_summary.csv"), index=False)

print(f"\n结果保存至: {OUTPUT_DIR}")
print(f"图表保存至: {FIGURE_DIR}")
print("\n========== scVelo RNA Velocity 分析完成 ==========")