# ==================== L1 表型锚定：CIRI-铜死亡急性期差异景观（24h）====================
# 输入数据：GSE174574（24h MCAO vs Sham，scRNA-seq）
# 主方案：Scanpy 标准流程 → Leiden 分群 → 人工注释 → 各亚群 Wilcoxon/MAST 差异分析
# 输出：24h 细胞类型特异性铜死亡差异表达谱

import scanpy as sc
import anndata as ann
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os
import warnings
warnings.filterwarnings('ignore')

# ==================== 0. 配置参数 ====================
SEED = 42
np.random.seed(SEED)

CUPROPTOSIS_CORE = ['FDX1', 'LIAS', 'LIPT1', 'DLD', 'DLAT', 'PDHA1', 'PDHB', 'MTF1', 'GLS', 'CDKN2A']
CUPROPTOSIS_EXTENDED = ['SIRT7', 'ATP7B', 'SLC31A1', 'COX17', 'ATOX1', 'CCS']
CUPROPTOSIS_ALL = CUPROPTOSIS_CORE + CUPROPTOSIS_EXTENDED

OUTPUT_DIR = '../results/L1_phenotype_anchoring'
FIGURE_DIR = '../figures/L1'
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(FIGURE_DIR, exist_ok=True)

# ==================== 1. 数据加载与预处理 ====================
def load_and_preprocess_scrna(h5ad_path):
    print("=== 步骤1: 加载 scRNA-seq 数据 ===")
    
    adata = sc.read_h5ad(h5ad_path)
    print(f"原始数据: {adata.n_obs} cells x {adata.n_vars} genes")
    
    # 质量控制
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    # 过滤低质量细胞
    adata = adata[adata.obs.n_genes_by_counts < 2500, :]
    adata = adata[adata.obs.pct_counts_mt < 20, :]
    adata = adata[adata.obs.n_counts > 500, :]
    
    print(f"质控后: {adata.n_obs} cells x {adata.n_vars} genes")
    
    # 标准化
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # 高变基因
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    
    # PCA
    sc.pp.pca(adata, random_state=SEED)
    sc.pp.neighbors(adata, random_state=SEED)
    
    # Leiden 分群
    sc.tl.leiden(adata, resolution=0.5, random_state=SEED)
    
    # UMAP
    sc.tl.umap(adata, random_state=SEED)
    
    return adata

# ==================== 2. 细胞类型注释 ====================
def annotate_cell_types(adata, marker_genes=None):
    print("\n=== 步骤2: 细胞类型注释 ===")
    
    if marker_genes is None:
        marker_genes = {
            'Microglia': ['PTPRC', 'AIF1', 'CX3CR1', 'TMEM119', 'P2RY12'],
            'Neuron': ['SNAP25', 'SYT1', 'NEFL', 'NEFM', 'RBFOX3'],
            'Astrocyte': ['GFAP', 'AQP4', 'SLC1A3', 'ALDH1L1', 'GFAP'],
            'Endothelial': ['PECAM1', 'VWF', 'CLDN5', 'CDH5', 'FLT1'],
            'Oligodendrocyte': ['MBP', 'PLP1', 'MOG', 'MAG', 'OPALIN'],
            'OPC': ['PDGFRA', 'VCAN', 'CSPG4', 'OLIG1', 'OLIG2'],
            'Pericyte': ['RGS5', 'PDGFRB', 'ACTA2', 'MCAM', 'CSPG4'],
        }
    
    # 计算各细胞类型的 marker 基因评分
    sc.tl.score_genes(adata, gene_list=marker_genes['Microglia'], score_name='Microglia_score')
    sc.tl.score_genes(adata, gene_list=marker_genes['Neuron'], score_name='Neuron_score')
    sc.tl.score_genes(adata, gene_list=marker_genes['Astrocyte'], score_name='Astrocyte_score')
    sc.tl.score_genes(adata, gene_list=marker_genes['Endothelial'], score_name='Endothelial_score')
    sc.tl.score_genes(adata, gene_list=marker_genes['Oligodendrocyte'], score_name='Oligodendrocyte_score')
    
    # 基于 marker 评分自动注释
    adata.obs['cell_type'] = 'Unknown'
    for ctype in marker_genes.keys():
        score_col = f'{ctype}_score'
        if score_col in adata.obs.columns:
            threshold = adata.obs[score_col].quantile(0.75)
            adata.obs.loc[adata.obs[score_col] > threshold, 'cell_type'] = ctype
    
    # 统计各细胞类型数量
    print("\n细胞类型分布:")
    print(adata.obs['cell_type'].value_counts())
    
    return adata

# ==================== 3. 差异表达分析 ====================
def differential_expression(adata, condition_key='condition', reference='Sham', groupby='cell_type'):
    print("\n=== 步骤3: 差异表达分析 ===")
    
    cell_types = adata.obs[groupby].unique()
    results = {}
    
    for ct in cell_types:
        print(f"\n--- {ct} ---")
        ct_data = adata[adata.obs[groupby] == ct].copy()
        
        if len(ct_data) < 10:
            print(f"  细胞数 < 10，跳过")
            continue
            
        # Wilcoxon 检验
        sc.tl.rank_genes_groups(ct_data, groupby=condition_key, reference=reference, method='wilcoxon')
        
        # 提取铜死亡基因差异
        de_results = sc.get.rank_genes_groups_df(ct_data, group='MCAO')
        de_results['cell_type'] = ct
        
        # 筛选铜死亡基因
        cupro_de = de_results[de_results['names'].isin(CUPROPTOSIS_ALL)]
        
        # 显著性筛选 |log2FC|>0.25, FDR<0.05
        sig_cupro = cupro_de[(abs(cupro_de['logfoldchanges']) > 0.25) & (cupro_de['pvals_adj'] < 0.05)]
        
        results[ct] = {
            'all_de': de_results,
            'cupro_de': cupro_de,
            'sig_cupro': sig_cupro
        }
        
        print(f"  铜死亡差异基因: {len(sig_cupro)} (上调: {sum(sig_cupro['logfoldchanges'] > 0)}, 下调: {sum(sig_cupro['logfoldchanges'] < 0)})")
    
    return results

# ==================== 4. 可视化 ====================
def plot_violin_plots(adata, genes, groupby='cell_type', condition_key='condition'):
    print("\n=== 步骤4: 生成可视化 ===")
    
    # 铜死亡基因小提琴图
    fig, axes = plt.subplots(2, 4, figsize=(20, 10))
    axes = axes.flatten()
    
    for i, gene in enumerate(genes[:8]):
        if gene in adata.var_names:
            sc.pl.violin(adata, keys=gene, groupby=groupby, color=condition_key, 
                        ax=axes[i], show=False, title=gene)
        else:
            axes[i].text(0.5, 0.5, f'{gene}\nNot detected', ha='center', va='center')
    
    plt.tight_layout()
    plt.savefig(f'{FIGURE_DIR}/cuproptosis_violin_plots.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    
    # UMAP 铜死亡基因表达
    sc.pl.umap(adata, color=CUPROPTOSIS_CORE[:6], vmin=0, vmax='p99', 
               save='_cuproptosis_core.pdf', show=False)
    
    print(f"  可视化已保存至 {FIGURE_DIR}")

# ==================== 5. 输出结果 ====================
def save_results(de_results, output_dir):
    print("\n=== 步骤5: 保存结果 ===")
    
    # 汇总所有细胞类型的铜死亡差异基因
    all_cupro_results = []
    for ct, res in de_results.items():
        if 'sig_cupro' in res and len(res['sig_cupro']) > 0:
            df = res['sig_cupro'].copy()
            df['cell_type'] = ct
            all_cupro_results.append(df)
    
    if all_cupro_results:
        combined = pd.concat(all_cupro_results)
        combined.to_csv(f'{output_dir}/cuproptosis_DEGs_24h.csv', index=False)
        print(f"  铜死亡差异基因汇总: {len(combined)} 条记录")
    
    # 各亚群差异基因详情
    for ct, res in de_results.items():
        if 'all_de' in res:
            res['all_de'].to_csv(f'{output_dir}/{ct}_DEGs.csv', index=False)
    
    # 铜死亡基因在各细胞类型的表达矩阵
    print("  结果已保存")

# ==================== 主流程 ====================
def main():
    print("=" * 60)
    print("L1 表型锚定：CIRI-铜死亡急性期差异景观（24h）")
    print("=" * 60)
    
    # 数据路径
    H5AD_PATH = '../data/GSE174574_24h_processed.h5ad'
    
    if not os.path.exists(H5AD_PATH):
        print(f"错误: 未找到数据文件 {H5AD_PATH}")
        print("请先下载 GSE174574 数据并转换为 h5ad 格式")
        return
    
    # 1. 加载与预处理
    adata = load_and_preprocess_scrna(H5AD_PATH)
    
    # 2. 细胞类型注释
    adata = annotate_cell_types(adata)
    
    # 3. 差异表达分析
    de_results = differential_expression(adata)
    
    # 4. 可视化
    plot_violin_plots(adata, CUPROPTOSIS_ALL)
    
    # 5. 保存结果
    save_results(de_results, OUTPUT_DIR)
    
    # 自检标准
    print("\n" + "=" * 60)
    print("自检标准检查")
    print("=" * 60)
    
    total_sig = sum(len(res['sig_cupro']) for res in de_results.values() if 'sig_cupro' in res)
    print(f"✓ scRNA-seq 铜死亡差异基因交集 N = {total_sig} (要求 ≥ 5)")
    
    if total_sig < 5:
        print("⚠ 警告: 差异基因数 < 5，建议扩展基因集或降低阈值")
    
    print("\nL1 分析完成！")

if __name__ == '__main__':
    main()