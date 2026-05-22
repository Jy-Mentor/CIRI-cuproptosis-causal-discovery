# L1 表型锚定：GSE174574 scRNA-seq 24h MCAO vs Sham 差异分析
import scanpy as sc
import anndata
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import tarfile
import gzip
import shutil
import os
from scipy.io import mmread

SEED = 42
np.random.seed(SEED)

DATA_DIR = r"D:\反向网络药理学\L1 数据集\RNA-seq"
RAW_TAR = os.path.join(DATA_DIR, "GSE174574_RAW.tar")
EXTRACT_DIR = os.path.join(DATA_DIR, "GSE174574_extracted")
OUTPUT_DIR = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\CIRI-cuproptosis-causal-discovery\results\L1_phenotype_anchoring"
FIGURE_DIR = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\CIRI-cuproptosis-causal-discovery\figures\L1"
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(FIGURE_DIR, exist_ok=True)

CUPROPTOSIS_CORE = ['Fdx1', 'Lias', 'Lipt1', 'Dld', 'Dlat', 'Pdha1', 'Pdhb', 'Mtf1', 'Gls', 'Cdkn2a']
CUPROPTOSIS_EXTENDED = ['Sirt7', 'Atp7b', 'Slc31a1', 'Cox17', 'Atox1', 'Ccs']
CUPROPTOSIS_ALL = CUPROPTOSIS_CORE + CUPROPTOSIS_EXTENDED

SAMPLES = [
    {"id": "sham1", "condition": "Sham", "files": ["GSM5319987_sham1"]},
    {"id": "sham2", "condition": "Sham", "files": ["GSM5319988_sham2"]},
    {"id": "sham3", "condition": "Sham", "files": ["GSM5319989_sham3"]},
    {"id": "mcao1", "condition": "MCAO", "files": ["GSM5319990_MCAO1"]},
    {"id": "mcao2", "condition": "MCAO", "files": ["GSM5319991_MCAO2"]},
    {"id": "mcao3", "condition": "MCAO", "files": ["GSM5319992_MCAO3"]},
]

def extract_tar():
    print("=== 步骤1: 解压 GSE174574_RAW.tar ===")
    os.makedirs(EXTRACT_DIR, exist_ok=True)
    with tarfile.open(RAW_TAR, "r") as tar:
        tar.extractall(EXTRACT_DIR)
    files = os.listdir(EXTRACT_DIR)
    print(f"  解压 {len(files)} 个文件")
    return files

def load_10x_sample(prefix):
    matrix_path = os.path.join(EXTRACT_DIR, f"{prefix}_matrix.mtx.gz")
    genes_path = os.path.join(EXTRACT_DIR, f"{prefix}_genes.tsv.gz")
    barcodes_path = os.path.join(EXTRACT_DIR, f"{prefix}_barcodes.tsv.gz")
    for p in [matrix_path, genes_path, barcodes_path]:
        if not os.path.exists(p):
            print(f"  警告: 未找到 {p}")
            return None
    with gzip.open(genes_path, 'rt') as f:
        genes = [line.strip().split('\t')[1] for line in f]
    with gzip.open(barcodes_path, 'rt') as f:
        barcodes = [line.strip() for line in f]
    with gzip.open(matrix_path, 'rb') as f:
        matrix = mmread(f).tocsr()
    adata = sc.AnnData(matrix.T)
    adata.var_names = genes
    adata.obs_names = [f"{bc}_{prefix}" for bc in barcodes]
    return adata

def merge_samples():
    print("\n=== 步骤2: 合并所有样本 ===")
    adatas = []
    for sample in SAMPLES:
        prefix = sample["files"][0]
        print(f"  加载 {sample['id']} ({sample['condition']})...")
        adata = load_10x_sample(prefix)
        if adata is not None:
            adata.obs["sample_id"] = sample["id"]
            adata.obs["condition"] = sample["condition"]
            adatas.append(adata)
    if not adatas:
        raise ValueError("未成功加载任何样本")
    for i, adata in enumerate(adatas):
        adata.obs_names_make_unique()
        adata.var_names_make_unique()
    merged = anndata.concat(adatas, axis=0, join="outer", label="batch", keys=[s["id"] for s in SAMPLES])
    merged.obs_names_make_unique()
    print(f"  合并后: {merged.n_obs} cells x {merged.n_vars} genes")
    return merged

def quality_control(adata):
    print("\n=== 步骤3: 质量控制 ===")
    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], inplace=True)
    print(f"  质控前: {adata.n_obs} cells")
    adata = adata[adata.obs.n_genes_by_counts < 5000, :]
    adata = adata[adata.obs.pct_counts_mt < 20, :]
    adata = adata[adata.obs.total_counts > 500, :]
    print(f"  质控后: {adata.n_obs} cells")
    return adata

def normalize_and_cluster(adata):
    print("\n=== 步骤4: 标准化与聚类 ===")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]
    sc.pp.scale(adata, max_value=10)
    sc.pp.pca(adata, random_state=SEED)
    sc.pp.neighbors(adata, random_state=SEED)
    sc.tl.leiden(adata, resolution=0.5, random_state=SEED)
    sc.tl.umap(adata, random_state=SEED)
    print(f"  Leiden clusters: {adata.obs['leiden'].nunique()}")
    return adata

def annotate_cell_types(adata):
    print("\n=== 步骤5: 细胞类型注释 ===")
    markers = {
        "Microglia": ["Ptprc", "Aif1", "Cx3cr1", "Tmem119", "P2ry12", "C1qa", "C1qb"],
        "Neuron": ["Snap25", "Syt1", "Nefl", "Rbfox3", "Syn1"],
        "Astrocyte": ["Gfap", "Aqp4", "Slc1a3", "Aldh1l1"],
        "Endothelial": ["Pecam1", "Vwf", "Cldn5", "Cdh5"],
        "Oligodendrocyte": ["Mbp", "Plp1", "Mog", "Mag"],
        "OPC": ["Pdgfra", "Vcan", "Cspg4", "Olig1", "Olig2"],
    }
    adata.obs["cell_type"] = "Unknown"
    for ctype, genes in markers.items():
        present = [g for g in genes if g in adata.var_names]
        if present:
            scores = adata[:, present].X.mean(axis=1)
            if hasattr(scores, "A1"):
                scores = scores.A1
            threshold = np.percentile(scores, 80)
            adata.obs.loc[np.array(scores).flatten() > threshold, "cell_type"] = ctype
    print("细胞类型分布:")
    print(adata.obs["cell_type"].value_counts())
    return adata

def differential_expression(adata):
    print("\n=== 步骤6: 差异表达分析 ===")
    sc.tl.rank_genes_groups(adata, groupby="condition", reference="Sham", method="wilcoxon")
    de_df = sc.get.rank_genes_groups_df(adata, group="MCAO")
    cupro_de = de_df[de_df["names"].isin(CUPROPTOSIS_ALL)].copy()
    sig_cupro = cupro_de[(abs(cupro_de["logfoldchanges"]) > 0.25) & (cupro_de["pvals_adj"] < 0.05)]
    all_sig = de_df[(abs(de_df["logfoldchanges"]) > 0.25) & (de_df["pvals_adj"] < 0.05)]
    print(f"  总差异基因: {len(all_sig)}")
    print(f"  铜死亡差异基因: {len(sig_cupro)}")
    return de_df, cupro_de, sig_cupro

def plot_results(adata, de_df, cupro_de):
    print("\n=== 步骤7: 可视化 ===")
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    sc.pl.umap(adata, color="condition", ax=axes[0], show=False, title="Condition")
    sc.pl.umap(adata, color="leiden", ax=axes[1], show=False, title="Clusters")
    sc.pl.umap(adata, color="cell_type", ax=axes[2], show=False, title="Cell Types")
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURE_DIR, "umap_overview.pdf"), dpi=300)
    plt.close()
    present_cupro = [g for g in CUPROPTOSIS_CORE if g in adata.var_names]
    if present_cupro:
        sc.pl.umap(adata, color=present_cupro[:6], vmin=0, vmax="p99",
                   save="_cuproptosis.pdf", show=False)
    sc.pl.rank_genes_groups(adata, n_genes=25, save="_deg.pdf", show=False)
    print("  可视化已保存")

def save_outputs(de_df, cupro_de, sig_cupro):
    print("\n=== 步骤8: 保存结果 ===")
    de_df.to_csv(os.path.join(OUTPUT_DIR, "GSE174574_all_DEGs.csv"), index=False)
    cupro_de.to_csv(os.path.join(OUTPUT_DIR, "GSE174574_cuproptosis_DEGs.csv"), index=False)
    sig_cupro.to_csv(os.path.join(OUTPUT_DIR, "GSE174574_sig_cuproptosis_DEGs.csv"), index=False)
    print(f"  结果已保存至 {OUTPUT_DIR}")

def main():
    print("=" * 60)
    print("L1 表型锚定：GSE174574 scRNA-seq 24h MCAO vs Sham")
    print("=" * 60)
    extract_tar()
    adata = merge_samples()
    adata = quality_control(adata)
    adata = normalize_and_cluster(adata)
    adata = annotate_cell_types(adata)
    de_df, cupro_de, sig_cupro = differential_expression(adata)
    plot_results(adata, de_df, cupro_de)
    save_outputs(de_df, cupro_de, sig_cupro)
    print("\n" + "=" * 60)
    print(f"铜死亡差异基因: {len(sig_cupro)} (要求 >= 5)")
    print("L1 分析完成！")
    print("=" * 60)

if __name__ == "__main__":
    main()
