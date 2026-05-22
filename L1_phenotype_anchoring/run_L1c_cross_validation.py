# L1c: scRNA-seq × Bulk 交叉验证
# 检查 Bulk 显著基因在 scRNA-seq 各细胞类型中的特异性表达
import scanpy as sc
import anndata
import numpy as np
import pandas as pd
import os
import tarfile
import gzip
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

# Bulk 显著基因（adj.P < 0.05，大写转小鼠首字母大写）
BULK_SIG_GENES = {
    'Cdkn2a': '上调', 'Cox17': '下调', 'Dbt': '上调', 'Dld': '下调',
    'Gls': '下调', 'Nfe2l2': '上调', 'Pdha1': '下调', 'Pdhb': '下调',
    'Slc31a1': '上调', 'Slc31a2': '上调', 'Atox1': '上调', 'Mt2': '上调', 'Cp': '上调'
}

# Bulk 不显著但检出的基因
BULK_NSIG_GENES = {
    'Atp7a': '下调', 'Atp7b': '下调', 'Lias': '下调', 'Fdx1': '上调',
    'Mtf1': '下调', 'Lipt2': '上调', 'Dlst': '下调', 'Gcsh': '下调',
    'Cox11': '下调', 'Sco1': '上调', 'Alb': '下调', 'Sod1': '下调',
    'Sod3': '下调', 'Commd1': '上调', 'Slc11a2': '下调', 'Steap3': '上调'
}

ALL_CHECK_GENES = list(BULK_SIG_GENES.keys()) + list(BULK_NSIG_GENES.keys())

SAMPLES = [
    {"id": "sham1", "condition": "Sham", "files": ["GSM5319987_sham1"]},
    {"id": "sham2", "condition": "Sham", "files": ["GSM5319988_sham2"]},
    {"id": "sham3", "condition": "Sham", "files": ["GSM5319989_sham3"]},
    {"id": "mcao1", "condition": "MCAO", "files": ["GSM5319990_MCAO1"]},
    {"id": "mcao2", "condition": "MCAO", "files": ["GSM5319991_MCAO2"]},
    {"id": "mcao3", "condition": "MCAO", "files": ["GSM5319992_MCAO3"]},
]

MARKERS = {
    "Microglia": ["Ptprc", "Aif1", "Cx3cr1", "Tmem119", "P2ry12", "C1qa", "C1qb"],
    "Neuron": ["Snap25", "Syt1", "Nefl", "Rbfox3", "Syn1"],
    "Astrocyte": ["Gfap", "Aqp4", "Slc1a3", "Aldh1l1"],
    "Endothelial": ["Pecam1", "Vwf", "Cldn5", "Cdh5"],
    "Oligodendrocyte": ["Mbp", "Plp1", "Mog", "Mag"],
    "OPC": ["Pdgfra", "Vcan", "Cspg4", "Olig1", "Olig2"],
}

def extract_tar():
    print("=== 解压 GSE174574_RAW.tar ===")
    os.makedirs(EXTRACT_DIR, exist_ok=True)
    with tarfile.open(RAW_TAR, "r") as tar:
        tar.extractall(EXTRACT_DIR)

def load_10x_sample(prefix):
    matrix_path = os.path.join(EXTRACT_DIR, f"{prefix}_matrix.mtx.gz")
    genes_path = os.path.join(EXTRACT_DIR, f"{prefix}_genes.tsv.gz")
    barcodes_path = os.path.join(EXTRACT_DIR, f"{prefix}_barcodes.tsv.gz")
    for p in [matrix_path, genes_path, barcodes_path]:
        if not os.path.exists(p):
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
    print("=== 合并样本 ===")
    adatas = []
    for sample in SAMPLES:
        prefix = sample["files"][0]
        adata = load_10x_sample(prefix)
        if adata is not None:
            adata.obs["sample_id"] = sample["id"]
            adata.obs["condition"] = sample["condition"]
            adatas.append(adata)
    if not adatas:
        raise ValueError("未加载任何样本")
    for adata in adatas:
        adata.obs_names_make_unique()
        adata.var_names_make_unique()
    merged = anndata.concat(adatas, axis=0, join="outer", label="batch", keys=[s["id"] for s in SAMPLES])
    merged.obs_names_make_unique()
    print(f"  合并后: {merged.n_obs} cells x {merged.n_vars} genes")
    return merged

def quality_control(adata):
    print("=== 质量控制 ===")
    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata = adata[adata.obs.n_genes_by_counts < 5000, :]
    adata = adata[adata.obs.pct_counts_mt < 20, :]
    adata = adata[adata.obs.total_counts > 500, :]
    print(f"  质控后: {adata.n_obs} cells, {adata.n_vars} genes")
    return adata

def normalize_and_annotate(adata):
    print("=== 标准化与细胞注释 ===")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    for ctype, genes in MARKERS.items():
        present = [g for g in genes if g in adata.var_names]
        if present:
            sc.tl.score_genes(adata, gene_list=present, score_name=f"score_{ctype}", use_raw=False)
    adata.obs["cell_type"] = "Unknown"
    score_cols = [f"score_{ct}" for ct in MARKERS if f"score_{ct}" in adata.obs.columns]
    if score_cols:
        best_scores = adata.obs[score_cols].idxmax(axis=1)
        adata.obs["cell_type"] = best_scores.str.replace("score_", "")
    print("细胞类型分布:")
    print(adata.obs["cell_type"].value_counts())
    return adata

def cell_type_de(adata):
    print("\n=== 各细胞类型差异分析 ===")
    results = []
    present_genes = [g for g in ALL_CHECK_GENES if g in adata.var_names]
    missing_genes = [g for g in ALL_CHECK_GENES if g not in adata.var_names]
    if missing_genes:
        print(f"  scRNA-seq 中缺失基因: {missing_genes}")

    for ctype in adata.obs["cell_type"].unique():
        if ctype == "Unknown":
            continue
        subset = adata[adata.obs["cell_type"] == ctype]
        if subset.n_obs < 100:
            print(f"  {ctype}: 细胞数 < 100，跳过")
            continue
        sc.tl.rank_genes_groups(subset, groupby="condition", reference="Sham", method="wilcoxon", use_raw=True)
        de_df = sc.get.rank_genes_groups_df(subset, group="MCAO")
        for gene in present_genes:
            row = de_df[de_df["names"] == gene]
            if len(row) > 0:
                row = row.iloc[0]
                log2fc = row["logfoldchanges"]
                pval = row["pvals_adj"]
                sig = pval < 0.05
                direction = "上调" if log2fc > 0 else "下调"
                bulk_dir = BULK_SIG_GENES.get(gene, BULK_NSIG_GENES.get(gene, "?"))
                consistent = (direction == bulk_dir) if sig else "N/A"
                results.append({
                    "gene": gene,
                    "cell_type": ctype,
                    "n_cells": subset.n_obs,
                    "log2FC": log2fc,
                    "adj.P": pval,
                    "sc_sig": sig,
                    "sc_direction": direction,
                    "bulk_direction": bulk_dir,
                    "bulk_sig": gene in BULK_SIG_GENES,
                    "direction_consistent": consistent
                })

    return pd.DataFrame(results)

def cross_validation_summary(de_results):
    print("\n=== 交叉验证总结 ===")
    bulk_sig = list(BULK_SIG_GENES.keys())
    bulk_nsig = list(BULK_NSIG_GENES.keys())

    sig_results = de_results[de_results["gene"].isin(bulk_sig)]
    nsig_results = de_results[de_results["gene"].isin(bulk_nsig)]

    print(f"\nBulk 显著基因 ({len(bulk_sig)} 个) 在 scRNA-seq 中的细胞特异性:")
    for gene in bulk_sig:
        gene_res = sig_results[sig_results["gene"] == gene]
        sig_ctypes = gene_res[gene_res["sc_sig"]]["cell_type"].tolist()
        if sig_ctypes:
            print(f"  {gene:12s}: 在 {', '.join(sig_ctypes)} 中显著 (Bulk={BULK_SIG_GENES[gene]})")
        else:
            print(f"  {gene:12s}: 在所有细胞类型中均不显著")

    print(f"\nBulk 不显著基因 ({len(bulk_nsig)} 个) 在 scRNA-seq 中的细胞特异性:")
    for gene in bulk_nsig:
        gene_res = nsig_results[nsig_results["gene"] == gene]
        sig_ctypes = gene_res[gene_res["sc_sig"]]["cell_type"].tolist()
        if sig_ctypes:
            print(f"  {gene:12s}: 在 {', '.join(sig_ctypes)} 中显著（细胞异质性掩盖）")

    print("\n方向一致性（Bulk 显著基因）:")
    consistent = sig_results[sig_results["sc_sig"] & (sig_results["direction_consistent"] == True)]
    inconsistent = sig_results[sig_results["sc_sig"] & (sig_results["direction_consistent"] == False)]
    print(f"  方向一致: {len(consistent)}/{len(sig_results[sig_results['sc_sig']])}")
    if len(inconsistent) > 0:
        for _, row in inconsistent.iterrows():
            print(f"  [X] {row['gene']} ({row['cell_type']}): scRNA={row['sc_direction']}, Bulk={row['bulk_direction']}")

    return de_results

def save_outputs(de_results):
    print("\n=== 保存结果 ===")
    de_results.to_csv(os.path.join(OUTPUT_DIR, "cross_validation_scRNA_bulk.csv"), index=False)

    summary = de_results.groupby("gene").agg({
        "sc_sig": "sum",
        "log2FC": "mean"
    }).reset_index()
    summary["bulk_sig"] = summary["gene"].isin(BULK_SIG_GENES.keys()).astype(int)
    summary.to_csv(os.path.join(OUTPUT_DIR, "cross_validation_summary.csv"), index=False)
    print(f"  结果已保存至 {OUTPUT_DIR}")

def main():
    print("=" * 60)
    print("L1c: scRNA-seq × Bulk 交叉验证")
    print("=" * 60)
    extract_tar()
    adata = merge_samples()
    adata = quality_control(adata)
    adata = normalize_and_annotate(adata)
    de_results = cell_type_de(adata)
    de_results = cross_validation_summary(de_results)
    save_outputs(de_results)
    print("\n交叉验证完成！")

if __name__ == "__main__":
    main()
