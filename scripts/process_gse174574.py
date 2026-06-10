# -*- coding: utf-8 -*-
"""
Process GSE174574: 10x scRNA-seq with biological replicates
============================================================
Optimized for speed: Cluster-based annotation instead of cell-level.
6 independent samples: 3 Sham (rep1/2/3) + 3 MCAO (rep1/2/3)
"""

import os
import sys
import tarfile
import gzip
import logging
import pandas as pd
import numpy as np
import scipy.sparse as sp
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    import scanpy as sc
except ImportError:
    print("Error: scanpy not installed")
    sys.exit(1)

# Paths
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAR_FILE = r"C:\Users\Jy-Mentor-7\Downloads\GSE174574_RAW (1).tar"
OUTPUT_DIR = os.path.join(BASE_DIR, "results", "stage7_gse174574")
os.makedirs(OUTPUT_DIR, exist_ok=True)

logger = logging.getLogger("gse174574")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler(sys.stdout)
handler.setFormatter(logging.Formatter("%(message)s"))
logger.addHandler(handler)

SAMPLES = [
    {"gsm": "GSM5319987", "sample_id": "sham_rep1", "fname": "sham1"},
    {"gsm": "GSM5319988", "sample_id": "sham_rep2", "fname": "sham2"},
    {"gsm": "GSM5319989", "sample_id": "sham_rep3", "fname": "sham3"},
    {"gsm": "GSM5319990", "sample_id": "MCAO_rep1", "fname": "MCAO1"},
    {"gsm": "GSM5319991", "sample_id": "MCAO_rep2", "fname": "MCAO2"},
    {"gsm": "GSM5319992", "sample_id": "MCAO_rep3", "fname": "MCAO3"},
]

def load_data():
    """Load and merge all 6 samples"""
    logger.info(f"Loading GSE174574 from {TAR_FILE}...")
    all_data = []
    all_barcodes = []
    all_genes = None

    with tarfile.open(TAR_FILE, 'r') as tf:
        for s in SAMPLES:
            logger.info(f"  Extracting {s['sample_id']}...")
            # Files
            bc_name = f"{s['gsm']}_{s['fname']}_barcodes.tsv.gz"
            gene_name = f"{s['gsm']}_{s['fname']}_genes.tsv.gz"
            mtx_name = f"{s['gsm']}_{s['fname']}_matrix.mtx.gz"
            
            # Extract barcodes
            for m in tf:
                if m.name == bc_name:
                    with tf.extractfile(m) as f:
                        barcodes = [line.decode('utf-8').strip() for line in gzip.GzipFile(fileobj=f)]
                    break
            # Prefix barcodes with sample ID to ensure uniqueness
            barcodes = [f"{s['sample_id']}_{bc}" for bc in barcodes]
            all_barcodes.extend(barcodes)

            # Extract genes (only need to do this once ideally, but check consistency)
            if all_genes is None:
                for m in tf:
                    if m.name == gene_name:
                        with tf.extractfile(m) as f:
                            df = pd.read_csv(gzip.GzipFile(fileobj=f), sep='\t', header=None)
                            # Gene symbol is usually col 1, fallback to col 0
                            all_genes = df.iloc[:, 1].tolist() if df.shape[1] > 1 else df.iloc[:, 0].tolist()
                        break
            
            # Extract matrix (genes x cells) -> transpose to (cells x genes)
            for m in tf:
                if m.name == mtx_name:
                    with tf.extractfile(m) as f:
                        from scipy.io import mmread
                        mat = mmread(gzip.GzipFile(fileobj=f)).T.tocsr()
                        all_data.append(mat)
                    break

    combined_matrix = sp.vstack(all_data, format='csr')
    logger.info(f"Loaded {combined_matrix.shape[0]} cells, {combined_matrix.shape[1]} genes")
    return combined_matrix, all_genes, all_barcodes

def process_and_save(matrix, genes, barcodes):
    """Process data, cluster, annotate, and save"""
    logger.info("Building AnnData...")
    
    # Create obs
    sample_ids = [bc.split("_")[0] + "_" + bc.split("_")[1] for bc in barcodes] # e.g. sham_rep1
    # Fix parsing: sample_id is first two parts of barcode
    sample_ids = []
    for bc in barcodes:
        parts = bc.split("_")
        # Format: sham_rep1_XXXX -> sham_rep1
        if len(parts) >= 3:
            sample_ids.append(f"{parts[0]}_{parts[1]}")
        else:
            sample_ids.append(parts[0])

    conditions = ["Sham" if "sham" in s else "MCAO" for s in sample_ids]

    # Make gene names unique
    gene_names = pd.Index([g.upper() for g in genes])
    gene_names = gene_names.str.replace("^MT-", "mt-", regex=True) # Standardize MT prefix
    # Handle duplicates by appending suffix
    seen = {}
    final_genes = []
    for g in gene_names:
        if g not in seen:
            seen[g] = 0
            final_genes.append(g)
        else:
            seen[g] += 1
            final_genes.append(f"{g}_{seen[g]}")
    gene_names = pd.Index(final_genes)

    adata = sc.AnnData(
        X=matrix,
        obs=pd.DataFrame({'sample_id': sample_ids, 'condition': conditions}, index=barcodes),
        var=pd.DataFrame(index=gene_names)
    )
    adata.var_names_make_unique()

    # QC
    adata.var['mt'] = adata.var_names.str.startswith('mt-') | adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
    logger.info(f"Raw cells: {adata.n_obs}")
    adata = adata[adata.obs.n_genes_by_counts > 200]
    adata = adata[adata.obs.n_genes_by_counts < 7500]
    adata = adata[adata.obs.pct_counts_mt < 20]
    logger.info(f"After QC: {adata.n_obs}")

    # Normalization & PCA
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True)
    sc.pp.pca(adata, n_comps=30)
    
    # Clustering
    logger.info("Clustering (Leiden)...")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    sc.tl.leiden(adata, resolution=0.5)
    logger.info(f"Found {adata.obs['leiden'].nunique()} clusters")

    # Cell Type Annotation based on Cluster Markers
    logger.info("Annotating cell types...")
    sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
    
    # Define markers for common types
    markers = {
        'Microglia': ['C1QA', 'C1QB', 'C1QC', 'CX3CR1', 'TYROBP', 'FCER1G', 'AIF1'],
        'Neuron': ['SNAP25', 'SYT1', 'RBFOX3', 'NEFL', 'NEFM', 'STMN2'],
        'Astrocyte': ['GFAP', 'AQP4', 'SLC1A2', 'ALDH1L1', 'FGFR3', 'GJA1'],
        'Oligodendrocyte': ['MBP', 'PLP1', 'MOG', 'MAG', 'CNP', 'OLIG1'],
        'OPC': ['PDGFRA', 'CSPG4', 'VCAN', 'OLIG2'],
        'Endothelial': ['PECAM1', 'CLDN5', 'VWF', 'CDH5', 'FLT1'],
        'Pericyte': ['RGS5', 'PDGFRB', 'ACTA2', 'MCAM'],
        'T_cells': ['CD3D', 'CD3E', 'CD3G', 'CD4', 'CD8A'],
    }

    # Score clusters based on markers
    cluster_types = {}
    # Get top genes per cluster
    top_genes = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
    
    # Simple heuristic: if top 50 genes of a cluster contain specific markers
    for cluster in adata.obs['leiden'].unique():
        cluster_genes = top_genes[str(cluster)].tolist()
        cluster_genes_upper = [g.upper() for g in cluster_genes]
        
        best_type = 'Unknown'
        max_score = 0
        
        for ctype, markers_list in markers.items():
            matches = sum(1 for m in markers_list if m in cluster_genes_upper)
            if matches > max_score:
                max_score = matches
                best_type = ctype
        
        # Refine: if no clear winner, check specific top markers
        if max_score == 0:
            # Check if Neuron (usually high UMI, high genes)
            if any(g in cluster_genes_upper for g in ['SNAP25', 'RBFOX3']):
                best_type = 'Neuron'
            elif any(g in cluster_genes_upper for g in ['C1QA', 'TYROBP']):
                best_type = 'Microglia'
            elif any(g in cluster_genes_upper for g in ['GFAP', 'AQP4']):
                best_type = 'Astrocyte'
        
        cluster_types[cluster] = best_type

    # Map clusters to cell types
    adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_types)
    logger.info(f"Cell types: {adata.obs['cell_type'].value_counts().to_dict()}")

    # Save Full AnnData
    adata.write(os.path.join(OUTPUT_DIR, "gse174574_adata.h5ad"))
    
    # Pseudo-bulk Aggregation (Sample x CellType)
    # FIX:[v7][use raw counts sum → CPM → log1p instead of mean of log1p]
    logger.info("Generating pseudo-bulk matrix (raw counts → CPM → log1p)...")
    pseudo_data = []
    pseudo_meta = []
    
    # Store raw counts before normalization
    raw_counts = matrix.copy()
    
    groups = adata.obs.groupby(['sample_id', 'cell_type']).groups
    for (sid, ct), indices in groups.items():
        if len(indices) < 3: # Skip if too few cells
            continue
        
        # Sum raw counts for these cells
        subset_raw = raw_counts[indices]
        sum_counts = np.asarray(subset_raw.sum(axis=0)).flatten().astype(np.float64)
        
        # CPM normalization (counts per million)
        total_counts = sum_counts.sum()
        if total_counts > 0:
            cpm = (sum_counts / total_counts) * 1e6
        else:
            cpm = np.zeros_like(sum_counts)
        
        # log1p transformation
        log_cpm = np.log1p(cpm)
        
        pseudo_data.append(log_cpm)
        pseudo_meta.append({
            'sample_id': sid,
            'cell_type': ct,
            'condition': 'Sham' if 'sham' in sid else 'MCAO',
            'n_cells': len(indices)
        })

    pseudo_adata = sc.AnnData(
        X=np.array(pseudo_data),
        obs=pd.DataFrame(pseudo_meta),
        var=pd.DataFrame(index=adata.var_names)
    )
    pseudo_adata.write(os.path.join(OUTPUT_DIR, "gse174574_pseudobulk.h5ad"))
    
    logger.info(f"Pseudo-bulk shape: {pseudo_adata.shape}")
    logger.info("Done! Saved to results/stage7_gse174574/")

if __name__ == "__main__":
    mat, genes, barcodes = load_data()
    process_and_save(mat, genes, barcodes)
