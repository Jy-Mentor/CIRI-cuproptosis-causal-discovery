# -*- coding: utf-8 -*-
"""
TaRGET II 多源特征提取器
- RNA FPKM: 32 samples → log2(FPKM+1) per gene → PCA → 16 dim
- ATAC narrowPeak: 26 samples → promoter openness → PCA → 16 dim
- ENSMUSG → human gene symbol via mygene
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd
from pathlib import Path

warnings.filterwarnings("ignore")

# Paths
TOXI_DIR = Path(r"D:\反向网络药理学\GAT拓展维度\Toxi")
FPKM_DIR = TOXI_DIR / "rna_fpkm"
ATAC_DIR = TOXI_DIR / "atac"
OUTPUT_DIR = Path(__file__).parent
CACHE_DIR = OUTPUT_DIR / "toxirna_cache"
CACHE_DIR.mkdir(exist_ok=True)

# PCA target dimensions
FPKM_PCA_DIM = 16
ATAC_PCA_DIM = 16

print("=" * 60)
print("  TaRGET II Multi-Source Feature Extractor")
print("=" * 60)


# ============================================================================
# 1. FPKM Expression Features
# ============================================================================

def extract_fpkm_features():
    """Parse all FPKM files, compute mean log2(FPKM+1) per ENSMUSG gene."""
    cache_file = CACHE_DIR / "fpkm_mean.npz"
    if cache_file.exists():
        print("[FPKM] Loading from cache...")
        data = np.load(cache_file, allow_pickle=True)
        return data["genes"], data["fpkm_mean"]

    fpkm_files = sorted(FPKM_DIR.glob("*.tsv"))
    print(f"[FPKM] Found {len(fpkm_files)} FPKM files")

    all_data = {}
    gene_ids_all = set()

    for i, f in enumerate(fpkm_files):
        try:
            df = pd.read_csv(f, sep="\t")
            # Strip version from ENSMUSG ID (e.g., ENSMUSG00000000001.4 → ENSMUSG00000000001)
            gene_ids = df["gene_id"].str.split(".").str[0].values
            fpkm_vals = df["FPKM"].fillna(0).values.astype(np.float32)
            for gid, val in zip(gene_ids, fpkm_vals):
                if gid not in all_data:
                    all_data[gid] = []
                all_data[gid].append(val)
            gene_ids_all.update(gene_ids)
            if (i + 1) % 8 == 0:
                print(f"  Processed {i+1}/{len(fpkm_files)} files...")
        except Exception as e:
            print(f"  [WARN] Failed {f.name}: {e}")

    genes = sorted(gene_ids_all)
    fpkm_mean = np.zeros(len(genes), dtype=np.float32)
    fpkm_max = np.zeros(len(genes), dtype=np.float32)
    fpkm_cv = np.zeros(len(genes), dtype=np.float32)

    for i, g in enumerate(genes):
        vals = np.array(all_data.get(g, [0]), dtype=np.float32)
        # log2(FPKM+1) transform to stabilize variance
        vals_log = np.log2(vals + 1.0)
        fpkm_mean[i] = vals_log.mean()
        fpkm_max[i] = vals_log.max()
        fpkm_cv[i] = vals_log.std() / (vals_log.mean() + 1e-8)

    # Combine into feature matrix
    fpkm_mat = np.column_stack([fpkm_mean, fpkm_max, fpkm_cv])
    print(f"[FPKM] Extracted {len(genes)} genes × {fpkm_mat.shape[1]} features")

    np.savez_compressed(cache_file, genes=np.array(genes), fpkm_mean=fpkm_mat)
    return np.array(genes), fpkm_mat


# ============================================================================
# 2. ATAC Promoter Openness
# ============================================================================

def download_mouse_tss():
    """Download mouse TSS coordinates from Ensembl BioMart."""
    cache_file = CACHE_DIR / "mouse_tss.csv"
    if cache_file.exists():
        print("[ATAC] Loading TSS from cache...")
        return pd.read_csv(cache_file)

    print("[ATAC] Downloading mouse gene coordinates from Ensembl...")
    try:
        import urllib.request
        import gzip
        from io import BytesIO

        # Use Ensembl REST API to get mouse gene info
        # Alternative: download from BioMart
        url = ("https://rest.ensembl.org/info/data/?content-type=application/json")
        # Actually use the lookup endpoint for each gene
        # For efficiency, use the POST endpoint for multiple IDs
        
        # Strategy: Use Ensembl BioMart via requests
        import requests
        
        # Get all mouse genes with chromosomes and TSS
        biomart_url = "http://www.ensembl.org/biomart/martservice"
        
        # Query for mouse genes with coordinates
        query = """<?xml version="1.0" encoding="UTF-8"?>
        <!DOCTYPE Query>
        <Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1"
               count="" datasetConfigVersion="0.6">
          <Dataset name="mmusculus_gene_ensembl" interface="default">
            <Attribute name="ensembl_gene_id"/>
            <Attribute name="chromosome_name"/>
            <Attribute name="start_position"/>
            <Attribute name="end_position"/>
            <Attribute name="strand"/>
          </Dataset>
        </Query>"""
        
        resp = requests.post(biomart_url, data={"query": query}, timeout=120)
        lines = resp.text.strip().split("\n")
        data = []
        for ln in lines:
            parts = ln.strip().split("\t")
            if len(parts) >= 5:
                try:
                    chrom = parts[1]
                    start = int(parts[2])
                    end = int(parts[3])
                    strand = int(parts[4])
                    data.append([parts[0], chrom, start, end, strand])
                except (ValueError, IndexError):
                    continue

        df = pd.DataFrame(data, columns=["ensembl_gene_id", "chr", "start", "end", "strand"])
        df.to_csv(cache_file, index=False)
        print(f"[ATAC] Downloaded {len(df)} mouse gene coordinates")
        return df

    except Exception as e:
        print(f"[WARN] BioMart download failed: {e}")
        print("[ATAC] Using fallback: chromosome-level peak density")
        return None


def extract_atac_features(ensmusg_ids):
    """Compute promoter openness from narrowPeak files."""
    cache_file = CACHE_DIR / "atac_promoter.npz"
    if cache_file.exists():
        print("[ATAC] Loading from cache...")
        data = np.load(cache_file, allow_pickle=True)
        cached_ids = data["genes"]
        cached_mat = data["atac_mat"]
        # Re-index to match requested gene order
        id_map = {g: i for i, g in enumerate(cached_ids)}
        result = np.zeros((len(ensmusg_ids), cached_mat.shape[1]), dtype=np.float32)
        for i, g in enumerate(ensmusg_ids):
            if g in id_map:
                result[i] = cached_mat[id_map[g]]
        return result

    tss_df = download_mouse_tss()
    atac_files = sorted(ATAC_DIR.glob("*.narrowPeak"))
    print(f"[ATAC] Found {len(atac_files)} narrowPeak files")

    if tss_df is not None and len(tss_df) > 0:
        # Build gene coordinate lookup
        gene_coords = {}
        for _, row in tss_df.iterrows():
            gid = row["ensembl_gene_id"]
            if isinstance(row["chr"], str) and row["chr"].isdigit() or \
               (isinstance(row["chr"], str) and row["chr"] in ("X", "Y", "MT")):
                chrom = f"chr{row['chr']}"
            else:
                continue
            strand = row["strand"]
            if strand == 1:
                tss = row["start"]
            else:
                tss = row["end"]
            # Promoter: TSS ± 2000bp
            promoter_start = max(0, tss - 2000)
            promoter_end = tss + 2000
            gene_coords[gid] = (chrom, promoter_start, promoter_end)

        print(f"[ATAC] Built coords for {len(gene_coords)} genes")

        import bisect

        # Compute promoter signal per sample (optimized with bisect)
        n_samples = len(atac_files)
        atac_mat = np.zeros((len(ensmusg_ids), n_samples), dtype=np.float32)

        for s_idx, atac_f in enumerate(atac_files):
            print(f"  ATAC sample {s_idx+1}/{n_samples}...")

            # Parse narrowPeak and group by chromosome
            peaks_by_chrom = {}
            with open(atac_f, "r") as fh:
                for ln in fh:
                    parts = ln.strip().split("\t")
                    if len(parts) >= 7:
                        try:
                            chrom = parts[0]
                            start = int(parts[1])
                            end = int(parts[2])
                            signal = float(parts[6])
                            if chrom not in peaks_by_chrom:
                                peaks_by_chrom[chrom] = []
                            peaks_by_chrom[chrom].append((start, end, signal))
                        except (ValueError, IndexError):
                            continue

            # Sort peaks by start position within each chromosome + precompute start arrays
            chrom_peak_starts = {}
            for chrom in peaks_by_chrom:
                peaks_by_chrom[chrom].sort(key=lambda x: x[0])
                chrom_peak_starts[chrom] = [p[0] for p in peaks_by_chrom[chrom]]

            # For each gene, use bisect to find overlapping peaks
            for g_idx, g in enumerate(ensmusg_ids):
                if g not in gene_coords:
                    continue
                g_chrom, g_start, g_end = gene_coords[g]
                if g_chrom not in peaks_by_chrom:
                    continue

                chrom_peaks = peaks_by_chrom[g_chrom]
                peak_starts = chrom_peak_starts[g_chrom]

                total_signal = 0.0
                # Use bisect to find the first peak that starts near gene_start
                idx = max(0, bisect.bisect_left(peak_starts, g_start) - 1)

                # Scan forward only while peak_start < gene_end
                while idx < len(chrom_peaks):
                    p_start, p_end, p_signal = chrom_peaks[idx]
                    if p_start >= g_end:
                        break
                    if p_end > g_start:
                        total_signal += p_signal
                    idx += 1

                atac_mat[g_idx, s_idx] = np.log2(total_signal + 1.0)

        # Impute zeros with column median
        for j in range(n_samples):
            col = atac_mat[:, j]
            nonzero = col[col > 0]
            if len(nonzero) > 0:
                col[col == 0] = np.median(nonzero)

        print(f"[ATAC] Extracted {len(ensmusg_ids)} genes x {n_samples} samples "
              f"(non-zero ratio: {(atac_mat > 0).sum() / atac_mat.size:.3f})")

        np.savez_compressed(cache_file, genes=np.array(ensmusg_ids), atac_mat=atac_mat)
        return atac_mat

    else:
        # Fallback: chromosome-level peak density
        print("[ATAC] Using fallback: chromosome-level peak density")
        # Simple: count peaks per chromosome, normalize
        n_samples = len(atac_files)
        atac_mat = np.zeros((len(ensmusg_ids), n_samples), dtype=np.float32)
        print(f"[ATAC] Fallback: placeholder (no TSS data available)")
        return atac_mat


# ============================================================================
# 3. ENSMUSG → Human Gene Symbol Mapping
# ============================================================================

def map_mouse_to_human(ensmusg_ids):
    """Map mouse ENSMUSG IDs to human gene symbols using mygene."""
    cache_file = CACHE_DIR / "mouse2human.csv"
    if cache_file.exists():
        print("[Map] Loading mouse→human mapping from cache...")
        df = pd.read_csv(cache_file)
        return dict(zip(df["ensmusg"], df["human_symbol"]))

    print(f"[Map] Querying mygene for {len(ensmusg_ids)} mouse genes...")
    try:
        import mygene
        mg = mygene.MyGeneInfo()

        # Query in batches
        batch_size = 1000
        mapping = {}

        for i in range(0, len(ensmusg_ids), batch_size):
            batch = ensmusg_ids[i:i + batch_size]
            if (i // batch_size + 1) % 5 == 0:
                print(f"  Querying batch {i//batch_size + 1}/{(len(ensmusg_ids)-1)//batch_size + 1}...")

            try:
                results = mg.querymany(
                    batch,
                    scopes="ensembl.gene",
                    species="mouse",
                    fields="ensembl.gene,symbol,homologene",
                    returnall=True,
                )
                for r in results.get("out", []):
                    mus_id = r.get("query", "")
                    # Try to get human ortholog
                    human_symbol = None
                    if "ensembl" in r:
                        ens_data = r["ensembl"]
                        if isinstance(ens_data, list):
                            for entry in ens_data:
                                if isinstance(entry, dict) and "gene" in entry:
                                    human_symbol = entry.get("gene", human_symbol)
                    if human_symbol is None:
                        human_symbol = r.get("symbol", None)
                    if human_symbol:
                        # Strip version
                        mus_base = mus_id.split(".")[0] if "." in mus_id else mus_id
                        mapping[mus_base] = human_symbol
            except Exception as e:
                print(f"  [WARN] Batch query failed: {e}")
                continue

        print(f"[Map] Mapped {len(mapping)}/{len(ensmusg_ids)} genes to human")

        # Save cache
        pd.DataFrame([
            {"ensmusg": k, "human_symbol": v} for k, v in mapping.items()
        ]).to_csv(cache_file, index=False)
        return mapping

    except Exception as e:
        print(f"[ERROR] mygene mapping failed: {e}")
        print("[Map] Using direct symbol matching (strip ENSMUSG prefix)")
        # Fallback: use the ENSMUSG ID as is, remove version
        return {g.split(".")[0]: g.split(".")[0] for g in ensmusg_ids}


# ============================================================================
# 4. PCA Dimensionality Reduction
# ============================================================================

def apply_pca(feature_mat, n_components, feature_name, gene_ids):
    """Apply PCA to reduce feature dimensions, save scaler for reproducibility."""
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    cache_dir = CACHE_DIR / "pca"
    cache_dir.mkdir(exist_ok=True)
    pca_cache = cache_dir / f"{feature_name}_pca.npz"

    if pca_cache.exists():
        print(f"[PCA:{feature_name}] Loading from cache...")
        data = np.load(pca_cache)
        return data["reduced"]

    # Handle NaN/Inf
    feature_mat = np.nan_to_num(feature_mat, nan=0.0, posinf=0.0, neginf=0.0)

    # Standardize
    scaler = StandardScaler()
    scaled = scaler.fit_transform(feature_mat)

    # PCA
    actual_dim = min(n_components, min(scaled.shape[0], scaled.shape[1]) - 1)
    if actual_dim < 2:
        print(f"[PCA:{feature_name}] Too few features ({actual_dim}), skipping PCA")
        # Pad or truncate
        result = np.zeros((scaled.shape[0], n_components), dtype=np.float32)
        result[:, :scaled.shape[1]] = scaled[:, :n_components]
        return result

    pca = PCA(n_components=actual_dim, random_state=42)
    reduced = pca.fit_transform(scaled)
    explained = pca.explained_variance_ratio_.sum()
    print(f"[PCA:{feature_name}] {scaled.shape[1]} → {actual_dim} dims "
          f"(explained variance: {explained:.3f})")

    # Pad to target n_components if needed
    if actual_dim < n_components:
        padded = np.zeros((reduced.shape[0], n_components), dtype=np.float32)
        padded[:, :actual_dim] = reduced
        reduced = padded

    np.savez_compressed(pca_cache, reduced=reduced.astype(np.float32))
    return reduced.astype(np.float32)


# ============================================================================
# 5. Main Assembly
# ============================================================================

def main():
    # --- Step 1: Extract FPKM ---
    ensmusg_ids, fpkm_mat = extract_fpkm_features()

    # --- Step 2: Extract ATAC ---
    atac_mat = extract_atac_features(ensmusg_ids)

    # --- Step 3: Map ENSMUSG → Human ---
    mouse2human = map_mouse_to_human(ensmusg_ids)

    # --- Step 4: Filter to genes with human mapping ---
    human_symbols = []
    valid_indices = []
    for i, g in enumerate(ensmusg_ids):
        g_base = g.split(".")[0] if "." in g else g
        h = mouse2human.get(g_base, None)
        if h is not None:
            human_symbols.append(h)
            valid_indices.append(i)

    print(f"\n[Assembly] {len(human_symbols)}/{len(ensmusg_ids)} genes mapped to human")

    fpkm_mat = fpkm_mat[valid_indices]
    if atac_mat.shape[0] == len(ensmusg_ids):
        atac_mat = atac_mat[valid_indices]
    else:
        atac_mat = np.zeros((len(human_symbols), max(1, atac_mat.shape[1])), dtype=np.float32)

    # --- Step 5: Handle duplicates (keep first occurrence) ---
    seen = {}
    unique_idx = []
    for i, sym in enumerate(human_symbols):
        if sym not in seen:
            seen[sym] = i
            unique_idx.append(i)

    if len(unique_idx) < len(human_symbols):
        print(f"[Assembly] Removed {len(human_symbols) - len(unique_idx)} duplicate mappings")
        human_symbols = [human_symbols[i] for i in unique_idx]
        fpkm_mat = fpkm_mat[unique_idx]
        atac_mat = atac_mat[unique_idx]

    # --- Step 6: PCA ---
    fpkm_pca = apply_pca(fpkm_mat, FPKM_PCA_DIM, "fpkm", human_symbols)
    atac_pca = apply_pca(atac_mat, ATAC_PCA_DIM, "atac", human_symbols)

    # --- Step 7: Combine and save ---
    combined = np.concatenate([fpkm_pca, atac_pca], axis=1)
    print(f"\n[Assembly] Final features: {combined.shape[0]} genes × {combined.shape[1]} dims")

    # Save
    out_csv = OUTPUT_DIR / "toxirna_enhanced_features.csv"
    df = pd.DataFrame(
        combined,
        index=human_symbols,
        columns=[f"toxirna_{i:03d}" for i in range(combined.shape[1])]
    )
    df.index.name = "gene_symbol"
    df.to_csv(out_csv)
    print(f"[Output] Saved to {out_csv}")

    out_npy = OUTPUT_DIR / "toxirna_enhanced_features.npy"
    np.save(out_npy, combined.astype(np.float32))
    print(f"[Output] Saved to {out_npy}")

    # Save gene list for reference
    gene_list_out = OUTPUT_DIR / "toxirna_gene_list.txt"
    with open(gene_list_out, "w") as f:
        f.write("gene_symbol\n")
        f.write("\n".join(human_symbols))
    print(f"[Output] Gene list saved to {gene_list_out}")

    print("\n" + "=" * 60)
    print("  Feature extraction complete!")
    print(f"  Genes: {len(human_symbols)}")
    print(f"  Dimensions: FPKM={fpkm_pca.shape[1]} + ATAC={atac_pca.shape[1]} "
          f"= {combined.shape[1]}")
    print("=" * 60)

    return combined, human_symbols


if __name__ == "__main__":
    main()