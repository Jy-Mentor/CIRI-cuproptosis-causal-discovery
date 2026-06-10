#!/usr/bin/env python3
"""
ARCHS4 Gene Co-expression Edge Query Script

Queries the ARCHS4 API (https://maayanlab.cloud/matrixapi/coltop) to retrieve
gene-gene co-expression correlation values for a given gene list, and outputs
edges with |correlation| > 0.7.

Input:  subgraph_genes.txt  (one gene symbol per line, ~15,000 genes)
Output: gene_coexp_edges.txt (gene1 gene2 correlation, tab-separated)

API: POST https://maayanlab.cloud/matrixapi/coltop
     Body: {"id": "GENE", "count": 200}
     Returns top-N correlated genes with Pearson correlation values.

Fallback: POST https://maayanlab.cloud/sigpy/data/correlation
          Body: {"gene": "GENE", "species": "human", "meta": ""}
"""

import argparse
import sys
import time
import random
from pathlib import Path

import requests


ARCHS4_MATRIX_API = "https://maayanlab.cloud/matrixapi/coltop"
ARCHS4_SIGPY_API = "https://maayanlab.cloud/sigpy/data/correlation"
REQUEST_DELAY_MIN = 0.3
REQUEST_DELAY_MAX = 0.5
MAX_RETRIES = 3
RETRY_BACKOFF = 1.5
TOP_N_CORRELATIONS = 200
CORRELATION_THRESHOLD = 0.7
REQUEST_TIMEOUT = 20


def query_matrix_api(gene, count=TOP_N_CORRELATIONS):
    payload = {"id": gene, "count": count}
    response = requests.post(ARCHS4_MATRIX_API, json=payload, timeout=REQUEST_TIMEOUT)
    response.raise_for_status()
    data = response.json()
    rowids = data.get("rowids", [])
    values = data.get("values", [])
    if len(rowids) == 0:
        return []
    results = []
    for i in range(len(rowids)):
        if rowids[i].upper() == gene.upper():
            continue
        results.append((rowids[i].upper(), values[i]))
    return results


def query_sigpy_api(gene):
    payload = {"gene": gene, "species": "human", "meta": ""}
    response = requests.post(ARCHS4_SIGPY_API, json=payload, timeout=REQUEST_TIMEOUT)
    response.raise_for_status()
    data = response.json()
    results = []
    for entry in data.get("positive_correlated_genes", []):
        results.append((entry["gene"].upper(), entry["correlation"]))
    for entry in data.get("negative_correlated_genes", []):
        results.append((entry["gene"].upper(), entry["correlation"]))
    return results


def query_gene_correlations(gene, use_fallback=False):
    last_exception = None
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            if not use_fallback:
                results = query_matrix_api(gene)
            else:
                results = query_sigpy_api(gene)
            return results, None
        except requests.exceptions.HTTPError as e:
            last_exception = e
            status = e.response.status_code if e.response is not None else None
            if status == 404:
                return [], None
            if attempt < MAX_RETRIES:
                wait = RETRY_BACKOFF ** attempt
                time.sleep(wait)
        except requests.exceptions.ConnectionError as e:
            last_exception = e
            if attempt < MAX_RETRIES:
                wait = RETRY_BACKOFF ** attempt
                time.sleep(wait)
        except requests.exceptions.Timeout as e:
            last_exception = e
            if attempt < MAX_RETRIES:
                wait = RETRY_BACKOFF ** attempt
                time.sleep(wait)
        except Exception as e:
            last_exception = e
            break

    if not use_fallback:
        return query_gene_correlations(gene, use_fallback=True)
    return [], str(last_exception)


def load_genes(filepath):
    genes = set()
    with open(filepath, "r") as f:
        first_line = True
        for line in f:
            if first_line:
                first_line = False
                stripped = line.strip().upper()
                if stripped in ("GENE", "GENES", "SYMBOL", "GENE_SYMBOL"):
                    continue
            gene = line.strip().upper()
            if gene:
                genes.add(gene)
    return sorted(genes)


def load_existing_pairs(output_path):
    pairs = set()
    last_gene = None
    if not output_path.exists():
        return pairs, last_gene
    with open(output_path, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                g1, g2 = parts[0].upper(), parts[1].upper()
                pairs.add((g1, g2))
                last_gene = g1
    return pairs, last_gene


def main():
    parser = argparse.ArgumentParser(
        description="Query ARCHS4 API for gene co-expression edges"
    )
    parser.add_argument(
        "-i", "--input",
        default="subgraph_genes.txt",
        help="Input file with gene symbols, one per line (default: subgraph_genes.txt)",
    )
    parser.add_argument(
        "-o", "--output",
        default="gene_coexp_edges.txt",
        help="Output file for co-expression edges (default: gene_coexp_edges.txt)",
    )
    parser.add_argument(
        "-t", "--threshold",
        type=float,
        default=CORRELATION_THRESHOLD,
        help=f"Correlation threshold (absolute value, default: {CORRELATION_THRESHOLD})",
    )
    parser.add_argument(
        "-n", "--top-n",
        type=int,
        default=TOP_N_CORRELATIONS,
        help=f"Number of top correlated genes to fetch per query (default: {TOP_N_CORRELATIONS})",
    )
    parser.add_argument(
        "--delay-min",
        type=float,
        default=REQUEST_DELAY_MIN,
        help=f"Minimum delay between requests in seconds (default: {REQUEST_DELAY_MIN})",
    )
    parser.add_argument(
        "--delay-max",
        type=float,
        default=REQUEST_DELAY_MAX,
        help=f"Maximum delay between requests in seconds (default: {REQUEST_DELAY_MAX})",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume from the last processed gene in the existing output file",
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)
    threshold = args.threshold
    delay_min = args.delay_min
    delay_max = args.delay_max

    if not input_path.exists():
        print(f"[ERROR] Input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    gene_list = load_genes(input_path)
    gene_set = set(gene_list)
    total_genes = len(gene_list)
    print(f"[INFO] Loaded {total_genes} unique genes from {input_path}")
    sys.stdout.flush()

    if total_genes == 0:
        print("[ERROR] No genes found in input file", file=sys.stderr)
        sys.exit(1)

    written_pairs = set()
    total_edges = 0
    errors = 0
    not_found = 0
    start_idx = 0

    if args.resume:
        written_pairs, last_gene = load_existing_pairs(output_path)
        total_edges = len(written_pairs)
        print(f"[RESUME] Loaded {total_edges} existing edges from {output_path}")
        if last_gene:
            for i, g in enumerate(gene_list):
                if g >= last_gene:
                    start_idx = i
                    break
            print(f"[RESUME] Skipping to gene index {start_idx} (gene: {gene_list[start_idx]})")
        sys.stdout.flush()

    if not args.resume:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as f_out:
            f_out.write(f"# gene1\tgene2\tcorrelation\n")
            f_out.write(f"# Source: ARCHS4 (https://maayanlab.cloud/archs4/)\n")
            f_out.write(f"# Threshold: |correlation| > {threshold}\n")
            f_out.write(f"# Input genes: {total_genes}\n")
            f_out.write(f"# Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")

    start_time = time.time()

    for idx in range(start_idx, total_genes):
        gene = gene_list[idx]
        if (idx + 1) % 100 == 0 or idx == start_idx:
            elapsed = time.time() - start_time
            rate = (idx + 1) / elapsed if elapsed > 0 else 0
            eta = (total_genes - idx - 1) / rate if rate > 0 else 0
            print(
                f"[PROGRESS] {idx + 1}/{total_genes} genes queried "
                f"({100 * (idx + 1) / total_genes:.1f}%) | "
                f"Edges: {total_edges} | Errors: {errors} | "
                f"NotFound: {not_found} | "
                f"ETA: {eta / 60:.1f} min"
            )
            sys.stdout.flush()

        t0 = time.time()
        results, error = query_gene_correlations(gene)
        t1 = time.time()
        if idx < start_idx + 5:
            print(f"  [DEBUG] Gene '{gene}' query took {t1 - t0:.1f}s, got {len(results)} results")

        if error:
            errors += 1
            if errors <= 5:
                print(f"[WARN] Gene '{gene}' failed after retries: {error}", file=sys.stderr)

        if not results:
            not_found += 1
        else:
            new_edges = []
            for target_gene, corr in results:
                if target_gene not in gene_set:
                    continue
                if abs(corr) <= threshold:
                    continue
                g1, g2 = (gene, target_gene) if gene < target_gene else (target_gene, gene)
                pair_key = (g1, g2)
                if pair_key in written_pairs:
                    continue
                written_pairs.add(pair_key)
                new_edges.append((g1, g2, corr))

            if new_edges:
                with open(output_path, "a") as f_out:
                    for g1, g2, corr in new_edges:
                        f_out.write(f"{g1}\t{g2}\t{corr:.6f}\n")
                total_edges += len(new_edges)

        if idx < total_genes - 1:
            sleep_time = random.uniform(delay_min, delay_max)
            time.sleep(sleep_time)

    elapsed = time.time() - start_time
    print(f"\n[DONE] Completed in {elapsed / 60:.1f} minutes")
    print(f"[RESULT] Total edges (|corr| > {threshold}): {total_edges}")
    print(f"[RESULT] Genes queried: {total_genes}")
    print(f"[RESULT] Errors: {errors}")
    print(f"[RESULT] Genes not found in ARCHS4: {not_found}")
    print(f"[RESULT] Output saved to: {output_path}")


if __name__ == "__main__":
    main()