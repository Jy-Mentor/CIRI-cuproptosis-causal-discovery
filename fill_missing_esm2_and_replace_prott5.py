# -*- coding: utf-8 -*-
"""
ESM-2 补全 + PCA 替换 ProtT5 (v2 - 统一维度策略)
  背景: 现有 subgraph_embeddings 为 1024-dim ProtT5，与 ESM-2 t6 (320-dim) 维度不兼容
  策略: 为 feature table 中所有基因重新计算 ESM-2 t6，获得统一 320-dim 矩阵
        已有 UniProt 缓存复用，避免重复 API 调用
  1. 读取 feature table，识别所有需要 ESM-2 的基因
  2. 并发获取/复用缓存的蛋白质序列
  3. ESM-2 t6 批量推理 (320-dim) → PCA 50-dim
  4. PCA → 50-dim 替换旧 ProtT5 列
  5. 输出完整 ESM-2 文件和更新后的特征表
"""
import os, sys, time, json, logging, warnings
from concurrent.futures import ThreadPoolExecutor, as_completed

os.environ.setdefault("HF_ENDPOINT", "https://hf-mirror.com")
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

import numpy as np
import pandas as pd
import requests
from sklearn.decomposition import PCA
from tqdm import tqdm

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s", datefmt="%H:%M:%S")
logger = logging.getLogger(__name__)

# ============ CONFIG ============
ESM2_MODEL = "facebook/esm2_t6_8M_UR50D"
ESM2_DIM = 320
BATCH_SIZE = 500
MAX_SEQ_LEN = 1022
REQUEST_TIMEOUT = 30
REQUEST_DELAY = 0.15
N_WORKERS = 20
N_PCA = 50

RF_DIR = r"C:\Users\Jy-Mentor-7\Desktop\随机森林"
PRJ_DIR = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙"
CACHE_FILE = os.path.join(PRJ_DIR, "uniprot_cache.json")
FT_PATH = os.path.join(RF_DIR, "gene_features_table_with_gat_emb.csv")

NEW_ESM_OUT = os.path.join(RF_DIR, "subgraph_embeddings_esm2_320d.csv")
FT_OUT = os.path.join(RF_DIR, "gene_features_table_with_esm2.csv")

# ============ CACHE ============
def load_cache():
    if os.path.exists(CACHE_FILE):
        with open(CACHE_FILE, "r", encoding="utf-8") as f:
            return json.load(f)
    return {}

def save_cache(cache):
    with open(CACHE_FILE, "w", encoding="utf-8") as f:
        json.dump(cache, f, ensure_ascii=False, indent=1)

# ============ API ============
def map_symbol_to_uniprot(gene):
    url = "https://mygene.info/v3/query"
    params = {"q": gene, "scopes": "symbol", "fields": "uniprot", "species": "human", "size": 3}
    try:
        resp = requests.get(url, params=params, timeout=REQUEST_TIMEOUT)
        resp.raise_for_status()
        hits = resp.json().get("hits", [])
        for hit in hits:
            ui = hit.get("uniprot", {})
            if isinstance(ui, dict):
                sp = ui.get("Swiss-Prot")
                if sp:
                    return sp[0] if isinstance(sp, list) else sp
        for hit in hits:
            ui = hit.get("uniprot", {})
            if isinstance(ui, dict):
                for k in ("TrEMBL",):
                    v = ui.get(k)
                    if v:
                        return v[0] if isinstance(v, list) else v
    except Exception:
        pass
    return None

def fetch_sequence(uid):
    try:
        resp = requests.get(f"https://rest.uniprot.org/uniprotkb/{uid}.fasta", timeout=REQUEST_TIMEOUT)
        resp.raise_for_status()
        lines = resp.text.strip().splitlines()
        return "".join(l.strip() for l in lines if not l.startswith(">") and l.strip()) or None
    except Exception:
        return None

def resolve_gene(gene):
    uniprot_id = map_symbol_to_uniprot(gene)
    seq = None
    if uniprot_id:
        seq = fetch_sequence(uniprot_id)
    time.sleep(REQUEST_DELAY)
    return gene, seq

# ============ AAC FALLBACK ============
AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"
def compute_aac(seq):
    if not seq:
        return np.zeros(len(AA_ORDER), dtype=np.float32)
    counts = {aa: 0 for aa in AA_ORDER}
    for aa in seq:
        if aa in counts:
            counts[aa] += 1
    total = sum(counts.values())
    if total == 0:
        return np.zeros(len(AA_ORDER), dtype=np.float32)
    return np.array([counts[aa] / total for aa in AA_ORDER], dtype=np.float32)

# ============ ESM-2 ============
def load_esm2():
    from transformers import AutoModel, AutoTokenizer
    logger.info("Loading %s ...", ESM2_MODEL)
    tok = AutoTokenizer.from_pretrained(ESM2_MODEL)
    model = AutoModel.from_pretrained(ESM2_MODEL)
    model.eval()
    return model, tok

def extract_batch(seqs, genes, model, tok):
    import torch
    truncated = [s[:MAX_SEQ_LEN] for s in seqs]
    inputs = tok(truncated, return_tensors="pt", padding=True, truncation=True,
                 max_length=MAX_SEQ_LEN + 2, add_special_tokens=True)
    with torch.no_grad():
        outputs = model(**inputs)
        last = outputs.last_hidden_state
    results = {}
    for i, gene in enumerate(genes):
        slen = (inputs["attention_mask"][i] == 1).sum()
        emb = last[i, :slen].mean(dim=0).cpu().numpy().astype(np.float32)
        results[gene] = emb
    return results

# ============ MAIN ============
def main():
    logger.info("=" * 60)
    logger.info("ESM-2 补全 + PCA 替换 ProtT5 (v2 统一维度)")
    logger.info("=" * 60)

    # ---- 1. Load feature table ----
    logger.info("[1/5] Loading feature table...")
    ft = pd.read_csv(FT_PATH)
    ft["gene_symbol"] = ft["gene_symbol"].str.upper()
    all_ft_genes = sorted(ft["gene_symbol"].unique())
    n_total = len(all_ft_genes)
    dt_count = int(ft["is_drug_target"].sum())
    dg_count = int(ft["is_disease_gene"].sum())
    logger.info("  %d genes (DT+ %d, DG+ %d)", n_total, dt_count, dg_count)

    # ---- 2. Get sequences for ALL feature table genes ----
    logger.info("[2/5] Resolving sequences for %d genes...", n_total)
    cache = load_cache()
    logger.info("  Cache: %d entries", len(cache))

    gene_to_seq = {}
    unresolved = []
    for g in all_ft_genes:
        entry = cache.get(g.upper())
        if entry and entry.get("sequence"):
            gene_to_seq[g] = entry["sequence"]
        else:
            unresolved.append(g)

    logger.info("  Cache hit: %d, need API: %d", len(gene_to_seq), len(unresolved))

    if unresolved:
        logger.info("  Fetching sequences via mygene.info + UniProt (workers=%d)...", N_WORKERS)
        with ThreadPoolExecutor(max_workers=N_WORKERS) as ex:
            futures = {ex.submit(resolve_gene, g): g for g in unresolved}
            n_done = 0
            for f in tqdm(as_completed(futures), total=len(unresolved), desc="UniProt API"):
                gene, seq = f.result()
                if seq:
                    gene_to_seq[gene] = seq
                    cache[gene.upper()] = {"sequence": seq}
                else:
                    cache[gene.upper()] = {"sequence": None}
                n_done += 1
                if n_done % 500 == 0:
                    save_cache(cache)
        save_cache(cache)

    n_with_seq = len(gene_to_seq)
    n_no_seq = n_total - n_with_seq
    logger.info("  Sequences obtained: %d/%d (%.1f%%), no-seq: %d", n_with_seq, n_total,
                100 * n_with_seq / n_total, n_no_seq)

    # ---- 3. ESM-2 inference for ALL genes ----
    logger.info("[3/5] ESM-2 inference (t6 8M, 320-dim)...")
    model, tok = load_esm2()

    all_embeddings = {}
    seq_items = [(g, gene_to_seq[g]) for g in all_ft_genes if g in gene_to_seq]
    batches = [seq_items[i:i + BATCH_SIZE] for i in range(0, len(seq_items), BATCH_SIZE)]
    logger.info("  %d batches (batch_size=%d), %d genes with sequence",
                len(batches), BATCH_SIZE, len(seq_items))

    for bi, batch in enumerate(batches):
        bg = [g for g, s in batch]
        bs = [s for g, s in batch]
        logger.info("  Batch %d/%d: %d genes", bi + 1, len(batches), len(bg))
        try:
            feats = extract_batch(bs, bg, model, tok)
            all_embeddings.update(feats)
        except Exception as e:
            logger.error("  Batch failed: %s — fallback one-by-one", e)
            for g, s in batch:
                try:
                    feats = extract_batch([s], [g], model, tok)
                    all_embeddings.update(feats)
                except Exception as e2:
                    logger.error("    %s failed: %s — using AAC", g, e2)
                    all_embeddings[g] = compute_aac(s)

    logger.info("  ESM-2 computed: %d/%d genes", len(all_embeddings), len(seq_items))

    # AAC fallback for no-sequence genes
    n_aac = 0
    for g in all_ft_genes:
        if g not in all_embeddings:
            seq = gene_to_seq.get(g)
            all_embeddings[g] = compute_aac(seq)
            n_aac += 1
    logger.info("  AAC fallback: %d genes (no sequence)", n_aac)

    # ---- 4. Build uniform matrix + PCA ----
    logger.info("[4/5] Building uniform embedding matrix + PCA → %d-dim...", N_PCA)

    all_genes = sorted(all_ft_genes)
    emb_matrix = np.array([all_embeddings[g] for g in all_genes], dtype=np.float32)
    logger.info("  Embedding matrix: %s", emb_matrix.shape)

    pca = PCA(n_components=N_PCA, random_state=42)
    emb_pca = pca.fit_transform(emb_matrix)
    var_sum = pca.explained_variance_ratio_.sum()
    logger.info("  PCA explained variance: %.4f (%d-dim)", var_sum, N_PCA)

    # Save full ESM-2 embeddings
    feat_cols = [f"feat_{i+1:04d}" for i in range(emb_matrix.shape[1])]
    df_full_esm = pd.DataFrame(emb_matrix, columns=feat_cols)
    df_full_esm.insert(0, "gene_symbol", all_genes)
    df_full_esm.to_csv(NEW_ESM_OUT, index=False, encoding="utf-8")
    logger.info("  Full ESM-2: %s (%d genes x %d-dim)", NEW_ESM_OUT, len(df_full_esm), emb_matrix.shape[1])

    # ---- 5. Replace ProtT5 in feature table ----
    logger.info("[5/5] Replacing ProtT5 PCA columns in feature table...")

    pca_cols_new = [f"pca_{i}" for i in range(N_PCA)]
    df_pca = pd.DataFrame(emb_pca, columns=pca_cols_new)
    df_pca.insert(0, "gene_symbol", all_genes)

    old_pca = [c for c in ft.columns if c.startswith("pca_")]
    ft_clean = ft.drop(columns=old_pca)
    ft_new = ft_clean.merge(df_pca, on="gene_symbol", how="left")
    for c in pca_cols_new:
        ft_new[c] = ft_new[c].fillna(0.0)

    # Verify coverage
    zero_mask = ft_new[pca_cols_new].abs().sum(axis=1) < 1e-8
    n_zero = int(zero_mask.sum())
    dt_zero = int((ft_new["is_drug_target"] == 1) & zero_mask).sum()
    dg_zero = int((ft_new["is_disease_gene"] == 1) & zero_mask).sum()
    logger.info("  Zero-PCA genes: %d/%d (%.1f%%), DT+ zero: %d, DG+ zero: %d",
                n_zero, len(ft_new), 100 * n_zero / len(ft_new), dt_zero, dg_zero)

    ft_new.to_csv(FT_OUT, index=False, encoding="utf-8-sig")
    logger.info("  Feature table saved: %s (%d x %d)", FT_OUT, ft_new.shape[0], ft_new.shape[1])

    # Summary
    print("\n" + "=" * 60)
    print("  SUMMARY")
    print("=" * 60)
    print(f"  Feature table genes: {n_total} (DT+ {dt_count}, DG+ {dg_count})")
    print(f"  ESM-2 computed via GPU: {len(all_embeddings) - n_aac} genes")
    print(f"  AAC fallback: {n_aac} genes")
    print(f"  PCA explained variance: {var_sum:.4f} ({N_PCA}-dim)")
    print(f"  Zero-PCA genes: {n_zero} → {100*n_zero/n_total:.1f}% of feature table")
    print(f"  Output ESM-2: {NEW_ESM_OUT}")
    print(f"  Output feature table: {FT_OUT}")
    print("  Next: re-run GAT v7 → RF v13")
    logger.info("Done!")


if __name__ == "__main__":
    main()