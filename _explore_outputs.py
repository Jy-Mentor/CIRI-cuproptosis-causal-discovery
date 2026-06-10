import pandas as pd
import os

base = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙"

# 1. CIRI candidates
ciri = pd.read_csv(os.path.join(base, "bridge_pathway_cirI_candidates.csv"))
print("=== bridge_pathway_cirI_candidates.csv ===")
print(f"Shape: {ciri.shape}")
print(f"Columns: {list(ciri.columns)}")
print(f"Unique genes: {ciri['gene_symbol'].nunique()}")
print(f"Unique pathways: {ciri['pathway_name'].nunique()}")
print(f"Score range: [{ciri['score'].min():.4f}, {ciri['score'].max():.4f}]")
print("\nFirst 10 rows:")
print(ciri.head(10).to_string())
print()

# 2. Ensemble scores
ensemble = pd.read_csv(os.path.join(base, "bridge_pathway_scores_ensemble.csv"))
print("=== bridge_pathway_scores_ensemble.csv ===")
print(f"Shape: {ensemble.shape}")
print(f"Columns: {list(ensemble.columns)}")
print(f"Unique genes: {ensemble['gene_symbol'].nunique()}")
print(f"Unique pathways: {ensemble['pathway_name'].nunique()}")
print(f"Score range: [{ensemble['score'].min():.4f}, {ensemble['score'].max():.4f}]")
print(f"Eval modes: {ensemble['eval_mode'].value_counts().to_dict()}")
print()

# 3. Bridge scores
bridge = pd.read_csv(os.path.join(base, "bridge_pathway_scores.csv"))
print("=== bridge_pathway_scores.csv ===")
print(f"Shape: {bridge.shape}")
print(f"Columns: {list(bridge.columns)}")
print(f"Unique genes: {bridge['gene_symbol'].nunique()}")
print(f"Unique pathways: {bridge['pathway_name'].nunique()}")
print(f"Score range: [{bridge['score'].min():.4f}, {bridge['score'].max():.4f}]")
print()

# Cuproptosis gene lists
CUPROPTOSIS_EXECUTOR = {"FDX1", "LIAS", "LIPT1", "DLAT", "PDHB", "PDHX", "SLC31A1"}
CUPROPTOSIS_REGULATOR = {"ATP7A", "ATP7B", "ATOX1", "NFE2L2", "HIF1A", "MTOR", "NFKB1", "GPX4"}
CUPROPTOSIS_GENES = CUPROPTOSIS_EXECUTOR | CUPROPTOSIS_REGULATOR
EXTRA_COPPER_GENES = {"LIPT2", "MTF1", "SLC31A2", "DLD", "PDHA1", "LIAS",
                      "GCSH", "NFS1", "ISCU", "COX11", "SCO1", "COA6",
                      "COX10", "SURF1", "LYRM4", "ACSF2", "LIPTC1"}

ALL_CU_GENES = CUPROPTOSIS_GENES | EXTRA_COPPER_GENES

print("=" * 60)
print("铜死亡基因在三个输出文件中的存在情况")
print("=" * 60)

for fname, df in [("bridge_pathway_cirI_candidates.csv", ciri),
                   ("bridge_pathway_scores_ensemble.csv", ensemble),
                   ("bridge_pathway_scores.csv", bridge)]:
    print(f"\n--- {fname} ---")
    found_genes = df[df['gene_symbol'].isin(ALL_CU_GENES)]
    print(f"铜死亡基因总数: {len(found_genes['gene_symbol'].unique())} / {len(ALL_CU_GENES)}")
    if not found_genes.empty:
        cu_gene_stats = found_genes.groupby('gene_symbol').agg(
            pathway_count=('pathway_name', 'nunique'),
            mean_score=('score', 'mean'),
            max_score=('score', 'max')
        ).sort_values('max_score', ascending=False)
        print(cu_gene_stats.to_string())
    else:
        print("  未找到任何铜死亡基因")

    # Search for copper-related pathways
    copper_kw = ["copper", "cuproptosis", "CU"]
    copper_pathways = df[df['pathway_name'].str.contains('|'.join(copper_kw), case=False, na=False)]
    if not copper_pathways.empty:
        print(f"\n  铜相关通路 ({copper_pathways['pathway_name'].nunique()} unique):")
        for _, row in copper_pathways.head(20).iterrows():
            print(f"    {row['gene_symbol']}: {row['pathway_name'][:70]} | score={row['score']:.4f}")

print("\n" + "=" * 60)
print("CIRI候选文件中铜死亡基因的详细通路预测")
print("=" * 60)
if 'gene_symbol' in ciri.columns:
    cu_in_ciri = ciri[ciri['gene_symbol'].isin(ALL_CU_GENES)]
    if not cu_in_ciri.empty:
        for gene in sorted(cu_in_ciri['gene_symbol'].unique()):
            rows = cu_in_ciri[cu_in_ciri['gene_symbol'] == gene].sort_values('score', ascending=False)
            print(f"\n  {gene} ({len(rows)} pathways):")
            for _, r in rows.head(5).iterrows():
                print(f"    #{int(r['rank'])} {r['pathway_name'][:70]:70s} score={r['score']:.4f}")

print("\nDone!")