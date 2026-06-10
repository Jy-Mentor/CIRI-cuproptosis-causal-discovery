# -*- coding: utf-8 -*-
import pandas as pd, numpy as np

esm = pd.read_csv(r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\subgraph_embeddings.csv")
esm["gene_symbol"] = esm["gene_symbol"].str.upper()
ft = pd.read_csv(r"C:\Users\Jy-Mentor-7\Desktop\随机森林\gene_features_table_with_gat_emb.csv")
ft["gene_symbol"] = ft["gene_symbol"].str.upper()

ft_genes = set(ft["gene_symbol"])
esm_genes = set(esm["gene_symbol"])
missing = sorted(ft_genes - esm_genes)
print(f"Missing ESM-2 genes: {len(missing)}")

for g in missing[:30]:
    row = ft[ft["gene_symbol"] == g].iloc[0]
    tags = []
    if row["is_drug_target"]:
        tags.append("DT+")
    if row["is_disease_gene"]:
        tags.append("DG+")
    print(f"  {g}  {' '.join(tags)}")
if len(missing) > 30:
    print(f"  ... and {len(missing) - 30} more")

with open(r"C:\Users\Jy-Mentor-7\Desktop\subgraph_genes.txt") as f:
    sg = set(l.strip().upper() for l in f if l.strip() and not l.startswith("#") and l.strip() != "gene_symbol")
in_sg = [g for g in missing if g in sg]
print(f"\nIn subgraph_genes.txt: {len(in_sg)}/{len(missing)}")

# DT/DG stats
dt_missing = [g for g in missing if ft.loc[ft["gene_symbol"] == g, "is_drug_target"].values[0] == 1]
dg_missing = [g for g in missing if ft.loc[ft["gene_symbol"] == g, "is_disease_gene"].values[0] == 1]
print(f"DT+ missing: {len(dt_missing)}")
print(f"DG+ missing: {len(dg_missing)}")

# Save missing gene list for ESM-2 script
with open("missing_esm2_genes.txt", "w") as f:
    f.write("gene_symbol\n")
    for g in missing:
        f.write(g + "\n")
print(f"\nSaved missing_esm2_genes.txt ({len(missing)} genes)")