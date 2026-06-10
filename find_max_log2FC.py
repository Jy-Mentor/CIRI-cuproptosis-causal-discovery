# -*- coding: utf-8 -*-
import pandas as pd
import sys

xlsx_path = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\deg_5datasets_summary.xlsx"
txt_path = r"C:\Users\Jy-Mentor-7\Desktop\GAT\all_genes.txt"
output_path = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\max_log2FC_result.xlsx"

xl = pd.ExcelFile(xlsx_path)

deg_sheets = []
for s in xl.sheet_names:
    df_check = pd.read_excel(xl, sheet_name=s, nrows=1)
    if "gene_symbol" in df_check.columns and "log2FoldChange" in df_check.columns:
        deg_sheets.append(s)
print(f"DEG sheets with log2FC: {deg_sheets}")

gene_log2fc_map = {}

for sheet in deg_sheets:
    df = pd.read_excel(xl, sheet_name=sheet)
    for _, row in df.iterrows():
        gene = str(row["gene_symbol"]).strip()
        lfc = row["log2FoldChange"]
        if pd.isna(lfc):
            continue
        abs_lfc = abs(lfc)
        if gene not in gene_log2fc_map or abs_lfc > gene_log2fc_map[gene]:
            gene_log2fc_map[gene] = abs_lfc

print(f"Total unique genes with log2FC in xlsx: {len(gene_log2fc_map)}")

target_genes = []
with open(txt_path, "r") as f:
    for line in f:
        gene = line.strip()
        if gene and gene != "gene":
            target_genes.append(gene)

print(f"Target genes from txt: {len(target_genes)}")

# Build case-insensitive lookup
gene_lower_map = {k.lower(): v for k, v in gene_log2fc_map.items()}

results = []
found_count = 0
for gene in target_genes:
    abs_lfc = gene_lower_map.get(gene.lower(), "NA")
    if abs_lfc != "NA":
        found_count += 1
    results.append({"gene": gene, "max_abs_log2FC": abs_lfc})

print(f"Matched: {found_count} / {len(target_genes)}")

result_df = pd.DataFrame(results)
result_df.to_excel(output_path, index=False)
print(f"Result saved to: {output_path}")