import pandas as pd
import numpy as np

csv = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\CIRI-cuproptosis-causal-discovery\L1_phenotype_anchoring\GSE97537_cuproptosis_DEGs.csv"
df = pd.read_csv(csv)

log2fc_threshold = 0.585
old_sig = df["Significant"].value_counts()
print("Old Significant:", old_sig.to_dict())

df["adj_num"] = pd.to_numeric(df["adj.P.Val"], errors="coerce")
df["Significant"] = np.where(
    (df["adj_num"] < 0.05) & (df["log2FC"].abs() >= log2fc_threshold),
    "是", "否"
)
df = df.drop(columns=["adj_num"])

new_sig = df["Significant"].value_counts()
print("New Significant:", new_sig.to_dict())

for _, row in df.iterrows():
    if row["log2FC"] != 0 and abs(row["log2FC"]) < log2fc_threshold and row["Significant"] == "否":
        print(f"  {row['Human_Gene']}: log2FC={row['log2FC']}, Significant={row['Significant']}")

df.to_csv(csv, index=False)
print("CSV updated")