"""
修复铜死亡基因标签：将标签1改为-1
"""
import csv
from pathlib import Path

PROCESSED_DIR = Path("C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/processed")

CUPROPTOSIS_GENES = {
    "FDX1", "LIAS", "LIPT1", "DLAT", "PDHB", "PDHX", "SLC31A1",
    "ATP7A", "ATP7B", "ATOX1", "NFE2L2", "HIF1A", "MTOR", "NFKB1", "GPX4"
}

# 读取labels
labels = []
with open(PROCESSED_DIR / "labels.csv", 'r', encoding='utf-8-sig') as f:
    reader = csv.DictReader(f)
    for row in reader:
        labels.append(row)

# 修改铜死亡基因标签
n_modified = 0
for row in labels:
    if row["GeneSymbol"] in CUPROPTOSIS_GENES and row["Label"] == "1":
        row["Label"] = "-1"
        n_modified += 1
        print(f"修改 {row['GeneSymbol']}: 1 → -1")

print(f"\n共修改 {n_modified} 个铜死亡基因标签")

# 保存修改后的labels
with open(PROCESSED_DIR / "labels.csv", 'w', newline='', encoding='utf-8-sig') as f:
    writer = csv.DictWriter(f, fieldnames=["GeneSymbol", "Label"])
    writer.writeheader()
    writer.writerows(labels)

print(f"已保存修改后的labels.csv")

# 验证
label_dist = {}
for row in labels:
    lbl = row["Label"]
    label_dist[lbl] = label_dist.get(lbl, 0) + 1

print(f"\n修改后标签分布:")
for lbl, count in sorted(label_dist.items()):
    print(f"  标签 {lbl}: {count}")

# 验证铜死亡基因标签
print(f"\n铜死亡基因标签验证:")
for row in labels:
    if row["GeneSymbol"] in CUPROPTOSIS_GENES:
        print(f"  {row['GeneSymbol']}: {row['Label']}")
