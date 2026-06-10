"""
检查预测逻辑：为什么DLAT、LIPT1、SLC31A1不在预测结果中
"""
import pickle
import csv
import torch
from pathlib import Path

BASE_DIR = Path("C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙")
PROCESSED_DIR = BASE_DIR / "processed"
RESULTS_DIR = BASE_DIR / "results"

MISSING_GENES = {"DLAT", "LIPT1", "SLC31A1"}

# 1. 读取gene_symbols.pkl
with open(PROCESSED_DIR / "gene_symbols.pkl", "rb") as f:
    gene_symbols = pickle.load(f)

# 2. 读取gene_to_idx.pkl
with open(PROCESSED_DIR / "gene_to_idx.pkl", "rb") as f:
    gene_to_idx = pickle.load(f)

# 3. 读取labels.csv
labels = {}
with open(PROCESSED_DIR / "labels.csv", 'r', encoding='utf-8-sig') as f:
    reader = csv.DictReader(f)
    for row in reader:
        labels[row["GeneSymbol"]] = int(row["Label"])

# 4. 构建y张量（与主脚本相同逻辑）
y = torch.full((len(gene_symbols),), -1, dtype=torch.long)
for i, gene in enumerate(gene_symbols):
    if gene in labels:
        y[i] = labels[gene]

print(f"y张量形状: {y.shape}")
print(f"y张量标签分布:")
unique, counts = y.unique(return_counts=True)
for lbl, count in zip(unique.tolist(), counts.tolist()):
    print(f"  标签 {lbl}: {count} 个基因")

# 5. 检查缺失基因的索引和标签
print(f"\n缺失的铜死亡基因:")
for gene in sorted(MISSING_GENES):
    if gene in gene_to_idx:
        idx = gene_to_idx[gene]
        lbl = y[idx].item()
        print(f"  {gene}: 索引={idx}, y[{idx}]={lbl}")
    else:
        print(f"  {gene}: 不在gene_to_idx中")

# 6. 模拟预测逻辑
prediction_labels = [-1, 2]  # 从prediction_config.json读取
prediction_mask = torch.zeros_like(y, dtype=torch.bool)
for lbl in prediction_labels:
    prediction_mask |= (y == lbl)

prediction_indices = torch.where(prediction_mask)[0]
print(f"\n预测掩码统计:")
print(f"  预测基因数: {prediction_mask.sum().item()}")

# 检查缺失基因是否在预测掩码中
print(f"\n缺失基因是否在预测掩码中:")
for gene in sorted(MISSING_GENES):
    if gene in gene_to_idx:
        idx = gene_to_idx[gene]
        in_mask = prediction_mask[idx].item()
        lbl = y[idx].item()
        print(f"  {gene}: 索引={idx}, 标签={lbl}, 在预测掩码中={'是' if in_mask else '否'}")

# 7. 检查预测结果文件
print(f"\n检查预测结果文件:")
pred_file = RESULTS_DIR / "all_unknown_predictions.csv"
if pred_file.exists():
    pred_genes = set()
    with open(pred_file, 'r', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        for row in reader:
            pred_genes.add(row["GeneSymbol"])
    
    print(f"  预测结果包含 {len(pred_genes)} 个基因")
    
    for gene in sorted(MISSING_GENES):
        if gene in pred_genes:
            print(f"  ✓ {gene}: 在预测结果中")
        else:
            print(f"  ✗ {gene}: 不在预测结果中")
            
            # 诊断原因
            if gene in gene_to_idx:
                idx = gene_to_idx[gene]
                lbl = y[idx].item()
                print(f"    - 索引: {idx}")
                print(f"    - 标签: {lbl}")
                print(f"    - 预测标签列表: {prediction_labels}")
                print(f"    - 标签匹配: {lbl in prediction_labels}")
else:
    print(f"  预测结果文件不存在")
