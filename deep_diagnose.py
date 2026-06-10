"""
深入诊断：为什么DLAT, LIPT1, SLC31A1不在预测结果中
"""
import pandas as pd
import pickle
import torch
from pathlib import Path

PROCESSED_DIR = Path("C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/processed")
RESULTS_DIR = Path("C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/results")

MISSING_GENES = {"DLAT", "LIPT1", "SLC31A1"}

print("="*70)
print("深度诊断：缺失的铜死亡基因")
print("="*70)

# 1. 检查所有相关文件
print("\n[1] 文件存在性检查:")
files_to_check = [
    "gene_symbols.pkl", "gene_to_idx.pkl", "labels.csv", 
    "node_features.csv", "predict_mask.txt", "edge_index.pt"
]
for f in files_to_check:
    exists = (PROCESSED_DIR / f).exists()
    print(f"  {'✓' if exists else '✗'} {f}")

# 2. 读取gene_symbols.pkl
print("\n[2] gene_symbols.pkl:")
with open(PROCESSED_DIR / "gene_symbols.pkl", "rb") as f:
    gene_symbols = pickle.load(f)
print(f"  基因数: {len(gene_symbols)}")
for gene in sorted(MISSING_GENES):
    if gene in gene_symbols:
        idx = gene_symbols.index(gene)
        print(f"  ✓ {gene}: 索引={idx}")
    else:
        print(f"  ✗ {gene}: 不在gene_symbols中")

# 3. 读取gene_to_idx.pkl
print("\n[3] gene_to_idx.pkl:")
with open(PROCESSED_DIR / "gene_to_idx.pkl", "rb") as f:
    gene_to_idx = pickle.load(f)
print(f"  基因数: {len(gene_to_idx)}")
for gene in sorted(MISSING_GENES):
    if gene in gene_to_idx:
        idx = gene_to_idx[gene]
        print(f"  ✓ {gene}: 索引={idx}")
    else:
        print(f"  ✗ {gene}: 不在gene_to_idx中")

# 4. 读取labels.csv
print("\n[4] labels.csv:")
labels_df = pd.read_csv(PROCESSED_DIR / "labels.csv")
print(f"  基因数: {len(labels_df)}")
print(f"  标签分布:")
print(f"    {labels_df['Label'].value_counts().sort_index().to_dict()}")
for gene in sorted(MISSING_GENES):
    row = labels_df[labels_df["GeneSymbol"] == gene]
    if len(row) > 0:
        lbl = row.iloc[0]["Label"]
        print(f"  ✓ {gene}: 标签={lbl}")
    else:
        print(f"  ✗ {gene}: 不在labels.csv中")

# 5. 读取node_features.csv
print("\n[5] node_features.csv:")
features_df = pd.read_csv(PROCESSED_DIR / "node_features.csv")
print(f"  基因数: {len(features_df)}")
for gene in sorted(MISSING_GENES):
    row = features_df[features_df["GeneSymbol"] == gene]
    if len(row) > 0:
        print(f"  ✓ {gene}: 在特征矩阵中")
    else:
        print(f"  ✗ {gene}: 不在特征矩阵中")

# 6. 读取predict_mask.txt
print("\n[6] predict_mask.txt:")
with open(PROCESSED_DIR / "predict_mask.txt", 'r') as f:
    predict_genes = set(line.strip() for line in f if line.strip())
print(f"  预测基因数: {len(predict_genes)}")
for gene in sorted(MISSING_GENES):
    if gene in predict_genes:
        print(f"  ✓ {gene}: 在predict_mask中")
    else:
        print(f"  ✗ {gene}: 不在predict_mask中")

# 7. 构建y张量（模拟主脚本逻辑）
print("\n[7] 构建y张量（模拟主脚本逻辑）:")
labels_dict = dict(zip(labels_df["GeneSymbol"], labels_df["Label"]))
y = torch.full((len(gene_symbols),), -1, dtype=torch.long)
for i, gene in enumerate(gene_symbols):
    if gene in labels_dict:
        y[i] = labels_dict[gene]

print(f"  y张量形状: {y.shape}")
print(f"  y张量标签分布:")
unique, counts = y.unique(return_counts=True)
for lbl, count in zip(unique.tolist(), counts.tolist()):
    print(f"    标签 {lbl}: {count}")

# 检查缺失基因的y值
print(f"\n  缺失基因的y值:")
for gene in sorted(MISSING_GENES):
    if gene in gene_to_idx:
        idx = gene_to_idx[gene]
        lbl = y[idx].item()
        print(f"    {gene}: 索引={idx}, y[{idx}]={lbl}")
    else:
        print(f"    {gene}: 不在gene_to_idx中")

# 8. 模拟预测掩码
print("\n[8] 模拟预测掩码:")
prediction_labels = [-1, 2]  # 从prediction_config.json读取
prediction_mask = torch.zeros_like(y, dtype=torch.bool)
for lbl in prediction_labels:
    prediction_mask |= (y == lbl)

prediction_indices = torch.where(prediction_mask)[0]
print(f"  预测基因数: {prediction_mask.sum().item()}")

for gene in sorted(MISSING_GENES):
    if gene in gene_to_idx:
        idx = gene_to_idx[gene]
        in_mask = prediction_mask[idx].item()
        lbl = y[idx].item()
        print(f"  {gene}: 索引={idx}, 标签={lbl}, 在预测掩码={'是' if in_mask else '否'}")

# 9. 检查预测结果
print("\n[9] 预测结果文件:")
pred_df = pd.read_csv(RESULTS_DIR / "all_unknown_predictions.csv")
print(f"  预测结果基因数: {len(pred_df)}")
pred_genes = set(pred_df["GeneSymbol"])

for gene in sorted(MISSING_GENES):
    if gene in pred_genes:
        row = pred_df[pred_df["GeneSymbol"] == gene].iloc[0]
        print(f"  ✓ {gene}: Rank={row['Rank']}, P={row['P_target']}")
    else:
        print(f"  ✗ {gene}: 不在预测结果中")

# 10. 交叉验证
print("\n[10] 交叉验证:")
for gene in sorted(MISSING_GENES):
    issues = []
    
    in_symbols = gene in gene_symbols
    in_idx = gene in gene_to_idx
    in_labels = gene in labels_dict
    in_features = gene in features_df["GeneSymbol"].values
    in_predict_mask = gene in predict_genes
    in_pred_result = gene in pred_genes
    
    if not in_symbols:
        issues.append("不在gene_symbols.pkl")
    if not in_idx:
        issues.append("不在gene_to_idx.pkl")
    if not in_labels:
        issues.append("不在labels.csv")
    if not in_features:
        issues.append("不在node_features.csv")
    if not in_predict_mask:
        issues.append("不在predict_mask.txt")
    
    if in_symbols and in_idx and in_labels and in_features and in_predict_mask and not in_pred_result:
        # 所有条件都满足，但不在结果中 - 检查标签
        idx = gene_to_idx[gene]
        lbl = y[idx].item()
        if lbl not in prediction_labels:
            issues.append(f"标签={lbl} 不在预测标签{prediction_labels}中")
    
    if issues:
        print(f"  {gene}:")
        for issue in issues:
            print(f"    - {issue}")
    else:
        print(f"  ✓ {gene}: 所有检查通过")

print("\n" + "="*70)
print("诊断完成")
print("="*70)
