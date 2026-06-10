"""深度诊断：检查预测时labels.csv和gene_symbols的顺序一致性"""
import pickle
from pathlib import Path

PROCESSED_DIR = Path("processed")
MISSING_GENES = {"DLAT", "LIPT1", "SLC31A1"}

def main():
    print("="*60)
    print("诊断：labels.csv与gene_symbols.pkl的顺序一致性")
    print("="*60)
    
    # 读取gene_symbols.pkl
    with open(PROCESSED_DIR / "gene_symbols.pkl", 'rb') as f:
        gene_symbols = pickle.load(f)
    
    print(f"\n[1] gene_symbols.pkl:")
    print(f"  基因数: {len(gene_symbols)}")
    print(f"  类型: {type(gene_symbols)}")
    
    # 创建索引映射
    gene_to_idx = {g: i for i, g in enumerate(gene_symbols)}
    
    # 检查缺失基因的索引
    print(f"\n[2] 缺失基因在gene_symbols中的索引:")
    for gene in sorted(MISSING_GENES):
        if gene in gene_to_idx:
            idx = gene_to_idx[gene]
            print(f"  {gene}: 索引={idx}, gene_symbols[{idx}]={gene_symbols[idx]}")
        else:
            print(f"  {gene}: 不在gene_symbols中")
    
    # 读取labels.csv
    with open(PROCESSED_DIR / "labels.csv", 'r', encoding='utf-8') as f:
        label_lines = f.readlines()
    
    # 解析labels.csv
    labels_list = []
    for line in label_lines[1:]:
        parts = line.strip().split(',')
        if len(parts) >= 2:
            gene = parts[0].replace('\ufeff', '')
            label = parts[1]
            labels_list.append((gene, label))
    
    print(f"\n[3] labels.csv:")
    print(f"  基因数: {len(labels_list)}")
    
    # 检查labels.csv中缺失基因的位置
    print(f"\n[4] 缺失基因在labels.csv中的位置:")
    for i, (gene, label) in enumerate(labels_list):
        if gene in MISSING_GENES:
            print(f"  {gene}: 行{i}, 标签={label}")
    
    # 关键检查：gene_symbols和labels.csv的基因是否一一对应
    print(f"\n[5] 顺序一致性检查:")
    
    # 检查前10个基因
    print(f"  前10个基因对比:")
    for i in range(10):
        gs_gene = gene_symbols[i] if i < len(gene_symbols) else "N/A"
        lbl_gene = labels_list[i][0] if i < len(labels_list) else "N/A"
        match = "✓" if gs_gene == lbl_gene else "✗"
        print(f"    索引{i}: gene_symbols='{gs_gene}' vs labels.csv='{lbl_gene}' {match}")
    
    # 检查缺失基因索引位置的基因
    print(f"\n[6] 缺失基因索引位置的基因:")
    for gene in sorted(MISSING_GENES):
        if gene in gene_to_idx:
            idx = gene_to_idx[gene]
            if idx < len(labels_list):
                lbl_gene, lbl_value = labels_list[idx]
                match = "✓" if lbl_gene == gene else "✗"
                print(f"  {gene} (索引{idx}): labels.csv[{idx}]='{lbl_gene}', 标签={lbl_value} {match}")
            else:
                print(f"  {gene} (索引{idx}): 超出labels.csv范围")
    
    # 创建labels.csv的基因到索引映射
    labels_to_idx = {gene: i for i, (gene, _) in enumerate(labels_list)}
    
    # 检查两个映射是否一致
    print(f"\n[7] 索引映射一致性:")
    consistent = 0
    inconsistent = 0
    for gene in gene_symbols:
        if gene in gene_to_idx and gene in labels_to_idx:
            if gene_to_idx[gene] == labels_to_idx[gene]:
                consistent += 1
            else:
                inconsistent += 1
    
    print(f"  一致: {consistent}")
    print(f"  不一致: {inconsistent}")
    
    # 关键问题：预测时使用的是哪个labels？
    # 根据脚本代码，预测时应该是：
    # 1. 从gene_symbols.pkl获取基因列表
    # 2. 从labels.csv获取标签
    # 3. 创建y = torch.tensor(labels_df["Label"].values)
    # 4. prediction_mask = (y == -1)
    # 5. unknown_indices = torch.where(prediction_mask)[0]
    # 6. 结果使用gene_symbols[i] for i in unknown_indices
    
    # 问题可能在于：labels_df["Label"].values的顺序与gene_symbols的顺序不一致
    
    # 检查：如果labels.csv的顺序与gene_symbols不同，预测会出错
    print(f"\n[8] 预测逻辑模拟:")
    
    # 模拟脚本中的预测逻辑
    # 脚本中使用pd.read_csv读取labels.csv，然后与gene_symbols配对
    
    # 检查labels.csv中每个基因的标签
    label_dict = {gene: label for gene, label in labels_list}
    
    # 模拟创建y张量
    y_labels = []
    for gene in gene_symbols:
        lbl = label_dict.get(gene, '-1')  # 默认-1
        y_labels.append(int(lbl))
    
    # 统计标签分布
    from collections import Counter
    label_counts = Counter(y_labels)
    print(f"  模拟y张量标签分布:")
    for lbl, count in sorted(label_counts.items()):
        print(f"    标签={lbl}: {count}")
    
    # 检查-1标签的索引
    minus1_indices = [i for i, lbl in enumerate(y_labels) if lbl == -1]
    print(f"\n  标签=-1的索引数: {len(minus1_indices)}")
    
    # 检查缺失基因是否在-1索引中
    print(f"\n[9] 缺失基因在模拟y中的状态:")
    for gene in sorted(MISSING_GENES):
        if gene in gene_to_idx:
            idx = gene_to_idx[gene]
            lbl = y_labels[idx]
            print(f"  {gene}: 索引={idx}, 标签={lbl}")
            if lbl == -1:
                print(f"    ✓ 应该出现在预测结果中")
            else:
                print(f"    ✗ 不会出现在预测结果中（标签不是-1）")

if __name__ == "__main__":
    main()
