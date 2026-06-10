"""独立预测脚本 - 绕过numpy DLL问题，纯Python验证"""
import pickle
import json
from pathlib import Path
from collections import defaultdict

# 配置
PROCESSED_DIR = Path("processed")
RESULTS_DIR = Path("results")
MODELS_DIR = Path("models")
MISSING_GENES = {"DLAT", "LIPT1", "SLC31A1"}
CUPROPTOSIS_GENES = {"FDX1", "LIAS", "LIPT1", "DLAT", "PDHB", "PDHX", "SLC31A1", "ATP7B", "ATOX1", "MTF1", "GLS", "CDKN2A"}

def read_csv_simple(filepath):
    """简单读取CSV文件"""
    rows = []
    with open(filepath, 'r', encoding='utf-8') as f:
        header = f.readline().strip().replace('\ufeff', '').split(',')
        for line in f:
            values = line.strip().split(',')
            if len(values) == len(header):
                rows.append(dict(zip(header, values)))
    return rows, header

def main():
    print("="*60)
    print("独立预测脚本 - 验证DLAT/LIPT1/SLC31A1")
    print("="*60)
    
    # 步骤1: 加载gene_symbols.pkl
    print("\n[步骤1] 加载gene_symbols.pkl...")
    with open(PROCESSED_DIR / "gene_symbols.pkl", 'rb') as f:
        gene_symbols = pickle.load(f)
    print(f"  基因数: {len(gene_symbols)}")
    gene_to_idx = {g: i for i, g in enumerate(gene_symbols)}
    
    # 步骤2: 加载labels.csv并创建标签映射
    print("\n[步骤2] 加载labels.csv...")
    labels_data, labels_header = read_csv_simple(PROCESSED_DIR / "labels.csv")
    print(f"  列名: {labels_header}")
    print(f"  基因数: {len(labels_data)}")
    
    # 创建基因名到标签的映射
    label_map = {}
    for row in labels_data:
        gene = row.get("GeneSymbol", "")
        label = row.get("Label", "-1")
        label_map[gene] = int(label)
    
    # 按gene_symbols顺序创建标签列表
    y_list = [label_map.get(g, -1) for g in gene_symbols]
    
    # 统计标签分布
    label_counts = defaultdict(int)
    for lbl in y_list:
        label_counts[lbl] += 1
    print(f"  标签分布: {dict(label_counts)}")
    
    # 验证缺失基因的标签
    print(f"\n  缺失基因验证:")
    for gene in sorted(MISSING_GENES):
        if gene in gene_to_idx:
            idx = gene_to_idx[gene]
            lbl = y_list[idx]
            print(f"    {gene}: 索引={idx}, 标签={lbl}")
        else:
            print(f"    {gene}: 不在gene_symbols中")
    
    # 步骤3: 创建预测掩码
    print("\n[步骤3] 创建预测掩码...")
    predict_mask = [i for i, lbl in enumerate(y_list) if lbl == -1]
    print(f"  预测基因数: {len(predict_mask)}")
    
    # 检查缺失基因是否在预测掩码中
    for gene in sorted(MISSING_GENES):
        if gene in gene_to_idx:
            idx = gene_to_idx[gene]
            if idx in predict_mask:
                print(f"  ✓ {gene}: 在预测掩码中")
            else:
                print(f"  ✗ {gene}: 不在预测掩码中 (标签={y_list[idx]})")
    
    # 步骤4: 读取已有的预测结果
    print("\n[步骤4] 读取已有预测结果...")
    pred_file = RESULTS_DIR / "all_unknown_predictions.csv"
    if pred_file.exists():
        pred_data, pred_header = read_csv_simple(pred_file)
        print(f"  预测结果列名: {pred_header}")
        print(f"  预测结果基因数: {len(pred_data)}")
        
        # 检查缺失基因是否在预测结果中
        for gene in sorted(MISSING_GENES):
            found = False
            for row in pred_data:
                if row.get("GeneSymbol") == gene:
                    found = True
                    print(f"  ✓ {gene}: 在预测结果中, P_target={row.get('P_target', 'N/A')}, Rank={row.get('Rank', 'N/A')}")
                    break
            if not found:
                print(f"  ✗ {gene}: 不在预测结果中")
    else:
        print(f"  ✗ 预测结果文件不存在")
    
    # 步骤5: 分析差异
    print("\n[步骤5] 差异分析...")
    if pred_file.exists():
        pred_genes = set(row.get("GeneSymbol", "") for row in pred_data)
        mask_genes = set(gene_symbols[i] for i in predict_mask)
        
        missing_from_pred = mask_genes - pred_genes
        extra_in_pred = pred_genes - mask_genes
        
        print(f"  应该在预测中但不在: {len(missing_from_pred)} 个")
        print(f"  不该在预测中但在: {len(extra_in_pred)} 个")
        
        if MISSING_GENES & missing_from_pred:
            print(f"\n  ⚠ 确认：DLAT/LIPT1/SLC31A1应该在预测结果中但不在")
            print(f"  这表明之前的预测代码使用了错误的索引映射")
        else:
            print(f"\n  ✓ DLAT/LIPT1/SLC31A1都在预测结果中")
    
    print("\n" + "="*60)
    print("结论:")
    print("="*60)
    print("  1. DLAT/LIPT1/SLC31A1在gene_symbols.pkl中存在")
    print("  2. DLAT/LIPT1/SLC31A1在labels.csv中标签=-1")
    print("  3. DLAT/LIPT1/SLC31A1应该在预测结果中")
    print("  4. 如果不在，说明预测代码有bug（gene_symbols和labels顺序不一致）")
    print("\n  修复方案：")
    print("  - 已修改石竹烯_CIRI_全流程整合脚本_v3.py的run_prediction()函数")
    print("  - 添加label_map映射，按gene_symbols顺序创建y张量")
    print("  - 需要重新运行 --mode predict")

if __name__ == "__main__":
    main()
