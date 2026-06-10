"""检查预测结果是否包含所有标签为-1的基因 - 修复版"""
from pathlib import Path

PROCESSED_DIR = Path("processed")
RESULTS_DIR = Path("results")
MISSING_GENES = {"DLAT", "LIPT1", "SLC31A1"}

def main():
    print("="*60)
    print("交叉验证：labels.csv vs 预测结果")
    print("="*60)
    
    # 读取labels.csv中所有标签为-1的基因
    with open(PROCESSED_DIR / "labels.csv", 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    header = lines[0].strip().split(',')
    print(f"labels.csv 列名: {header}")
    
    # 找到列索引
    gene_col = None
    label_col = None
    for i, col in enumerate(header):
        if 'gene' in col.lower() or 'symbol' in col.lower():
            gene_col = i
        if 'label' in col.lower():
            label_col = i
    
    if gene_col is None:
        gene_col = 0
    if label_col is None:
        label_col = 1
    
    print(f"使用列: GeneSymbol=列{gene_col}, Label=列{label_col}")
    
    label_minus1_genes = []
    for line in lines[1:]:
        parts = line.strip().split(',')
        if len(parts) > max(gene_col, label_col):
            if parts[label_col] == "-1":
                label_minus1_genes.append(parts[gene_col])
    
    label_minus1_set = set(label_minus1_genes)
    print(f"\n[1] labels.csv中标签=-1的基因: {len(label_minus1_genes)} 个")
    
    # 检查缺失基因是否在标签-1集合中
    for gene in sorted(MISSING_GENES):
        if gene in label_minus1_set:
            print(f"  ✓ {gene}: 标签=-1")
        else:
            print(f"  ✗ {gene}: 标签不是-1")
    
    # 读取预测结果
    pred_file = RESULTS_DIR / "all_unknown_predictions.csv"
    if not pred_file.exists():
        print(f"\n✗ 预测结果文件不存在: {pred_file}")
        return
    
    with open(pred_file, 'r', encoding='utf-8') as f:
        pred_lines = f.readlines()
    
    pred_header = pred_lines[0].strip().split(',')
    print(f"\n预测结果列名: {pred_header}")
    
    # 找到GeneSymbol列
    gene_col_idx = None
    for i, col in enumerate(pred_header):
        if 'gene' in col.lower() or 'symbol' in col.lower():
            gene_col_idx = i
            break
    
    if gene_col_idx is None:
        gene_col_idx = 1  # 默认第二列
    
    pred_genes = []
    for line in pred_lines[1:]:
        parts = line.strip().split(',')
        if len(parts) > gene_col_idx:
            pred_genes.append(parts[gene_col_idx])
    
    pred_gene_set = set(pred_genes)
    print(f"\n[2] 预测结果中的基因: {len(pred_genes)} 个")
    
    # 检查缺失基因是否在预测结果中
    for gene in sorted(MISSING_GENES):
        if gene in pred_gene_set:
            print(f"  ✓ {gene}: 在预测结果中")
        else:
            print(f"  ✗ {gene}: 不在预测结果中")
    
    # 找出差异
    missing_from_pred = label_minus1_set - pred_gene_set
    extra_in_pred = pred_gene_set - label_minus1_set
    
    print(f"\n[3] 差异分析:")
    print(f"  在labels.csv中标签=-1但不在预测结果中: {len(missing_from_pred)} 个")
    print(f"  在预测结果中但不在labels.csv标签=-1中: {len(extra_in_pred)} 个")
    
    if missing_from_pred:
        print(f"\n  缺失的基因(前30个):")
        for gene in sorted(missing_from_pred)[:30]:
            print(f"    {gene}")
        
        missing_three = MISSING_GENES & missing_from_pred
        if missing_three:
            print(f"\n  ⚠ 确认：以下基因在labels.csv中标签=-1但不在预测结果中:")
            for gene in sorted(missing_three):
                print(f"    {gene}")
    
    # 检查predict_mask.txt
    print(f"\n[4] 检查predict_mask.txt:")
    mask_file = PROCESSED_DIR / "predict_mask.txt"
    if mask_file.exists():
        with open(mask_file, 'r', encoding='utf-8') as f:
            mask_genes = [line.strip() for line in f if line.strip()]
        mask_gene_set = set(mask_genes)
        
        print(f"  predict_mask.txt中的基因: {len(mask_genes)} 个")
        
        for gene in sorted(MISSING_GENES):
            if gene in mask_gene_set:
                print(f"  ✓ {gene}: 在predict_mask.txt中")
            else:
                print(f"  ✗ {gene}: 不在predict_mask.txt中")
        
        mask_minus_pred = mask_gene_set - pred_gene_set
        print(f"\n  在predict_mask.txt中但不在预测结果中: {len(mask_minus_pred)} 个")
        
        if mask_minus_pred:
            missing_three_in_mask = MISSING_GENES & mask_minus_pred
            if missing_three_in_mask:
                print(f"\n  ⚠ 关键发现：以下基因在predict_mask.txt中但不在预测结果中:")
                for gene in sorted(missing_three_in_mask):
                    print(f"    {gene}")
    else:
        print(f"  ✗ predict_mask.txt不存在")

if __name__ == "__main__":
    main()
