"""深度诊断：检查labels.csv中这3个基因的标签，以及预测逻辑中的问题"""
from pathlib import Path

PROCESSED_DIR = Path("processed")
MISSING_GENES = {"DLAT", "LIPT1", "SLC31A1"}

def main():
    print("="*60)
    print("深度诊断：标签值和预测逻辑")
    print("="*60)
    
    # 读取labels.csv
    labels_file = PROCESSED_DIR / "labels.csv"
    if not labels_file.exists():
        print("✗ labels.csv 不存在")
        return
    
    with open(labels_file, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    header = lines[0].strip()
    print(f"\n[1] labels.csv 信息:")
    print(f"  列名: {header}")
    print(f"  总行数: {len(lines) - 1}")
    
    # 查找这3个基因
    print(f"\n[2] 缺失基因的标签:")
    for line in lines[1:]:
        parts = line.strip().split(',')
        if len(parts) >= 2 and parts[0] in MISSING_GENES:
            gene = parts[0]
            label_str = parts[1]
            print(f"  {gene}: 原始标签='{label_str}' (长度={len(label_str)}, 编码={[hex(ord(c)) for c in label_str]})")
    
    # 统计不同标签的数量
    label_counts = {}
    for line in lines[1:]:
        parts = line.strip().split(',')
        if len(parts) >= 2:
            label = parts[1]
            label_counts[label] = label_counts.get(label, 0) + 1
    
    print(f"\n[3] 标签分布:")
    for label, count in sorted(label_counts.items()):
        print(f"  标签='{label}': {count} 个基因")
    
    # 检查ogb_config.json
    import json
    ogb_config_file = PROCESSED_DIR / "ogb_config.json"
    if ogb_config_file.exists():
        with open(ogb_config_file, 'r', encoding='utf-8') as f:
            ogb_config = json.load(f)
        
        print(f"\n[4] ogb_config.json:")
        print(f"  predict_labels: {ogb_config.get('predict_labels', [])}")
        
        # 检查预测标签是否与labels.csv中的标签匹配
        predict_labels = ogb_config.get('predict_labels', [])
        for pred_label in predict_labels:
            pred_label_str = str(pred_label)
            if pred_label_str in label_counts:
                print(f"  ✓ 预测标签 '{pred_label_str}' 在labels.csv中存在 ({label_counts[pred_label_str]} 个基因)")
            else:
                print(f"  ✗ 预测标签 '{pred_label_str}' 在labels.csv中不存在！")
                print(f"    可能的标签值: {list(label_counts.keys())}")
    
    # 分析可能的编码问题
    print(f"\n[5] 编码问题检查:")
    print(f"  正常'-1'的UTF-8编码: {hex(ord('-'))}, {hex(ord('1'))}")
    for line in lines[1:]:
        parts = line.strip().split(',')
        if len(parts) >= 2 and parts[0] in MISSING_GENES:
            gene = parts[0]
            label_str = parts[1]
            # 检查是否有不可见字符
            has_invisible = any(ord(c) < 32 or (ord(c) > 126 and ord(c) < 160) for c in label_str)
            if has_invisible:
                print(f"  ⚠ {gene}: 标签包含不可见字符！")
            else:
                print(f"  ✓ {gene}: 标签无不可见字符")
    
    print(f"\n[6] 结论:")
    # 检查-1是否存在
    if '-1' in label_counts:
        print(f"  标签'-1'存在，共{label_counts['-1']}个基因")
        missing_count = sum(1 for g in MISSING_GENES if g in [line.strip().split(',')[0] for line in lines[1:] if len(line.strip().split(',')) >= 2 and line.strip().split(',')[1] == '-1'])
        print(f"  其中缺失基因: {missing_count}/{len(MISSING_GENES)}")
        if missing_count == 0:
            print(f"  ⚠ 问题：这3个基因的标签不是'-1'！")
            # 找出它们的实际标签
            for line in lines[1:]:
                parts = line.strip().split(',')
                if len(parts) >= 2 and parts[0] in MISSING_GENES:
                    print(f"    {parts[0]}: 实际标签='{parts[1]}'")
    else:
        print(f"  ✗ 标签'-1'不存在！")

if __name__ == "__main__":
    main()
