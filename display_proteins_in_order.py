import pandas as pd

# 读取之前生成的蛋白度值数据
data_file = 'protein_degrees.csv'
print(f"读取蛋白度值数据: {data_file}")

try:
    protein_data = pd.read_csv(data_file)
    print(f"成功读取文件，包含 {len(protein_data)} 个蛋白")
except Exception as e:
    print(f"读取文件时出错: {e}")
    exit()

# 用户指定的蛋白顺序
user_order = [
    "NR3C1", "STAT3", "MAPK14", "SMARCA4", "CASP9",
    "PTPRC", "MDM2", "CTSB", "ARG1", "NFE2L2",
    "ESR1", "APP", "CDC42", "STAT5A", "PTGS2",
    "NFKB1", "JAK2", "LYN", "PPARD", "NOS3",
    "HMGCR", "IL1B", "PPARA", "MMP9", "PPARG",
    "IDO1", "HMOX1", "CASP8", "PTGS1", "EGR1"
]

print(f"\n用户指定的蛋白顺序数量: {len(user_order)}")

# 按照用户指定的顺序重新排列数据
ordered_data = []
for protein in user_order:
    # 查找该蛋白的数据
    protein_row = protein_data[protein_data['protein'] == protein]
    if not protein_row.empty:
        ordered_data.append(protein_row.iloc[0])
    else:
        # 如果找不到，创建一个空行
        empty_row = pd.Series({'protein': protein, 'degree_centrality': None, 'found': False})
        ordered_data.append(empty_row)

# 转换为DataFrame
ordered_df = pd.DataFrame(ordered_data)

# 显示结果
print("\n按照用户指定顺序呈现的蛋白度值:")
print(ordered_df[['protein', 'degree_centrality']])

# 保存排序后的结果
ordered_df.to_csv('protein_degrees_ordered.csv', index=False)
print(f"\n排序后的结果已保存到 protein_degrees_ordered.csv")
