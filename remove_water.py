# -*- coding: utf-8 -*-
import os

# 定义要处理的PDB文件列表（使用相对路径）
pdb_files = [
    "output/processed_receptors/processed_adora1.pdb",
    "output/processed_receptors/processed_ctsd_1lyw.pdb",
    "output/processed_receptors/processed_ctsd_6qcb.pdb",
    "output/processed_receptors/processed_fasn_6c7p.pdb",
    "output/processed_receptors/processed_fasn_8vmc.pdb",
    "output/processed_receptors/processed_fdx1.pdb",
    "output/processed_receptors/processed_lias.pdb",
    "output/processed_receptors/processed_nfe2l2.pdb",
    "output/processed_receptors/processed_rage_3o3u.pdb",
    "output/processed_receptors/processed_rage_5cjy.pdb",
    "output/processed_receptors/processed_timp1.pdb",
    "output/processed_receptors/processed_tspo.pdb"
]

# 定义要删除的水分子残基名称
water_residues = {'HOH', 'WAT', 'DOD', 'H2O', 'SOL', 'TIP'}

def remove_water(input_file):
    """删除PDB文件中的水分子"""
    try:
        with open(input_file, 'r') as f:
            lines = f.readlines()
    except Exception as e:
        print("Error reading file:", input_file)
        print("Error:", str(e))
        return -1
    
    cleaned_lines = []
    removed_count = 0
    
    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            res_name = line[17:20].strip().upper()
            if res_name in water_residues:
                removed_count += 1
                continue
        cleaned_lines.append(line)
    
    # 写回原文件
    try:
        with open(input_file, 'w') as f:
            f.writelines(cleaned_lines)
        print("%25s | Removed %d water molecules" % (os.path.basename(input_file), removed_count))
        return removed_count
    except Exception as e:
        print("Error writing file:", input_file)
        print("Error:", str(e))
        return -1

# 批量处理
print("=" * 60)
print("Water Removal Report")
print("=" * 60)

total_removed = 0
for pdb_file in pdb_files:
    count = remove_water(pdb_file)
    if count > 0:
        total_removed += count

print("=" * 60)
print("Total water molecules removed:", total_removed)
print("=" * 60)