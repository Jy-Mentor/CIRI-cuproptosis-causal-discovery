#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os

# 使用当前目录
current_dir = os.getcwd()
output_dir = os.path.join(current_dir, "docking_results_merged")

# 创建输出目录
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# 定义文件路径
receptor_dir = os.path.join(current_dir, "output", "receptors")
ligand_dir = os.path.join(current_dir, "output", "docking_results")

# 要处理的PDB ID列表
pdb_ids = ["1zmc", "2flu", "3o3u", "3p1m", "3rzf", "5n2s", "6dvh", "6qcb", "8vmc", "9soq"]

for pdb_id in pdb_ids:
    # 构建文件路径
    receptor_file = os.path.join(receptor_dir, "{0}.pdb".format(pdb_id))
    ligand_file = os.path.join(ligand_dir, "{0}_complex.pdb".format(pdb_id))
    output_file = os.path.join(output_dir, "{0}_complex.pdb".format(pdb_id))
    
    # 检查文件是否存在
    if not os.path.exists(receptor_file):
        print("Error: {0}: Receptor file not found".format(pdb_id))
        continue
    if not os.path.exists(ligand_file):
        print("Error: {0}: Ligand file not found".format(pdb_id))
        continue
    
    # 读取受体文件
    receptor_lines = []
    try:
        with open(receptor_file, 'r') as f:
            for line in f:
                if line.startswith(('ATOM', 'HETATM', 'TER')):
                    receptor_lines.append(line)
    except Exception as e:
        print("Error reading receptor file {0}: {1}".format(receptor_file, e))
        continue
    
    # 读取配体文件
    ligand_lines = []
    try:
        with open(ligand_file, 'r') as f:
            for line in f:
                if line.startswith(('ATOM', 'HETATM')) and 'UNL' in line[17:20]:
                    ligand_lines.append(line)
    except Exception as e:
        print("Error reading ligand file {0}: {1}".format(ligand_file, e))
        continue
    
    # 如果没有找到UNL残基，尝试读取所有原子
    if not ligand_lines:
        print("Warning: {0}: No UNL residue found in ligand file, trying to copy all ATOM...".format(pdb_id))
        try:
            with open(ligand_file, 'r') as f:
                for line in f:
                    if line.startswith(('ATOM', 'HETATM')):
                        ligand_lines.append(line)
        except Exception as e:
            print("Error reading ligand file {0}: {1}".format(ligand_file, e))
            continue
    
    # 合并文件
    merged_lines = receptor_lines + ['TER\n'] + ligand_lines + ['END\n']
    
    # 写入输出文件
    try:
        with open(output_file, 'w') as f:
            f.writelines(merged_lines)
    except Exception as e:
        print("Error writing output file {0}: {1}".format(output_file, e))
        continue
    
    # 统计原子数
    rec_atoms = len([l for l in receptor_lines if l.startswith('ATOM')])
    lig_atoms = len(ligand_lines)
    print("Success: {0}: {1} receptor atoms + {2} ligand atoms -> {3}".format(pdb_id, rec_atoms, lig_atoms, output_file))

print("\n===============================")
print("Merge completed! Please use files in docking_results_merged directory to re-run analysis")
print("===============================")
