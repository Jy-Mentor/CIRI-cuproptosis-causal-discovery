#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os 
import shutil 

BASE = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\output" 
RECEPTOR_DIR = os.path.join(BASE, "receptors") 
LIGAND_DIR = os.path.join(BASE, "docking_results")  # 你下载的"假"complex文件（实际只有配体） 
OUTPUT_DIR = os.path.join(BASE, "docking_results_merged") 
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR) 

pdb_ids = ["1zmc", "2flu", "3o3u", "3p1m", "3rzf", "5n2s", "6dvh", "6qcb", "8vmc", "9soq"] 

for pdb_id in pdb_ids: 
    receptor_file = os.path.join(RECEPTOR_DIR, "{0}.pdb".format(pdb_id)) 
    ligand_file = os.path.join(LIGAND_DIR, "{0}_complex.pdb".format(pdb_id)) 
    output_file = os.path.join(OUTPUT_DIR, "{0}_complex.pdb".format(pdb_id)) 
    
    # 检查文件是否存在 
    if not os.path.exists(receptor_file): 
        print("❌ {0}: 受体文件不存在".format(pdb_id)) 
        continue 
    if not os.path.exists(ligand_file): 
        print("❌ {0}: 配体文件不存在".format(pdb_id)) 
        continue 
    
    # 读取受体（只保留ATOM/HETATM/TER记录，去除CONNECT等） 
    with open(receptor_file, 'r') as f: 
        receptor_lines = [l for l in f if l.startswith(('ATOM', 'HETATM', 'TER'))] 
    
    # 读取配体（从下载的文件中提取UNL部分） 
    with open(ligand_file, 'r') as f: 
        lines = f.readlines() 
        # 提取UNL残基（配体） 
        ligand_lines = [l for l in lines if l.startswith(('ATOM', 'HETATM')) and 'UNL' in l[17:20]] 
    
    if not ligand_lines: 
        print("⚠️ {0}: 配体文件中没有UNL残基，尝试直接复制所有ATOM...".format(pdb_id)) 
        ligand_lines = [l for l in lines if l.startswith(('ATOM', 'HETATM'))] 
    
    # 合并：受体 + 配体 + END 
    merged = receptor_lines + ['TER\n'] + ligand_lines + ['END\n'] 
    
    with open(output_file, 'w') as f: 
        f.writelines(merged) 
    
    rec_atoms = len([l for l in receptor_lines if l.startswith('ATOM')]) 
    lig_atoms = len(ligand_lines) 
    print("✅ {0}: 受体{1}原子 + 配体{2}原子 → {3}".format(pdb_id, rec_atoms, lig_atoms, output_file)) 

print("\n===============================") 
print("合并完成！请使用 docking_results_merged 目录中的文件重新运行分析") 
print("===============================")