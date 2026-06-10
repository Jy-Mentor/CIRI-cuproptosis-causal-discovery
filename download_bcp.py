#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
from rdkit import Chem
from rdkit.Chem import AllChem

# 定义输入和输出目录
output_dir = "output/new_receptors"
os.makedirs(output_dir, exist_ok=True)

# BCP的SMILES式
bcp_smiles = "[H][C@@]12CC\\C(C)=C\\CCC(=C)[C@H]1CC2(C)C"

print("Generating BCP 3D structure from SMILES...")
print(f"SMILES: {bcp_smiles}")
print(f"Output directory: {os.path.abspath(output_dir)}")
print()

try:
    # 从SMILES创建分子
    mol = Chem.MolFromSmiles(bcp_smiles)
    if mol is None:
        raise Exception("Failed to create molecule from SMILES")
    
    # 添加氢（仅用于3D构象生成）
    mol = Chem.AddHs(mol)
    
    # 生成3D构象
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)
    
    # 移除氢（按照用户要求）
    mol = Chem.RemoveHs(mol)
    
    # 保存为PDB文件
    output_file = os.path.join(output_dir, "bcp.pdb")
    with open(output_file, "w") as f:
        f.write(Chem.MolToPDBBlock(mol))
    
    print(f"Successfully generated BCP PDB file: {output_file}")
    
except Exception as e:
    print(f"Error generating BCP structure: {e}")

print("\nBCP structure generation complete!")
