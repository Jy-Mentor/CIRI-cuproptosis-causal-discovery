#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from rdkit import Chem

# 读取SDF文件（含多个构象/重复属性）
input_file = 'output/new_receptors/BCP.sdf'
output_file = 'output/new_receptors/BCP_clean.sdf'

print(f"Reading SDF file: {input_file}")

suppl = Chem.SDMolSupplier(input_file, removeHs=False)

# 取第一个有效构象
mol = None
for m in suppl:
    if m is not None:
        mol = m
        break

if mol:
    # 清理属性（删除PubChem特有的大量属性字段）
    for prop in list(mol.GetPropNames()):
        mol.ClearProp(prop)
    
    # 保存为干净SDF
    writer = Chem.SDWriter(output_file)
    writer.write(mol)
    writer.close()
    print(f"已生成干净文件：{output_file}")
else:
    print("读取失败")
