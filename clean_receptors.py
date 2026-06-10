import os
from pathlib import Path
from rdkit import Chem

BASE_DIR = Path(r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\output")
RECEPTOR_DIR = BASE_DIR / "receptors"

# 受体列表
receptors = ["1zmc", "3o3u", "3p1m", "3rzf", "5n2s", "8vmc", "9soq"]

for pdb_id in receptors:
    print(f"\n{'='*50}")
    print(f"清洗受体: {pdb_id}")
    print(f"{'='*50}")
    
    # 输入文件
    input_file = RECEPTOR_DIR / f"{pdb_id}.pdb"
    if not input_file.exists():
        print(f"❌ 找不到输入文件")
        continue
    
    # 输出文件
    output_file = RECEPTOR_DIR / f"{pdb_id}_clean.pdb"
    
    try:
        # 加载PDB文件
        mol = Chem.MolFromPDBFile(str(input_file), removeHs=False)
        if not mol:
            print(f"❌ 无法加载PDB文件")
            continue
        
        print(f"原始原子数: {mol.GetNumAtoms()}")
        
        # 保留蛋白质原子（标准氨基酸）
        standard_amino_acids = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'}
        
        # 创建新分子，只包含标准氨基酸
        clean_mol = Chem.Mol()
        builder = Chem.RWMol(clean_mol)
        
        for atom in mol.GetAtoms():
            resinfo = atom.GetPDBResidueInfo()
            if resinfo:
                resname = resinfo.GetResidueName()
                if resname in standard_amino_acids:
                    builder.AddAtom(atom)
        
        # 转换回Mol对象
        clean_mol = builder.GetMol()
        
        print(f"清洗后原子数: {clean_mol.GetNumAtoms()}")
        
        # 保存清洗后的文件
        Chem.MolToPDBFile(clean_mol, str(output_file))
        print(f"✅ 清洗完成，保存到: {output_file}")
        
    except Exception as e:
        print(f"❌ 清洗时出错: {e}")
        continue

print("\n所有受体清洗完成！")