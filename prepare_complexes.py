import os
import argparse
from pathlib import Path

# 成功受体列表
SUCCESSFUL_RECEPTORS = ["1zmc", "3o3u", "3p1m", "3rzf", "5n2s", "8vmc", "9soq"]

def parse_args():
    parser = argparse.ArgumentParser(description='Prepare complex files by combining receptor and ligand')
    parser.add_argument('--receptor_dir', 
                        default="C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/output/receptors/", 
                        help='Directory containing receptor PDB files')
    parser.add_argument('--docking_dir', 
                        default="C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/output/docking_results/", 
                        help='Directory containing docking results')
    parser.add_argument('--output_dir', 
                        default="C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/output/complexes/", 
                        help='Output directory for combined complex files')
    return parser.parse_args()

def pdbqt_to_pdb(pdbqt_path, pdb_path):
    """将 PDBQT 文件转换为 PDB 文件"""
    with open(pdbqt_path, 'r') as f:
        lines = f.readlines()
    
    pdb_lines = []
    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            # 截断到 66 列，移除电荷信息
            pdb_line = line[:66] + '\n'
            pdb_lines.append(pdb_line)
    
    with open(pdb_path, 'w') as f:
        f.writelines(pdb_lines)
    
    return len(pdb_lines) > 0

def extract_ligand_from_complex(complex_pdb, output_ligand_pdb):
    """从复合物文件中提取配体"""
    with open(complex_pdb, 'r') as f:
        lines = f.readlines()
    
    ligand_lines = []
    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            resname = line[17:20].strip()
            if resname in ['UNL', 'LIG', 'UNK', 'BCP']:
                ligand_lines.append(line)
    
    with open(output_ligand_pdb, 'w') as f:
        f.writelines(ligand_lines)
    
    return len(ligand_lines) > 0

def combine_receptor_ligand(receptor_pdb, ligand_pdb, output_complex_pdb):
    """将受体和配体合并为复合物文件"""
    # 读取受体
    with open(receptor_pdb, 'r') as f:
        receptor_lines = f.readlines()
    
    # 读取配体
    with open(ligand_pdb, 'r') as f:
        ligand_lines = f.readlines()
    
    # 合并
    combined_lines = receptor_lines + ligand_lines
    
    # 保存
    with open(output_complex_pdb, 'w') as f:
        f.writelines(combined_lines)
    
    return len(combined_lines) > 0

def main():
    args = parse_args()
    
    # 确保输出目录存在
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    
    for pdb_id in SUCCESSFUL_RECEPTORS:
        print(f"Processing {pdb_id}...")
        
        # 1. 找到受体文件
        receptor_pdb = Path(args.receptor_dir) / f"{pdb_id}.pdb"
        if not receptor_pdb.exists():
            print(f"Receptor file not found: {receptor_pdb}")
            continue
        
        # 2. 找到配体文件（复合物文件）
        complex_pdb = Path(args.docking_dir) / f"{pdb_id}_complex.pdb"
        if not complex_pdb.exists():
            print(f"Complex file not found: {complex_pdb}")
            continue
        
        # 3. 提取配体
        temp_ligand = output_dir / f"{pdb_id}_ligand.pdb"
        if not extract_ligand_from_complex(complex_pdb, temp_ligand):
            print(f"Failed to extract ligand from {complex_pdb}")
            continue
        
        # 4. 合并受体和配体
        output_complex = output_dir / f"{pdb_id}_full_complex.pdb"
        if combine_receptor_ligand(receptor_pdb, temp_ligand, output_complex):
            print(f"Created full complex: {output_complex}")
        else:
            print(f"Failed to create full complex for {pdb_id}")
        
        # 清理临时文件
        if temp_ligand.exists():
            temp_ligand.unlink()
    
    print("Processing completed!")

if __name__ == "__main__":
    main()