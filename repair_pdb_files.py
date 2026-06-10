import os
import re
import argparse
from pathlib import Path
import pandas as pd
from Bio.PDB import PDBParser, PDBIO, Select

# 失败受体列表
FAILED_RECEPTORS = ["2flu", "6dvh", "6qcb"]

def parse_args():
    parser = argparse.ArgumentParser(description='Diagnose and repair PDB files')
    parser.add_argument('--input_dir', 
                        default="C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/output/receptors/", 
                        help='Directory containing input PDB files')
    parser.add_argument('--output_dir', 
                        default="C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/output/receptors_repaired/", 
                        help='Directory for repaired PDB files')
    return parser.parse_args()

class StandardSelect(Select):
    """选择标准ATOM记录，去除水分子"""
    def accept_residue(self, residue):
        resname = residue.get_resname()
        if resname == 'HOH':
            return False
        return True

def diagnose_pdb(pdb_path):
    """诊断PDB文件问题"""
    errors = []
    
    try:
        # 尝试使用Bio.PDB解析
        parser = PDBParser(QUIET=False)
        structure = parser.get_structure('structure', str(pdb_path))
        print(f"{pdb_path.name} parsed successfully")
    except Exception as e:
        print(f"{pdb_path.name} parsing failed: {e}")
        errors.append(f"General parsing error: {e}")
    
    # 逐行检查
    line_number = 0
    atom_serial = set()
    model_count = 0
    has_chain_id = False
    
    with open(pdb_path, 'r') as f:
        for line in f:
            line_number += 1
            
            # 检查模型数量
            if line.startswith('MODEL'):
                model_count += 1
            
            # 检查ATOM记录
            if line.startswith('ATOM') or line.startswith('HETATM'):
                # 检查链ID
                chain_id = line[21:22].strip()
                if chain_id:
                    has_chain_id = True
                
                # 检查原子序列号
                try:
                    serial = int(line[6:11].strip())
                    if serial in atom_serial:
                        errors.append(f"Line {line_number}: Duplicate atom serial number {serial}")
                    atom_serial.add(serial)
                except:
                    errors.append(f"Line {line_number}: Invalid atom serial number")
                
                # 检查坐标
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    if x == 0.0 and y == 0.0 and z == 0.0:
                        errors.append(f"Line {line_number}: Zero coordinates")
                except:
                    errors.append(f"Line {line_number}: Invalid coordinates")
                
                # 检查元素符号
                element = line[76:78].strip()
                if not element:
                    errors.append(f"Line {line_number}: Missing element symbol")
    
    # 检查链ID
    if not has_chain_id:
        errors.append("Missing chain ID")
    
    # 检查模型数量
    if model_count > 1:
        errors.append(f"Multiple models ({model_count}), only first will be kept")
    
    return errors

def repair_pdb(input_path, output_path):
    """修复PDB文件"""
    repairs = []
    
    # 读取原始文件
    with open(input_path, 'r') as f:
        lines = f.readlines()
    
    # 修复过程
    new_lines = []
    atom_counter = 1
    model_kept = False
    
    for line in lines:
        # 只保留第一个模型
        if line.startswith('MODEL'):
            if not model_kept:
                new_lines.append(line)
                model_kept = True
            continue
        if line.startswith('ENDMDL'):
            if model_kept:
                new_lines.append(line)
            continue
        
        # 处理ATOM记录
        if line.startswith('ATOM') or line.startswith('HETATM'):
            # 修复链ID
            if line[21:22].strip() == '':
                new_line = line[:21] + 'A' + line[22:]
                repairs.append("Added chain ID 'A'")
            else:
                new_line = line
            
            # 修复原子序列号
            new_line = new_line[:6] + f"{atom_counter:5d}" + new_line[11:]
            atom_counter += 1
            
            # 检查并修复坐标
            try:
                x = float(new_line[30:38].strip())
                y = float(new_line[38:46].strip())
                z = float(new_line[46:54].strip())
                if x == 0.0 and y == 0.0 and z == 0.0:
                    repairs.append(f"Removed atom with zero coordinates")
                    continue
            except:
                repairs.append(f"Removed atom with invalid coordinates")
                continue
            
            # 修复B因子
            try:
                b_factor = float(new_line[60:66].strip())
            except:
                new_line = new_line[:60] + "  1.00" + new_line[66:]
                repairs.append("Fixed B-factor")
            
            # 确保元素符号
            element = new_line[76:78].strip()
            if not element:
                # 从原子名推断元素符号
                atom_name = new_line[12:16].strip()
                if atom_name.startswith('C'):
                    element = 'C'
                elif atom_name.startswith('N'):
                    element = 'N'
                elif atom_name.startswith('O'):
                    element = 'O'
                elif atom_name.startswith('S'):
                    element = 'S'
                else:
                    element = 'C'
                new_line = new_line[:76] + element.rjust(2) + new_line[78:]
                repairs.append("Added element symbol")
            
            new_lines.append(new_line)
        
        # 保留其他重要记录
        elif line.startswith('HEADER') or line.startswith('TITLE') or line.startswith('COMPND') or \
             line.startswith('AUTHOR') or line.startswith('RESOLUTION') or line.startswith('REMARK') or \
             line.startswith('ATOM') or line.startswith('HETATM') or line.startswith('CONECT') or \
             line.startswith('END'):
            new_lines.append(line)
    
    # 保存修复后的文件
    with open(output_path, 'w') as f:
        f.writelines(new_lines)
    
    # 使用Bio.PDB进行最终清洗
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('structure', str(output_path))
        
        io = PDBIO()
        io.set_structure(structure)
        io.save(str(output_path), StandardSelect())
        
        repairs.append("Removed water molecules")
    except Exception as e:
        repairs.append(f"Bio.PDB cleaning failed: {e}")
    
    return repairs

def verify_pdb(pdb_path):
    """验证修复后的PDB文件"""
    try:
        from rdkit import Chem
        mol = Chem.MolFromPDBFile(str(pdb_path), removeHs=False)
        if mol:
            return True, "RDKit load successful"
        else:
            return False, "RDKit load failed"
    except Exception as e:
        return False, f"RDKit error: {e}"

def main():
    args = parse_args()
    
    # 确保输出目录存在
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    
    # 生成诊断报告
    diagnostic_report = []
    repair_report = []
    
    for pdb_id in FAILED_RECEPTORS:
        input_path = Path(args.input_dir) / f"{pdb_id}.pdb"
        output_path = output_dir / f"{pdb_id}_repaired.pdb"
        
        print(f"\nProcessing {pdb_id}...")
        
        # 诊断
        errors = diagnose_pdb(input_path)
        diagnostic_report.append({"PDB_ID": pdb_id, "Errors": errors})
        
        # 修复
        repairs = repair_pdb(input_path, output_path)
        repair_report.append({"PDB_ID": pdb_id, "Repairs": repairs})
        
        # 验证
        valid, message = verify_pdb(output_path)
        print(f"Verification: {message}")
    
    # 保存诊断报告
    diagnostic_file = output_dir / "diagnostic_report.txt"
    with open(diagnostic_file, 'w') as f:
        f.write("PDB Diagnostic Report\n")
        f.write("====================\n\n")
        for entry in diagnostic_report:
            f.write(f"PDB: {entry['PDB_ID']}\n")
            if entry['Errors']:
                for error in entry['Errors']:
                    f.write(f"  - {error}\n")
            else:
                f.write("  No errors found\n")
            f.write("\n")
    
    # 保存修复报告
    repair_df = pd.DataFrame(repair_report)
    repair_df['Repair_Description'] = repair_df['Repairs'].apply(lambda x: ', '.join(x) if x else 'No repairs needed')
    repair_csv = output_dir / "repair_report.csv"
    repair_df[['PDB_ID', 'Repair_Description']].to_csv(repair_csv, index=False)
    
    print(f"\nDiagnostic report saved to: {diagnostic_file}")
    print(f"Repair report saved to: {repair_csv}")
    print("\nRepair process completed!")

if __name__ == "__main__":
    main()