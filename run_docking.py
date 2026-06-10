import os
import subprocess
import argparse
import sys
from pathlib import Path

# 添加 utils 目录到路径
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from utils.config_loader import get_config

# 获取配置
config = get_config()

# 受体列表
RECEPTORS = ["1zmc", "2flu", "3o3u", "3p1m", "3rzf", "5n2s", "6dvh", "6qcb", "8vmc", "9soq"]

def parse_args():
    parser = argparse.ArgumentParser(description='Run molecular docking with AutoDock Vina')
    parser.add_argument('--receptor_dir', 
                        default="C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/output/receptors/", 
                        help='Directory containing receptor PDB files')
    parser.add_argument('--ligand_path', 
                        default="C:/Users/Jy-Mentor-7/AppData/Local/Temp/BCP.pdbqt", 
                        help='Path to ligand PDBQT file')
    parser.add_argument('--output_dir', 
                        default="C:/Users/Jy-Mentor-7/AppData/Local/Temp/docking_results/", 
                        help='Output directory for docking results')
    parser.add_argument('--vina_path', 
                        default=config.get_vina_exe(), 
                        help='Path to AutoDock Vina executable')
    return parser.parse_args()

def prepare_receptor(pdb_path, output_pdbqt):
    """使用prepare_receptor4.py准备受体"""
    try:
        # 检查是否存在prepare_receptor4.py
        prepare_script = config.get_mgltools_prepare_receptor()
        if not os.path.exists(prepare_script):
            print("prepare_receptor4.py not found at: " + prepare_script)
            return False
        
        # 运行prepare_receptor4.py
        cmd = [
            config.get_mgltools_python(),
            prepare_script,
            "-r", str(pdb_path),
            "-o", str(output_pdbqt),
            "-A", "hydrogens"
        ]
        
        print("Preparing receptor: " + str(pdb_path))
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode == 0:
            print("Receptor prepared successfully: " + str(output_pdbqt))
            return True
        else:
            print("Error preparing receptor: " + result.stderr)
            return False
    except Exception as e:
        print("Exception preparing receptor: " + str(e))
        return False

def run_docking(vina_path, receptor_pdbqt, ligand_pdbqt, output_pdb, log_file):
    """运行AutoDock Vina对接"""
    try:
        # 设置对接参数
        cmd = [
            vina_path,
            "--receptor", str(receptor_pdbqt),
            "--ligand", str(ligand_pdbqt),
            "--out", str(output_pdb),
            "--log", str(log_file),
            "--center_x", "0",
            "--center_y", "0",
            "--center_z", "0",
            "--size_x", "20",
            "--size_y", "20",
            "--size_z", "20",
            "--num_modes", "9",
            "--energy_range", "4"
        ]
        
        print(f"Running docking for: {os.path.basename(receptor_pdbqt)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode == 0:
            print(f"Docking completed successfully: {output_pdb}")
            return True
        else:
            print(f"Error running docking: {result.stderr}")
            return False
    except Exception as e:
        print(f"Exception running docking: {e}")
        return False

def main():
    args = parse_args()
    
    # 确保输出目录存在
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    
    # 对每个受体进行对接
    for pdb_id in RECEPTORS:
        print(f"\nProcessing {pdb_id}...")
        
        # 构建路径
        receptor_pdb = Path(args.receptor_dir) / f"{pdb_id}.pdb"
        receptor_pdbqt = output_dir / f"{pdb_id}_receptor.pdbqt"
        output_pdb = output_dir / f"{pdb_id}_complex.pdb"
        log_file = output_dir / f"{pdb_id}_docking.log"
        
        # 检查受体文件是否存在
        if not receptor_pdb.exists():
            print(f"Receptor file not found: {receptor_pdb}")
            continue
        
        # 准备受体
        if not prepare_receptor(receptor_pdb, receptor_pdbqt):
            continue
        
        # 运行对接
        if not run_docking(args.vina_path, receptor_pdbqt, args.ligand_path, output_pdb, log_file):
            continue
    
    print("\nDocking process completed!")

if __name__ == "__main__":
    main()