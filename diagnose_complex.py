#!/usr/bin/env python3
"""
诊断复合物文件中配体是否在受体口袋里
"""

from pymol import cmd
from pathlib import Path
import math

def diagnose_complex():
    """诊断复合物文件"""
    BASE_DIR = Path(r"C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/output/docking_results")

    for pdb_id in ["1zmc", "3o3u", "3rzf", "5n2s", "8vmc"]:
        print(f"\n{'='*40}")
        print(f"诊断: {pdb_id}")
        
        complex_file = BASE_DIR / f"{pdb_id}_complex.pdb"
        if not complex_file.exists():
            print("❌ 文件不存在")
            continue
        
        cmd.reinitialize()
        cmd.load(str(complex_file), "complex")
        
        # 检测配体
        resnames = []
        cmd.iterate("complex", "resnames.append(resn)", space={'resnames': resnames})
        unique_resnames = list(set(resnames))
        print(f"所有残基名: {unique_resnames}")
        
        # 分离假设配体是 UNL
        cmd.create("receptor", "complex and not resn UNL")
        cmd.create("ligand", "complex and resn UNL")
        
        # 计算原子数量
        rec_atoms = cmd.count_atoms("receptor")
        lig_atoms = cmd.count_atoms("ligand")
        print(f"受体原子数: {rec_atoms}")
        print(f"配体原子数: {lig_atoms}")
        
        if rec_atoms == 0:
            print("❌ 错误：复合物文件中没有受体原子！")
            print("   这解释了为什么相互作用为0 - 因为只有配体，没有受体")
        else:
            # 计算质心距离
            rec_center = cmd.centerofmass("receptor")
            lig_center = cmd.centerofmass("ligand")
            
            distance = math.sqrt(sum([(a-b)**2 for a, b in zip(rec_center, lig_center)]))
            print(f"受体-配体质心距离: {distance:.2f} Å")
            
            if distance > 20:
                print("⚠️ 警告：距离太远！配体不在口袋里")
            elif distance < 10:
                print("✅ 距离正常，应在口袋内")
            else:
                print("⚠️ 距离中等，可能在口袋边缘")

if __name__ == "__main__":
    diagnose_complex()
