import os
from pathlib import Path
from collections import Counter

BASE_DIR = Path(r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\output")
DOCKING_DIR = BASE_DIR / "docking_results"

pdb_ids = ["1zmc", "3o3u", "3rzf", "5n2s", "8vmc"]  # 出问题的5个

for pdb_id in pdb_ids:
    print(f"\n{'='*50}")
    print(f"诊断: {pdb_id}")
    print(f"{'='*50}")
    
    # 找文件
    complex_file = None
    for ext in ["_complex.pdb", "_docked.pdb"]:
        candidate = DOCKING_DIR / f"{pdb_id}{ext}"
        if candidate.exists():
            complex_file = candidate
            break
    
    if not complex_file:
        print(f"❌ 找不到文件")
        continue
    
    print(f"📁 文件: {complex_file}")
    
    # 检查内容
    with open(complex_file, 'r') as f:
        lines = f.readlines()
    
    # 统计残基名
    het_resnames = [line[17:20].strip() for line in lines if line.startswith('HETATM')]
    model_count = sum(1 for line in lines if line.startswith('MODEL'))
    
    print(f"📊 MODEL数量: {model_count}")
    print(f"📊 HETATM残基名统计: {Counter(het_resnames)}")
    
    # 检查配体坐标（取前5个原子）
    print(f"📍 配体前5个原子坐标:")
    count = 0
    for line in lines:
        if line.startswith('HETATM') or (line.startswith('ATOM') and line[17:20].strip() in het_resnames):
            x, y, z = line[30:38], line[38:46], line[46:54]
            print(f"   {line[17:20].strip()} {line[12:16].strip()}: ({x}, {y}, {z})")
            count += 1
            if count >= 5:
                break
    
    # 检查受体坐标（用于对比）
    print(f"📍 受体前5个原子坐标:")
    count = 0
    for line in lines:
        if line.startswith('ATOM') and line[17:20].strip() not in het_resnames:
            x, y, z = line[30:38], line[38:46], line[46:54]
            print(f"   {line[17:20].strip()} {line[12:16].strip()}: ({x}, {y}, {z})")
            count += 1
            if count >= 5:
                break