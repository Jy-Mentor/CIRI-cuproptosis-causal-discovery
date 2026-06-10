from pathlib import Path
from collections import Counter

DOCKING_DIR = Path(r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\output\docking_results")

for pdb_id in ["1zmc", "3o3u", "3rzf", "5n2s", "8vmc"]:  # 之前失败的5个
    print(f"\n{pdb_id}:")
    
    # 找文件
    file_path = None
    for ext in ["_complex.pdb", "_docked.pdb"]:
        if (DOCKING_DIR / f"{pdb_id}{ext}").exists():
            file_path = DOCKING_DIR / f"{pdb_id}{ext}"
            break
    
    if not file_path:
        print("  文件不存在")
        continue
    
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    # 统计所有非氨基酸、非水残基
    common_aa = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
                 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
                 'THR', 'TRP', 'TYR', 'VAL'}
    
    resnames = []
    for line in lines:
        if line.startswith(('ATOM', 'HETATM')):
            resname = line[17:20].strip()
            if resname not in common_aa and resname not in ['HOH', 'WAT']:
                resnames.append(resname)
    
    if resnames:
        print(f"  检测到非氨基酸残基: {Counter(resnames)}")
    else:
        print("  未检测到任何非氨基酸残基（配体可能被标记为氨基酸名或缺失）")