import requests
import os
import time

# 创建保存PDB文件的目录
pdb_dir = "pdb_files"
os.makedirs(pdb_dir, exist_ok=True)

# 蛋白到PDB ID的映射（使用已知的PDB ID）
protein_pdb_map = {
    "NR3C1": "1R27",  # 糖皮质激素受体
    "STAT3": "1BG1",  # STAT3 DNA结合域
    "MAPK14": "2ERK",  # p38 MAP激酶
    "SMARCA4": "5Z3U",  # BRG1 ATPase域
    "CASP9": "1NW9",  # 半胱天冬酶-9
    "PTPRC": "1YFO",  # CD45磷酸酶
    "MDM2": "1T4E",  # MDM2-p53相互作用
    "CTSB": "1CSB",  # 组织蛋白酶B
    "ARG1": "1N11",  # 精氨酸酶1
    "NFE2L2": "4CBT",  # Nrf2
    "ESR1": "1ERE",  # 雌激素受体α
    "APP": "1Z7S",  # β-淀粉样前体蛋白
    "CDC42": "1AN0",  # CDC42 GTPase
    "STAT5A": "1Y1U",  # STAT5A
    "PTGS2": "5F19",  # COX-2
    "NFKB1": "1VKB",  # NF-κB p50
    "JAK2": "4Z16",  # JAK2激酶
    "LYN": "1FGN",  # LYN激酶
    "PPARD": "3DZU",  # PPARδ
    "NOS3": "1NQD",  # eNOS
    "HMGCR": "1HWK",  # 羟甲基戊二酰辅酶A还原酶
    "IL1B": "1ITB",  # 白细胞介素-1β
    "PPARA": "2P54",  # PPARα
    "MMP9": "1GKC",  # 基质金属蛋白酶-9
    "PPARG": "2PRG",  # PPARγ
    "IDO1": "2D0T",  # 吲哚胺2,3-双加氧酶1
    "HMOX1": "1N45",  # 血红素加氧酶-1
    "CASP8": "1QTN",  # 半胱天冬酶-8
    "PTGS1": "3N8Y",  # COX-1
    "EGR1": "1X4Y"   # 早期生长反应蛋白1
}

# PDB数据库API URL
pdb_api_url = "https://files.rcsb.org/download/"

print(f"开始下载 {len(protein_pdb_map)} 个蛋白的PDB文件...")
print(f"PDB文件将保存到: {os.path.abspath(pdb_dir)}")
print()

success_count = 0
fail_count = 0
failed_proteins = []

for protein, pdb_id in protein_pdb_map.items():
    print(f"正在下载 {protein} (PDB ID: {pdb_id})...")
    
    try:
        # 下载PDB文件
        pdb_file_url = f"{pdb_api_url}{pdb_id}.pdb"
        pdb_response = requests.get(pdb_file_url, timeout=10)
        pdb_response.raise_for_status()
        
        # 保存PDB文件
        pdb_file_path = os.path.join(pdb_dir, f"{protein}_{pdb_id}.pdb")
        with open(pdb_file_path, "wb") as f:
            f.write(pdb_response.content)
        
        print(f"成功下载 {protein} 的PDB文件: {pdb_file_path}")
        success_count += 1
        
    except Exception as e:
        print(f"下载 {protein} 时出错: {e}")
        fail_count += 1
        failed_proteins.append(protein)
    
    # 避免请求过于频繁
    time.sleep(1)
    print()

print("\n下载完成！")
print(f"成功下载: {success_count} 个蛋白")
print(f"下载失败: {fail_count} 个蛋白")

if failed_proteins:
    print("\n下载失败的蛋白:")
    for protein in failed_proteins:
        print(f"- {protein}")

print(f"\nPDB文件已保存到: {os.path.abspath(pdb_dir)}")
