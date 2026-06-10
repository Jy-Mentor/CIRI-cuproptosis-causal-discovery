import os
import urllib.request
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

# 蛋白质列表
proteins = [
    "TP53", "IL1B", "IL6", "TNF", "STAT3", "BCL2", "NFKB1", "PTGS2", "TLR4", "SRC",
    "STAT1", "RELA", "ICAM1", "CCL2", "CCL5", "CASP8", "VCAM1", "TGFB1", "PTPRC", "IKBKB",
    "STAT5A", "CCND1", "HMOX1", "TIMP1", "NLRP3", "CDK4", "PARP1", "CCR5", "FAS", "MAPK9",
    "NFE2L2", "SREBF1", "IRF1", "IL10RA", "CXCR3", "PGR", "BID", "EGR1", "F3", "AIF1",
    "CTSS", "PTGS1", "IRAK4", "LYN", "SREBF2", "TOP2A", "GFAP", "CCNA2", "PTGES", "PTPN2",
    "ERBB4", "CTSD", "CTSB", "C3", "SQLE", "HMGCR", "LSS", "CYP51A1"
]

# 保存目录
output_dir = "pdb_files_batch"
os.makedirs(output_dir, exist_ok=True)

# 已知的PDB ID映射（用于一些常见蛋白质）
pdb_id_mapping = {
    "TP53": "1TUP",
    "IL1B": "1ITB",
    "IL6": "1ALU",
    "TNF": "2AZ5",
    "STAT3": "1BG1",
    "BCL2": "1G5M",
    "NFKB1": "1VKB",
    "PTGS2": "5F19",
    "TLR4": "3FXI",
    "SRC": "2SRC",
    "STAT1": "1YVL",
    "RELA": "1RAM",
    "ICAM1": "1IC1",
    "CCL2": "2M2N",
    "CCL5": "2M2N",  # 与CCL2相同结构
    "CASP8": "1QTN",
    "VCAM1": "1G9V",
    "TGFB1": "1KLC",
    "PTPRC": "1YFO",
    "IKBKB": "4KIK",
    "STAT5A": "1Y1U",
    "CCND1": "2GDK",
    "HMOX1": "1N45",
    "TIMP1": "1UEA",
    "NLRP3": "6NPY",
    "CDK4": "2W96",
    "PARP1": "4DQY",
    "CCR5": "4MBS",
    "FAS": "3EZQ",
    "MAPK9": "1EEJ",
    "NFE2L2": "4CBT",
    "SREBF1": "2Q21",
    "IRF1": "1T2K",
    "IL10RA": "1ILK",
    "CXCR3": "4RWS",
    "PGR": "3U4E",
    "BID": "2BID",
    "EGR1": "1F2E",
    "F3": "1DQV",
    "AIF1": "1SRL",
    "CTSS": "1PPJ",
    "PTGS1": "3PGH",
    "IRAK4": "2NRU",
    "LYN": "1FGN",
    "SREBF2": "2QMJ",
    "TOP2A": "1ZXM",
    "GFAP": "3KLT",
    "CCNA2": "1E32",
    "PTGES": "1CI4",
    "PTPN2": "2HNP",
    "ERBB4": "3U7S",
    "CTSD": "1IDG",
    "CTSB": "1CSB",
    "C3": "2I6Q",
    "SQLE": "4P4Q",
    "HMGCR": "1HWK",
    "LSS": "5N8E",
    "CYP51A1": "3LD6"
}

def download_pdb(protein_name):
    """下载单个蛋白质的PDB文件"""
    try:
        # 检查是否有已知的PDB ID映射
        if protein_name in pdb_id_mapping:
            pdb_id = pdb_id_mapping[protein_name]
            print(f"正在下载 {protein_name} (PDB ID: {pdb_id})...")
        else:
            # 尝试直接使用蛋白质名称搜索（可能不成功）
            pdb_id = protein_name
            print(f"正在尝试下载 {protein_name}...")
        
        # 构建PDB文件URL
        url = f"https://files.rcsb.org/view/{pdb_id}.pdb"
        output_file = os.path.join(output_dir, f"{protein_name}_{pdb_id}.pdb")
        
        # 下载文件
        urllib.request.urlretrieve(url, output_file)
        print(f"成功下载: {protein_name} -> {output_file}")
        return (protein_name, True, output_file)
    except Exception as e:
        print(f"下载失败: {protein_name} - {str(e)}")
        return (protein_name, False, str(e))

if __name__ == "__main__":
    print(f"开始下载 {len(proteins)} 个蛋白质的PDB文件...")
    print(f"保存目录: {output_dir}")
    print("=" * 80)
    
    # 使用线程池并行下载
    results = []
    with ThreadPoolExecutor(max_workers=5) as executor:
        # 提交所有下载任务
        future_to_protein = {executor.submit(download_pdb, protein): protein for protein in proteins}
        
        # 收集结果
        for future in as_completed(future_to_protein):
            protein = future_to_protein[future]
            try:
                result = future.result()
                results.append(result)
            except Exception as e:
                print(f"处理 {protein} 时出错: {str(e)}")
            
            # 添加小延迟，避免请求过于频繁
            time.sleep(0.5)
    
    print("=" * 80)
    print("下载完成！")
    
    # 统计结果
    success_count = sum(1 for r in results if r[1])
    failed_count = len(results) - success_count
    
    print(f"成功: {success_count} 个")
    print(f"失败: {failed_count} 个")
    
    if failed_count > 0:
        print("\n失败的蛋白质:")
        for r in results:
            if not r[1]:
                print(f"- {r[0]}: {r[2]}")
    
    print("\n注意: 对于失败的蛋白质，可能需要通过其他方式获取其PDB结构，")
    print("例如通过UniProt ID映射或使用AlphaFold预测结构。")
