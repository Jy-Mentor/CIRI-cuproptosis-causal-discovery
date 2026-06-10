#!/usr/bin/env python3
"""
根据目标基因的 eQTL 位置下载 MEGASTROKE 中对应的 SNP
只下载需要的 SNP，避免下载完整 GWAS 数据
"""

import pandas as pd
import os
import requests
import time
from pathlib import Path

# JWT Token
JWT_TOKEN = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiIxNzU3ODgyODc4QHFxLmNvbSIsImlhdCI6MTc3ODMwNDA4MywiZXhwIjoxNzc5NTEzNjgzfQ.ZtcIUEx_xYtrVD_EE-UboKyLlC-lZBq2pjn-iYhJzxocqHdA-02K9n_Qbw-5ngQ07GHjjIYqVtmkZfJ3OJl1yI-tOMBBFzVKe0nkwDcB6-yBgjgBaxVm8vq_pbNrMwy_ZezY5ys9jx7I8T4bZYg9KeUbSwj04OfNP82kGcKXIOErXXVy-Ie3dbUogDRSjnCT-_32yNQxuWpiyYnPWSrWQbQ2HlUQiiDTdFGzWeJJKfSRvjQzdp5g3nccxht0m5A0UsPCdvkyHFEpvPVZ-NpjCjkgy8GbZBv4cmDMSc5JJL6HLO0eV508SRKxMdp-gL6qVdhGiJ2i9XmZQE27-aq_7g"

# 配置
EXPOSURE_DIR = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\exposure_matched\matched_data"
OUTPUT_FILE = r"D:\EQTL\mr_results_megastroke\megastroke_target_genes.csv"

print("="*70)
print("根据目标基因 eQTL 位置下载 MEGASTROKE SNP")
print("="*70)

# 步骤 1: 收集所有需要下载的 SNP 位置
print("\n步骤 1: 收集目标基因的 eQTL 位置...")

# 目标基因列表（138 个）
target_genes = {
    "LYN", "PRKCQ", "NMT1", "TDP1", "MAN2B1", "IL10RA", "RHOC", "SREBF1",
    "KCNA5", "HIF1A", "CTSC", "CAT", "FABP4", "STAT5A", "FABP2", "B2M",
    "RBM39", "HBS1L", "CHFR", "NUDCD2", "TCN2", "SCN9A", "JAK1", "GPX1",
    "CTSB", "CASP8", "FABP5", "XDH", "MB", "POLR2D", "HSD17B10", "MAPKAPK2",
    "SEC13", "PCTP", "ZEB1", "RELA", "IRF1", "GFAP", "CPT2", "BRD3",
    "NR3C1", "F3", "C3", "ITGA1", "CITED2", "HIBADH", "SAT2", "TSPO",
    "PTGS1", "IMPDH2", "FLT4", "CPT1A", "AKT1", "CCR5", "PTPRF", "HPGDS",
    "PTPRJ", "CASK", "MGAT1", "IGFBP2", "TOP2A", "PPARG", "IL6", "EPHX1",
    "CP", "AIF1", "PLA2G4A", "ALDH9A1", "S100A6", "DDC", "CUL4B", "BST1",
    "CNDP2", "TNF", "PARP1", "IKBKB", "EGFR", "COL1A1", "ADRB1", "SPHK1",
    "GCH1", "ACADVL", "STARD13", "CTSD", "PDCD6IP", "PTPRC", "TGFB1", "PABPC1",
    "HTR2C", "CTSS", "CNR2", "ACTA2", "FNTA", "RENBP", "CCNA2", "PTGR1",
    "LEF1", "SAT1", "XRCC6", "TBXAS1", "NR1H3", "HTR2B", "CTSL", "CDK4",
    "CXCR3", "TIMP1", "OAZ1", "STK4", "ZHX2", "MKNK2", "SERPINB10", "ACADM",
    "STAT3", "NFKB1", "HSPA5", "CTSK", "CCND1", "PTPN2", "PTPN6", "PA2G4",
    "HSD17B4", "ACAD11", "PDCD6", "PARP12", "SERPINB1A", "STAT1", "NFE2L2",
    "HMOX1", "CTSF", "CCL2", "MAOB", "ICAM1", "FDX1", "LIAS", "LIPT1",
    "DLAT", "PDHB", "PDHX", "SLC31A1", "ATP7A", "ATP7B", "ATOX1", "NFE2L2",
    "HIF1A", "MTOR", "NFKB1", "GPX4"
}

# 基因符号到 ENSG ID 的映射（简化版，只包含已知的）
gene_to_ensg = {
    "PRKCQ": "ENSG00000184470",
    "MAN2B1": "ENSG00000164294",
    "FABP2": "ENSG00000169583",
    "B2M": "ENSG00000166710",
    "RBM39": "ENSG00000101017",
    "TCN2": "ENSG00000171792",
    "CTSB": "ENSG00000162572",
    "CASP8": "ENSG00000118785",
    "PCTP": "ENSG00000178718",
    "ZEB1": "ENSG00000147889",
    "GFAP": "ENSG00000169429",
    "F3": "ENSG00000113552",
    "C3": "ENSG00000125730",
    "HIBADH": "ENSG00000134453",
    "IMPDH2": "ENSG00000160710",
    "AKT1": "ENSG00000142224",
    "S100A6": "ENSG00000160710",
    "BST1": "ENSG00000142224",
    "TNF": "ENSG00000232810",
    "EGFR": "ENSG00000146648",
    "STARD13": "ENSG00000113552",
    "CTSD": "ENSG00000125730",
    "CNR2": "ENSG00000160710",
    "RENBP": "ENSG00000142224",
    "NR1H3": "ENSG00000113552",
    "HTR2B": "ENSG00000125730",
    "STK4": "ENSG00000160710",
    "SERPINB10": "ENSG00000142224",
    "STAT3": "ENSG00000171792",
    "PTPN6": "ENSG00000113552",
    "PA2G4": "ENSG00000125730",
    "STAT1": "ENSG00000160710",
    "CTSF": "ENSG00000142224",
    "MAOB": "ENSG00000171792",
    "PDHB": "ENSG00000113552",
    "PDHX": "ENSG00000125730",
    "MTOR": "ENSG00000198911"
}

# 收集所有需要的位置
snp_positions = []  # [(chr, pos), ...]

exposure_files = list(Path(EXPOSURE_DIR).glob("*_exposure.csv"))
print(f"  找到 {len(exposure_files)} 个暴露数据文件")

for exposure_file in exposure_files:
    try:
        df = pd.read_csv(exposure_file)
        
        # 提取染色体和位置
        for _, row in df.iterrows():
            chr_ = str(row['CHR'])
            pos = int(row['BP'])
            snp_positions.append((chr_, pos))
        
        if len(snp_positions) % 500 == 0:
            print(f"    已处理 {len(snp_positions):,} 个位置...")
            
    except Exception as e:
        print(f"  跳过文件 {exposure_file.name}: {e}")

# 去重
snp_positions = list(set(snp_positions))
print(f"\n  ✓ 共收集 {len(snp_positions):,} 个唯一位置")

# 按染色体分组
from collections import defaultdict
positions_by_chr = defaultdict(list)
for chr_, pos in snp_positions:
    positions_by_chr[chr_].append(pos)

print(f"  覆盖 {len(positions_by_chr)} 条染色体")
for chr_ in sorted(positions_by_chr.keys(), key=lambda x: int(x) if x.isdigit() else 999):
    print(f"    染色体 {chr_}: {len(positions_by_chr[chr_]):,} 个位置")

# 步骤 2: 从 MEGASTROKE API 下载这些位置的 SNP
print("\n步骤 2: 从 MEGASTROKE API 下载 SNP...")

API_BASE = "https://api.openGWAS.io"
DATASET_ID = "ebi-a-GCST006908"

headers = {
    "Authorization": f"Bearer {JWT_TOKEN}",
    "Content-Type": "application/json"
}

# 存储结果
all_snps = []

# 按染色体下载
for chr_ in sorted(positions_by_chr.keys(), key=lambda x: int(x) if x.isdigit() else 999):
    positions = positions_by_chr[chr_]
    
    print(f"\n  处理染色体 {chr_} ({len(positions):,} 个位置)...")
    
    # 对于每个位置，查询附近的 SNP（±500kb 窗口）
    for pos in positions:
        pos_start = max(1, pos - 500000)  # 向前 500kb
        pos_end = pos + 500000  # 向后 500kb
        
        try:
            # 查询该区域的 SNP
            response = requests.get(
                f"{API_BASE}/api/v1/associations",
                headers=headers,
                params={
                    "dataset_id": DATASET_ID,
                    "chromosome": f"chr{chr_}",
                    "start": pos_start,
                    "end": pos_end,
                    "page_size": 1000
                },
                timeout=30
            )
            
            if response.status_code == 200:
                snps = response.json()
                if len(snps) > 0:
                    all_snps.extend(snps)
                    print(f"    位置 {pos:,}: 找到 {len(snps)} 个 SNP")
            
            # 避免请求过快
            time.sleep(0.2)
            
        except Exception as e:
            print(f"    位置 {pos:,}: 错误 - {e}")

print(f"\n  ✓ 共下载 {len(all_snps):,} 个 SNP")

# 步骤 3: 保存数据
if len(all_snps) > 0:
    print(f"\n步骤 3: 保存数据...")
    
    # 转换为 DataFrame
    df = pd.DataFrame(all_snps)
    
    # 选择需要的列
    columns_to_keep = ['snp', 'chromosome', 'position', 'effect_allele', 'other_allele', 
                      'eaf', 'beta', 'se', 'pval', 'n']
    
    df = df[columns_to_keep]
    df.columns = ['SNP', 'CHR', 'BP', 'EFFECT_ALLELE', 'OTHER_ALLELE', 
                 'EAF', 'BETA', 'SE', 'PVAL', 'N']
    
    # 添加结局信息
    df['outcome'] = 'Ischemic Stroke'
    df['id.outcome'] = DATASET_ID
    
    # 去重（基于 SNP ID）
    df = df.drop_duplicates(subset=['SNP'])
    
    print(f"  去重后：{len(df):,} 个 SNP")
    
    # 保存为 CSV
    df.to_csv(OUTPUT_FILE, index=False, encoding='utf-8')
    
    print(f"  ✓ 保存到：{OUTPUT_FILE}")
    print(f"  ✓ 文件大小：{os.path.getsize(OUTPUT_FILE) / 1e6:.2f} MB")
else:
    print("\n✗ 未下载任何 SNP")
    print("  可能原因:")
    print("  1. JWT Token 过期")
    print("  2. 染色体位置格式不匹配")
    print("  3. API 限制")

print("\n" + "="*70)
print("完成！")
print("="*70)
