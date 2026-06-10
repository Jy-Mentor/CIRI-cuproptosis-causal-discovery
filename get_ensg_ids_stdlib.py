#!/usr/bin/env python3
"""
通过 Ensembl REST API 批量查询 138 个基因的 ENSG ID
纯标准库版本：使用 urllib，单个基因查询
"""

import urllib.request
import urllib.error
import json
from typing import Dict, List, Optional

# 目标基因列表（138 个，去重后）
target_genes = [
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
]

# 去重
target_genes = sorted(list(set(target_genes)))
print(f"目标基因数：{len(target_genes)}\n")

def get_ensg_id_single(gene_symbol: str, species: str = "homo_sapiens") -> Optional[str]:
    """
    通过 Ensembl REST API 查询单个基因的 ENSG ID
    使用 GET 端点
    """
    server = "https://rest.ensembl.org"
    ext = f"/xrefs/symbol/{species}/{gene_symbol}?external_db=HGNC"
    url = server + ext
    
    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json"
    }
    
    try:
        req = urllib.request.Request(url, headers=headers)
        
        with urllib.request.urlopen(req, timeout=30) as response:
            result = response.read().decode('utf-8')
            data = json.loads(result)
            
            # 查找 ENSG ID（type 为 gene 且 id 以 ENSG 开头）
            for xref in data:
                if xref.get('type') == 'gene':
                    ensg_id = xref.get('id')
                    if ensg_id and ensg_id.startswith('ENSG'):
                        return ensg_id
            
            return None
            
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return None  # 基因不存在
        print(f"HTTP 错误 {gene_symbol}: {e.code}")
        return None
    except Exception as e:
        print(f"查询失败 {gene_symbol}: {e}")
        return None

# 逐个查询 ENSG ID
ensg_map = {}
print("正在查询 Ensembl API...")
for i, gene in enumerate(target_genes):
    ensg = get_ensg_id_single(gene)
    if ensg:
        ensg_map[gene] = ensg
        print(f"  [{i+1:3d}/{len(target_genes)}] {gene:15s} -> {ensg}")
    else:
        print(f"  [{i+1:3d}/{len(target_genes)}] {gene:15s} -> NOT FOUND")

# 统计结果
print(f"\n查询完成！")
print(f"成功映射：{len(ensg_map)}/{len(target_genes)} 个基因\n")

# 找出未映射的基因
missing_genes = [gene for gene in target_genes if gene not in ensg_map]
if missing_genes:
    print(f"未映射的基因 ({len(missing_genes)} 个):")
    print(", ".join(missing_genes[:20]) + ("..." if len(missing_genes) > 20 else ""))
    print()

# 保存为 CSV 格式（手动构建）
output_file = "gene_ensg_mapping.csv"
with open(output_file, 'w', encoding='utf-8-sig') as f:
    f.write("gene_symbol,ensg_id\n")
    for gene in sorted(ensg_map.keys()):
        f.write(f"{gene},{ensg_map[gene]}\n")

print(f"✓ 保存映射文件：{output_file}")
print(f"  共 {len(ensg_map)} 行")

# 显示示例
print(f"\n映射示例 (前 20 个):")
for i, (gene, ensg) in enumerate(sorted(ensg_map.items())[:20]):
    print(f"  {gene:15s} -> {ensg}")

# 生成 R 语言格式的映射列表
print(f"\nR 语言格式映射列表:")
print("gene_to_ensg <- list(")
for gene in sorted(ensg_map.keys()):
    print(f'  "{gene}" = "{ensg_map[gene]}",')
print(")")

print(f"\n完成！")
