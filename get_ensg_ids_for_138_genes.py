#!/usr/bin/env python3
"""
通过 Ensembl REST API 批量查询 138 个基因的 ENSG ID
"""

import requests
import pandas as pd
import json
from typing import Dict, List, Optional

# 目标基因列表（138 个，去重后 144 个）
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

def get_ensg_ids_batch(gene_symbols: List[str], species: str = "homo_sapiens") -> Dict[str, Optional[str]]:
    """
    通过 Ensembl REST API 批量查询 ENSG ID
    使用 POST 端点一次性查询所有基因
    """
    server = "https://rest.ensembl.org"
    ext = f"/xrefs/symbol/{species}"
    
    # 准备 POST 数据
    post_data = {
        "external_db": "HGNC",
        "symbols": gene_symbols
    }
    
    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json"
    }
    
    print("正在查询 Ensembl API...")
    print(f"查询 {len(gene_symbols)} 个基因...\n")
    
    try:
        # POST 请求批量查询
        r = requests.post(
            server + ext,
            headers=headers,
            json=post_data,
            timeout=120
        )
        
        if r.status_code == 200:
            all_results = r.json()
            
            # 构建映射字典
            ensg_map = {}
            for result in all_results:
                symbol = result.get('symbol')
                ensg_id = result.get('primary_id')
                
                if symbol and ensg_id and ensg_id.startswith('ENSG'):
                    # 如果已有映射，保留第一个（主转录本）
                    if symbol not in ensg_map:
                        ensg_map[symbol] = ensg_id
            
            return ensg_map
        else:
            print(f"API 错误：{r.status_code}")
            print(f"响应：{r.text[:200]}")
            return {}
            
    except Exception as e:
        print(f"查询失败：{e}")
        return {}

# 查询 ENSG ID
ensg_map = get_ensg_ids_batch(target_genes)

# 统计结果
print(f"\n查询完成！")
print(f"成功映射：{len(ensg_map)}/{len(target_genes)} 个基因\n")

# 找出未映射的基因
missing_genes = [gene for gene in target_genes if gene not in ensg_map]
if missing_genes:
    print(f"未映射的基因 ({len(missing_genes)} 个):")
    print(", ".join(missing_genes[:20]) + ("..." if len(missing_genes) > 20 else ""))
    print()

# 创建 DataFrame 并保存
df = pd.DataFrame({
    'gene_symbol': list(ensg_map.keys()),
    'ensg_id': list(ensg_map.values())
})

# 按基因符号排序
df = df.sort_values('gene_symbol').reset_index(drop=True)

# 保存为 CSV
output_file = "gene_ensg_mapping.csv"
df.to_csv(output_file, index=False, encoding='utf-8-sig')
print(f"✓ 保存映射文件：{output_file}")
print(f"  共 {len(df)} 行")

# 显示示例
print(f"\n映射示例:")
print(df.head(20).to_string(index=False))

# 生成 R 语言格式的映射列表
print(f"\nR 语言格式映射列表:")
print("gene_to_ensg <- list(")
for _, row in df.iterrows():
    print(f'  "{row["gene_symbol"]}" = "{row["ensg_id"]}",')
print(")")

print(f"\n完成！")
