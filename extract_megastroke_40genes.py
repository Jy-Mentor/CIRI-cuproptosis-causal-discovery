#!/usr/bin/env python3
"""
重新提取MEGASTROKE中40个基因的所有SNP数据
使用IEU OpenGWAS API
"""

import requests
import pandas as pd
import time
import json

# 配置
TOKEN = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiIxNzU3ODgyODc4QHFxLmNvbSIsImlhdCI6MTc3NzA5MzY5OSwiZXhwIjoxNzc4MzAzMjk5fQ.c9QuwWsF2C_XtupG3juu9OVw3PjVvJx2NrI9tU8xTgtI5E9Jay-DN1exaCDpY2r0j52QeVMiBLQMiAxsVAdQF9VWjCsbZdPO_78fqthP9zshGe-g0DviJOM9W4FdKkaePOG3zA1Qus9e8cZlhJSQvjT6fr8IES7gViZ5e5M8zljCjQ3nPi2FKGX_zgrrpTXlM78Kmev8TWFD6ZQ29IQb7wFriEgcVWbNRWCebcWZpTVPTTSwLOJ1JsegyuFOjkdyNjKg85YJCPw89OHUYwWUjc0fRbH_sFW5RQEGOV-UvEjo9K094TwIVCuTv8wHkkfF4Fro9QNa7ElWrd9sM55_qg"
GWAS_ID = "ebi-a-GCST006908"  # MEGASTROKE
HEADERS = {"Authorization": f"Bearer {TOKEN}"}

# 40个基因列表
genes = ["NFKB1", "RELA", "IKBKB", "STAT3", "STAT1", "JAK1", "IRF1", "TNF", "IL6", "TGFB1",
         "DLAT", "ATP7B", "CP", "NFE2L2", "HMOX1",
         "CAT", "GPX1", "XDH",
         "CASP8", "CTSB", "CTSL", "AIF1",
         "CPT1A", "FABP4", "PPARG", "SREBF1",
         "F3", "ICAM1", "CCL2", "TBXAS1", "PTGS1", "SPHK1",
         "HIF1A", "HSPA5", "MTOR", "EGFR", "AKT1", "MAPKAPK2", "TSPO", "ADRB1", "PARP1", "PTPRC"]

print("=" * 60)
print("提取MEGASTROKE中40个基因的所有SNP")
print(f"GWAS ID: {GWAS_ID}")
print(f"基因数: {len(genes)}")
print("=" * 60)

# 步骤1: 读取所有暴露数据，获取完整的SNP列表
print("\n步骤1: 读取暴露数据获取SNP列表...")

all_snps = set()

for pval in ["1e-05", "5e-08"]:
    exposure_file = f"D:/EQTL/clump/eQTLgen_allgene_p_{pval}_kb_1000_r2_0.01.xlsx"
    try:
        df = pd.read_excel(exposure_file)
        gene_df = df[df['gene'].str.upper().isin([g.upper() for g in genes])]
        snps = gene_df['SNP'].dropna().unique()
        all_snps.update(snps)
        print(f"  p={pval}: 找到 {len(snps)} 个SNP")
    except Exception as e:
        print(f"  p={pval}: 错误 - {e}")

rsid_list = list(all_snps)
print(f"\n总SNP数: {len(rsid_list)}")

if len(rsid_list) == 0:
    print("错误：未找到任何SNP")
    exit(1)

# 步骤2: 分批提取结局数据
print("\n步骤2: 提取MEGASTROKE数据...")

batch_size = 60  # API限制: N(id) * N(variant) <= 64
n_batches = (len(rsid_list) + batch_size - 1) // batch_size

all_results = []

for i in range(n_batches):
    start = i * batch_size
    end = min((i + 1) * batch_size, len(rsid_list))
    batch_snps = rsid_list[start:end]
    
    print(f"\n批次 {i+1}/{n_batches} ({len(batch_snps)} 个SNP)...")
    
    url = "https://api.opengwas.io/api/associations"
    data = {
        "id": [GWAS_ID],
        "variant": batch_snps
    }
    
    try:
        response = requests.post(url, headers=HEADERS, json=data, timeout=30)
        
        if response.status_code == 200:
            result_data = response.json()
            if isinstance(result_data, list) and len(result_data) > 0:
                df = pd.DataFrame(result_data)
                all_results.append(df)
                print(f"  ✓ 成功提取 {len(df)} 行")
            else:
                print(f"  ℹ️ 无数据返回")
        else:
            print(f"  ✗ 失败: {response.status_code}")
            print(f"     {response.text[:200]}")
            
    except Exception as e:
        print(f"  ✗ 错误: {str(e)[:100]}")
    
    if i < n_batches - 1:
        time.sleep(0.5)

# 步骤3: 合并并保存
print("\n" + "=" * 60)

if len(all_results) == 0:
    print("错误：未提取到任何结局数据！")
    exit(1)

outcome_data = pd.concat(all_results, ignore_index=True)
print(f"总匹配 SNP 数: {len(outcome_data)}")
print(f"唯一 SNP 数: {outcome_data['rsid'].nunique()}")

# 重命名列以匹配TwoSampleMR格式
column_mapping = {
    'rsid': 'SNP',
    'chromosome': 'chr',
    'position': 'pos.outcome',
    'ea': 'effect_allele.outcome',
    'nea': 'other_allele.outcome',
    'beta': 'beta.outcome',
    'se': 'se.outcome',
    'p': 'pval.outcome',
    'eaf': 'eaf.outcome',
    'n': 'samplesize.outcome'
}

existing_cols = {k: v for k, v in column_mapping.items() if k in outcome_data.columns}
outcome_data = outcome_data.rename(columns=existing_cols)

# 添加其他必要列
outcome_data['outcome'] = 'Ischemic Stroke'
outcome_data['id.outcome'] = GWAS_ID

# 保存
output_dir = "D:/EQTL/mr_results_megastroke"
import os
os.makedirs(output_dir, exist_ok=True)

output_file = os.path.join(output_dir, "megastroke_outcome_40genes.csv")
outcome_data.to_csv(output_file, index=False)
print(f"\n已保存到: {output_file}")

# 验证
print("\n" + "=" * 60)
print("数据验证:")
print(f"id.outcome: {outcome_data['id.outcome'].iloc[0]}")
print(f"outcome: {outcome_data['outcome'].iloc[0]}")

if 'samplesize.outcome' in outcome_data.columns:
    ss_max = outcome_data['samplesize.outcome'].max()
    print(f"samplesize.outcome 最大值: {ss_max}")
    if ss_max > 100000:
        print("✅ 验证通过：GWAS大样本数据")

print("\n下一步：重新运行MR分析，使用新的结局文件")
print(f'outcome_file <- "{output_file}"')
