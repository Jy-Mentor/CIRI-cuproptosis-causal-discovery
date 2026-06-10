#!/usr/bin/env python3
"""
从IEU OpenGWAS提取MEGASTROKE结局数据
使用正确的API端点
"""

import requests
import pandas as pd
import time
import json

# ============================================
# 配置
# ============================================
TOKEN = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiIxNzU3ODgyODc4QHFxLmNvbSIsImlhdCI6MTc3NzA5MzY5OSwiZXhwIjoxNzc4MzAzMjk5fQ.c9QuwWsF2C_XtupG3juu9OVw3PjVvJx2NrI9tU8xTgtI5E9Jay-DN1exaCDpY2r0j52QeVMiBLQMiAxsVAdQF9VWjCsbZdPO_78fqthP9zshGe-g0DviJOM9W4FdKkaePOG3zA1Qus9e8cZlhJSQvjT6fr8IES7gViZ5e5M8zljCjQ3nPi2FKGX_zgrrpTXlM78Kmev8TWFD6ZQ29IQb7wFriEgcVWbNRWCebcWZpTVPTTSwLOJ1JsegyuFOjkdyNjKg85YJCPw89OHUYwWUjc0fRbH_sFW5RQEGOV-UvEjo9K094TwIVCuTv8wHkkfF4Fro9QNa7ElWrd9sM55_qg"
GWAS_ID = "ebi-a-GCST006908"  # MEGASTROKE
HEADERS = {"Authorization": f"Bearer {TOKEN}"}

print("=" * 60)
print("提取MEGASTROKE结局数据")
print(f"GWAS ID: {GWAS_ID}")
print("=" * 60)

# ============================================
# 步骤1: 读取暴露数据，获取 SNP 列表
# ============================================
print("\n步骤1: 读取暴露数据...")

exposure_file = "D:/EQTL/clump/eQTLgen_allgene_p_1e-05_kb_1000_r2_0.01.xlsx"
exposure_data = pd.read_excel(exposure_file)

# 15个基因
genes = ["NFKB1", "STAT3", "HIF1A", "HSPA5", "HMOX1",
         "RELA", "NFE2L2", "CP", "LIAS", "IKBKB",
         "JAK1", "PARP1", "CASP8", "MTOR", "PTPRC"]

# 筛选目标基因
exp_snps = exposure_data[exposure_data['gene'].str.upper().isin([g.upper() for g in genes])]
print(f"目标基因SNP数: {len(exp_snps)}")

# 获取唯一 rsID
rsid_list = exp_snps['SNP'].dropna().unique().tolist()
print(f"暴露 SNP 总数: {len(rsid_list)}")

# ============================================
# 步骤2: 分批提取结局数据
# ============================================
print("\n步骤2: 提取MEGASTROKE数据...")

batch_size = 100
n_batches = (len(rsid_list) + batch_size - 1) // batch_size

all_results = []

for i in range(n_batches):
    start = i * batch_size
    end = min((i + 1) * batch_size, len(rsid_list))
    batch_snps = rsid_list[start:end]
    
    print(f"\n批次 {i+1}/{n_batches} ({len(batch_snps)} 个SNP)...")
    
    # 使用associations端点 (POST方法)
    url = f"https://api.opengwas.io/api/associations"
    data = {
        "id": [GWAS_ID],
        "variant": batch_snps
    }
    
    try:
        response = requests.post(url, headers=HEADERS, json=data, timeout=30)
        
        if response.status_code == 200:
            data = response.json()
            if isinstance(data, list) and len(data) > 0:
                df = pd.DataFrame(data)
                all_results.append(df)
                print(f"  ✓ 成功提取 {len(df)} 行")
            else:
                print(f"  ℹ️ 无数据返回")
        else:
            print(f"  ✗ 失败: {response.status_code}")
            print(f"     {response.text[:200]}")
            
    except Exception as e:
        print(f"  ✗ 错误: {str(e)[:100]}")
    
    # 避免限流
    if i < n_batches - 1:
        time.sleep(0.5)

# ============================================
# 步骤3: 合并并保存
# ============================================
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
    'chromosome': 'chr.outcome',
    'position': 'pos.outcome',
    'ea': 'effect_allele.outcome',
    'nea': 'other_allele.outcome',
    'beta': 'beta.outcome',
    'se': 'se.outcome',
    'p': 'pval.outcome',
    'eaf': 'eaf.outcome',
    'n': 'samplesize.outcome'
}

# 只保留存在的列
existing_cols = {k: v for k, v in column_mapping.items() if k in outcome_data.columns}
outcome_data = outcome_data.rename(columns=existing_cols)

# 添加其他必要列
outcome_data['outcome'] = 'Ischemic Stroke'
outcome_data['id.outcome'] = GWAS_ID

# 保存
output_dir = "D:/EQTL/mr_results_megastroke"
import os
os.makedirs(output_dir, exist_ok=True)

output_file = os.path.join(output_dir, "megastroke_outcome.csv")
outcome_data.to_csv(output_file, index=False)
print(f"\n已保存到: {output_file}")

# ============================================
# 步骤4: 验证
# ============================================
print("\n" + "=" * 60)
print("数据验证:")
print(f"id.outcome: {outcome_data['id.outcome'].iloc[0]}")
print(f"outcome: {outcome_data['outcome'].iloc[0]}")

if 'samplesize.outcome' in outcome_data.columns:
    ss_max = outcome_data['samplesize.outcome'].max()
    print(f"samplesize.outcome 最大值: {ss_max}")
    if ss_max > 100000:
        print("✅ 验证通过：GWAS大样本数据")
    else:
        print("⚠️ 样本量较小")

print("\n下一步：修改MR脚本，使用新的结局文件")
print(f'outcome_file <- "{output_file}"')
