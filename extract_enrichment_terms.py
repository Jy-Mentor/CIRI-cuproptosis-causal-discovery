# -*- coding: utf-8 -*-
import csv

# 定义文件路径
bp_file = "大创\GO_BP_enrichment_results.tsv"
cc_file = "大创\GO_CC_enrichment_results.tsv"
mf_file = "大创\GO_MF_enrichment_results.tsv"
kegg_file = "大创\KEGG_enrichment_results.tsv"

# 读取TSV文件并解析为字典列表
def read_tsv_file(file_path):
    data = []
    with open(file_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            # 转换数值字段
            if 'pvalue' in row:
                row['pvalue'] = float(row['pvalue'])
            if 'p.adjust' in row:
                row['p.adjust'] = float(row['p.adjust'])
            data.append(row)
    return data

# 读取文件
data_bp = read_tsv_file(bp_file)
data_cc = read_tsv_file(cc_file)
data_mf = read_tsv_file(mf_file)
data_kegg = read_tsv_file(kegg_file)

# 过滤显著富集的条目 (p<0.05 且 FDR<0.05)
data_bp_sig = [row for row in data_bp if row['pvalue'] < 0.05 and row['p.adjust'] < 0.05]
data_cc_sig = [row for row in data_cc if row['pvalue'] < 0.05 and row['p.adjust'] < 0.05]
data_mf_sig = [row for row in data_mf if row['pvalue'] < 0.05 and row['p.adjust'] < 0.05]
data_kegg_sig = [row for row in data_kegg if row['pvalue'] < 0.05 and row['p.adjust'] < 0.05]

# 按 Adjusted P-value 排序
data_bp_sig_sorted = sorted(data_bp_sig, key=lambda x: x['p.adjust'])
data_cc_sig_sorted = sorted(data_cc_sig, key=lambda x: x['p.adjust'])
data_mf_sig_sorted = sorted(data_mf_sig, key=lambda x: x['p.adjust'])
data_kegg_sig_sorted = sorted(data_kegg_sig, key=lambda x: x['p.adjust'])

# 提取 Top 条目
top3_bp = [row['Description'] for row in data_bp_sig_sorted[:3]]
top4_cc = [row['Description'] for row in data_cc_sig_sorted[:4]]
top4_mf = [row['Description'] for row in data_mf_sig_sorted[:4]]
top4_kegg = [row['Description'] for row in data_kegg_sig_sorted[:4]]

# 输出结果
print("Top 3 GO BP 显著富集条目:")
for i, term in enumerate(top3_bp, 1):
    print("{0}. {1}".format(i, term))

print("\nTop 4 GO CC 显著富集条目:")
for i, term in enumerate(top4_cc, 1):
    print("{0}. {1}".format(i, term))

print("\nTop 4 GO MF 显著富集条目:")
for i, term in enumerate(top4_mf, 1):
    print("{0}. {1}".format(i, term))

print("\nTop 4 KEGG 显著富集通路:")
for i, term in enumerate(top4_kegg, 1):
    print("{0}. {1}".format(i, term))

# 保存结果到文件
with open('enrichment_terms_summary.txt', 'w') as f:
    f.write("Top 3 GO BP 显著富集条目:\n")
    for term in top3_bp:
        f.write("{0}\n".format(term))
    
    f.write("\nTop 4 GO CC 显著富集条目:\n")
    for term in top4_cc:
        f.write("{0}\n".format(term))
    
    f.write("\nTop 4 GO MF 显著富集条目:\n")
    for term in top4_mf:
        f.write("{0}\n".format(term))
    
    f.write("\nTop 4 KEGG 显著富集通路:\n")
    for term in top4_kegg:
        f.write("{0}\n".format(term))

print("\n结果已保存到 enrichment_terms_summary.txt 文件")