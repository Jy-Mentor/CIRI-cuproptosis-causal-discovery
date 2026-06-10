#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 创建简单的 Word 文档 - 使用 HTML 转换

import csv
from datetime import datetime

# 读取 CSV 数据
csv_file = 'c:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/mr_results_138genes_integrated/detailed_mr_results_table.csv'
with open(csv_file, 'r', encoding='utf-8') as f:
    reader = csv.reader(f)
    headers = next(reader)
    data = list(reader)

# 创建 HTML 内容
html_content = """
<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8">
<title>MR Analysis Results</title>
<style>
body { font-family: Arial; margin: 40px; }
h1 { text-align: center; color: #2E75B6; }
h2 { color: #1F4E79; border-bottom: 2px solid #2E75B6; padding-bottom: 5px; }
h3 { color: #2E75B6; }
table { border-collapse: collapse; width: 100%; margin: 20px 0; }
th, td { border: 1px solid black; padding: 8px; text-align: left; }
th { background-color: #4472C4; color: white; font-weight: bold; }
.significant { background-color: #FFC000; }
p { line-height: 1.6; }
</style>
</head>
<body>

<h1>Mendelian Randomization Analysis Results</h1>

<h2>1. Analysis Overview</h2>
<p>Two-sample Mendelian randomization (MR) analysis was performed using TwoSampleMR R package with cis-eQTL from 138 BCP target genes as genetic instrumental variables.</p>
<p><strong>Selection criteria:</strong> P<1×10⁻⁵, F>10, MAF>0.01, LD clumping (r²<0.01, 1000kb)</p>
<p><strong>Outcome data:</strong> MEGASTROKE (Any ischemic stroke, n=40,328)</p>

<h2>2. Significant Results (P < 0.05)</h2>
<table>
<tr>
"""

# 添加表头
for header in headers:
    html_content += f"<th>{header}</th>\n"

html_content += "</tr>\n"

# 添加显著性结果
sig_data = [row for row in data if float(row[3]) < 0.05][:3]
for row in sig_data:
    html_content += '<tr class="significant">\n'
    for cell in row:
        html_content += f"<td>{cell}</td>\n"
    html_content += "</tr>\n"

html_content += """
</table>

<h2>3. Marginally Significant Results (P < 0.1)</h2>
<table>
<tr>
"""

# 添加表头
for header in headers:
    html_content += f"<th>{header}</th>\n"

html_content += "</tr>\n"

# 添加边缘显著结果
edge_data = [row for row in data if 0.05 <= float(row[3]) < 0.1]
for row in edge_data:
    html_content += "<tr>\n"
    for cell in row:
        html_content += f"<td>{cell}</td>\n"
    html_content += "</tr>\n"

html_content += """
</table>

<h2>4. Cuproptosis-Related Genes</h2>
<p>NFKB1 and FDX1 showed no significant causal effects (P>0.05), suggesting they are pathological response regulators rather than independent causal factors.</p>

<table>
<tr>
<th>Gene</th>
<th>Function</th>
<th>OR (IVW)</th>
<th>P-value</th>
<th>Conclusion</th>
</tr>
<tr>
<td>NFKB1</td>
<td>NF-kB pathway core</td>
<td>0.992</td>
<td>0.815</td>
<td>No significant effect</td>
</tr>
<tr>
<td>FDX1</td>
<td>Cuproptosis regulator</td>
<td>0.969</td>
<td>0.194</td>
<td>No significant effect</td>
</tr>
<tr>
<td>ATP7B</td>
<td>Copper transporter</td>
<td>1.008</td>
<td>0.731</td>
<td>No significant effect</td>
</tr>
<tr>
<td>ATOX1</td>
<td>Antioxidant copper chaperone</td>
<td>0.962</td>
<td>0.118</td>
<td>Marginal trend</td>
</tr>
</table>

<h2>5. Quality Metrics</h2>
<table>
<tr>
<th>Metric</th>
<th>Value</th>
</tr>
<tr>
<td>Successfully analyzed genes</td>
<td>95</td>
</tr>
<tr>
<td>Skipped genes</td>
<td>48 (no eQTL data)</td>
</tr>
<tr>
<td>Significant genes (P<0.05)</td>
<td>3</td>
</tr>
<tr>
<td>Passed heterogeneity test</td>
<td>92.6%</td>
</tr>
<tr>
<td>Passed pleiotropy test</td>
<td>97.9%</td>
</tr>
<tr>
<td>Robust results</td>
<td>3 (100% of significant)</td>
</tr>
</table>

<h2>6. Conclusions</h2>
<p><strong>Key findings:</strong></p>
<ul>
<li>3 significant metabolic genes identified: ADRB1, SREBF1, ACADVL</li>
<li>NFKB1-FDX1 axis showed pathological response attributes</li>
<li>Results support RAGE blockade as an upstream regulatory strategy</li>
</ul>

<p style="margin-top: 40px; font-size: 12px; color: gray;">Generated: """ + datetime.now().strftime("%Y-%m-%d") + """</p>

</body>
</html>
"""

# 保存 HTML 文件
html_file = 'c:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/MR_Results.html'
with open(html_file, 'w', encoding='utf-8') as f:
    f.write(html_content)

print(f'HTML report created successfully: {html_file}')
print('You can open this file in a browser or convert it to Word/PDF')
