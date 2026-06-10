#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
5数据集并集 → 人类同源基因转换 (v2 - 修复 mygene.info API)

两步转换:
  Step 1: 查询小鼠/大鼠 → 从 homologene 提取人类 Entrez Gene ID
  Step 2: 按 Entrez ID → 查询人类 gene symbol
"""

import pandas as pd
import requests
import time
import openpyxl

INPUT = r'c:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\deg_5datasets_summary.xlsx'
OUTPUT = r'c:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\human_targets_5datasets.xlsx'
UNMAPPED_OUTPUT = r'c:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\unmapped_5datasets.xlsx'

BATCH_SIZE = 100
SLEEP = 0.3

# ============================================================
# 1. 读取各数据集显著DEG
# ============================================================
print('=' * 60)
print('[1/5] 读取各数据集显著DEG...')
print('=' * 60)

wb = openpyxl.load_workbook(INPUT)

SHEET_LABEL_MAP = {
    'GSE16561': 'GSE16561', 'GSE37587': 'GSE37587',
    'GSE61616': 'GSE61616', 'GSE97537': 'GSE97537',
    '104036_all': 'GSE104036_all', '104036_3hr': 'GSE104036_3hr',
    '104036_6hr': 'GSE104036_6hr', '104036_12hr': 'GSE104036_12hr',
    '104036_24hr': 'GSE104036_24hr',
}

sheet_genes = {}
for sn, label in SHEET_LABEL_MAP.items():
    if sn not in wb.sheetnames:
        continue
    ws = wb[sn]
    rows = list(ws.values)
    hdr = rows[0]
    ci = {h: i for i, h in enumerate(hdr) if h}
    if 'gene_symbol' not in ci or 'significant' not in ci:
        continue
    genes = set()
    for r in rows[1:]:
        if len(r) > max(ci['gene_symbol'], ci['significant']):
            if r[ci['significant']] == True or str(r[ci['significant']]).lower() == 'true':
                genes.add(str(r[ci['gene_symbol']]).strip())
    sheet_genes[label] = genes
    print(f'  {label}: {len(genes)}')

# ============================================================
# 2. 按物种分离
# ============================================================
print('\n[2/5] 按物种分离...')

human_genes = set()
rat_genes = set()
mouse_genes = set()

for label, genes in sheet_genes.items():
    if label in ('GSE16561', 'GSE37587'):
        human_genes.update(genes)
    elif label in ('GSE61616', 'GSE97537'):
        rat_genes.update(genes)
    elif 'GSE104036' in label:
        mouse_genes.update(genes)

print(f'  人类: {len(human_genes)} | 大鼠: {len(rat_genes)} | 小鼠: {len(mouse_genes)}')

# ============================================================
# 3. Step 1: 查询同源基因 → 获取人类 Entrez ID
# ============================================================
print('\n[3/5] Step 1 — 查询同源基因获取人类 Entrez ID...')

def query_homologene(gene_list, species):
    """查询物种基因，从 homologene 提取人类 Entrez ID"""
    url = 'https://mygene.info/v3/query'
    entrez_map = {}  # {orig_gene: human_entrez_id}
    failed = []

    for i in range(0, len(gene_list), BATCH_SIZE):
        batch = gene_list[i:i + BATCH_SIZE]
        payload = {
            'q': batch,
            'scopes': 'symbol',
            'species': species,
            'fields': 'symbol,homologene',
        }
        try:
            r = requests.post(url, json=payload, timeout=60)
            if r.status_code == 200:
                hits = r.json()  # POST returns list directly
                for hit in hits:
                    orig = hit.get('symbol', '')
                    hg = hit.get('homologene')
                    if hg and 'genes' in hg:
                        for entry in hg['genes']:
                            if len(entry) >= 2 and entry[0] == 9606:
                                entrez_map[orig] = entry[1]
                                break
                hit_symbols = {h.get('symbol', '') for h in hits}
                for gene in batch:
                    if gene not in hit_symbols and gene not in entrez_map:
                        failed.append(gene)
            else:
                print(f'  HTTP {r.status_code}')
                failed.extend(batch)
        except Exception as e:
            print(f'  Error: {e}')
            failed.extend(batch)

        if (i // BATCH_SIZE + 1) % 5 == 0:
            print(f'  进度 {species}: {min(i + BATCH_SIZE, len(gene_list))}/{len(gene_list)}')
        time.sleep(SLEEP)

    return entrez_map, failed


def entrez_to_symbol(entrez_ids):
    """批量 Entrez ID → gene symbol"""
    url = 'https://mygene.info/v3/query'
    mapping = {}
    ids = list(set(entrez_ids.values()))

    for i in range(0, len(ids), BATCH_SIZE):
        batch = ids[i:i + BATCH_SIZE]
        payload = {
            'q': [str(e) for e in batch],
            'scopes': 'entrezgene',
            'species': 'human',
            'fields': 'symbol',
        }
        try:
            r = requests.post(url, json=payload, timeout=60)
            if r.status_code == 200:
                hits = r.json()
                for hit in hits:
                    eid = hit.get('_id', '')
                    sym = hit.get('symbol', '')
                    if eid and sym:
                        mapping[int(eid)] = sym
        except Exception as e:
            print(f'  Error: {e}')

        if (i // BATCH_SIZE + 1) % 5 == 0:
            print(f'  Symbol lookup: {min(i + BATCH_SIZE, len(ids))}/{len(ids)}')
        time.sleep(SLEEP)

    return mapping


# 小鼠
mouse_list = sorted(mouse_genes)
mouse_entrez, mouse_failed = query_homologene(mouse_list, 'mouse')
print(f'  小鼠→人类Entrez: {len(mouse_entrez)}/{len(mouse_list)}')

# 大鼠
rat_list = sorted(rat_genes)
rat_entrez, rat_failed = query_homologene(rat_list, 'rat')
print(f'  大鼠→人类Entrez: {len(rat_entrez)}/{len(rat_list)}')

# ============================================================
# 4. Step 2: Entrez ID → Human Symbol
# ============================================================
print('\n[4/5] Step 2 — 人类 Entrez ID → Gene Symbol...')

all_entrez = {}
for k, v in {**mouse_entrez, **rat_entrez}.items():
    all_entrez[k] = v

eid_to_symbol = entrez_to_symbol(all_entrez)
print(f'  Entrez→Symbol: {len(eid_to_symbol)}/{len(set(all_entrez.values()))}')

# 构建最终映射
final_map = {}  # {orig_mouse/rat_gene: human_symbol}
for orig_gene, eid in mouse_entrez.items():
    final_map[orig_gene] = eid_to_symbol.get(eid, f'ENTREZ_{eid}')
for orig_gene, eid in rat_entrez.items():
    final_map[orig_gene] = eid_to_symbol.get(eid, f'ENTREZ_{eid}')

# ============================================================
# 5. 合并输出
# ============================================================
print('\n[5/5] 合并输出...')

human_targets = set(human_genes)
human_targets.update(final_map.values())
human_targets = {g for g in human_targets if g and not g.startswith('ENTREZ_')}

# 来源追溯
def find_sources(hg):
    src = set()
    if hg in human_genes:
        for label in ('GSE16561', 'GSE37587'):
            if hg in sheet_genes.get(label, set()):
                src.add(label)
    for orig, mapped in final_map.items():
        if mapped == hg:
            for label, genes in sheet_genes.items():
                if orig in genes:
                    src.add(label)
    return '; '.join(sorted(src))

rows = []
for gene in sorted(human_targets):
    rows.append({'human_gene_symbol': gene, 'source_datasets': find_sources(gene)})

pd.DataFrame(rows).to_excel(OUTPUT, index=False, sheet_name='human_targets')
print(f'  输出: {OUTPUT}')
print(f'  人类靶基因数: {len(rows)}')

# 未映射
unmapped = []
for g in sorted(set(mouse_failed)):
    unmapped.append({'species': 'mouse', 'gene_symbol': g})
for g in sorted(set(rat_failed)):
    unmapped.append({'species': 'rat', 'gene_symbol': g})
if unmapped:
    pd.DataFrame(unmapped).to_excel(UNMAPPED_OUTPUT, index=False, sheet_name='unmapped')
    print(f'  未映射: {len(unmapped)} → {UNMAPPED_OUTPUT}')

print('\n' + '=' * 60)
print('完成!')
print('=' * 60)