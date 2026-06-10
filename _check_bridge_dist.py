import csv

rows = []
with open(r'C:\Users\Jy-Mentor-7\Desktop\随机森林\bridge_predictions_v8_final.csv', 'r', encoding='utf-8-sig') as f:
    reader = csv.DictReader(f)
    for r in reader:
        rows.append(r)

total = len(rows)
keys = list(rows[0].keys())
print(f'列名: {keys}')
print(f'总基因数: {total}')

# 直接用索引
GENE = keys[0]
IS_DT = keys[1]
IS_DG = keys[2] 
IS_BG = keys[3]
P_DT = keys[4]
P_DG = keys[5]
P_BG = keys[6]

bridges = [r for r in rows if r[IS_BG] == '1']
drug_only = [r for r in rows if r[IS_DT] == '1' and r[IS_BG] == '0']
disease_only = [r for r in rows if r[IS_DG] == '1' and r[IS_BG] == '0']
neither = [r for r in rows if r[IS_DT] == '0' and r[IS_DG] == '0']

print(f'\n=== 类别分布 ===')
print(f'桥基因(is_bridge=1):      {len(bridges):6d} ({len(bridges)/total*100:.1f}%)')
print(f'仅DrugTarget:             {len(drug_only):6d} ({len(drug_only)/total*100:.1f}%)')
print(f'仅DiseaseGene:            {len(disease_only):6d} ({len(disease_only)/total*100:.1f}%)')
print(f'两者皆非:                 {len(neither):6d} ({len(neither)/total*100:.1f}%)')

print(f'\n=== prob_drug_target 分段 ===')
pdts = [float(r[P_DT]) for r in rows]
for t in [0.999, 0.99, 0.9, 0.7, 0.5, 0.3, 0.1, 0.01, 0.001]:
    c = sum(1 for p in pdts if p > t)
    print(f'  >{t:6.3f}: {c:6d} ({c/total*100:5.1f}%)')

print(f'\n=== prob_disease_gene 分段 ===')
pds = [float(r[P_DG]) for r in rows]
for t in [0.999, 0.99, 0.9, 0.7, 0.5, 0.3, 0.1, 0.01, 0.001]:
    c = sum(1 for p in pds if p > t)
    print(f'  >{t:6.3f}: {c:6d} ({c/total*100:5.1f}%)')

print(f'\n=== bridge_prob 分段 ===')
bps = [float(r[P_BG]) for r in rows]
for t in [0.99, 0.9, 0.5, 0.1, 0.01, 0.001, 0.0001]:
    c = sum(1 for b in bps if b > t)
    print(f'  >{t}: {c:6d} ({c/total*100:5.1f}%)')
print(f'  median: {sorted(bps)[len(bps)//2]:.6f}')

print(f'\n=== Top 20 bridge_prob (is_bridge_gene=1) ===')
bridges_sorted = sorted(bridges, key=lambda x: float(x[P_BG]), reverse=True)
for i, r in enumerate(bridges_sorted[:20]):
    g = r[GENE]
    print(f'  {i+1:2d}. {g:12s}  bg={float(r[P_BG]):.4f}  dt={float(r[P_DT]):.4f}  dg={float(r[P_DG]):.4f}')

print(f'\n=== NLRP3 / CP / NFE2L2 ===')
for r in rows:
    g = r[GENE]
    if g in ['NLRP3', 'CP', 'NFE2L2']:
        print(f'  {g}: is_dt={r[IS_DT]} is_dg={r[IS_DG]} is_bg={r[IS_BG]}')
        print(f'       prob_dt={float(r[P_DT]):.4f}  prob_dg={float(r[P_DG]):.4f}  bridge_prob={float(r[P_BG]):.4f}')

print(f'\n=== 在{len(bridges)}个桥基因中的排名 ===')
for rank, r in enumerate(bridges_sorted):
    g = r[GENE]
    if g in ['NLRP3', 'CP', 'NFE2L2']:
        print(f'  #{rank+1} {g}: bridge_prob={float(r[P_BG]):.4f}')

high_dt = [r for r in rows if float(r[P_DT]) > 0.99]
print(f'\n=== 高置信DrugTarget (>0.99) ===')
print(f'  总数: {len(high_dt)}')
print(f'  其中 prob_disease > 0.5: {sum(1 for r in high_dt if float(r[P_DG]) > 0.5)} 个')
print(f'  其中 prob_disease > 0.9: {sum(1 for r in high_dt if float(r[P_DG]) > 0.9)} 个')
print(f'  其中 is_bridge=1: {sum(1 for r in high_dt if r[IS_BG]=="1")} 个')
