"""
铜死亡核心基因在 HGT 预测结果中的检索分析
"""

import pandas as pd
import warnings
warnings.filterwarnings('ignore')

cuproptosis_genes = [
    'FDX1', 'LIAS', 'LIPT1', 'DLD', 'DLAT', 'PDHA1', 'PDHB', 'MTF1', 'GLS',
    'CDKN2A', 'COX11', 'MFN2', 'TOMM20', 'NDUFB9', 'ATP6V1E1', 'NFE2L2',
    'NLRP3', 'ATP7B', 'ATP7A', 'SLC31A1', 'LIPT2', 'DBT', 'GCSH', 'DLST',
    'SURF1', 'NDUFB2', 'NDUFB6', 'NDUFA8', 'NDUFA1', 'NDUFC1'
]

print('=' * 80)
print('铜死亡核心基因在 HGT 预测结果中的检索')
print('=' * 80)

found_genes = []
all_results = []

for i, chunk in enumerate(pd.read_csv('bridge_pathway_scores.csv', chunksize=1000000)):
    mask = chunk['gene_symbol'].isin(cuproptosis_genes)
    if mask.any():
        matches = chunk[mask]
        found_genes.extend(matches['gene_symbol'].unique().tolist())
        high = matches[matches['score'] >= 0.8]
        if len(high) > 0:
            all_results.append(high)
    
    if (i + 1) % 10 == 0:
        print(f'  已扫描 {i+1} 个数据块...')

found_genes = sorted(list(set(found_genes)))
print(f'\n  在 30 个铜死亡基因中，{len(found_genes)} 个出现在预测结果中:')
print(f'  {found_genes}')

not_found = [g for g in cuproptosis_genes if g not in found_genes]
print(f'\n  未出现在预测结果中的基因 ({len(not_found)} 个):')
print(f'  {not_found}')

df_high = pd.concat(all_results, ignore_index=True) if all_results else pd.DataFrame()
print(f'\n  高分预测 (score >= 0.8): {len(df_high):,} 条')

if len(df_high) > 0:
    print(f'\n  每个铜死亡基因的 Top-5 预测通路:')
    print(f'  {"基因":<12} {"通路":<75} {"分数":<10} {"模式":<12} {"排名"}')
    print('  ' + '-' * 115)
    
    for gene in cuproptosis_genes:
        gene_data = df_high[df_high['gene_symbol'] == gene].sort_values('rank').head(5)
        if len(gene_data) > 0:
            for _, row in gene_data.iterrows():
                print(f'  {row["gene_symbol"]:<12} {row["pathway_name"][:73]:<75} {row["score"]:<10.4f} {row["eval_mode"]:<12} {int(row["rank"])}')
            print()
    
    print(f'\n  铜死亡基因相关通路汇总 (score >= 0.8, 按最高分降序):')
    pathway_summary = df_high.groupby('pathway_name').agg({
        'gene_symbol': lambda x: ', '.join(sorted(set(x))),
        'score': 'max'
    }).reset_index()
    pathway_summary.columns = ['pathway', 'genes', 'max_score']
    pathway_summary = pathway_summary.sort_values('max_score', ascending=False)
    
    for _, row in pathway_summary.head(40).iterrows():
        print(f'  {row["pathway"][:70]:<70} {row["max_score"]:.4f}  [{row["genes"]}]')
    
    pathway_summary.to_csv('cuproptosis_genes_pathway_results.csv', index=False, encoding='utf-8-sig')
    print(f'\n  结果已保存至: cuproptosis_genes_pathway_results.csv')
