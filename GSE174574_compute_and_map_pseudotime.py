"""
计算伪时间并映射到真实时间
"""

import scanpy as sc
import numpy as np
import pandas as pd
from pathlib import Path
import json
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

DATA_FILE = Path(r'C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\GSE174574_analysis\results\GSE174574_processed.h5ad')
OUTPUT_DIR = Path(r'C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\GSE174574_analysis')

BULK_TIME_POINTS = {
    'Homeostatic': 0,
    'M2': 2,
    'M1': 8,
    'DAM': 24
}

adata = sc.read_h5ad(DATA_FILE)
print(f"数据: {adata.n_obs} cells × {adata.n_vars} genes")
print(f"细胞类型: {adata.obs['cell_type'].value_counts().to_dict()}")

homeo_mask = adata.obs['cell_type'] == 'Homeostatic'
homeo_cells = np.where(homeo_mask)[0]
root_cell_idx = homeo_cells[0]

adata.uns['iroot'] = root_cell_idx
print(f"根细胞索引: {root_cell_idx} (Homeostatic)")

sc.tl.dpt(adata, n_dcs=10)

adata.obs['pseudotime'] = adata.obs['dpt_pseudotime']

print("\n伪时间统计:")
for ct in adata.obs['cell_type'].unique():
    mask = adata.obs['cell_type'] == ct
    pt = adata.obs.loc[mask, 'pseudotime'].dropna()
    if len(pt) > 0:
        print(f"  {ct}: mean={pt.mean():.4f}, median={pt.median():.4f}")

cell_type_median = {}
for ct in adata.obs['cell_type'].unique():
    mask = adata.obs['cell_type'] == ct
    cell_type_median[ct] = float(adata.obs.loc[mask, 'pseudotime'].median())

sorted_ct = sorted(cell_type_median.items(), key=lambda x: x[1])
print(f"\n伪时间排序: {sorted_ct}")

mapping_table = {}
pseudo_times = []
real_times = []

for ct, real_time in BULK_TIME_POINTS.items():
    if ct in cell_type_median:
        mapping_table[ct] = {
            'pseudotime_median': cell_type_median[ct],
            'real_time_hours': real_time,
            'cell_count': int((adata.obs['cell_type'] == ct).sum())
        }
        pseudo_times.append(cell_type_median[ct])
        real_times.append(real_time)
        print(f"  {ct}: pseudo={cell_type_median[ct]:.4f} → {real_time}h")

if len(pseudo_times) >= 2:
    mapping_func = interp1d(pseudo_times, real_times, kind='linear', fill_value='extrapolate')
    
    pseudo_times_all = adata.obs['pseudotime'].values
    real_times_all = mapping_func(pseudo_times_all)
    real_times_all = np.clip(real_times_all, 0, 24)
    
    adata.obs['real_time_hours'] = real_times_all
    
    print(f"\n映射后:")
    print(f"  范围: {adata.obs['real_time_hours'].min():.2f} - {adata.obs['real_time_hours'].max():.2f} h")
    
    for ct in adata.obs['cell_type'].unique():
        mask = adata.obs['cell_type'] == ct
        ct_times = adata.obs.loc[mask, 'real_time_hours']
        print(f"  {ct}: {ct_times.mean():.2f} ± {ct_times.std():.2f} h")
    
    l2c_dir = OUTPUT_DIR / 'l2c_interface'
    l2c_dir.mkdir(parents=True, exist_ok=True)
    
    l2c_config = {
        'time_mapping': {
            'method': 'DPT pseudotime → real hours',
            'reference': 'GSE23160 Bulk (0h, 2h, 8h, 24h)',
            'mapping_table': mapping_table,
            'cell_ordering': [ct for ct, _ in sorted_ct]
        },
        'l2c_usage': {
            'time_column': 'real_time_hours',
            'time_unit': 'hours',
            'time_range': [0, 24],
            'velocity_field_init': 'Use GSE23160 temporal logFC',
            'cell_state_path': ['Homeostatic', 'M2', 'M1', 'DAM']
        }
    }
    
    with open(l2c_dir / 'l2c_real_time_config.json', 'w', encoding='utf-8') as f:
        json.dump(l2c_config, f, indent=2, ensure_ascii=False)
    
    adata.write_h5ad(l2c_dir / 'GSE174574_with_real_time.h5ad')
    
    cell_data = adata.obs[['cell_type', 'pseudotime', 'real_time_hours']].copy()
    cell_data.to_csv(l2c_dir / 'cell_time_mapping.csv')
    
    print("\n✅ 完成! 输出文件:")
    print(f"  - l2c_interface/l2c_real_time_config.json")
    print(f"  - l2c_interface/GSE174574_with_real_time.h5ad")
    print(f"  - l2c_interface/cell_time_mapping.csv")
else:
    print("\n⚠️ 映射点不足")
