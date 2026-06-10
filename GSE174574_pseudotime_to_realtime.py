"""
伪时间→真实时间映射校准
将scRNA伪时间转换为L2c可用的真实时间变量
"""

import scanpy as sc
import numpy as np
import pandas as pd
from pathlib import Path
import json
from scipy.interpolate import interp1d

DATA_FILE = Path(r'C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\GSE174574_analysis\results\GSE174574_processed.h5ad')
OUTPUT_DIR = Path(r'C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\GSE174574_analysis')

# 真实时间点 (小时)
BULK_TIME_POINTS = {
    'Homeostatic': 0,   # 基线
    'M2': 2,           # 2h
    'M1': 8,           # 8h
    'DAM': 24          # 24h
}

def map_pseudotime_to_real_time():
    print("="*70)
    print("伪时间 → 真实时间映射校准")
    print("="*70)
    
    adata = sc.read_h5ad(DATA_FILE)
    print(f"\n数据: {adata.n_obs} cells × {adata.n_vars} genes")
    
    if 'pseudotime' not in adata.obs.columns:
        print("  ❌ 伪时间列不存在!")
        return None
    
    cell_types = adata.obs['cell_type'].unique()
    print(f"\n细胞类型: {list(cell_types)}")
    
    pseudotime_stats = {}
    for ct in cell_types:
        mask = adata.obs['cell_type'] == ct
        pt_values = adata.obs.loc[mask, 'pseudotime']
        pseudotime_stats[ct] = {
            'mean': float(pt_values.mean()),
            'median': float(pt_values.median()),
            'min': float(pt_values.min()),
            'max': float(pt_values.max()),
            'std': float(pt_values.std())
        }
    
    print("\n伪时间统计:")
    for ct, stats in pseudotime_stats.items():
        print(f"  {ct}: mean={stats['mean']:.2f}, median={stats['median']:.2f}")
    
    mapping_table = {}
    for ct, real_time in BULK_TIME_POINTS.items():
        if ct in pseudotime_stats:
            mapping_table[ct] = {
                'pseudotime_median': pseudotime_stats[ct]['median'],
                'real_time_hours': real_time,
                'cell_count': int((adata.obs['cell_type'] == ct).sum())
            }
    
    print("\n伪时间→真实时间映射:")
    for ct, mapping in mapping_table.items():
        print(f"  {ct}: pseudo={mapping['pseudotime_median']:.2f} → {mapping['real_time_hours']}h")
    
    return adata, mapping_table, pseudotime_stats

def create_pseudotime_mapping_function(mapping_table):
    """创建伪时间→真实时间的插值函数"""
    print("\n" + "="*70)
    print("创建时间映射函数")
    print("="*70)
    
    pseudo_times = [v['pseudotime_median'] for v in mapping_table.values()]
    real_times = [v['real_time_hours'] for v in mapping_table.values()]
    
    print(f"\n映射点:")
    for pt, rt in zip(pseudo_times, real_times):
        print(f"  pseudo={pt:.2f} → real={rt}h")
    
    mapping_func = interp1d(pseudo_times, real_times, 
                           kind='linear', fill_value='extrapolate')
    
    return mapping_func

def apply_mapping_to_cells(adata, mapping_func):
    """将映射应用到所有细胞"""
    print("\n" + "="*70)
    print("应用时间映射到所有细胞")
    print("="*70)
    
    adata_mapped = adata.copy()
    
    pseudo_times = adata_mapped.obs['pseudotime'].values
    real_times = mapping_func(pseudo_times)
    
    adata_mapped.obs['real_time_hours'] = real_times
    adata_mapped.obs['real_time_hours'] = adata_mapped.obs['real_time_hours'].clip(lower=0)
    
    print(f"\n映射后真实时间统计:")
    print(f"  范围: {adata_mapped.obs['real_time_hours'].min():.2f} - {adata_mapped.obs['real_time_hours'].max():.2f} h")
    print(f"  均值: {adata_mapped.obs['real_time_hours'].mean():.2f} h")
    print(f"  中位数: {adata_mapped.obs['real_time_hours'].median():.2f} h")
    
    for ct in adata_mapped.obs['cell_type'].unique():
        mask = adata_mapped.obs['cell_type'] == ct
        ct_times = adata_mapped.obs.loc[mask, 'real_time_hours']
        print(f"  {ct}: {ct_times.mean():.2f} ± {ct_times.std():.2f} h")
    
    return adata_mapped

def prepare_l2c_with_real_time(adata_mapped, mapping_table):
    """准备带真实时间的L2c数据"""
    print("\n" + "="*70)
    print("准备L2c数据 (含真实时间)")
    print("="*70)
    
    l2c_dir = OUTPUT_DIR / 'l2c_interface'
    l2c_dir.mkdir(parents=True, exist_ok=True)
    
    l2c_config = {
        'time_mapping': {
            'method': 'Pseudotime to real-time interpolation',
            'reference': 'GSE23160 Bulk time points (2h, 8h, 24h)',
            'mapping_table': mapping_table,
            'notes': [
                'Pseudotime mapped to real hours using linear interpolation',
                'Homeostatic=0h, M2=2h, M1=8h, DAM=24h',
                'Use real_time_hours column in L2c NeuralODE'
            ]
        },
        'l2c_usage': {
            'time_column': 'real_time_hours',
            'velocity_field_init': 'Use GSE23160 temporal logFC as initial condition',
            'cell_state_order': ['Homeostatic', 'M2', 'M1', 'DAM'],
            'time_range': [0, 24],
            'time_unit': 'hours'
        }
    }
    
    with open(l2c_dir / 'l2c_real_time_config.json', 'w', encoding='utf-8') as f:
        json.dump(l2c_config, f, indent=2, ensure_ascii=False)
    print(f"  ✓ l2c_real_time_config.json")
    
    adata_mapped.write_h5ad(l2c_dir / 'GSE174574_with_real_time.h5ad')
    print(f"  ✓ GSE174574_with_real_time.h5ad")
    
    cell_data = adata_mapped.obs[['cell_type', 'pseudotime', 'real_time_hours', 'condition']].copy()
    cell_data.to_csv(l2c_dir / 'cell_time_mapping.csv')
    print(f"  ✓ cell_time_mapping.csv")
    
    return l2c_config

if __name__ == '__main__':
    result = map_pseudotime_to_real_time()
    if result:
        adata, mapping_table, pseudotime_stats = result
        
        mapping_func = create_pseudotime_mapping_function(mapping_table)
        
        adata_mapped = apply_mapping_to_cells(adata, mapping_func)
        
        l2c_config = prepare_l2c_with_real_time(adata_mapped, mapping_table)
        
        print("\n" + "="*70)
        print("✅ 伪时间→真实时间映射完成!")
        print("="*70)
        print(f"\n输出文件:")
        print(f"  - {OUTPUT_DIR / 'l2c_interface' / 'l2c_real_time_config.json'}")
        print(f"  - {OUTPUT_DIR / 'l2c_interface' / 'GSE174574_with_real_time.h5ad'}")
        print(f"  - {OUTPUT_DIR / 'l2c_interface' / 'cell_time_mapping.csv'}")
        print(f"\nL2c使用方式:")
        print(f"  1. 加载 GSE174574_with_real_time.h5ad")
        print(f"  2. 使用 real_time_hours 列作为时间变量 t")
        print(f"  3. NeuralODE求解: dy/dt = f(y, t), t ∈ [0, 24] hours")
