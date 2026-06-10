"""
重新计算伪时间并映射到真实时间
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

def compute_pseudotime_and_map():
    print("="*70)
    print("重新计算伪时间并映射到真实时间")
    print("="*70)
    
    adata = sc.read_h5ad(DATA_FILE)
    print(f"\n数据: {adata.n_obs} cells × {adata.n_vars} genes")
    print(f"obs列: {adata.obs.columns.tolist()}")
    print(f"细胞类型: {adata.obs['cell_type'].value_counts().to_dict()}")
    
    print(f"\n伪时间统计:")
    for ct in adata.obs['cell_type'].unique():
        mask = adata.obs['cell_type'] == ct
        pt = adata.obs.loc[mask, 'pseudotime'].dropna()
        if len(pt) > 0:
            print(f"  {ct}: mean={pt.mean():.4f}, median={pt.median():.4f}, range=[{pt.min():.4f}, {pt.max():.4f}]")
    
    cell_type_median = {}
    for ct in adata.obs['cell_type'].unique():
        mask = adata.obs['cell_type'] == ct
        cell_type_median[ct] = float(adata.obs.loc[mask, 'pseudotime'].median())
    
    sorted_ct = sorted(cell_type_median.items(), key=lambda x: x[1])
    print(f"\n伪时间排序: {sorted_ct}")
    
    mapping_table = {}
    for ct, real_time in BULK_TIME_POINTS.items():
        if ct in cell_type_median:
            mapping_table[ct] = {
                'pseudotime_median': cell_type_median[ct],
                'real_time_hours': real_time,
                'cell_count': int((adata.obs['cell_type'] == ct).sum())
            }
    
    print(f"\n伪时间→真实时间映射:")
    pseudo_times = []
    real_times = []
    for ct, mapping in mapping_table.items():
        print(f"  {ct}: pseudo={mapping['pseudotime_median']:.4f} → {mapping['real_time_hours']}h")
        pseudo_times.append(mapping['pseudotime_median'])
        real_times.append(mapping['real_time_hours'])
    
    if len(pseudo_times) >= 2:
        mapping_func = interp1d(pseudo_times, real_times, kind='linear', fill_value='extrapolate')
        
        pseudo_times_all = adata.obs['pseudotime'].values
        real_times_all = mapping_func(pseudo_times_all)
        real_times_all = np.clip(real_times_all, 0, 24)
        
        adata.obs['real_time_hours'] = real_times_all
        
        print(f"\n映射后真实时间:")
        print(f"  范围: {adata.obs['real_time_hours'].min():.2f} - {adata.obs['real_time_hours'].max():.2f} h")
        print(f"  均值: {adata.obs['real_time_hours'].mean():.2f} h")
        
        for ct in adata.obs['cell_type'].unique():
            mask = adata.obs['cell_type'] == ct
            ct_times = adata.obs.loc[mask, 'real_time_hours']
            print(f"  {ct}: {ct_times.mean():.2f} ± {ct_times.std():.2f} h")
    else:
        print("\n⚠️ 映射点不足，无法创建插值函数")
        return
    
    l2c_dir = OUTPUT_DIR / 'l2c_interface'
    l2c_dir.mkdir(parents=True, exist_ok=True)
    
    l2c_config = {
        'time_mapping': {
            'method': 'DPT pseudotime → real hours (linear interpolation)',
            'reference': 'GSE23160 Bulk (0h, 2h, 8h, 24h)',
            'mapping_table': mapping_table,
            'pseudocyte_ordering': [ct for ct, _ in sorted_ct]
        },
        'l2c_usage': {
            'time_column': 'real_time_hours',
            'time_unit': 'hours',
            'time_range': [0, 24],
            'velocity_field_init': 'Use GSE23160 temporal logFC as initial condition',
            'cell_state_path': ['Homeostatic', 'M2', 'M1', 'DAM']
        }
    }
    
    with open(l2c_dir / 'l2c_real_time_config.json', 'w', encoding='utf-8') as f:
        json.dump(l2c_config, f, indent=2, ensure_ascii=False)
    print(f"\n✓ l2c_real_time_config.json")
    
    adata.write_h5ad(l2c_dir / 'GSE174574_with_real_time.h5ad')
    print(f"✓ GSE174574_with_real_time.h5ad")
    
    cell_data = adata.obs[['cell_type', 'pseudotime', 'real_time_hours', 'condition']].copy()
    cell_data.to_csv(l2c_dir / 'cell_time_mapping.csv')
    print(f"✓ cell_time_mapping.csv")
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    sc.pl.umap(adata, color='pseudotime', ax=axes[0], show=False,
               title='DPT Pseudotime', cmap='viridis', size=5)
    
    sc.pl.umap(adata, color='real_time_hours', ax=axes[1], show=False,
               title='Mapped Real Time (hours)', cmap='plasma', size=5)
    
    for ct in adata.obs['cell_type'].unique():
        mask = adata.obs['cell_type'] == ct
        axes[2].hist(adata.obs.loc[mask, 'real_time_hours'], alpha=0.6, label=ct, bins=30)
    axes[2].set_xlabel('Real Time (hours)')
    axes[2].set_ylabel('Cell Count')
    axes[2].set_title('Time Distribution by Cell Type')
    axes[2].legend()
    axes[2].axvline(x=2, color='gray', linestyle='--', alpha=0.5)
    axes[2].axvline(x=8, color='gray', linestyle='--', alpha=0.5)
    axes[2].axvline(x=24, color='gray', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    fig_path = OUTPUT_DIR / 'figures' / 'Figure_Pseudotime_to_RealTime.png'
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✓ Figure_Pseudotime_to_RealTime.png")
    
    print("\n" + "="*70)
    print("✅ 伪时间→真实时间映射完成!")
    print("="*70)
    print(f"\nL2c使用方式:")
    print(f"  1. 加载 GSE174574_with_real_time.h5ad")
    print(f"  2. 使用 real_time_hours 列作为时间变量 t")
    print(f"  3. NeuralODE: dy/dt = f(y, t), t ∈ [0, 24] hours")

if __name__ == '__main__':
    compute_pseudotime_and_map()
