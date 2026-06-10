"""
GSE174574 - L2c 模块数据接口
准备用于 L2c NeuralODE 模块的数据和参数
"""

import scanpy as sc
import numpy as np
import pandas as pd
from pathlib import Path
import json

DATA_FILE = Path(r'C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\GSE174574_analysis\results\GSE174574_processed.h5ad')
OUTPUT_DIR = Path(r'C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\GSE174574_analysis')

def prepare_l2c_data():
    print("="*70)
    print("GSE174574 - L2c NeuralODE 模块数据准备")
    print("="*70)
    
    adata = sc.read_h5ad(DATA_FILE)
    print(f"\n数据加载: {adata.n_obs} cells × {adata.n_vars} genes")
    
    cell_types = adata.obs['cell_type'].value_counts()
    print(f"\n细胞类型:")
    for ct, count in cell_types.items():
        print(f"  {ct}: {count} ({count/adata.n_obs*100:.1f}%)")
    
    l2c_dir = OUTPUT_DIR / 'l2c_interface'
    l2c_dir.mkdir(parents=True, exist_ok=True)
    
    cell_type_markers = {
        'DAM': ['Trem2', 'Apoe', 'Cst7', 'Lpl', 'Spp1', 'Ccl12'],
        'Homeostatic': ['P2ry12', 'Tmem119', 'Hexb', 'Sall1', 'Olfml3'],
        'M1': ['Nos2', 'Il1b', 'Tnf', 'Cxcl10', 'Il6'],
        'M2': ['Arg1', 'Cd163', 'Mrc1', 'Il10', 'Tgfb1']
    }
    
    detected_markers = {}
    for ct, genes in cell_type_markers.items():
        found = [g for g in genes if g in adata.var_names]
        detected_markers[ct] = found
    
    print(f"\n细胞标记基因检测:")
    for ct, genes in detected_markers.items():
        status = f"{len(genes)}/{len(cell_type_markers[ct])} detected"
        print(f"  {ct}: {status} - {genes}")
    
    print(f"\n生成L2c数据文件...")
    
    adata_subset = adata[adata.obs['cell_type'].isin(['DAM', 'Homeostatic'])].copy()
    
    transition_genes = []
    for gene in adata_subset.var_names:
        dam_mask = adata_subset.obs['cell_type'] == 'DAM'
        homeo_mask = adata_subset.obs['cell_type'] == 'Homeostatic'
        
        if dam_mask.sum() > 0 and homeo_mask.sum() > 0:
            dam_expr = adata_subset[dam_mask, gene].X.mean()
            homeo_expr = adata_subset[homeo_mask, gene].X.mean()
            
            if hasattr(dam_expr, 'toarray'):
                dam_expr = dam_expr.toarray().flatten().mean()
            if hasattr(homeo_expr, 'toarray'):
                homeo_expr = homeo_expr.toarray().flatten().mean()
            
            fold_change = abs(dam_expr - homeo_expr) / (homeo_expr + 1e-10)
            if fold_change > 0.5:
                transition_genes.append({
                    'gene': gene,
                    'dam_expression': float(dam_expr),
                    'homeostatic_expression': float(homeo_expr),
                    'fold_change': float(fold_change)
                })
    
    transition_genes = sorted(transition_genes, key=lambda x: x['fold_change'], reverse=True)
    top_50_genes = transition_genes[:50]
    
    l2c_config = {
        'metadata': {
            'source_dataset': 'GSE174574',
            'species': 'Mouse',
            'platform': 'GPL21103 (10x Genomics)',
            'total_cells': int(adata.n_obs),
            'cell_types': {k: int(v) for k, v in adata.obs['cell_type'].value_counts().to_dict().items()},
            'analysis_date': '2026-05-22'
        },
        'temporal_calibration': {
            'method': 'GSE23160_Bulk_timeseries',
            'reason': 'Pseudotime represents state ordering, not actual reperfusion time',
            'cell_states': ['Homeostatic', 'M2', 'M1', 'DAM'],
            'transition_path': 'Homeostatic -> M2/M1 -> DAM'
        },
        'transition_genes': top_50_genes,
        'cuproptosis_status': {
            'genes_in_hvg': 0,
            'total_genes': 16,
            'note': 'Cuproptosis genes not in top 2000 HVGs. Use full expression matrix for cuproptosis analysis.'
        },
        'pseudotime_info': {
            'available': 'pseudotime' in adata.obs.columns,
            'interpretation': 'State ordering only - DO NOT use as actual time variable',
            'usage': 'Visualization and trajectory validation only'
        }
    }
    
    with open(l2c_dir / 'l2c_config.json', 'w', encoding='utf-8') as f:
        json.dump(l2c_config, f, indent=2, ensure_ascii=False)
    print(f"  ✓ l2c_config.json")
    
    expression_matrix = adata.to_df()
    expression_matrix.to_csv(l2c_dir / 'expression_matrix.csv')
    print(f"  ✓ expression_matrix.csv ({expression_matrix.shape[0]} × {expression_matrix.shape[1]})")
    
    cell_annotations = adata.obs.copy()
    cell_annotations.to_csv(l2c_dir / 'cell_annotations.csv')
    print(f"  ✓ cell_annotations.csv")
    
    umap_coords = adata.obsm['X_umap']
    umap_df = pd.DataFrame(umap_coords, columns=['UMAP1', 'UMAP2'], index=adata.obs_names)
    umap_df['cell_type'] = adata.obs['cell_type']
    umap_df['condition'] = adata.obs['condition'] if 'condition' in adata.obs.columns else 'unknown'
    umap_df.to_csv(l2c_dir / 'umap_coordinates.csv')
    print(f"  ✓ umap_coordinates.csv")
    
    print(f"\n{'='*70}")
    print(f"✅ L2c 数据准备完成!")
    print(f"{'='*70}")
    print(f"\n输出文件 (位于 {l2c_dir}):")
    print(f"  1. l2c_config.json     - L2c模块配置参数")
    print(f"  2. expression_matrix.csv - 基因表达矩阵")
    print(f"  3. cell_annotations.csv  - 细胞注释")
    print(f"  4. umap_coordinates.csv  - UMAP坐标")
    
    print(f"\nL2c 使用建议:")
    print(f"  - 使用 GSE23160 Bulk 时序数据进行时间校准")
    print(f"  - 伪时间仅用于可视化，不作为实际时间变量")
    print(f"  - 关注 Homeostatic -> DAM 状态转换")
    print(f"  - 前50个转换基因已提取")
    
    return l2c_config

if __name__ == '__main__':
    prepare_l2c_data()
