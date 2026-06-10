"""
GSE174574 从loom文件加载数据 (用于scVelo分析)
"""

import scanpy as sc
import scvelo as scv
import numpy as np
import pandas as pd
from pathlib import Path

LOOM_FILE = Path(r'D:\反向网络药理学\L1 数据集\RNA-seq\GSE174574_merged.velocyto.loom')
OUTPUT_DIR = Path(r'C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\GSE174574_analysis')

def load_loom_and_run_velocity():
    """从loom文件加载并运行scVelo"""
    
    print("从loom文件加载数据...")
    adata = scv.read(str(LOOM_FILE), cache=True)
    
    print(f"数据维度: {adata.n_obs} cells × {adata.n_vars} genes")
    print(f"包含layers: {list(adata.layers.keys())}")
    
    print("\n预处理...")
    scv.pp.filter_and_normalize(adata)
    scv.pp.moments(adata, n_pcs=30, n_cells=30)
    
    print("计算RNA velocity...")
    scv.tl.velocity(adata, mode='stochastic')
    scv.tl.velocity_graph(adata)
    
    print("计算velocity embedding...")
    scv.tl.velocity_embedding(adata, basis='umap')
    
    print("\n保存velocity matrix...")
    output_dir = OUTPUT_DIR / 'results'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    velocity_matrix_file = output_dir / 'velocity_matrix.npy'
    np.save(velocity_matrix_file, adata.obsm['velocity_umap'])
    print(f"Velocity matrix保存: {velocity_matrix_file}")
    
    velocity_adata_file = output_dir / 'adata_velocity.h5ad'
    adata.write(velocity_adata_file)
    print(f"Velocity AnnData保存: {velocity_adata_file}")
    
    print("\n生成Figure 1: Velocity Flow Field...")
    fig_dir = OUTPUT_DIR / 'figures'
    fig_dir.mkdir(parents=True, exist_ok=True)
    
    scv.pl.velocity_embedding_stream(
        adata, basis='umap', color='clusters',
        save=str(fig_dir / 'Figure_1_Velocity_Flow_Field.png'),
        title='RNA Velocity Flow Field'
    )
    
    print("\n✅ 全部完成!")
    print(f"Velocity matrix可用于L2c NeuralODE模块")
    
    return adata

if __name__ == '__main__':
    adata = load_loom_and_run_velocity()
