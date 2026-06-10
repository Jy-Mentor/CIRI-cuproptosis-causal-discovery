import pandas as pd
import numpy as np
from pathlib import Path

# Read KO data
ko_df = pd.read_csv('C:/Users/Jy-Mentor-7/Desktop/生物信息学/KO/data/summary_all_KO.csv')
print(f'KO data shape: {ko_df.shape}')
print(f'Columns: {list(ko_df.columns)}')
print(f'Unique genes: {ko_df["gene"].nunique()}')
print(f'Unique cell types: {ko_df["cell_type"].nunique()}')
print(f'Cell types: {sorted(ko_df["cell_type"].unique())}')
print(f'Status values: {ko_df["status"].value_counts().to_dict()}')
print(f'n_sig stats: count={ko_df["n_sig"].notna().sum()}, max={ko_df["n_sig"].max()}')
print(f'n_corr stats: count={ko_df["n_corr"].notna().sum()}, max={ko_df["n_corr"].max()}')
print(f'First 20 rows:')
print(ko_df.head(20).to_string())
