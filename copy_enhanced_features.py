import os
import shutil
from datetime import datetime

gat_dir = r"C:\Users\Jy-Mentor-7\Desktop\GAT"
src = r"D:\反向网络药理学\GAT拓展维度\enhanced_gene_features.csv"
dst = os.path.join(gat_dir, "subgraph_embeddings.csv")
ts = datetime.now().strftime("%Y%m%d_%H%M%S")
backup = os.path.join(gat_dir, f"subgraph_embeddings_backup_{ts}.csv")

if os.path.exists(dst):
    shutil.copy2(dst, backup)
    print(f"Backed up to: {backup}")

shutil.copy2(src, dst)
print(f"Copied enhanced features to: {dst}")

import pandas as pd
df = pd.read_csv(dst, nrows=2)
print(f"New file cols: {len(df.columns)} (expect 1073 = 1 + 1072)")
print(f"First gene: {df.iloc[0, 0]}")