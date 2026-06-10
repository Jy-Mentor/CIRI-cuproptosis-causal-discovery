#!/usr/bin/env python3
import scanpy as sc

# 检查多个 h5ad 文件的基因数
h5ad_files = [
    r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\GSE174574_analysis\results\GSE174574_processed.h5ad",
    r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\GSE174574_analysis\l2c_interface\GSE174574_with_gpr_time.h5ad",
    r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\GSE174574_analysis\l2c_interface\GSE174574_with_real_time.h5ad",
    r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\results\cuproptosis_singlecell\sc_adata_cuproptosis.h5ad",
]

for f in h5ad_files:
    try:
        adata = sc.read_h5ad(f)
        print(f"{f.split('/')[-1]}: {adata.shape}")
        # 检查几个铜死亡基因
        test_genes = ["Fdx1", "Lias", "Dld", "Gls", "Cdkn2a"]
        found = [g for g in test_genes if g in adata.var_names]
        print(f"  Found test genes: {found}")
    except Exception as e:
        print(f"{f.split('/')[-1]}: ERROR - {e}")
    print()
