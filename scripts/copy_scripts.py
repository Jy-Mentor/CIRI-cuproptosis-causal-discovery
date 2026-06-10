import shutil, os

src = r"c:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙"
dst = os.path.join(src, "pipeline_scripts")

files = [
    os.path.join(src, "config.py"),
    os.path.join(src, "scripts", "stage1_rma_degs.R"),
    os.path.join(src, "scripts", "stage2_single_cell.py"),
    os.path.join(src, "scripts", "stage3_enrichment.py"),
    os.path.join(src, "scripts", "stage4_seed_wgcna.py"),
    os.path.join(src, "scripts", "stage4_wgcna.R"),
    os.path.join(src, "scripts", "stage5_string_ppi.py"),
    os.path.join(src, "scripts", "stage6_grn_knockout.py"),
    os.path.join(src, "scripts", "stage7_ml_shap.py"),
    os.path.join(src, "scripts", "stage8_final_targets.py"),
    os.path.join(src, "scripts", "stage9_ppi_gat.py"),
    os.path.join(src, "scripts", "utils.py"),
]

for f in files:
    if os.path.exists(f):
        shutil.copy2(f, dst)
        print(f"✓ {os.path.basename(f)}")
    else:
        print(f"✗ 不存在: {f}")

print(f"\n已复制 {len([f for f in files if os.path.exists(f)])}/{len(files)} 文件到 {dst}")
