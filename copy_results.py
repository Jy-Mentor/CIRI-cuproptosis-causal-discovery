import os
import shutil

src = os.path.dirname(os.path.abspath(__file__))
dest = os.path.join(src, "铜死亡CIRI_最终交付")

scripts_dir = os.path.join(dest, "scripts")
results_dir = os.path.join(dest, "results")

os.makedirs(scripts_dir, exist_ok=True)
os.makedirs(results_dir, exist_ok=True)

print("Copying scripts...")
script_files = [
    "cuproptosis_gsva.py",
    "cuproptosis_gsea.py",
    "cuproptosis_singlecell.py",
    "cuproptosis_wgcna.py",
    "cuproptosis_immunology.py",
    "cuproptosis_ppi_neighbors.py",
    "cuproptosis_hallmark_gsva.py",
    "run_cuproptosis_modules.py",
    "export_results_to_excel.py",
    "utils.py",
    "prepare_geo_data.py"
]

for f in script_files:
    src_path = os.path.join(src, "scripts", f)
    dest_path = os.path.join(scripts_dir, f)
    if os.path.exists(src_path):
        shutil.copy2(src_path, dest_path)
        print(f"  OK: {f}")
    else:
        print(f"  MISSING: {f}")

print("Copying config.py...")
shutil.copy2(os.path.join(src, "config.py"), os.path.join(dest, "config.py"))

print("Copying results...")
result_dirs = [
    "cuproptosis_gsva",
    "cuproptosis_gsea",
    "cuproptosis_singlecell",
    "cuproptosis_wgcna",
    "cuproptosis_immunology",
    "cuproptosis_ppi",
    "cuproptosis_hallmark_gsva",
    "cuproptosis_analysis",
    "stage1_rma_degs"
]

total_files = 0
for d in result_dirs:
    src_path = os.path.join(src, "results", d)
    dest_path = os.path.join(results_dir, d)
    if os.path.exists(src_path):
        if os.path.exists(dest_path):
            shutil.rmtree(dest_path)
        shutil.copytree(src_path, dest_path)
        fc = sum([len(files) for _, _, files in os.walk(dest_path)])
        total_files += fc
        print(f"  OK: {d} ({fc} files)")
    else:
        print(f"  MISSING: {d}")

print(f"\nDone! Total files: {total_files}")
print(f"Path: {dest}")
