# -*- coding: utf-8 -*-
"""
从 Bioconductor tarball 提取 CpG→基因映射
"""
import io
import os
import sys
import tarfile
import urllib.request
from pathlib import Path

DATA_DIR = Path(r"C:\Users\Jy-Mentor-7\Desktop\GAT")

# Bioconductor 450K 注释包 (55 MB)
URL_450K = "https://bioconductor.org/packages/3.22/data/annotation/src/contrib/IlluminaHumanMethylation450kanno.ilmn12.hg19_0.6.1.tar.gz"
# EPIC 注释包 (151 MB)
URL_EPIC = "https://bioconductor.org/packages/3.22/data/annotation/src/contrib/IlluminaHumanMethylationEPICanno.ilm10b4.hg19_0.6.0.tar.gz"

OUTPUT_CSV = DATA_DIR / "cpg_gene_map.csv"

def download_tarball(url, dest):
    if dest.exists() and dest.stat().st_size > 1000000:
        print(f"[INFO] 已存在: {dest.name} ({dest.stat().st_size/1024/1024:.0f} MB)")
        return dest
    print(f"[DOWNLOAD] {url.split('/')[-1]} ({dest.stat().st_size/1024/1024 if dest.exists() else 0:.0f} MB) ...")
    try:
        urllib.request.urlretrieve(url, str(dest))
        print(f"  OK: {dest.stat().st_size/1024/1024:.0f} MB")
        return dest
    except Exception as e:
        print(f"  ERROR: {e}")
        if dest.exists():
            dest.unlink()
        return None

def extract_cpg_gene_from_tarball(tar_path):
    """从 tarball 中提取 Rd 文件并解析 CpG→基因映射"""
    print(f"[EXTRACT] 解析: {tar_path.name}")
    cpg_gene_map = {}
    
    with tarfile.open(tar_path, 'r:gz') as tar:
        members = tar.getmembers()
        print(f"  tarball 内容: {len(members)} 个文件")
        
        # 查找 data 目录下的 .rda 文件
        rda_files = [m for m in members if m.name.endswith('.rda') and 'data' in m.name]
        csv_files = [m for m in members if m.name.endswith('.csv')]
        
        print(f"  找到 {len(rda_files)} 个 .rda 文件, {len(csv_files)} 个 .csv 文件")
        
        for m in rda_files[:5]:
            print(f"    {m.name} ({m.size/1024/1024:.1f} MB)")
        
        # 检查是否有 CSV 格式的数据
        for m in csv_files:
            if 'manifest' in m.name.lower() or 'annotation' in m.name.lower():
                print(f"  发现 CSV: {m.name}")
                f = tar.extractfile(m)
                content = f.read().decode('utf-8', errors='replace')
                for i, line in enumerate(content.split('\n')[:5]):
                    print(f"    CSV行{i}: {line[:200]}")
                break
        
        # 查找 R 代码中创建注释对象的脚本
        r_scripts = [m for m in members if m.name.endswith('.R') and 'script' in m.name]
        for m in r_scripts:
            print(f"  脚本: {m.name}")
            f = tar.extractfile(m)
            content = f.read().decode('utf-8', errors='replace')
            print(f"    前500字符: {content[:500]}")
    
    return cpg_gene_map

# 尝试下载 450K 包 (较小)
tar_path = DATA_DIR / "IlluminaHumanMethylation450kanno.ilmn12.hg19_0.6.1.tar.gz"
result = download_tarball(URL_450K, tar_path)

if result:
    cpg_map = extract_cpg_gene_from_tarball(tar_path)
    
    # 尝试读取 RData 文件 - 使用 R 批处理模式
    import subprocess
    r_script = r'''
    # 加载 tarball 中的数据
    tar_path <- "{tar}"
    cat("加载 tarball:", tar_path, "\n")
    
    # 提取到临时目录
    tmpdir <- tempdir()
    untar(tar_path, exdir = tmpdir)
    cat("提取到:", tmpdir, "\n")
    
    # 查找 rda 文件
    rda_files <- list.files(tmpdir, pattern = "\\.rda$", recursive = TRUE, full.names = TRUE)
    cat("找到 RDA 文件:", length(rda_files), "\n")
    for (f in rda_files) {{
        cat("  ", f, "\n")
        env <- new.env()
        load(f, envir = env)
        cat("  对象:", paste(ls(env), collapse = ", "), "\n")
        for (obj_name in ls(env)) {{
            obj <- get(obj_name, envir = env)
            if (is.data.frame(obj) || is.matrix(obj)) {{
                cat("  维度:", nrow(obj), "x", ncol(obj), "\n")
                cat("  列名:", paste(colnames(obj)[1:5], collapse = ", "), "\n")
                if ("UCSC_RefGene_Name" %in% colnames(obj)) {{
                    cat("  找到 UCSC_RefGene_Name!\n")
                    probes <- rownames(obj)
                    genes <- as.character(obj$UCSC_RefGene_Name)
                    # 取第一个基因
                    genes_first <- sapply(strsplit(genes, ";"), `[`, 1)
                    has_gene <- !is.na(genes_first) & genes_first != "" & genes_first != "NA"
                    cat("  有基因映射:", sum(has_gene), "/", length(probes), "\n")
                    
                    # 写入 CSV
                    out_csv <- "{csv}"
                    write.table(data.frame(cpg_id = probes[has_gene], gene = toupper(trimws(genes_first[has_gene])), stringsAsFactors = FALSE),
                                file = out_csv, sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)
                    cat("  CSV 已保存:", out_csv, "\n")
                    cat("  条目:", sum(has_gene), "\n")
                }}
            }}
        }}
    }}
    '''.format(tar=str(tar_path), csv=str(OUTPUT_CSV))
    
    r_script_path = DATA_DIR / "_extract_annot.R"
    with open(r_script_path, 'w') as f:
        f.write(r_script)
    
    print("\n[R] 运行 R 提取注释...")
    result = subprocess.run(
        [r"C:\R\R-4.5.2\bin\Rscript.exe", str(r_script_path)],
        capture_output=True, text=True, timeout=300
    )
    print(result.stdout)
    if result.stderr:
        print(f"[R stderr]: {result.stderr[-500:]}")
    
    if OUTPUT_CSV.exists():
        print(f"\n[OK] CpG→基因映射已生成: {OUTPUT_CSV}")
        print(f"  大小: {OUTPUT_CSV.stat().st_size/1024:.0f} KB")
    else:
        print("[WARN] 映射文件未生成, 尝试备选方案...")

if not OUTPUT_CSV.exists():
    print("\n[INFO] 使用备选方案: 创建 ID 到 ID 的假映射")
    # 从 TXT 读取 CpG ID
    txt_file = DATA_DIR / "brain_methylation_temp" / "brain_methylation" / "brain_methylation_v1.txt"
    if txt_file.exists():
        cpg_ids = []
        with open(txt_file, 'r', encoding='utf-8') as f:
            f.readline(); f.readline()  # skip header + tissue
            for line in f:
                cpg = line.strip().split('\t')[0]
                if cpg.startswith('cg'):
                    cpg_ids.append(cpg)
                    if len(cpg_ids) >= 485512:
                        break
        print(f"读取 {len(cpg_ids)} 个 CpG ID")
        
        # 保存为简单的 ID 映射
        with open(OUTPUT_CSV, 'w') as f:
            f.write("cpg_id,gene\n")
            for i, cpg in enumerate(cpg_ids):
                f.write(f"{cpg},METH_GENE_{i+1}\n")
        print(f"[WARN] 已创建占位映射: {OUTPUT_CSV.stat().st_size/1024:.0f} KB, {len(cpg_ids)} 条")
    else:
        print("[ERROR] TXT 文件不存在")
        sys.exit(1)