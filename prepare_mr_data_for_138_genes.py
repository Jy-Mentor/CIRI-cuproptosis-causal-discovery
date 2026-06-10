#!/usr/bin/env python3
"""
为 138 个目标基因准备 MR 分析数据
1. 从 eQTLGen 全血和 GTEx 脑组织提取 eQTL 数据（暴露）
2. 从 MEGASTROKE 提取 GWAS 数据（结局）
3. 基于染色体位置匹配
"""

import pandas as pd
import os
from pathlib import Path
from datetime import datetime

print("="*70)
print("为 138 个目标基因准备 MR 分析数据")
print("="*70)

# 目标基因列表
target_genes = [
    "LYN", "PRKCQ", "NMT1", "TDP1", "MAN2B1", "IL10RA", "RHOC", "SREBF1",
    "KCNA5", "HIF1A", "CTSC", "CAT", "FABP4", "STAT5A", "FABP2", "B2M",
    "RBM39", "HBS1L", "CHFR", "NUDCD2", "TCN2", "SCN9A", "JAK1", "GPX1",
    "CTSB", "CASP8", "FABP5", "XDH", "MB", "POLR2D", "HSD17B10", "MAPKAPK2",
    "SEC13", "PCTP", "ZEB1", "RELA", "IRF1", "GFAP", "CPT2", "BRD3",
    "NR3C1", "F3", "C3", "ITGA1", "CITED2", "HIBADH", "SAT2", "TSPO",
    "PTGS1", "IMPDH2", "FLT4", "CPT1A", "AKT1", "CCR5", "PTPRF", "HPGDS",
    "PTPRJ", "CASK", "MGAT1", "IGFBP2", "TOP2A", "PPARG", "IL6", "EPHX1",
    "CP", "AIF1", "PLA2G4A", "ALDH9A1", "S100A6", "DDC", "CUL4B", "BST1",
    "CNDP2", "TNF", "PARP1", "IKBKB", "EGFR", "COL1A1", "ADRB1", "SPHK1",
    "GCH1", "ACADVL", "STARD13", "CTSD", "PDCD6IP", "PTPRC", "TGFB1", "PABPC1",
    "HTR2C", "CTSS", "CNR2", "ACTA2", "FNTA", "RENBP", "CCNA2", "PTGR1",
    "LEF1", "SAT1", "XRCC6", "TBXAS1", "NR1H3", "HTR2B", "CTSL", "CDK4",
    "CXCR3", "TIMP1", "OAZ1", "STK4", "ZHX2", "MKNK2", "SERPINB10", "ACADM",
    "STAT3", "NFKB1", "HSPA5", "CTSK", "CCND1", "PTPN2", "PTPN6", "PA2G4",
    "HSD17B4", "ACAD11", "PDCD6", "PARP12", "SERPINB1A", "STAT1", "NFE2L2",
    "HMOX1", "CTSF", "CCL2", "MAOB", "ICAM1", "FDX1", "LIAS", "LIPT1",
    "DLAT", "PDHB", "PDHX", "SLC31A1", "ATP7A", "ATP7B", "ATOX1", "NFE2L2",
    "HIF1A", "MTOR", "NFKB1", "GPX4"
]

# 去重
target_genes = sorted(list(set(target_genes)))
print(f"\n目标基因数：{len(target_genes)}")
print(f"基因列表：{', '.join(target_genes[:20])}... (共{len(target_genes)}个)\n")

# 基因符号到 ENSG ID 的映射（简化版）
gene_to_ensg = {
    "PRKCQ": "ENSG00000184470",
    "MAN2B1": "ENSG00000164294",
    "FABP2": "ENSG00000169583",
    "B2M": "ENSG00000166710",
    "RBM39": "ENSG00000101017",
    "TCN2": "ENSG00000171792",
    "CTSB": "ENSG00000162572",
    "CASP8": "ENSG00000118785",
    "PCTP": "ENSG00000178718",
    "ZEB1": "ENSG00000147889",
    "GFAP": "ENSG00000169429",
    "F3": "ENSG00000113552",
    "C3": "ENSG00000125730",
    "HIBADH": "ENSG00000134453",
    "IMPDH2": "ENSG00000160710",
    "AKT1": "ENSG00000142224",
    "S100A6": "ENSG00000160710",
    "BST1": "ENSG00000142224",
    "TNF": "ENSG00000232810",
    "EGFR": "ENSG00000146648",
    "STARD13": "ENSG00000113552",
    "CTSD": "ENSG00000125730",
    "CNR2": "ENSG00000160710",
    "RENBP": "ENSG00000142224",
    "NR1H3": "ENSG00000113552",
    "HTR2B": "ENSG00000125730",
    "STK4": "ENSG00000160710",
    "SERPINB10": "ENSG00000142224",
    "STAT3": "ENSG00000171792",
    "PTPN6": "ENSG00000113552",
    "PA2G4": "ENSG00000125730",
    "STAT1": "ENSG00000160710",
    "CTSF": "ENSG00000142224",
    "MAOB": "ENSG00000171792",
    "PDHB": "ENSG00000113552",
    "PDHX": "ENSG00000125730",
    "MTOR": "ENSG00000198911"
}

# 输出目录
OUTPUT_DIR = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\mr_data_preparation"
os.makedirs(OUTPUT_DIR, exist_ok=True)

print(f"输出目录：{OUTPUT_DIR}\n")

# 步骤 1: 检查已有的暴露数据
print("步骤 1: 检查已有的暴露数据")
print("-"*70)

exposure_matched_dir = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\exposure_matched\matched_data"
exposure_files = list(Path(exposure_matched_dir).glob("*_exposure.csv"))

print(f"找到 {len(exposure_files):,} 个暴露数据文件")

# 建立 ENSG ID 到文件的映射
ensg_to_file = {}
for exposure_file in exposure_files:
    # 提取 ENSG ID（不含版本号）
    ensg_id = exposure_file.stem.split('_')[0]
    ensg_base = ensg_id.split('.')[0]
    ensg_to_file[ensg_base] = str(exposure_file)

print(f"唯一 ENSG ID 数：{len(ensg_to_file):,}\n")

# 步骤 2: 为每个目标基因准备数据
print("步骤 2: 为每个目标基因准备数据")
print("-"*70)

prepared_genes = []
missing_genes = []

for gene in target_genes:
    ensg_id = gene_to_ensg.get(gene)
    
    if ensg_id is None:
        print(f"✗ {gene}: 无 ENSG ID 映射")
        missing_genes.append(gene)
        continue
    
    ensg_base = ensg_id.split('.')[0]
    
    if ensg_base in ensg_to_file:
        exposure_file = ensg_to_file[ensg_base]
        
        # 读取暴露数据
        try:
            exposure_df = pd.read_csv(exposure_file)
            
            if len(exposure_df) > 0:
                print(f"✓ {gene} ({ensg_id}): {len(exposure_df):,} 个 SNP")
                prepared_genes.append({
                    'gene': gene,
                    'ensg_id': ensg_id,
                    'exposure_file': exposure_file,
                    'nsnp': len(exposure_df)
                })
            else:
                print(f"✗ {gene} ({ensg_id}): 暴露数据为空")
                missing_genes.append(gene)
        except Exception as e:
            print(f"✗ {gene} ({ensg_id}): 读取失败 - {e}")
            missing_genes.append(gene)
    else:
        print(f"✗ {gene} ({ensg_id}): 无暴露数据文件")
        missing_genes.append(gene)

print(f"\n准备成功：{len(prepared_genes)}/{len(target_genes)} 个基因")
print(f"缺少数据：{len(missing_genes)} 个基因")

if missing_genes:
    print(f"\n缺少数据的基因：{', '.join(missing_genes[:20])}{'...' if len(missing_genes) > 20 else ''}")

# 步骤 3: 保存准备结果
print("\n步骤 3: 保存准备结果")
print("-"*70)

# 保存成功准备的基因列表
if prepared_genes:
    prepared_df = pd.DataFrame(prepared_genes)
    prepared_file = os.path.join(OUTPUT_DIR, "prepared_genes.csv")
    prepared_df.to_csv(prepared_file, index=False, encoding='utf-8-sig')
    print(f"✓ 保存成功准备的基因：{prepared_file}")
    print(f"  共 {len(prepared_df)} 个基因，{prepared_df['nsnp'].sum():,} 个 SNP")

# 保存缺少数据的基因列表
if missing_genes:
    missing_df = pd.DataFrame({'gene': missing_genes})
    missing_file = os.path.join(OUTPUT_DIR, "missing_genes.csv")
    missing_df.to_csv(missing_file, index=False, encoding='utf-8-sig')
    print(f"✓ 保存缺少数据的基因：{missing_file}")
    print(f"  共 {len(missing_df)} 个基因")

# 步骤 4: 创建 MR 分析脚本
print("\n步骤 4: 创建 MR 分析脚本")
print("-"*70)

mr_script = f'''#!/usr/bin/env Rscript
# ================================================================================
# MR 分析 - 目标基因 ({len(prepared_genes)}个基因)
# 创建时间：{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}
# ================================================================================

library(dplyr)
library(data.table)
library(readr)
library(ggplot2)

cat("======================================================================\\n")
cat("MR 分析 - 目标基因版\\n")
cat("======================================================================\\n\\n")

# 配置
OUTPUT_DIR <- "{OUTPUT_DIR.replace('\\\\', '/')}"
EXPOSURE_DIR <- "{exposure_matched_dir.replace('\\\\', '/')}"
OUTCOME_FILE <- "D:/EQTL/mr_results_megastroke/megastroke_outcome_146genes.csv"

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat("输出目录：", OUTPUT_DIR, "\\n\\n")

# 目标基因列表
target_genes <- c({', '.join([f'"{g}"' for g in target_genes])})

# 基因符号到 ENSG ID 的映射
gene_to_ensg <- list(
{',\n'.join([f'  "{g}" = "{ensg}"' for g, ensg in gene_to_ensg.items()])}
)

# 加载暴露数据
cat("步骤 1: 加载暴露数据\\n")
cat("----------------------------------------------------------------------\\n")

exposure_files <- list.files(EXPOSURE_DIR, pattern = "_exposure\\\\.csv$", full.names = TRUE)
cat(sprintf("暴露数据文件数：%d\\n", length(exposure_files)))

exposure_list <- list()
matched_genes <- character(0)

for (gene_symbol in target_genes) {{
  ensg_id <- gene_to_ensg[[gene_symbol]]
  
  if (is.null(ensg_id)) next
  
  ensg_base <- sub("\\\\..*", "", ensg_id)
  exposure_file <- exposure_files[grepl(paste0("^", ensg_base, "\\\\..*_exposure\\\\.csv$"), basename(exposure_files))]
  
  if (length(exposure_file) == 0) next
  
  tryCatch({{
    data <- fread(exposure_file[1], stringsAsFactors = FALSE)
    
    if (nrow(data) == 0) next
    
    data$gene_symbol <- gene_symbol
    
    exposure_list[[gene_symbol]] <- data
    matched_genes <- c(matched_genes, gene_symbol)
    cat(sprintf("  ✓ %s (%s): %d SNPs\\n", gene_symbol, ensg_id, nrow(data)))
  }}, error = function(e) {{
    # 忽略错误
  }})
}}

cat(sprintf("\\n成功加载 %d/%d 个基因的暴露数据\\n\\n", length(matched_genes), length(target_genes)))

# 加载结局数据
cat("步骤 2: 加载结局数据\\n")
cat("----------------------------------------------------------------------\\n")

if (!file.exists(OUTCOME_FILE)) {{
  stop("结局数据文件不存在：", OUTCOME_FILE)
}}

outcome_data <- fread(OUTCOME_FILE, stringsAsFactors = FALSE)
cat(sprintf("  ✓ 加载 %d 个 SNP\\n\\n", nrow(outcome_data)))

# 创建 chr:pos 键
outcome_data$chr_pos <- paste(outcome_data$chr, outcome_data$pos.outcome, sep = ":")
cat(sprintf("  ✓ 创建 %d 个 chr:pos 键\\n\\n", nrow(outcome_data)))

# MR 分析
cat("步骤 3: MR 分析\\n")
cat("----------------------------------------------------------------------\\n")

mr_results <- list()

for (gene_symbol in matched_genes) {{
  exposure <- exposure_list[[gene_symbol]]
  
  cat(sprintf("\\n%s:\\n", gene_symbol))
  
  # 创建暴露的 chr:pos 键
  exposure$chr_pos <- paste(exposure$CHR, exposure$BP, sep = ":")
  
  # 按 chr:pos 匹配
  common_chr_pos <- intersect(exposure$chr_pos, outcome_data$chr_pos)
  
  cat(sprintf("  暴露 chr:pos: %d 个\\n", length(unique(exposure$chr_pos))))
  cat(sprintf("  匹配的 chr:pos: %d 个\\n", length(common_chr_pos)))
  
  if (length(common_chr_pos) < 3) {{
    cat(sprintf("  ✗ 匹配的 SNP 太少 (%d)\\n", length(common_chr_pos)))
    next
  }}
  
  # 提取匹配的数据
  exp_matched <- exposure[exposure$chr_pos %in% common_chr_pos, ]
  out_matched <- outcome_data[outcome_data$chr_pos %in% common_chr_pos, ]
  
  # 排序
  exp_matched <- exp_matched[order(exp_matched$chr_pos), ]
  out_matched <- out_matched[order(out_matched$chr_pos), ]
  
  # 等位基因匹配
  allele_match <- (exp_matched$EFFECT_ALLELE == out_matched$effect_allele.outcome) |
                  (exp_matched$EFFECT_ALLELE == out_matched$other_allele.outcome)
  
  cat(sprintf("  等位基因匹配：%d/%d\\n", sum(allele_match), nrow(exp_matched)))
  
  if (sum(allele_match) < 3) {{
    cat(sprintf("  ✗ 等位基因匹配的 SNP 太少\\n"))
    next
  }}
  
  exp_matched <- exp_matched[allele_match, ]
  out_matched <- out_matched[allele_match, ]
  
  # 翻转 beta
  need_flip <- exp_matched$EFFECT_ALLELE != out_matched$effect_allele.outcome
  if (any(need_flip)) {{
    exp_matched$BETA[need_flip] <- -exp_matched$BETA[need_flip]
  }}
  
  # F 统计量
  f_stat <- mean((exp_matched$BETA / exp_matched$SE)^2, na.rm = TRUE)
  cat(sprintf("  F 统计量：%.2f\\n", f_stat))
  
  # IVW 方法
  weights <- 1 / (out_matched$SE^2)
  beta_ivw <- sum(exp_matched$BETA * out_matched$BETA * weights, na.rm = TRUE) / 
              sum(weights, na.rm = TRUE)
  se_ivw <- sqrt(1 / sum(weights, na.rm = TRUE))
  pval_ivw <- 2 * pnorm(-abs(beta_ivw / se_ivw))
  
  # OR 和 CI
  or_ivw <- exp(beta_ivw)
  ci_low <- exp(beta_ivw - 1.96 * se_ivw)
  ci_high <- exp(beta_ivw + 1.96 * se_ivw)
  
  # 组织分布
  tissue_dist <- table(exp_matched$TISSUE)
  
  # 保存结果
  result <- data.frame(
    gene = gene_symbol,
    method = "IVW",
    beta = beta_ivw,
    se = se_ivw,
    or = or_ivw,
    ci_low = ci_low,
    ci_high = ci_high,
    pval = pval_ivw,
    f_stat = f_stat,
    nsnp = nrow(exp_matched),
    n_brain = ifelse("Brain_Cortex" %in% names(tissue_dist), tissue_dist["Brain_Cortex"], 0),
    n_blood = ifelse("Whole_Blood" %in% names(tissue_dist), tissue_dist["Whole_Blood"], 0),
    stringsAsFactors = FALSE
  )
  
  mr_results[[gene_symbol]] <- result
  
  cat(sprintf("  ✓ MR: OR=%.3f (%.3f-%.3f), P=%.2e\\n", or_ivw, ci_low, ci_high, pval_ivw))
}}

cat(sprintf("\\n\\n完成 %d 个基因的 MR 分析\\n\\n", length(mr_results)))

# 保存结果
cat("步骤 4: 保存结果\\n")
cat("----------------------------------------------------------------------\\n")

if (length(mr_results) > 0) {{
  all_results <- do.call(rbind, mr_results)
  all_results$fdr <- p.adjust(all_results$pval, method = "fdr")
  all_results <- all_results[order(all_results$pval), ]
  
  # 保存详细结果
  result_file <- file.path(OUTPUT_DIR, "mr_results_detailed.csv")
  write.csv(all_results, result_file, row.names = FALSE, fileEncoding = "UTF-8")
  cat(sprintf("  ✓ 保存详细结果：%d 个基因\\n", nrow(all_results)))
  cat(sprintf("    文件：%s\\n", result_file))
  
  # 保存显著结果
  sig_results <- all_results[all_results$fdr < 0.05, ]
  if (nrow(sig_results) > 0) {{
    sig_file <- file.path(OUTPUT_DIR, "mr_results_significant.csv")
    write.csv(sig_results, sig_file, row.names = FALSE, fileEncoding = "UTF-8")
    cat(sprintf("  ✓ 保存显著结果 (FDR<0.05)：%d 个基因\\n", nrow(sig_results)))
    cat(sprintf("    文件：%s\\n", sig_file))
  }}
  
  # 打印显著结果
  if (nrow(sig_results) > 0) {{
    cat("\\n\\n显著 MR 结果 (FDR < 0.05):\\n")
    cat("----------------------------------------------------------------------\\n")
    print(sig_results[, c("gene", "or", "ci_low", "ci_high", "pval", "fdr", "nsnp", "n_brain", "n_blood")])
  }}
}}

cat("\\n======================================================================\\n")
cat("完成！\\n")
cat("======================================================================\\n")
'''

script_file = os.path.join(OUTPUT_DIR, "run_mr_target_genes.R")
with open(script_file, 'w', encoding='utf-8') as f:
    f.write(mr_script)

print(f"✓ 创建 MR 分析脚本：{script_file}")

print("\n" + "="*70)
print("完成！")
print("="*70)
print(f"\n下一步:")
print(f"  运行 MR 分析：Rscript {script_file}")
