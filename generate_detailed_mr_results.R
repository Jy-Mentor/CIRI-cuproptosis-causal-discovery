#!/usr/bin/env Rscript
# ================================================================================
# 生成详细的 MR 结果表格（包含功能注释、异质性、多效性检验）
# ================================================================================

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
  library(dplyr)
})

cat("======================================================================\n")
cat("生成详细的 MR 结果表格\n")
cat("======================================================================\n\n")

# 配置
RESULTS_DIR <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/mr_results_138genes_integrated"
OUTPUT_FILE <- file.path(RESULTS_DIR, "detailed_mr_results_table.csv")

# 加载汇总结果
cat("加载 MR 结果...\n")
summary_results <- fread(file.path(RESULTS_DIR, "summary_results.csv"))
cat(sprintf("  加载了 %d 个基因的结果\n\n", nrow(summary_results)))

# 初始化结果表格
detailed_results <- data.frame(
  gene = character(),
  function_annotation = character(),
  or_ivw = numeric(),
  pval_ivw = numeric(),
  nsnp = integer(),
  heterogeneity_test = character(),
  pleiotropy_test = character(),
  robustness = character(),
  stringsAsFactors = FALSE
)

cat("处理每个基因的敏感性分析...\n\n")

for (i in 1:nrow(summary_results)) {
  gene <- summary_results$gene[i]
  cat(sprintf("[%d/%d] %s: ", i, nrow(summary_results), gene))
  
  # 读取 harmonised 数据
  harmonised_file <- file.path(RESULTS_DIR, paste0(gene, "_harmonised.csv"))
  if (!file.exists(harmonised_file)) {
    cat("无 harmonised 数据\n")
    next
  }
  
  dat <- fread(harmonised_file)
  
  # 读取 MR 结果
  mr_file <- file.path(RESULTS_DIR, paste0(gene, "_mr_results.csv"))
  mr_res <- fread(mr_file)
  
  # 提取 IVW 结果
  ivw_res <- mr_res[mr_res$method == "Inverse variance weighted", ]
  if (nrow(ivw_res) == 0) {
    cat("无 IVW 结果\n")
    next
  }
  
  # 1. OR 和 P 值
  or_ivw <- exp(ivw_res$b[1])
  pval_ivw <- ivw_res$pval[1]
  nsnp <- ivw_res$nsnp
  
  # 2. 异质性检验 (Cochran's Q)
  if (nrow(dat) >= 3) {
    het_res <- tryCatch({
      mr_heterogeneity(dat)
    }, error = function(e) {
      NULL
    })
    
    if (!is.null(het_res) && nrow(het_res) > 0) {
      q_pval <- het_res$Q_pval[1]
      if (is.na(q_pval) || q_pval > 0.05) {
        heterogeneity <- "无显著异质性"
      } else {
        heterogeneity <- "存在异质性"
      }
    } else {
      heterogeneity <- "NA"
    }
  } else {
    heterogeneity <- "SNP 不足"
  }
  
  # 3. 多效性检验 (MR-Egger 截距)
  if (nrow(dat) >= 3) {
    egger_res <- tryCatch({
      mr_pleiotropy_test(dat)
    }, error = function(e) {
      NULL
    })
    
    if (!is.null(egger_res) && nrow(egger_res) > 0) {
      egger_pval <- egger_res$pval[1]
      if (is.na(egger_pval) || egger_pval > 0.05) {
        pleiotropy <- "无显著多效性"
      } else {
        pleiotropy <- "存在多效性"
      }
    } else {
      pleiotropy <- "NA"
    }
  } else {
    pleiotropy <- "SNP 不足"
  }
  
  # 4. 结果稳健性评估
  # 检查不同方法的一致性
  methods <- c("Inverse variance weighted", "Weighted median", "MR-Egger")
  available_methods <- mr_res[mr_res$method %in% methods, ]
  
  if (nrow(available_methods) >= 2) {
    # 检查方向一致性
    signs <- sign(available_methods$b)
    if (all(signs == signs[1])) {
      direction_consistent <- TRUE
    } else {
      direction_consistent <- FALSE
    }
    
    # 检查显著性一致性
    if (pval_ivw < 0.05) {
      if (direction_consistent && heterogeneity == "无显著异质性" && pleiotropy == "无显著多效性") {
        robustness <- "稳健"
      } else if (direction_consistent) {
        robustness <- "中等"
      } else {
        robustness <- "不稳健"
      }
    } else {
      robustness <- "不显著"
    }
  } else {
    robustness <- "NA"
  }
  
  # 5. 功能注释（基于基因名的简单注释）
  function_map <- list(
    "ADRB1" = "β-1 肾上腺素能受体，调节心率和心肌收缩力",
    "SREBF1" = "固醇调节元件结合转录因子，调控脂质代谢",
    "ACADVL" = "极长链酰基-CoA 脱氢酶，脂肪酸β氧化",
    "PABPC1" = "多聚腺苷酸结合蛋白，mRNA 稳定性调控",
    "PTPRJ" = "蛋白酪氨酸磷酸酶，细胞信号传导",
    "RHOC" = "Rho C GTP 酶，细胞骨架重组",
    "AIF1" = "同种移植炎症因子 1，炎症反应",
    "CNDP2" = "细胞质二肽酶 2，氨基酸代谢",
    "HSD17B4" = "17β-羟类固醇脱氢酶 4，类固醇代谢",
    "ATOX1" = "抗氧化蛋白 1，铜离子转运",
    "CAT" = "过氧化氢酶，抗氧化防御",
    "PLA2G4A" = "磷脂酶 A2 组 IVA，花生四烯酸释放",
    "LIPT1" = "脂酰转移酶 1，线粒体代谢",
    "FDX1" = "铁氧还蛋白 1，氧化还原反应",
    "PRKCQ" = "蛋白激酶 C θ，T 细胞激活",
    "BRD3" = "溴结构域蛋白 3，染色质调控",
    "PDHB" = "丙酮酸脱氢酶 E1β，糖代谢",
    "TNF" = "肿瘤坏死因子，炎症细胞因子",
    "CTSF" = "组织蛋白酶 F，蛋白质降解",
    "GCH1" = "GTP 环化水解酶 1，四氢生物蝶呤合成",
    "ACTA2" = "α-平滑肌肌动蛋白，细胞收缩",
    "NMT1" = "N-肉豆蔻酰转移酶 1，蛋白修饰",
    "SPHK1" = "鞘氨醇激酶 1，鞘脂代谢",
    "CTSC" = "组织蛋白酶 C，免疫调节",
    "ITGA1" = "整合素α1，细胞粘附",
    "RBM39" = "RNA 结合蛋白 39，RNA 剪接",
    "STAT1" = "信号转导和转录激活因子 1，免疫应答",
    "LIAS" = "硫辛酸合成酶，线粒体代谢",
    "HMOX1" = "血红素加氧酶 1，抗氧化应激",
    "MAPKAPK2" = "MAPK 激活蛋白激酶 2，应激反应",
    "PDCD6IP" = "程序性细胞死亡 6 相互作用蛋白，凋亡调控",
    "CITED2" = "Cbp/p300 相互作用转录激活因子，转录调控",
    "PARP12" = "聚 ADP-核糖聚合酶 12，DNA 修复",
    "ZHX2" = "锌指同源框 2，转录抑制",
    "AKT1" = "蛋白激酶 B，细胞存活和增殖",
    "PDCD6" = "程序性细胞死亡 6，凋亡调控",
    "SAT2" = "多胺乙酰转移酶 2，多胺代谢",
    "SLC31A1" = "溶质载体家族 31 成员 1，铜离子转运",
    "BST1" = "骨髓/淋巴细胞标记物，ADP-核糖环化酶",
    "PTGR1" = "前列腺素还原酶 1，脂质代谢",
    "NFE2L2" = "核因子 E2 相关因子 2，抗氧化反应",
    "CCR5" = "C-C 趋化因子受体 5，免疫细胞迁移",
    "PCTP" = "磷脂酰胆碱转移蛋白，脂质转运",
    "HSPA5" = "热休克蛋白 70 kDa 成员 5，蛋白折叠",
    "PPARG" = "过氧化物酶体增殖物激活受体γ，脂质代谢",
    "EPHX1" = "环氧化物水解酶 1，解毒代谢",
    "GPX4" = "谷胱甘肽过氧化物酶 4，抗氧化防御",
    "TDP1" = "酪氨酰-DNA 磷酸二酯酶 1，DNA 修复",
    "POLR2D" = "RNA 聚合酶 II 亚基 D，转录",
    "SEC13" = "SEC13 同源物，囊泡运输",
    "PARP1" = "聚 ADP-核糖聚合酶 1，DNA 修复",
    "PDHX" = "丙酮酸脱氢酶复合物 X 组分，糖代谢",
    "CNR2" = "大麻素受体 2，免疫调节",
    "TBXAS1" = "血栓烷 A 合酶，血小板聚集",
    "CTSB" = "组织蛋白酶 B，蛋白质降解",
    "STAT3" = "信号转导和转录激活因子 3，免疫应答",
    "FLT4" = "fms 相关酪氨酸激酶 4，血管生成",
    "IL10RA" = "白细胞介素 10 受体α，抗炎信号",
    "SERPINB10" = "丝氨酸蛋白酶抑制剂 B10，蛋白酶抑制",
    "CHFR" = "检查点蛋白 FR，细胞周期调控",
    "HBS1L" = "HBS1 样翻译因子，核糖体循环",
    "TCN2" = "转钴胺素 II，维生素 B12 转运",
    "NR1H3" = "核受体 1H3，胆固醇代谢",
    "PTGS1" = "前列腺素内过氧化物合酶 1，炎症",
    "CCL2" = "C-C 趋化因子配体 2，单核细胞趋化",
    "C3" = "补体成分 3，免疫应答",
    "CTSD" = "组织蛋白酶 D，蛋白质降解",
    "ATP7B" = "ATP 酶铜转运型 7B，铜离子转运",
    "IKBKB" = "IκB 激酶β，NF-κB 激活",
    "CASP8" = "caspase 8，凋亡执行",
    "ICAM1" = "细胞间粘附分子 1，免疫细胞粘附",
    "PTPN2" = "蛋白酪氨酸磷酸酶非受体型 2，信号调控",
    "MTOR" = "哺乳动物雷帕霉素靶蛋白，细胞生长",
    "MGAT1" = "α-1,3-甘露糖基糖蛋白β-1,2-N-乙酰葡糖胺转移酶，糖基化",
    "ALDH9A1" = "醛脱氢酶 9 家族成员 A1，神经递质代谢",
    "NFKB1" = "核因子κB 亚基 1，炎症和免疫",
    "CPT1A" = "肉碱棕榈酰转移酶 1A，脂肪酸氧化",
    "HIF1A" = "缺氧诱导因子 1α，缺氧应答",
    "CTSL" = "组织蛋白酶 L，蛋白质降解",
    "CCND1" = "细胞周期蛋白 D1，细胞周期进程",
    "IL6" = "白细胞介素 6，炎症细胞因子",
    "MKNK2" = "MAPK 相互作用激酶 2，翻译调控",
    "MAN2B1" = "甘露糖苷酶α2B1，糖蛋白降解",
    "SERPINB1" = "丝氨酸蛋白酶抑制剂 B1，蛋白酶抑制",
    "PTPRF" = "蛋白酪氨酸磷酸酶受体型 F，细胞信号",
    "HIBADH" = "3-羟基异丁酸脱氢酶，缬氨酸代谢",
    "SCN9A" = "电压门控钠通道 9，神经兴奋性",
    "ACADM" = "中链酰基-CoA 脱氢酶，脂肪酸氧化",
    "IGFBP2" = "胰岛素样生长因子结合蛋白 2，生长调控",
    "XRCC6" = "X 射线修复交叉互补蛋白 6，DNA 修复",
    "LYN" = "Lyn 酪氨酸激酶，免疫受体信号",
    "TSPO" = "18 kDa 转运蛋白，胆固醇转运",
    "CTSK" = "组织蛋白酶 K，骨吸收",
    "STAT5A" = "信号转导和转录激活因子 5A，细胞因子信号",
    "CPT2" = "肉碱棕榈酰转移酶 2，脂肪酸氧化"
  )
  
  function_annotation <- ifelse(gene %in% names(function_map), 
                                 function_map[[gene]], 
                                 "功能未知")
  
  # 添加到结果表格
  result_row <- data.frame(
    gene = gene,
    function_annotation = function_annotation,
    or_ivw = or_ivw,
    pval_ivw = pval_ivw,
    nsnp = nsnp,
    heterogeneity_test = heterogeneity,
    pleiotropy_test = pleiotropy,
    robustness = robustness,
    stringsAsFactors = FALSE
  )
  
  detailed_results <- rbind(detailed_results, result_row)
  
  cat(sprintf("完成\n"))
}

# 排序（按 P 值）
detailed_results <- detailed_results[order(detailed_results$pval_ivw), ]

# 保存结果
cat("\n保存详细结果...\n")
write.csv(detailed_results, OUTPUT_FILE, row.names = FALSE, fileEncoding = "UTF-8")
cat(sprintf("  保存到：%s\n\n", OUTPUT_FILE))

# 打印摘要
cat("======================================================================\n")
cat("详细结果摘要\n")
cat("======================================================================\n\n")

cat("显著性基因 (P < 0.05):\n")
sig_genes <- detailed_results[detailed_results$pval_ivw < 0.05, ]
if (nrow(sig_genes) > 0) {
  print(sig_genes[, c("gene", "function_annotation", "or_ivw", "pval_ivw", "nsnp", "heterogeneity_test", "pleiotropy_test", "robustness")])
} else {
  cat("  无显著基因\n")
}

cat("\n\n稳健的显著结果 (P < 0.05 且 稳健=稳健):\n")
robust_sig <- detailed_results[detailed_results$pval_ivw < 0.05 & detailed_results$robustness == "稳健", ]
if (nrow(robust_sig) > 0) {
  print(robust_sig[, c("gene", "function_annotation", "or_ivw", "pval_ivw", "nsnp", "heterogeneity_test", "pleiotropy_test", "robustness")])
} else {
  cat("  无稳健的显著结果\n")
}

cat("\n\n完成！\n")
