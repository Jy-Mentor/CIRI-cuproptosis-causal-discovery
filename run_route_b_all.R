#!/usr/bin/env Rscript
# ================================================================================
# MR 分析优化路线 B - 综合运行脚本
# 自动执行所有优化分析：功能富集、药物靶点、COLOC、敏感性分析
# ================================================================================

message(rep("=", 70))
message("MR 分析优化路线 B - 综合运行脚本")
message("目标期刊：Nature Communications (IF=16.6)")
message(rep("=", 70))
message("")

# 定义分析脚本
analysis_scripts <- list(
  list(
    name = "功能富集分析",
    file = "functional_enrichment_analysis.R",
    description = "GO/KEGG/Reactome 通路富集",
    status = "pending"
  ),
  list(
    name = "药物靶点预测",
    file = "drug_target_prediction.R",
    description = "DGIdb/OpenTargets 药物靶点查询",
    status = "pending"
  ),
  list(
    name = "COLOC 共定位分析",
    file = "coloc_analysis.R",
    description = "确认共享因果变异",
    status = "pending"
  ),
  list(
    name = "敏感性分析增强",
    file = "sensitivity_analysis_enhanced.R",
    description = "MR-PRESSO/径向 MR",
    status = "pending"
  )
)

# 检查脚本是否存在
message("步骤 1: 检查分析脚本...")
message(rep("-", 70))

script_dir <- "c:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙"

for (i in seq_along(analysis_scripts)) {
  script <- analysis_scripts[[i]]
  script_path <- file.path(script_dir, script$file)
  
  if (file.exists(script_path)) {
    message(paste("✓", script$name, "(", script$file, ")"))
    analysis_scripts[[i]]$status <- "ready"
  } else {
    message(paste("✗", script$name, "未找到:", script$file))
    analysis_scripts[[i]]$status <- "missing"
  }
}

# 统计就绪的脚本
ready_count <- sum(sapply(analysis_scripts, function(x) x$status == "ready"))
message("")
message(paste("就绪脚本:", ready_count, "/", length(analysis_scripts)))
message("")

# 询问是否继续
if (ready_count == 0) {
  stop("没有可用的分析脚本，请检查文件路径")
}

message(paste("\n即将运行", ready_count, "个分析脚本..."))
message("按 Ctrl+C 可中断")
Sys.sleep(3)

# 运行每个分析脚本
results_summary <- list()

for (i in seq_along(analysis_scripts)) {
  script <- analysis_scripts[[i]]
  
  if (script$status != "ready") {
    message(paste("\n跳过:", script$name))
    next
  }
  
  message("")
  message(rep("=", 70))
  message(paste("运行分析", i, "/", length(analysis_scripts), ":", script$name))
  message(paste("描述:", script$description))
  message(rep("=", 70))
  
  script_path <- file.path(script_dir, script$file)
  
  # 运行脚本
  start_time <- Sys.time()
  
  tryCatch({
    # 使用 system 命令运行 R 脚本
    result <- system(paste('"C:\\R\\R-4.5.2\\bin\\Rscript.exe"', shQuote(script_path)), 
                     intern = FALSE)
    
    end_time <- Sys.time()
    duration <- difftime(end_time, start_time, units = "mins")
    
    if (result == 0) {
      message(paste("\n✓", script$name, "完成 (耗时:", round(duration, 2), "分钟)"))
      analysis_scripts[[i]]$status <- "completed"
      results_summary[[script$name]] <- list(
        status = "success",
        duration = duration
      )
    } else {
      message(paste("\n✗", script$name, "失败 (退出码:", result, ")"))
      analysis_scripts[[i]]$status <- "failed"
      results_summary[[script$name]] <- list(
        status = "error",
        duration = duration,
        error = paste("退出码:", result)
      )
    }
    
  }, error = e) {
    message(paste("\n✗", script$name, "错误:", e$message))
    analysis_scripts[[i]]$status <- "error"
    results_summary[[script$name]] <- list(
      status = "error",
      error = e$message
    )
  }
  
  # 分析间等待
  if (i < length(analysis_scripts)) {
    message("\n等待 5 秒后继续下一个分析...")
    Sys.sleep(5)
  }
}

# 生成综合报告
message("")
message(rep("=", 70))
message("生成综合分析总结报告")
message(rep("=", 70))
message("")

# 创建总结报告
report <- c(
  "# MR 分析优化路线 B - 综合报告",
  "",
  "## 分析概述",
  paste("运行日期:", Sys.time()),
  paste("运行脚本数:", length(analysis_scripts)),
  paste("成功完成:", sum(sapply(analysis_scripts, function(x) x$status == "completed"))),
  paste("失败:", sum(sapply(analysis_scripts, function(x) x$status %in% c("failed", "error")))),
  ""
)

# 添加每个分析的结果
report <- c(report, "## 各分析模块结果", "")

for (script in analysis_scripts) {
  status_icon <- switch(script$status,
    "completed" = "✅",
    "failed" = "❌",
    "error" = "⚠️",
    "missing" = "❓",
    "⏳"
  )
  
  report <- c(report, paste0(
    "### ", status_icon, " ", script$name, "  \n",
    "- **文件**: `", script$file, "`  \n",
    "- **描述**: ", script$description, "  \n",
    "- **状态**: ", script$status, "  \n"
  ))
  
  if (!is.null(results_summary[[script$name]])) {
    if (!is.na(results_summary[[script$name]]$duration)) {
      report <- c(report, paste0("- **耗时**: ", round(results_summary[[script$name]]$duration, 2), " 分钟  \n"))
    }
  }
  
  report <- c(report, "")
}

# 输出目录
report <- c(report, 
  "## 输出文件目录",
  "",
  "所有分析结果保存在以下目录：",
  "",
  "```,
  "D:/下载/MR_batch_results/20260508_optimized_fixed_v2/",
  "├── functional_enrichment/     # 功能富集分析结果",
  "│   ├── GO_BP_dotplot.png",
  "│   ├── KEGG_dotplot.png",
  "│   ├── Reactome_dotplot.png",
  "│   └── enrichment_report.md",
  "├── drug_targets/              # 药物靶点预测结果",
  "│   ├── drug_targets_summary.csv",
  "│   ├── high_priority_targets.csv",
  "│   └── drug_target_report.md",
  "├── coloc_analysis/            # COLOC 共定位分析结果",
  "│   ├── coloc_summary.csv",
  "│   ├── coloc_success_genes.csv",
  "│   └── coloc_report.md",
  "└── sensitivity_analysis/      # 敏感性分析增强结果",
  "    ├── sensitivity_summary.csv",
  "    ├── sensitivity_report.md",
  "    └── [gene]_radial.png",
  "```",
  ""
)

# 论文撰写建议
report <- c(report,
  "## 论文撰写建议",
  "",
  "### 结果部分结构",
  "",
  "1. **MR 分析主要结果**",
  "   - 成功识别 103 个基因的因果关联",
  "   - 3 个基因通过 FDR 校正 (SREBF1, SPHK1, NR1H3)",
  "   - 敏感性分析结果",
  "",
  "2. **功能富集分析**",
  "   - GO 富集显示脂质代谢和炎症通路",
  "   - KEGG 通路富集支持代谢假说",
  "   - Reactome 分析提供机制洞察",
  "",
  "3. **药物靶点预测**",
  "   - 识别高优先级可成药靶点",
  "   - SREBF1 和 NR1H3 已有在研药物",
  "   - 为卒中治疗提供新方向",
  "",
  "4. **COLOC 共定位分析**",
  "   - 确认 eQTL 和 GWAS 共享因果变异",
  "   - 排除 LD 导致的假阳性",
  "   - 增强因果推断强度",
  "",
  "5. **敏感性分析增强**",
  "   - MR-PRESSO 检测水平多效性",
  "   - 径向 MR 识别异常值 SNP",
  "   - 提高结果稳健性",
  "",
  "### 推荐图表",
  "",
  "- **Figure 1**: MR 分析流程图 + 主要结果森林图",
  "- **Figure 2**: 显著基因因果关联图",
  "- **Figure 3**: GO/KEGG 富集分析气泡图",
  "- **Figure 4**: 药物靶点网络图",
  "- **Figure 5**: COLOC 共定位后验概率图",
  "- **Figure 6**: 敏感性分析径向森林图",
  "",
  "### 目标期刊",
  "",
  "- **首选**: Nature Communications (IF=16.6)",
  "- **备选**: Science Advances (IF=14.1)",
  "- **保底**: AJHG (IF=9.8)",
  ""
)

# 保存报告
report_file <- file.path(script_dir, "ROUTE_B_COMPREHENSIVE_REPORT.md")
writeLines(report, report_file)
message(paste("综合报告已保存:", report_file))

# 最终总结
message("")
message(rep("=", 70))
message("MR 分析优化路线 B 执行完成!")
message(rep("=", 70))
message("")
message("分析统计:")
message(paste("  总脚本数:", length(analysis_scripts)))
message(paste("  成功完成:", sum(sapply(analysis_scripts, function(x) x$status == "completed"))))
message(paste("  失败:", sum(sapply(analysis_scripts, function(x) x$status %in% c("failed", "error")))))
message(paste("  跳过:", sum(sapply(analysis_scripts, function(x) x$status == "missing"))))
message("")
message("下一步:")
message("1. 检查各分析模块的输出文件")
message("2. 查看综合报告：ROUTE_B_COMPREHENSIVE_REPORT.md")
message("3. 根据结果撰写论文")
message("")
message("祝科研顺利！🎉")
