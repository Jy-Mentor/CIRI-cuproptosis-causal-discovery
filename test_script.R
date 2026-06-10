#!/usr/bin/env Rscript

# 设置全局选项
options(stringsAsFactors = FALSE)
set.seed(42)

# 参数区
COR_THRESHOLD <- 0.5
ALPHA_BASE <- 0.05
SAVE_ADJACENCY <- FALSE  # 是否保存邻接矩阵
N_CORES <- 4  # 并行计算核心数

# 设置输出目录
output_dir <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

# 缓存文件路径
cache_file <- file.path(output_dir, "analysis_cache.rds")

# 加载缓存
cat("发现缓存文件，尝试加载...\n")
tryCatch({
  # 确保Seurat包已加载
  if (!require("Seurat", character.only = TRUE)) {
    install.packages("Seurat", dependencies = TRUE)
    library(Seurat)
  }
  cached_data <- readRDS(cache_file)
  # 检查缓存是否包含必要的对象
  if (all(c("neurons", "target_genes", "output_dir") %in% names(cached_data))) {
    cat("缓存加载成功，跳过步骤1-4\n")
    # 从缓存中加载对象
    neurons <- cached_data$neurons
    target_genes <- cached_data$target_genes
    output_dir <- cached_data$output_dir
  } else {
    cat("缓存文件不完整，重新运行步骤1-4\n")
    quit()
  }
}, error = function(e) {
  cat(sprintf("加载缓存失败: %s\n", e$message))
  quit()
})

# 步骤5：93基因提取（Windows兼容版）
cat(">>>> 步骤5: 93基因提取（Windows兼容版）\n") 

# 首先诊断：显示Seurat对象中的实际基因名样本
cat("Seurat对象中的基因名样本（前20个）:\n")
print(head(rownames(neurons), 20))

# 铜死亡核心基因检查
cat("铜死亡核心基因检查:\n")
core_candidates <- c("fdx1", "lias", "slc31a1", "dlat", "Fdx1", "Lias", "Slc31a1", "Dlat", 
                     "FDX1", "LIAS", "SLC31A1", "DLAT")
for(g in core_candidates) {
  if(g %in% rownames(neurons)) {
    cat(sprintf("  ✓ 找到: %s\n", g))
  }
}

# 加载基因映射表，获取小鼠基因符号
cat("加载基因映射表...\n")
gene_mapping_file <- file.path(output_dir, "gene_mapping_93.csv")

# 检查文件存在性
if (!file.exists(gene_mapping_file)) {
  cat("错误: 基因映射文件不存在，请先完成步骤1\n")
  flush.console()
  quit()
}

gene_mapping <- read.csv(gene_mapping_file, stringsAsFactors = FALSE)

# 最简匹配方案（忽略大小写）
target_lower <- tolower(gene_mapping$mouse_symbol[!is.na(gene_mapping$mouse_symbol)])
available_lower <- tolower(rownames(neurons))
matched_lower <- intersect(target_lower, available_lower)
# 转回原始格式
matched_genes <- rownames(neurons)[available_lower %in% matched_lower]

# 强制确保铜死亡核心基因被包含（即使表达低）
core_priority <- c("fdx1", "lias", "slc31a1", "dlat")  # 小写格式
core_found <- intersect(tolower(core_priority), matched_lower)
if(length(core_found) < 4) {
  cat(sprintf("警告: 仅找到%d/4个铜死亡核心基因，尝试查找其他格式...\n", length(core_found)))
  # 尝试查找任何包含这些基因名片段的基因
  for(core in c("fdx", "lias", "slc31a", "dlat")) {
    partial_match <- grep(core, rownames(neurons), value=TRUE, ignore.case=TRUE)
    if(length(partial_match) > 0) {
      cat(sprintf("  发现部分匹配 '%s': %s\n", core, paste(partial_match, collapse=", ")))
      matched_genes <- unique(c(matched_genes, partial_match))
    }
  }
}

cat(sprintf("目标基因匹配: %d/%d\n", length(matched_genes), length(target_lower)))

# 检查未匹配的基因
unmatched_genes <- setdiff(target_lower, matched_lower)
if (length(unmatched_genes) > 0) {
  cat("未匹配的基因: ", paste(unmatched_genes, collapse=", "), "\n")
}

# 2. 使用Seurat v5 LayerData接口直接提取（保持稀疏格式） 
cat("提取表达矩阵...\n")
expr_sparse <- LayerData(neurons, assay="RNA", layer="counts", features=matched_genes) 
# 保持稀疏矩阵格式，避免内存爆炸
exp_matrix <- t(expr_sparse)  # 转置为细胞×基因，保持稀疏格式

# 3. 智能过滤：保留高表达基因 + 强制保留核心基因（即使低表达）
gene_counts <- colSums(exp_matrix > 0) 
high_genes <- names(gene_counts[gene_counts >= 5]) 

# 动态检测核心基因实际格式
core_candidates <- c("fdx1", "lias", "slc31a1", "dlat", "mt2a", "atox1", "nfkb1") 
core_found <- c() 

cat("动态检测核心基因格式...\n")
for(core in core_candidates) { 
  # 尝试精确匹配（忽略大小写） 
  match_idx <- grep(paste0("^", core, "$"), colnames(exp_matrix), ignore.case=TRUE) 
  if(length(match_idx) > 0) { 
    actual_name <- colnames(exp_matrix)[match_idx[1]] 
    core_found <- c(core_found, actual_name) 
    cat(sprintf("  ✓ 核心基因匹配: %s → %s\n", core, actual_name)) 
  } else { 
    # 尝试模糊匹配（包含子串） 
    partial_idx <- grep(core, colnames(exp_matrix), ignore.case=TRUE) 
    if(length(partial_idx) > 0) { 
      actual_name <- colnames(exp_matrix)[partial_idx[1]] 
      core_found <- c(core_found, actual_name) 
      cat(sprintf("  ✓ 核心基因模糊匹配: %s → %s\n", core, actual_name)) 
    } 
  } 
}

# 强制合并（使用实际检测到的格式） 
keep_genes <- unique(c(high_genes, core_found)) 
exp_matrix <- exp_matrix[, keep_genes, drop=FALSE] 

cat(sprintf("高表达基因: %d个 | 强制保留核心基因: %d个 | 最终: %d基因\n",  
             length(high_genes), length(core_found), ncol(exp_matrix)))

# 检查铜死亡核心基因是否存在
core_genes <- c("fdx1", "lias", "slc31a1", "dlat")
# 转换为实际格式进行检查
missing_core <- c()
for(core in core_genes) {
  match_idx <- grep(paste0("^", core, "$"), colnames(exp_matrix), ignore.case=TRUE)
  if(length(match_idx) == 0) {
    missing_core <- c(missing_core, core)
  }
}
if(length(missing_core) > 0) {
  cat("警告: 铜死亡核心基因缺失: ", paste(missing_core, collapse=", "), "\n")
  flush.console()
}

# 保存key_genes变量，避免后续代码出错
key_genes <- c()

cat(">>>> 93基因提取与预处理完成 | 关键统计: 最终基因数 =", ncol(exp_matrix), "\n\n")

# 步骤6：分块PC网络（一举两得版）
cat(">>>> 步骤6: 分块PC网络分析（全部93基因保留）\n")

# 安装并加载必要的包
cat("安装并加载pcalg及其依赖包...\n")
flush.console()
tryCatch({
  # 设置CRAN镜像
  options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
  
  if (!require("pcalg", character.only = TRUE)) {
    # 安装BiocManager（如果未安装）
    if (!require("BiocManager", character.only = TRUE)) {
      install.packages("BiocManager")
      library(BiocManager)
    }
    # 安装RBGL（pcalg的依赖）
    BiocManager::install("RBGL")
    # 安装pcalg
    install.packages("pcalg")
    library(pcalg)
  }
  cat("pcalg包加载成功\n")
  flush.console()
}, error = function(e) {
  cat(sprintf("包安装失败: %s\n", e$message))
  flush.console()
  stop("无法安装必要的包，请手动安装pcalg及其依赖")
})

# 创建与exp_matrix行数匹配的逻辑向量
exp_cells <- rownames(exp_matrix)
cell_indices <- match(exp_cells, colnames(neurons))
cell_groups <- neurons$group[cell_indices]

# 自适应调整模块基因名格式以匹配exp_matrix
adjust_module_format <- function(module_genes, available_genes) {
  adjusted <- c()
  for(g in module_genes) {
    idx <- grep(paste0("^", g, "$"), available_genes, ignore.case=TRUE)
    if(length(idx) > 0) adjusted <- c(adjusted, available_genes[idx[1]])
  }
  return(unique(adjusted))
}

# 1. 定义功能模块（93个基因的科学分类）
gene_modules <- list(
  cuproptosis_core = adjust_module_format(
    c("fdx1", "lias", "dlat", "dld", "lipt1", "pdhx", "pdhb", "slc31a1"),
    colnames(exp_matrix)),
  copper_homeostasis = adjust_module_format(
    c("atp7b", "atp7a", "atox1", "commd1", "mt2a", "cp"),
    colnames(exp_matrix)),
  inflammation = adjust_module_format(
    c("nfkb1", "rela", "tlr4", "il6", "stat1", "stat3", "ccl2", "ptgs2"),
    colnames(exp_matrix)),
  others = setdiff(colnames(exp_matrix),
                   unlist(lapply(list(
                     c("fdx1", "lias", "dlat", "dld", "lipt1", "pdhx", "pdhb", "slc31a1"),
                     c("atp7b", "atp7a", "atox1", "commd1", "mt2a", "cp"),
                     c("nfkb1", "rela", "tlr4", "il6", "stat1", "stat3", "ccl2", "ptgs2")
                   ), function(x) adjust_module_format(x, colnames(exp_matrix))))
)

cat("模块基因格式调整完成:\n")
for(mod_name in names(gene_modules)) {
  cat(sprintf("  %s: %d个基因 (%s...)\n", mod_name, length(gene_modules[[mod_name]]), 
              paste(head(gene_modules[[mod_name]], 3), collapse=", ")))
}
flush.console()

# 2. 分块构建网络（每块基因数<25，计算极快）
build_module_network <- function(module_genes, module_name, data) {
  tryCatch({
    # 确保基因存在，尝试不同的大小写格式
    present_genes <- c()
    for (gene in module_genes) {
      # 尝试原始格式
      if (gene %in% colnames(data)) {
        present_genes <- c(present_genes, gene)
      } else {
        # 尝试小写
        gene_lower <- tolower(gene)
        if (gene_lower %in% colnames(data)) {
          present_genes <- c(present_genes, gene_lower)
        } else {
          # 尝试大写
          gene_upper <- toupper(gene)
          if (gene_upper %in% colnames(data)) {
            present_genes <- c(present_genes, gene_upper)
          } else {
            # 尝试首字母大写
            gene_title <- tools::toTitleCase(gene)
            if (gene_title %in% colnames(data)) {
              present_genes <- c(present_genes, gene_title)
            }
          }
        }
      }
    }
    
    # 去重
    present_genes <- unique(present_genes)
    
    if (length(present_genes) < 3) {
      cat(sprintf("%s模块: 基因数不足（%d个），跳过\n", module_name, length(present_genes)))
      flush.console()
      return(NULL)
    }
    
    module_data <- data[, present_genes, drop=FALSE]
    
    # 根据样本量动态调整alpha
    n_cells <- nrow(module_data)
    alpha <- ifelse(n_cells < 100, 0.1, ifelse(n_cells > 500, 0.01, ALPHA_BASE))
    cat(sprintf("%s模块: %d细胞，使用alpha=%.3f\n", module_name, n_cells, alpha))
    flush.console()
    
    # 临时转为dense矩阵用于计算，若细胞数过多则采样
    if (n_cells > 10000) {
      cat("细胞数过多，采样10000个细胞计算相关系数...\n")
      flush.console()
      sample_idx <- sample(n_cells, 10000)
      module_data_sample <- module_data[sample_idx, , drop=FALSE]
      module_data_dense <- as.matrix(module_data_sample)
    } else {
      module_data_dense <- as.matrix(module_data)
    }
    
    # PC算法（小数据量极快）
    library(pcalg)
    suff <- list(C=cor(module_data_dense), n=nrow(module_data_dense))
    pc.fit <- pc(suff, indepTest=gaussCItest, alpha=alpha, 
                 labels=colnames(module_data_dense), verbose=FALSE)
    
    # 提取边
    adj <- as(pc.fit@graph, "matrix")
    edges <- which(adj[upper.tri(adj)] == 1, arr.ind=TRUE)
    
    cat(sprintf("%s模块: %d基因, %d条边\n", module_name, length(present_genes), length(edges)))
    flush.console()
    return(list(adj=adj, genes=present_genes))
  }, error = function(e) {
    cat(sprintf("%s模块计算失败: %s\n", module_name, e$message))
    flush.console()
    return(NULL)
  })
}

# 3. 为Control和Stroke分别构建
all_edges_list <- list()
network_statistics <- data.frame()

# 加载必要的包
if (!require("data.table", character.only = TRUE)) {
  install.packages("data.table")
  library(data.table)
}

if (!require("parallel", character.only = TRUE)) {
  install.packages("parallel")
  library(parallel)
}

# Stroke组诊断
cat("\n==== Stroke组诊断 ====\n")
stroke_idx <- cell_groups == "stroke"
cat(sprintf("Stroke组细胞数: %d\n", sum(stroke_idx)))
if(sum(stroke_idx) > 0) {
  stroke_data <- exp_matrix[stroke_idx, , drop=FALSE]
  cat(sprintf("Stroke组表达矩阵维度: %d细胞 × %d基因\n", nrow(stroke_data), ncol(stroke_data)))
  cat("核心基因在Stroke组的表达情况:\n")
  for(g in core_found) {
    if(g %in% colnames(stroke_data)) {
      expr_cells <- sum(stroke_data[, g] > 0)
      cat(sprintf("  %s: %d/%d细胞表达 (%.1f%%)\n", 
                  g, expr_cells, nrow(stroke_data), 100*expr_cells/nrow(stroke_data)))
    }
  }
  if(sum(stroke_idx) < 30) {
    cat("警告: Stroke组细胞数<30，PC算法无法运行（样本量不足）\n")
  }
}
cat("====================\n\n")

for (group in c("control", "stroke")) {
  tryCatch({
    idx <- cell_groups == group
    group_data <- exp_matrix[idx, , drop=FALSE]
    
    cat(sprintf("\n构建%s组网络:\n", group))
    flush.console()  # 强制刷新控制台输出
    
    # 并行计算4个模块
    cat("并行计算模块网络...\n")
    flush.console()
    
    # 检测操作系统并选择并行方案
    is_windows <- Sys.info()["sysname"] == "Windows"
    
    if (!is_windows && N_CORES > 1) {
      # Mac/Linux使用mclapply
      library(parallel)
      cat(sprintf("使用mclapply进行%d核并行计算\n", N_CORES))
      networks <- mclapply(names(gene_modules), function(mod_name) {
        cat(sprintf("处理模块: %s\n", mod_name))
        build_module_network(gene_modules[[mod_name]], mod_name, group_data)
      }, mc.cores=N_CORES)
    } else if (is_windows && N_CORES > 1) {
      # Windows使用parLapply（支持多核）
      library(parallel)
      cat(sprintf("Windows系统：使用parLapply进行%d核并行计算\n", N_CORES))
      # 创建集群
      cl <- makeCluster(N_CORES)
      # 导出必要变量到集群
      clusterExport(cl, c("build_module_network", "gene_modules", "group_data", 
                          "ALPHA_BASE", "COR_THRESHOLD"), 
                    envir=environment())
      # 加载必要包到集群节点
      clusterEvalQ(cl, {
        library(pcalg)
        library(Matrix)
      })
      # 并行计算
      networks <- parLapply(cl, names(gene_modules), function(mod_name) {
        cat(sprintf("处理模块: %s\n", mod_name))
        build_module_network(gene_modules[[mod_name]], mod_name, group_data)
      })
      # 关闭集群（重要！）
      stopCluster(cl)
    } else {
      # 单核模式（Windows/Mac/Linux通用）
      cat("使用单核lapply计算\n")
      networks <- lapply(names(gene_modules), function(mod_name) {
        cat(sprintf("处理模块: %s\n", mod_name))
        build_module_network(gene_modules[[mod_name]], mod_name, group_data)
      })
    }
    
    names(networks) <- names(gene_modules)
    
    # 4. 合并所有模块的边（跨模块边用简单相关性补充）
    cat("合并模块边...\n")
    flush.console()  # 强制刷新控制台输出
    
    all_edges <- data.frame(from=character(), to=character(), weight=numeric())
    
    # 收集模块内边
    for (mod in networks) {
      if (is.null(mod)) next
      adj <- mod$adj
      genes <- mod$genes
      for (i in 1:nrow(adj)) {
        for (j in i:ncol(adj)) {
          if (adj[i,j] == 1) {
            all_edges <- rbind(all_edges, 
                             data.frame(from=genes[i], to=genes[j], weight=1))
          }
        }
      }
    }
    
    # 去重（使用data.table提高效率）
    cat("去重边列表...\n")
    flush.console()
    
    if (nrow(all_edges) > 0) {
      all_edges <- as.data.table(all_edges)
      all_edges[, edge_id := paste(pmin(from, to), pmax(from, to), sep="_")]
      all_edges <- unique(all_edges, by="edge_id")
      all_edges <- as.data.frame(all_edges[, c("from", "to", "weight")])
    }
    
    # 保存边列表（默认输出格式）
    cat("保存网络边列表...\n")
    flush.console()  # 强制刷新控制台输出
    write.csv(all_edges, file.path(output_dir, sprintf("network_%s_edges.csv", group)), row.names=FALSE)
    
    # 保存边列表（用于差异分析）
    all_edges_list[[group]] <- all_edges
    
    cat(sprintf("%s组网络构建完成\n", group))
    flush.console()
    
  }, error = function(e) {
    cat(sprintf("%s组网络构建失败: %s\n", group, e$message))
    flush.console()
    # 保存空的边列表
    all_edges <- data.frame(from=character(), to=character(), weight=numeric())
    write.csv(all_edges, file.path(output_dir, sprintf("network_%s_edges.csv", group)), row.names=FALSE)
    all_edges_list[[group]] <- all_edges
  })
}

cat(">>>> PC因果网络分析完成\n\n")

# 步骤7：差异网络分析
cat(">>>> 步骤7: 差异网络分析\n")

# 检查两组网络是否都构建成功
if (!is.null(all_edges_list$control) && !is.null(all_edges_list$stroke)) {
  # 提取两组的边
  control_edges <- all_edges_list$control
  stroke_edges <- all_edges_list$stroke
  
  # 创建边的唯一标识符
  control_edges$edge_id <- apply(control_edges[, c("from", "to")], 1, function(x) paste(sort(x), collapse="_"))
  stroke_edges$edge_id <- apply(stroke_edges[, c("from", "to")], 1, function(x) paste(sort(x), collapse="_"))
  
  # 识别Stroke特异性新边（在stroke中存在，在control中不存在）
  stroke_specific_edges <- stroke_edges[!stroke_edges$edge_id %in% control_edges$edge_id, ]
  
  if (nrow(stroke_specific_edges) > 0) {
    # 移除edge_id列，重命名列名
    stroke_specific_edges <- stroke_specific_edges[, c("from", "to", "weight")]
    colnames(stroke_specific_edges) <- c("From", "To", "weight")
    stroke_specific_edges$Type <- "Stroke_Specific"
    
    # 输出novel_stroke_edges.csv
    edges_output <- file.path(output_dir, "novel_stroke_edges.csv")
    write.csv(stroke_specific_edges, edges_output, row.names=FALSE)
    cat(sprintf("Stroke特异性新边已保存到: %s\n", edges_output))
    flush.console()
    
    # 统计：stroke特异性边数
    cat(sprintf("Stroke特异性边数: %d\n", nrow(stroke_specific_edges)))
    flush.console()
  } else {
    cat("未发现Stroke特异性新边\n")
    flush.console()
    # 创建空文件
    edges_output <- file.path(output_dir, "novel_stroke_edges.csv")
    write.csv(data.frame(From = character(), To = character(), Type = character(), weight = numeric()), 
              edges_output, row.names=FALSE)
  }
} else {
  cat("至少一组网络构建失败，无法进行差异分析\n")
  flush.console()
  # 创建空文件
  edges_output <- file.path(output_dir, "novel_stroke_edges.csv")
  write.csv(data.frame(From = character(), To = character(), Type = character(), weight = numeric()), 
            edges_output, row.names=FALSE)
}

cat(">>>> 差异网络分析完成\n\n")

# 验证所有输出文件是否生成
cat(">>>> 验证输出文件\n")
output_files <- c(
  file.path(output_dir, "gene_mapping_93.csv"),
  file.path(output_dir, "network_control_edges.csv"),
  file.path(output_dir, "network_stroke_edges.csv"),
  file.path(output_dir, "novel_stroke_edges.csv")
)

for (file in output_files) {
  if (file.exists(file)) {
    cat(sprintf("✓ %s 已生成\n", basename(file)))
  } else {
    cat(sprintf("✗ %s 未生成\n", basename(file)))
  }
  flush.console()
}

cat("\n分析完成！\n")
flush.console()