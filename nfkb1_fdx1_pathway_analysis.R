# ============================================
# NFKB1到FDX1最短路径分析
# 在g_k3网络（133节点主网络）上执行
# ============================================

cat("正在加载必要的R包...\n")

if (!require("igraph", quietly = TRUE)) {
  install.packages("igraph", repos = "https://cloud.r-project.org/")
  library(igraph)
}

# 设置路径
work_dir <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙"
result_dir <- file.path(work_dir, "String_Network_Systematic_Analysis")

# 读取g_k3网络对象
cat("\n读取g_k3网络对象...\n")
g_k3_file <- file.path(result_dir, "02_g_k3_network_object.rds")

if (!file.exists(g_k3_file)) {
  stop("g_k3网络对象文件不存在，请先运行系统性分析脚本")
}

g_k3 <- readRDS(g_k3_file)
cat(paste0("g_k3网络 - 节点数: ", vcount(g_k3), ", 边数: ", ecount(g_k3), "\n"))

# 检查NFKB1和FDX1是否在网络中
network_nodes <- V(g_k3)$name
cat(paste0("\n检查目标节点...\n"))
cat(paste0("NFKB1在网络中: ", "NFKB1" %in% network_nodes, "\n"))
cat(paste0("FDX1在网络中: ", "FDX1" %in% network_nodes, "\n"))

# 如果节点不在网络中，显示警告
if (!("NFKB1" %in% network_nodes)) {
  cat("警告: NFKB1不在g_k3网络中\n")
  cat("网络中的NFKB1类似节点:\n")
  print(grep("NFKB", network_nodes, value = TRUE))
}

if (!("FDX1" %in% network_nodes)) {
  cat("警告: FDX1不在g_k3网络中\n")
  cat("网络中的FDX1类似节点:\n")
  print(grep("FDX", network_nodes, value = TRUE))
}

# ==================== NFKB1到FDX1最短路径分析 ====================
cat("\n========================================\n")
cat("NFKB1到FDX1最短路径分析\n")
cat("========================================\n")

# 只有在两个节点都在网络中时才执行分析
if (("NFKB1" %in% network_nodes) && ("FDX1" %in% network_nodes)) {
  
  # 计算所有最短路径
  paths <- all_shortest_paths(g_k3, from="NFKB1", to="FDX1", mode="all")
  
  cat(paste0("\nNFKB1到FDX1的最短路径数: ", length(paths$res), "\n"))
  
  if (length(paths$res) > 0) {
    cat(paste0("最短路径长度: ", length(paths$res[[1]]) - 1, "\n\n"))
    
    # 打印每条路径
    for(i in 1:length(paths$res)) {
      path_nodes <- names(paths$res[[i]])
      cat(paste0("路径 ", i, ": ", paste(path_nodes, collapse=" → "), "\n"))
    }
    
    # 提取所有桥接节点
    all_path_nodes <- unique(unlist(lapply(paths$res, function(p) names(p))))
    bridge_candidates <- setdiff(all_path_nodes, c("NFKB1", "FDX1"))
    cat(paste0("\n可能的桥接节点 (", length(bridge_candidates), " 个):\n"))
    cat(paste(bridge_candidates, collapse=", "))
    cat("\n")
    
    # ==================== 保存结果 ====================
    cat("\n\n保存分析结果...\n")
    
    # 创建路径数据框
    path_data <- data.frame()
    for(i in 1:length(paths$res)) {
      path_nodes <- names(paths$res[[i]])
      path_data <- rbind(path_data, data.frame(
        Path_ID = i,
        Path_Length = length(path_nodes) - 1,
        Path = paste(path_nodes, collapse=" → "),
        Nodes = paste(path_nodes, collapse=","),
        stringsAsFactors = FALSE
      ))
    }
    
    write.table(path_data, 
                file = file.path(result_dir, "05_nfkb1_fdx1_shortest_paths.txt"), 
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    # 保存桥接节点
    bridge_df <- data.frame(
      Bridge_Node = bridge_candidates,
      stringsAsFactors = FALSE
    )
    write.table(bridge_df, 
                file = file.path(result_dir, "05_nfkb1_fdx1_bridge_nodes.txt"), 
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    # 创建详细的节点路径表
    detailed_paths <- data.frame()
    for(i in 1:length(paths$res)) {
      path_nodes <- names(paths$res[[i]])
      for(j in 1:length(path_nodes)) {
        detailed_paths <- rbind(detailed_paths, data.frame(
          Path_ID = i,
          Step = j - 1,
          Node = path_nodes[j],
          Is_Source = (path_nodes[j] == "NFKB1"),
          Is_Target = (path_nodes[j] == "FDX1"),
          Is_Bridge = !(path_nodes[j] %in% c("NFKB1", "FDX1")),
          stringsAsFactors = FALSE
        ))
      }
    }
    
    write.table(detailed_paths, 
                file = file.path(result_dir, "05_nfkb1_fdx1_detailed_paths.txt"), 
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    # 计算每个桥接节点在多少条路径中出现
    bridge_frequency <- table(unlist(lapply(paths$res, function(p) {
      nodes <- names(p)
      setdiff(nodes, c("NFKB1", "FDX1"))
    })))
    
    bridge_freq_df <- data.frame(
      Bridge_Node = names(bridge_frequency),
      Path_Count = as.numeric(bridge_frequency),
      Path_Frequency = as.numeric(bridge_frequency) / length(paths$res) * 100,
      stringsAsFactors = FALSE
    )
    bridge_freq_df <- bridge_freq_df[order(-bridge_freq_df$Path_Count), ]
    
    write.table(bridge_freq_df, 
                file = file.path(result_dir, "05_nfkb1_fdx1_bridge_frequency.txt"), 
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    cat("\n桥接节点频率统计:\n")
    print(bridge_freq_df)
    
    # ==================== 路径可视化准备 ====================
    cat("\n========================================\n")
    cat("路径可视化数据准备\n")
    cat("========================================\n")
    
    # 提取路径子图
    path_subgraph_nodes <- unique(unlist(lapply(paths$res, names)))
    path_subgraph <- induced_subgraph(g_k3, path_subgraph_nodes)
    
    cat(paste0("路径子图 - 节点数: ", vcount(path_subgraph), ", 边数: ", ecount(path_subgraph), "\n"))
    
    # 保存路径子图的边列表（用于Cytoscape可视化）
    path_edges <- as.data.frame(get.edgelist(path_subgraph))
    path_edges$interaction <- "pp"
    path_edges_file <- file.path(result_dir, "05_nfkb1_fdx1_pathway_edges.sif")
    write.table(path_edges[, c(1, 3, 2)], file = path_edges_file, 
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    cat(paste0("路径边列表已保存: ", path_edges_file, "\n"))
    
    # 创建节点属性文件（标记节点类型）
    path_node_attr <- data.frame(
      Node = V(path_subgraph)$name,
      Node_Type = ifelse(V(path_subgraph)$name == "NFKB1", "Source",
                        ifelse(V(path_subgraph)$name == "FDX1", "Target", "Bridge")),
      In_Path_Count = sapply(V(path_subgraph)$name, function(n) {
        sum(sapply(paths$res, function(p) n %in% names(p)))
      }),
      stringsAsFactors = FALSE
    )
    
    # 添加中心性信息
    path_node_attr$DC <- centrality_df$DC[match(path_node_attr$Node, centrality_df$Node)]
    path_node_attr$RRA_Rank <- rra_df$Rank[match(path_node_attr$Node, rra_df$Name)]
    
    path_node_file <- file.path(result_dir, "05_nfkb1_fdx1_pathway_nodes.txt")
    write.table(path_node_attr, file = path_node_file, 
                sep = "\t", quote = FALSE, row.names = FALSE)
    cat(paste0("路径节点属性已保存: ", path_node_file, "\n"))
    
    cat("\n========================================\n")
    cat("         最短路径分析完成\n")
    cat("========================================\n")
    
  } else {
    cat("\n错误: 未找到NFKB1到FDX1的路径\n")
    
    # 检查连通性
    if (are.connected(g_k3, "NFKB1", "FDX1")) {
      cat("节点是连通的，但可能存在问题\n")
    } else {
      cat("NFKB1和FDX1在g_k3网络中不连通\n")
      
      # 找到它们所属的连通分量
      components <- components(g_k3)
      nfkb1_comp <- components$membership["NFKB1"]
      fdx1_comp <- components$membership["FDX1"]
      cat(paste0("NFKB1属于连通分量: ", nfkb1_comp, "\n"))
      cat(paste0("FDX1属于连通分量: ", fdx1_comp, "\n"))
    }
  }
  
} else {
  cat("\n错误: NFKB1或FDX1不在g_k3网络中，无法执行路径分析\n")
  cat("\n网络中的节点示例:\n")
  print(head(network_nodes, 20))
}
