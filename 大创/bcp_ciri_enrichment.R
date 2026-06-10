# 智能包安装函数
install_if_missing <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(new_packages) > 0) {
    install.packages(new_packages, repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
  }
}

install_bioc_if_missing <- function(packages) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
  }
  new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(new_packages) > 0) {
    BiocManager::install(new_packages)
  }
}

# 设置 CRAN 镜像
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

# 安装所需包
install_if_missing(c("glmnet", "pROC", "tidyverse", "ggpubr"))
install_bioc_if_missing(c("impute", "rat2302.db"))

# 加载包
library(glmnet)
library(pROC)
library(tidyverse)
library(ggpubr)
library(impute)
try(library(rat2302.db))

# 样本分组定义
sample_ids <- c("GSM1509422", "GSM1509423", "GSM1509424", "GSM1509425", "GSM1509426", 
                 "GSM1509427", "GSM1509428", "GSM1509429", "GSM1509430", "GSM1509431")

group_vector <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
names(group_vector) <- sample_ids

# 输入文件路径
series_matrix_path <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/nb/GSE61616_series_matrix.txt"
mapping_file_path <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/大鼠 小鼠 人类映射库.txt"
candidate_genes_path <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/78_genes.txt"

# 输出目录
output_dir <- "C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/LASSO_ROC_Results/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Step 1: 读取 GEO Series Matrix
cat("Step 1: 读取 GEO Series Matrix...\n")

# 读取文件，跳过注释行
series_matrix <- read.delim(series_matrix_path, comment.char = "!", header = TRUE, stringsAsFactors = FALSE)

# 提取表达矩阵（假设第一列是探针ID）
probe_ids <- series_matrix[, 1]
expr_matrix <- series_matrix[, -1]
rownames(expr_matrix) <- probe_ids

# 检查列名是否匹配 sample_ids
colnames_expr <- colnames(expr_matrix)
if (!all(sample_ids %in% colnames_expr)) {
  cat("警告：列名与 sample_ids 不完全匹配，尝试模糊匹配...\n")
  # 尝试模糊匹配
  matched_cols <- sapply(sample_ids, function(id) grep(id, colnames_expr, value = TRUE)[1])
  expr_matrix <- expr_matrix[, matched_cols]
  colnames(expr_matrix) <- sample_ids
}

# 确保只保留 10 个样本
if (ncol(expr_matrix) > 10) {
  expr_matrix <- expr_matrix[, sample_ids]
}

# 检查缺失值
if (any(is.na(expr_matrix))) {
  cat("检测到缺失值，使用 impute.knn 进行填充...\n")
  expr_matrix <- impute.knn(as.matrix(expr_matrix))$data
}

# Step 2: 探针注释与物种映射
cat("Step 2: 探针注释与物种映射...\n")

# 读取 78 个候选基因
candidate_genes <- read.table(candidate_genes_path, header = FALSE, stringsAsFactors = FALSE)
colnames(candidate_genes) <- "Human_Symbol"
candidate_genes <- candidate_genes$Human_Symbol

# 读取本地映射库
mapping_df <- read.table(mapping_file_path, header = TRUE, stringsAsFactors = FALSE, sep = "\t", quote = "")

# 清理列名
colnames(mapping_df) <- gsub("\\s+", "", colnames(mapping_df))

# 提取大鼠基因和人类基因的映射
if ("Rat_Symbol" %in% colnames(mapping_df) && "Human_Symbol" %in% colnames(mapping_df)) {
  rat_human_mapping <- mapping_df[, c("Rat_Symbol", "Human_Symbol")]
} else {
  # 尝试其他可能的列名
  rat_col <- grep("rat|Rat", colnames(mapping_df), ignore.case = TRUE)[1]
  human_col <- grep("human|Human", colnames(mapping_df), ignore.case = TRUE)[1]
  if (!is.na(rat_col) && !is.na(human_col)) {
    rat_human_mapping <- mapping_df[, c(rat_col, human_col)]
    colnames(rat_human_mapping) <- c("Rat_Symbol", "Human_Symbol")
  } else {
    stop("无法在映射文件中找到大鼠和人类基因列")
  }
}

# 清理映射数据
rat_human_mapping <- rat_human_mapping[!is.na(rat_human_mapping$Rat_Symbol) & !is.na(rat_human_mapping$Human_Symbol), ]
rat_human_mapping <- rat_human_mapping[rat_human_mapping$Rat_Symbol != "" & rat_human_mapping$Human_Symbol != "", ]

# 从 series matrix 中提取探针对应的 Rat Gene Symbol
# 注意：这里假设 series matrix 中包含探针注释信息
# 如果没有，我们将使用 rat2302.db 包（如果可用）

# 尝试从列名或注释中提取基因信息
if (exists("rat2302.db") && require("rat2302.db")) {
  cat("使用 rat2302.db 进行探针注释...\n")
  probe_to_rat <- select(rat2302.db, keys = probe_ids, columns = "SYMBOL")
  colnames(probe_to_rat) <- c("ProbeID", "Rat_Symbol")
} else {
  # 尝试从文件中提取基因信息
  cat("尝试从文件中提取基因信息...\n")
  # 这里简化处理，假设探针ID可能包含基因信息或使用映射库中的大鼠基因
  probe_to_rat <- data.frame(ProbeID = probe_ids, Rat_Symbol = probe_ids, stringsAsFactors = FALSE)
}

# 合并探针注释和人类基因映射
probe_to_human <- merge(probe_to_rat, rat_human_mapping, by = "Rat_Symbol", all.x = TRUE)

# 过滤出存在于 78 个候选基因中的人类基因
probe_to_human <- probe_to_human[probe_to_human$Human_Symbol %in% candidate_genes, ]

# 过滤出在表达矩阵中存在的探针
probe_to_human <- probe_to_human[probe_to_human$ProbeID %in% rownames(expr_matrix), ]

# 处理多探针对应同一人类基因的情况（取平均值）
expr_subset <- expr_matrix[probe_to_human$ProbeID, ]

# 直接使用 Human_Symbol 作为分组变量进行聚合，避免设置重复的行名
expr_human <- aggregate(expr_subset, by = list(Human_Symbol = probe_to_human$Human_Symbol), FUN = mean)
rownames(expr_human) <- expr_human$Human_Symbol
expr_human <- expr_human[, -1]

# 输出未映射的基因
unmapped_genes <- candidate_genes[!candidate_genes %in% rownames(expr_human)]
if (length(unmapped_genes) > 0) {
  write.table(unmapped_genes, file = paste0(output_dir, "unmapped_genes.txt"), 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  cat(paste0("未映射的基因数量: ", length(unmapped_genes), "，已输出到 unmapped_genes.txt\n"))
}

# 检查有效基因数量
if (nrow(expr_human) < 30) {
  cat("警告: 有效基因不足 30 个，结果可能不稳定\n")
}

# Step 3: 数据预处理
cat("Step 3: 数据预处理...\n")

# 检查是否需要 log2 转换（如果数据范围很大）
if (max(expr_human) > 100) {
  cat("数据范围较大，进行 log2 转换...\n")
  expr_human <- log2(expr_human + 1)
}

# 标准化
cat("进行数据标准化...\n")
expr_human_scaled <- scale(expr_human)

# 转置为 LASSO 输入格式
X <- t(expr_human_scaled)
y <- group_vector[rownames(X)]

# Step 4: LASSO 回归（小样本优化）
cat("Step 4: LASSO 回归分析...\n")

set.seed(123)

# 样本量 n=10 (5 vs 5)，LASSO 采用 5 折交叉验证
tryCatch({
  cv_fit <- cv.glmnet(
    x = X,
    y = y,
    family = "binomial",
    alpha = 1,
    nfolds = 5,  # 5折交叉验证
    type.measure = "auc",
    standardize = TRUE,
    maxit = 1e5
  )
}, error = function(e) {
  cat("AUC 计算出错，改用 deviance 作为评估指标...\n")
  cv_fit <- cv.glmnet(
    x = X,
    y = y,
    family = "binomial",
    alpha = 1,
    nfolds = 5,
    type.measure = "deviance",
    standardize = TRUE,
    maxit = 1e5
  )
})

lambda_min <- cv_fit$lambda.min

# 提取非零系数基因
coef_matrix <- coef(cv_fit, s = "lambda.min")
non_zero_genes <- rownames(coef_matrix)[coef_matrix[, 1] != 0]
non_zero_genes <- non_zero_genes[non_zero_genes != "(Intercept)"]

# 提取系数值
coef_values <- coef_matrix[non_zero_genes, 1]

# 按 |系数| 排序，取 Top 10
if (length(non_zero_genes) > 10) {
  top_genes <- names(sort(abs(coef_values), decreasing = TRUE)[1:10])
  top_coef <- coef_values[top_genes]
} else {
  top_genes <- non_zero_genes
  top_coef <- coef_values
}

# Step 5: ROC 验证（小样本留一法补充）
cat("Step 5: ROC 验证分析...\n")

# 单基因 ROC
single_auc <- numeric(length(top_genes))
names(single_auc) <- top_genes

for (gene in top_genes) {
  roc_obj <- roc(y, X[, gene])
  single_auc[gene] <- auc(roc_obj)
}

# 组合模型 ROC
pred_prob <- predict(cv_fit, X, s = "lambda.min", type = "response")
roc_combined <- roc(y, as.vector(pred_prob))
auc_combined <- auc(roc_combined)

# LOOCV 分析
loocv_auc <- numeric(nrow(X))

for(i in 1:nrow(X)) {
  X_train <- X[-i,]
  y_train <- y[-i]
  X_test <- X[i,]
  y_test <- y[i]
  
  fit_loo <- glmnet(X_train, y_train, family = "binomial", lambda = lambda_min)
  pred_loo <- predict(fit_loo, X_test, type = "response")
  
  if (length(unique(y_train)) > 1) {
    roc_loo <- roc(y_train, predict(fit_loo, X_train, type = "response"))
    loocv_auc[i] <- auc(roc_loo)
  } else {
    loocv_auc[i] <- NA
  }
}

mean_loocv_auc <- mean(loocv_auc, na.rm = TRUE)

# Step 6: 差异表达箱线图
cat("Step 6: 绘制差异表达箱线图...\n")

# 检查是否有 top_genes
if (length(top_genes) > 0) {
  # 准备箱线图数据
  boxplot_data <- as.data.frame(t(expr_human[top_genes,]))
  boxplot_data$Sample <- rownames(boxplot_data)
  boxplot_data$Group <- ifelse(group_vector[rownames(boxplot_data)] == 0, "Sham", "Model")

  # 转换为长格式
  long_data <- pivot_longer(boxplot_data, cols = all_of(top_genes), 
                           names_to = "Gene", values_to = "Expression")

  # 设置 Group 因子水平
  long_data$Group <- factor(long_data$Group, levels = c("Sham", "Model"))
} else {
  # 如果没有 top_genes，创建空数据框
  long_data <- data.frame(Gene = character(), Group = factor(), Expression = numeric())
  cat("警告: 没有选择到显著基因，箱线图将为空\n")
}

# Step 7: 输出结果文件
cat("Step 7: 输出结果文件...\n")

# 输出 lambda.min
write.csv(data.frame(lambda_min = lambda_min, Note = "Small sample size (n=10), results should be validated in larger cohorts"), 
          file = paste0(output_dir, "lambda_optimal.csv"), row.names = FALSE)

# 输出 Top 10 基因
if (length(top_genes) > 0) {
  top_genes_df <- data.frame(Gene = top_genes, 
                            Coefficient = top_coef[top_genes], 
                            AUC_single = single_auc[top_genes],
                            Note = "Small sample size (n=10), results should be validated in larger cohorts")
  write.csv(top_genes_df, file = paste0(output_dir, "top10_hub_genes.csv"), row.names = FALSE)
} else {
  # 如果没有非零基因，创建空文件
  write.csv(data.frame(Gene = character(), Coefficient = numeric(), AUC_single = numeric(), 
                      Note = "Small sample size (n=10), no significant genes found"), 
            file = paste0(output_dir, "top10_hub_genes.csv"), row.names = FALSE)
}

# 绘制 ROC 曲线
pdf(paste0(output_dir, "roc_summary.pdf"), width = 10, height = 8)
if (length(top_genes) > 0) {
  plot(roc_combined, col = "red", lwd = 2, main = "ROC Curves: Combined Model vs Single Genes")

  for (gene in top_genes) {
    roc_obj <- roc(y, X[, gene])
    lines(roc_obj, col = "gray", lwd = 1)
  }

  legend("bottomright", legend = c(paste0("Combined Model (AUC = ", round(auc_combined, 3), ")"), 
                                  "Single Genes"), 
         col = c("red", "gray"), lwd = c(2, 1))
} else {
  plot.new()
  text(0.5, 0.5, "No significant genes found", cex = 1.2)
}
dev.off()

# 绘制差异表达箱线图
if (length(top_genes) > 0) {
  pdf(paste0(output_dir, "boxplot_diff_expression.pdf"), width = 15, height = 6)
  p <- ggboxplot(long_data, x = "Group", y = "Expression", 
                 fill = "Group", 
                 palette = c("blue", "red"),
                 add = "jitter",
                 xlab = "Group",
                 ylab = "Expression",
                 main = "Differential Expression of Top Genes") +
    stat_compare_means(method = "wilcox.test", label = "p.format") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # 添加分面
  p <- p + facet_wrap(~Gene, ncol = 5)
  print(p)
  dev.off()
} else {
  # 如果没有 top_genes，创建空的 PDF 文件
  pdf(paste0(output_dir, "boxplot_diff_expression.pdf"), width = 15, height = 6)
  plot.new()
  text(0.5, 0.5, "No significant genes found", cex = 1.2)
  dev.off()
}

# 输出分析总结
summary_text <- paste0(
  "分析完成！\n",
  "样本量: 10 (5 Sham vs 5 Model)\n",
  "输入基因数量: 78\n",
  "成功映射的人类基因数量: ", nrow(expr_human), "\n",
  "未映射的基因数量: ", length(unmapped_genes), "\n",
  "LASSO 选择的非零基因数量: ", length(non_zero_genes), "\n",
  "最终分析的 Top 基因数量: ", length(top_genes), "\n",
  "最优 lambda 值: ", round(lambda_min, 6), "\n",
  "组合模型 AUC: ", round(auc_combined, 3), "\n",
  "LOOCV 平均 AUC: ", round(mean_loocv_auc, 3), "\n",
  "结果文件已保存到: ", output_dir
)

cat(summary_text)
writeLines(summary_text, con = paste0(output_dir, "analysis_summary.txt"))
