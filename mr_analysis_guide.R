# ================================================================================
# 孟德尔随机化 (MR) 分析脚本使用指南
# ================================================================================

# 1. 数据格式要求
# --------------------------------------------------------------------------------

# 暴露因素数据框必需列:
# - SNP: SNP 标识符 (如 rs123456)
# - beta.exposure: 暴露因素的效应值
# - se.exposure: 暴露效应的标准误
# - effect_allele.exposure: 效应等位基因
# - other_allele.exposure: 其他等位基因
# - eaf.exposure: 效应等位基因频率
# - pval.exposure: P 值

# 结局因素数据框必需列:
# - SNP: SNP 标识符 (与暴露因素一致)
# - beta.outcome: 结局的效应值
# - se.outcome: 结局效应的标准误
# - effect_allele.outcome: 效应等位基因
# - other_allele.outcome: 其他等位基因
# - eaf.outcome: 效应等位基因频率
# - pval.outcome: P 值


# 2. 使用方法
# --------------------------------------------------------------------------------

# 方法 1: 从文件加载数据
source("mr_analysis.R")

results <- run_mr_pipeline(
  exposure_file = "exposure_data.xlsx",  # 暴露因素文件
  outcome_file = "outcome_data.xlsx",    # 结局因素文件
  pval_threshold = 5e-8,                 # 工具变量 P 值阈值
  clump_kb = 10000,                      # LD clump 距离 (kb)
  clump_r2 = 0.001,                      # LD clump 的 r2 阈值
  output_dir = "./mr_results",           # 结果输出目录
  save_plots = TRUE                      # 是否保存图表
)


# 方法 2: 直接使用数据框
source("mr_analysis.R")

# 准备数据
exposure_dat <- readxl::read_excel("exposure.xlsx")
outcome_dat <- readxl::read_excel("outcome.xlsx")

# 运行分析
results <- run_mr_pipeline(
  exposure_data = exposure_dat,
  outcome_data = outcome_dat,
  output_dir = "./results"
)


# 方法 3: 从 IEU GWAS 数据库获取数据
source("mr_analysis.R")

# 获取 GWAS 数据
gwas_data <- fetch_gwas_data(
  exposure_id = "ieu-a-2",    # 暴露因素 GWAS ID
  outcome_id = "ieu-b-4"      # 结局因素 GWAS ID
)

# 运行分析
results <- run_mr_pipeline(
  exposure_data = gwas_data$exposure,
  outcome_data = gwas_data$outcome,
  output_dir = "./ieu_gwas_results"
)


# 方法 4: 命令行运行
# Rscript mr_analysis.R exposure_data.xlsx outcome_data.xlsx ./results


# 3. 结果解读
# --------------------------------------------------------------------------------

# MR 分析结果包含以下方法:
# - Inverse variance weighted (IVW): 主要分析方法
# - MR Egger: 检测并校正水平多效性
# - Weighted median: 中位数估计
# - Simple mode: 简单模式
# - Weighted mode: 加权模式

# 关键指标:
# - b: 因果效应估计值
# - se: 标准误
# - pval: P 值
# - ci_lower/ci_upper: 95% 置信区间

# 敏感性分析:
# - 异质性检验 (Cochran's Q): 评估工具变量间异质性
# - 水平多效性检验 (MR Egger intercept): 检测定向多效性
# - 留一法 (Leave-one-out): 评估单个 SNP 影响


# 4. 输出文件
# --------------------------------------------------------------------------------

# 结果目录包含以下文件:
# - mr_results.csv: 主要 MR 分析结果
# - heterogeneity_results.csv: 异质性检验结果
# - pleiotropy_results.csv: 多效性检验结果
# - leave_one_out_results.csv: 留一法分析结果
# - harmonised_data.csv: 协调后数据
# - mr_summary.txt: 分析总结报告
# - mr_scatter_plot.png: 散点图
# - mr_forest_plot.png: 森林图
# - mr_funnel_plot.png: 漏斗图
# - mr_leave_one_out_plot.png: 留一法图


# 5. 常见问题
# --------------------------------------------------------------------------------

# Q1: 工具变量太少怎么办？
# A: 放宽 P 值阈值 (如 1e-5 或 1e-4)

# Q2: 数据协调后 SNP 数量为 0？
# A: 检查 SNP 命名是否一致，等位基因是否匹配

# Q3: 异质性显著怎么办？
# A: 使用 MR Egger 或 Weighted median 等稳健方法

# Q4: 多效性检验显著？
# A: 结果需谨慎解释，可能存在混杂因素


# 6. 数据示例
# --------------------------------------------------------------------------------

# 暴露因素数据示例:
# SNP        beta.exposure  se.exposure  effect_allele  other_allele  eaf  pval
# rs123456   0.05           0.01         A              G             0.3  1e-10
# rs789012   -0.03          0.008        T              C             0.45 5e-9

# 结局因素数据示例:
# SNP        beta.outcome  se.outcome  effect_allele  other_allele  eaf  pval
# rs123456   0.02          0.005       A              G             0.3  3e-5
# rs789012   -0.015        0.004       T              C             0.45 2e-4


# 7. 高级用法
# --------------------------------------------------------------------------------

# 自定义分析流程:
source("mr_analysis.R")

# 1. 加载数据
exposure_dat <- load_exposure_data("exposure.xlsx")
outcome_dat <- load_outcome_data("outcome.xlsx")

# 2. 选择工具变量
instruments <- select_instruments(
  exposure_dat,
  pval_threshold = 1e-5,
  clump_kb = 5000,
  clump_r2 = 0.01
)

# 3. 提取结局数据中的对应 SNP
outcome_snps <- outcome_dat %>%
  dplyr::filter(SNP %in% instruments$SNP)

# 4. 数据协调
harmonised <- harmonise_data(instruments, outcome_snps)

# 5. 运行 MR 分析
mr_results <- run_mr_analysis(harmonised)

# 6. 敏感性分析
sensitivity <- run_sensitivity_analysis(harmonised)

# 7. 可视化
create_mr_plots(harmonised, mr_results, sensitivity, "./custom_results")

# 8. 导出结果
export_results(mr_results, sensitivity, harmonised, "./custom_results")
