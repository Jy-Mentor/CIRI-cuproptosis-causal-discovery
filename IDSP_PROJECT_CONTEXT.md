# IDSP 项目上下文提示词

> 将此文件内容粘贴到新的 AI 会话中，即可让 AI 获得完整项目上下文。
> 需要查询原始代码时，告诉 AI 工作区路径为当前文件夹。

---

## 1. 项目概述

**项目名**: CIRI-cuproptosis-causal-discovery  
**科学目标**: 在脑缺血再灌注损伤（CIRI）中识别和验证「铁驱动的衰老程序（Iron-Driven Senescence Program, IDSP）」。  
**核心假说**: 铁死亡（急性、早期峰）和细胞衰老（慢性、持续）在 CIRI 中既有区分又有协同，IDSP Index 可量化二者的共激活程度。

**GitHub 仓库**: `https://github.com/Jy-Mentor/CIRI-cuproptosis-causal-discovery`

---

## 2. 工作区路径

```
c:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\
```

需要查看原始代码细节时，告诉 AI 从此路径读取文件。

---

## 3. 文件结构

| 文件 | 用途 |
|------|------|
| `l1_dual_analysis.py` | **核心文件**：L1 双评分分析流水线（~1800行） |
| `idsp_gene_sets.py` | 基因集定义：PURE_FERROPTOSIS (73), PURE_SENESCENCE (139), SHARED_GENES (12) |
| `ferro_aging_ciri_analysis.py` | 原始铁死亡-衰老分析脚本（数据加载函数来源） |
| `GPL6883.annot.gz` | 人类 Illumina 平台注释 |
| `l1_results/` | 输出目录（CSV + 图表 + 报告） |

---

## 4. 数据依赖

5 个 GEO 数据集，位于 `D:\反向网络药理学\L1 数据集\bulk\`：

| 数据集 | 物种 | 实验设计 | 样本量 |
|--------|------|---------|--------|
| GSE16561 | 人 | 全血 Stroke vs Control | 39 vs 24 |
| GSE37587 | 人 | 全血 Follow-Up vs Baseline（配对） | 34 vs 34 |
| GSE61616 | 大鼠 | MCAO 7d vs Sham | 5 vs 5 |
| GSE97537 | 大鼠 | MCAO 24h vs Sham | 7 vs 5 |
| GSE104036 | 小鼠 | RNA-seq 多时间点 (3/6/12/24h + Sham) | 3×4 vs 3 |

注释文件：
- `GPL6883.annot.gz`（人类，位于工作区根目录）
- `GPL1355-10794 (1).txt`（大鼠，位于 GSE61616 数据目录）

---

## 5. 核心方法论

### 5.1 L1 双评分分析流水线

```
数据加载 → 探针→基因折叠 → 秩和富集评分(铁死亡+衰老) → IDSP Index → 
Meta分析 → 高级验证(Bootstrap/置换/ROC/I²/LODO/RRA/JSD) → 报告+图表
```

### 5.2 IDSP Index 公式

```
IDSP Index = z(ferr) + z(sene) − |z(ferr) − z(sene)|
```
含义：铁死亡和衰老评分都高且差异小时，IDSP Index 最大。

### 5.3 基因集

- **PURE_FERROPTOSIS**: 73 基因（FerrDb Validated，排除衰老交集）
- **PURE_SENESCENCE**: 139 基因（CellAge + SenMayo，排除铁死亡交集）
- **SHARED_GENES**: 12 基因（铁死亡∩衰老）

### 5.4 验证标准

1. **双评分相关性**: mean_r < 0.6 → PASS（铁死亡和衰老不完全共线）
2. **GPX4 验证**: 高 IDSP 样本中 GPX4 不显著下降 → 排除典型铁死亡
3. **时间动态**: 铁死亡早期峰 vs 衰老持续激活（GSE104036 时间序列）

---

## 6. 完整函数清单（l1_dual_analysis.py）

### 数据加载

| 函数 | 签名 | 用途 |
|------|------|------|
| `find_file` | `(dir_path, keywords) -> Optional[str]` | 按关键词查找数据文件 |
| `parse_series_matrix` | `(filepath) -> pd.DataFrame` | 解析 GEO Series Matrix 文件 |
| `parse_gpl6883_annotation` | `(annot_path) -> Dict[str, str]` | 解析人类 GPL6883 探针注释 |
| `parse_gpl1355_annotation` | `(filepath) -> Dict[str, str]` | 解析大鼠 GPL1355 探针注释 |
| `collapse_probes` | `(expr_df, probe_map) -> pd.DataFrame` | 探针→基因折叠（最大表达值） |
| `_load_expr_gse16561` | `() -> Tuple[DataFrame, List, List]` | 加载 GSE16561 |
| `_load_expr_gse37587` | `() -> Tuple[DataFrame, List, List]` | 加载 GSE37587 |
| `_load_expr_gse61616` | `() -> Tuple[DataFrame, List, List]` | 加载 GSE61616（大鼠） |
| `_load_expr_gse97537` | `() -> Tuple[DataFrame, List, List]` | 加载 GSE97537（大鼠） |
| `_load_expr_gse104036` | `() -> Tuple[DataFrame, dict, List]` | 加载 GSE104036（时间序列） |

### 核心评分

| 函数 | 签名 | 用途 |
|------|------|------|
| `rank_sum_enrichment_score` | `(expr, gene_mask) -> float` | 单样本秩和富集评分 |
| `compute_enrichment_score_matrix` | `(expr_df, gene_set) -> pd.Series` | 对所有样本计算富集评分 |
| `calc_idsp_index` | `(ferr_score, sene_score) -> pd.Series` | 计算 IDSP Index（含零除保护） |
| `dual_enrichment_analysis` | `(expr_df, dataset_name, case_cols, control_cols) -> Tuple[DataFrame, dict]` | 双评分分析（核心） |

### 统计检验

| 函数 | 签名 | 用途 |
|------|------|------|
| `cohens_d` | `(case, control) -> float` | Cohen's d 效应量 |
| `fisher_meta_analysis` | `(p_values) -> Tuple[float, float]` | Fisher 合并 p 值 |
| `stouffer_meta` | `(p_values, weights=None, directions=None) -> float` | Stouffer Z-score Meta（含效应方向） |
| `random_effects_meta_analysis` | `(effect_sizes, variances) -> dict` | 随机效应 Meta（DerSimonian-Laird） |

### 高级分析

| 函数 | 签名 | 用途 |
|------|------|------|
| `bootstrap_idsp_ci` | `(scores_df, n_bootstrap=2000, ci=0.95, seed=42) -> dict` | Bootstrap IDSP Index 置信区间 |
| `permutation_enrichment_test` | `(scores, case_cols, control_cols, n_perm=2000, seed=42) -> dict` | 置换检验（直接接受评分 Series） |
| `dual_score_roc_auc` | `(scores_df, case_cols, control_cols) -> dict` | ROC/AUC 判别能力评估 |
| `i_squared_heterogeneity` | `(comparisons, effect_key, var_key) -> float` | I² 异质性（标准方差加权 Q） |
| `lodo_cross_validation` | `(comparisons, meta_func) -> pd.DataFrame` | 留一数据集交叉验证 |

### 前沿方法

| 函数 | 签名 | 用途 |
|------|------|------|
| `robust_rank_aggregation` | `(rank_matrix) -> pd.DataFrame` | RRA 跨数据集基因一致性（KS 检验） |
| `jsd_and_ks_comparison` | `(scores_df, case_cols, control_cols, dataset_name) -> dict` | JSD + KS 分布差异量化 |

### 验证与可视化

| 函数 | 签名 | 用途 |
|------|------|------|
| `gpx4_validation` | `(expr_df, scores_df, case_cols, control_cols, dataset_name) -> dict` | GPX4 分层验证（含 NaN 四分位保护） |
| `analyze_signature_genes` | `(expr_df, case_cols, control_cols, gene_set, dataset_name) -> pd.DataFrame` | 单基因差异分析 |
| `temporal_dual_analysis` | `(expr_df, timepoint_dict, sham_cols, dataset_name) -> pd.DataFrame` | 时间动态双评分 |
| `plot_forest_dual` | `(comparisons, save_path)` | 效应量森林图 |
| `plot_temporal_dual` | `(temporal_df, save_path)` | 时间动态双轴图 |
| `plot_scatter_dual` | `(all_scores_df, save_path)` | 双评分散点图 |
| `plot_gene_heatmap` | `(all_gene_dfs, save_path)` | 核心基因热图 |

---

## 7. 关键边界保护（已实现）

| 保护 | 位置 |
|------|------|
| `calc_idsp_index` 标准差为 0 返回全 0 Series | L399-402 |
| `gpx4_validation` NaN 四分位数保护 + n<4 检查 | L431-435 |
| `dual_enrichment_analysis` NaN 日志保护（`pd.notna(r_all)` 检查） | L385-389 |
| `cohens_d` 合并标准差为 0 返回 0 | L472 |
| `lodo_cross_validation` p 值与方向联合过滤避免长度不匹配 | L1032-1053 |
| `robust_rank_aggregation` rankdata 长度对齐（`pd.Series(..., index=series.index)`） | L1105-1106 |
| 报告写入使用 `safe_fmt` 防御 NaN | L1753-1754 |
| 时间动态 Sham 基线为水平虚线 | L1262-1269 |
| raw counts 检测使用浮点容差 | L778 |

---

## 8. 历史 Bug 修复清单

| 顺序 | Bug | 修复 |
|------|-----|------|
| 1 | `parse_gpl1355_annotation` 未定义 | 补充大鼠注释解析函数 |
| 2 | `calc_idsp_index` 可能除以零 | 标准差为 0 时返回全 0 Series |
| 3 | `gpx4_validation` 四分位数可能为 NaN | dropna + n<4 提前返回 |
| 4 | 置换检验传入全 1 伪矩阵 | 重写为直接接受评分 Series |
| 5 | `stouffer_meta` z 方向错误 | 增加 `directions` 参数 |
| 6 | I² 使用错误近似 | 改为标准方差加权 Q 公式 |
| 7 | `i_squared_heterogeneity_both` 未定义 | 统一使用 `i_squared_heterogeneity`，返回 float |
| 8 | `dual_enrichment_analysis` 未存储方差 | 补充 `var_ferroptosis`/`var_senescence` 字段 |
| 9 | LODO meta_func 签名不匹配 | 传双参数 `(pvals, dirs)` |
| 10 | LODO p 值与方向未对齐 | 联合过滤 `(p, d)` 元组 |
| 11 | 置换检验结果未保存 | 导出 `L1_permutation_tests.csv` |
| 12 | 高级 Meta 结果未写入报告/CSV | 导出 `L1_meta_analysis_summary.csv` + 写入验证报告 |
| 13 | `dual_enrichment_analysis` 日志 NaN 格式化崩溃 | `pd.notna(r_all)` 条件判断 |
| 14 | RRA rankdata 长度不匹配 | `pd.Series(rankdata, index=series.index)` |

---

## 9. 输出文件

所有输出位于 `l1_results/` 目录：

| 文件 | 内容 |
|------|------|
| `L1_dual_scores_all_datasets.csv` | 每个样本的铁死亡/衰老/IDSP 评分 |
| `L1_dual_comparison_summary.csv` | 各数据集 Cohen's d / p 值 / 方差 |
| `L1_temporal_dual_scores.csv` | GSE104036 时间动态评分 |
| `L1_gpx4_validation.csv` | GPX4 分层验证结果 |
| `L1_gene_level_analysis.csv` | 单基因差异分析 |
| `L1_idsp_index_all.csv` | IDSP Index |
| `L1_permutation_tests.csv` | 置换检验 p 值 |
| `L1_lodo_cross_validation.csv` | LODO 结果 |
| `L1_meta_analysis_summary.csv` | Fisher/Stouffer/随机效应 Meta 汇总 |
| `L1_rra_gene_consistency.csv` | RRA 基因一致性 |
| `L1_jsd_ks_distribution.csv` | JSD/KS 分布差异 |
| `L1_validation_report.txt` | 文本验证报告 |
| `figures/Fig1A_forest_dual.png` | 效应量森林图 |
| `figures/Fig1B_temporal_dual.png` | 时间动态双轴图 |
| `figures/Fig1C_scatter_dual.png` | 双评分散点图 |
| `figures/Fig1D_gene_heatmap.png` | 核心基因热图 |

---

## 10. 使用方法

```bash
cd "c:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙"
python l1_dual_analysis.py
```

依赖：`numpy`, `pandas`, `scipy`, `matplotlib`, `scikit-learn`, `statsmodels`

---

## 11. AI 使用说明

将此文件粘贴到新 AI 会话开头，并告知：
1. 工作区路径为此文件夹
2. 需要深入查看代码时，直接 read 对应文件
3. 当前 L1 阶段已完成，下一步为 L2 因果发现分析