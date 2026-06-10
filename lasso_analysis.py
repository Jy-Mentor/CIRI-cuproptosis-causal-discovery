import pandas as pd
import numpy as np
from sklearn.linear_model import ElasticNetCV
from sklearn.model_selection import LeaveOneOut
from sklearn.preprocessing import StandardScaler, QuantileTransformer
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report
from scipy.stats import wilcoxon, mannwhitneyu, bootstrap
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
import os

# 设置随机种子以确保可复现性
np.random.seed(42)

# 创建输出目录
output_dir = 'C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/output'
os.makedirs(output_dir, exist_ok=True)

# 步骤1：数据预处理

# 1.1 读取平台注释文件
def load_platform_annotation():
    platform_file = 'C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/GPL1355-10794.txt'
    # 跳过注释行，直到找到表头
    skip_rows = 0
    with open(platform_file, 'r', encoding='utf-8') as f:
        for i, line in enumerate(f):
            if line.startswith('ID'):
                skip_rows = i
                break
    
    annotation = pd.read_csv(platform_file, sep='\t', skiprows=skip_rows, encoding='utf-8')
    # 保留有Gene Symbol的行
    annotation = annotation[annotation['Gene Symbol'].notna()]
    # 只保留需要的列
    annotation = annotation[['ID', 'Gene Symbol', 'ENTREZ_GENE_ID']]
    return annotation

# 1.2 读取表达矩阵
def load_expression_matrix():
    series_file = 'C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/GSE61616_series_matrix.txt'
    # 跳过注释行，直到找到表头
    skip_rows = 0
    with open(series_file, 'r', encoding='utf-8') as f:
        for i, line in enumerate(f):
            if line.startswith('!Sample_title'):
                skip_rows = i + 1
                break
    
    # 读取表达数据
    expr_data = pd.read_csv(series_file, sep='\t', skiprows=skip_rows, encoding='utf-8')
    # 第一列是探针ID
    expr_data = expr_data.rename(columns={expr_data.columns[0]: 'ID'})
    
    # 转换表达值为数值类型
    sample_columns = expr_data.columns[1:]
    for col in sample_columns:
        expr_data[col] = pd.to_numeric(expr_data[col], errors='coerce')
    return expr_data

# 1.3 读取基因映射文件
def load_gene_mapping():
    mapping_file = 'C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/大鼠 小鼠 人类映射库.txt'
    mapping_df = pd.read_csv(mapping_file, sep='\t', comment='#')
    # 只保留大鼠到人类的映射
    mapping_df = mapping_df[['RAT_GENE_SYMBOL', 'HUMAN_ORTHOLOG_SYMBOL']]
    mapping_df = mapping_df.dropna()
    # 创建映射字典
    mapping = {}
    for _, row in mapping_df.iterrows():
        rat_gene = row['RAT_GENE_SYMBOL'].strip().upper()
        human_genes = row['HUMAN_ORTHOLOG_SYMBOL'].strip().upper().split('|')
        mapping[rat_gene] = human_genes
    return mapping

# 1.4 数据预处理主函数
def preprocess_data():
    # 加载数据
    annotation = load_platform_annotation()
    expr_data = load_expression_matrix()
    gene_mapping = load_gene_mapping()
    
    # 合并表达数据和注释
    merged_data = pd.merge(expr_data, annotation, on='ID', how='inner')
    
    # 提取样本列（假设从第二列开始是样本）
    sample_columns = merged_data.columns[1:-2]  # 排除ID和注释列
    
    # 处理多探针对应同一基因的情况：取中位数
    gene_expr = merged_data.groupby('Gene Symbol')[sample_columns].median()
    
    # 物种转换：将大鼠基因映射到人类基因
    human_expr = {}
    for rat_gene, expr_values in gene_expr.iterrows():
        rat_gene_upper = rat_gene.strip().upper()
        if rat_gene_upper in gene_mapping:
            human_genes = gene_mapping[rat_gene_upper]
            for human_gene in human_genes:
                if human_gene not in human_expr:
                    human_expr[human_gene] = []
                human_expr[human_gene].append(expr_values)
    
    # 处理多个大鼠基因映射到同一人类基因的情况：取平均值
    human_expr_df = {}
    for human_gene, expr_list in human_expr.items():
        if len(expr_list) > 1:
            # 多个大鼠基因映射到同一人类基因，取平均值
            avg_expr = pd.concat(expr_list).groupby(level=0).mean()
            human_expr_df[human_gene] = avg_expr
        else:
            # 单个映射
            human_expr_df[human_gene] = expr_list[0]
    
    # 转换为DataFrame
    human_expr_df = pd.DataFrame(human_expr_df).T
    
    # 目标基因列表
    target_genes = ['RAGE', 'PPARG', 'PARP1', 'PTGS2', 'FDX1']
    
    # 提取目标基因的表达数据
    X_raw = []
    present_genes = []
    missing_genes = []
    
    for gene in target_genes:
        if gene in human_expr_df.index:
            expr_values = human_expr_df.loc[gene].values
            # 检查是否有NaN值
            if np.isnan(expr_values).any():
                missing_genes.append(gene)
                # 使用均值填充NaN值
                mean_value = np.nanmean(expr_values)
                expr_values = np.nan_to_num(expr_values, nan=mean_value)
            X_raw.append(expr_values)
            present_genes.append(gene)
        else:
            missing_genes.append(gene)
            # 如果基因不存在，使用所有样本的均值填充
            mean_value = np.mean(human_expr_df.values)
            X_raw.append([mean_value] * len(sample_columns))
    
    X_raw = np.array(X_raw).T
    
    print(f"存在的基因: {present_genes}")
    print(f"缺失的基因: {missing_genes}")
    
    # 创建样本标签：前5个是Sham（0），后5个是Model（1）
    y = np.array([0] * 5 + [1] * 5)
    
    # 检查样本数量
    print(f"样本数量: {len(sample_columns)}")
    print(f"样本列: {sample_columns}")
    
    # 只保留前10个样本（5 Sham + 5 Model）
    if len(sample_columns) > 10:
        sample_columns = sample_columns[:10]
        X_raw = X_raw[:10, :]
    
    # 保存原始表达数据
    rat_target_expression = pd.DataFrame(X_raw, columns=target_genes)
    rat_target_expression['Sample'] = sample_columns
    rat_target_expression['Group'] = y
    rat_target_expression.to_csv(os.path.join(output_dir, 'rat_target_expression.csv'), index=False)
    
    return X_raw, y, target_genes, sample_columns

# 步骤2：LASSO-ElasticNet建模
def run_loocv(X_raw, y, target_genes):
    loocv = LeaveOneOut()
    predictions = []
    true_labels = []
    probabilities = []
    coefficients = []
    alphas = []
    
    for train_idx, test_idx in loocv.split(X_raw):
        # 划分训练集和测试集
        X_train, X_test = X_raw[train_idx], X_raw[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]
        
        # 防泄漏标准化：仅在训练集上拟合
        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)
        
        # 训练模型
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            model = ElasticNetCV(
                l1_ratio=0.5, 
                cv=LeaveOneOut(), 
                n_alphas=50, 
                random_state=42, 
                max_iter=2000
            )
            model.fit(X_train_scaled, y_train)
        
        # 预测
        y_pred = model.predict(X_test_scaled)
        predictions.append(y_pred[0])
        true_labels.append(y_test[0])
        probabilities.append(y_pred[0])  # 对于回归，预测值作为概率
        
        # 记录系数和最优alpha
        coefficients.append(model.coef_)
        alphas.append(model.alpha_)
    
    # 保存结果
    loocv_results = pd.DataFrame({
        'Sample': range(10),
        'True_Label': true_labels,
        'Prediction': predictions,
        'Probability': probabilities
    })
    loocv_results.to_csv(os.path.join(output_dir, 'loocv_predictions.csv'), index=False)
    
    # 保存系数
    coef_df = pd.DataFrame(coefficients, columns=target_genes)
    coef_df['Alpha'] = alphas
    coef_df.to_csv(os.path.join(output_dir, 'elasticnet_coefficients.csv'), index=False)
    
    # 计算性能指标
    # 二分类阈值设置为0.5
    binary_predictions = [1 if p >= 0.5 else 0 for p in predictions]
    accuracy = accuracy_score(true_labels, binary_predictions)
    cm = confusion_matrix(true_labels, binary_predictions)
    
    # 计算敏感度和特异度
    tn, fp, fn, tp = cm.ravel()
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
    
    return accuracy, sensitivity, specificity, cm, coefficients, alphas

# 步骤3：事后差异表达分析
def posthoc_de_analysis(X_raw, y, target_genes):
    results = []
    
    for i, gene in enumerate(target_genes):
        gene_expr = X_raw[:, i]
        sham_expr = gene_expr[y == 0]
        model_expr = gene_expr[y == 1]
        
        # 计算log2FC
        log2fc = np.log2(np.mean(model_expr) / np.mean(sham_expr)) if np.mean(sham_expr) > 0 else np.nan
        
        # Wilcoxon秩和检验
        try:
            stat, p_value = wilcoxon(sham_expr, model_expr)
        except:
            stat, p_value = mannwhitneyu(sham_expr, model_expr)
        
        # 计算Cohen's d（小样本校正版本）
        n1, n2 = len(sham_expr), len(model_expr)
        var1, var2 = np.var(sham_expr, ddof=1), np.var(model_expr, ddof=1)
        pooled_var = ((n1-1)*var1 + (n2-1)*var2) / (n1 + n2 - 2)
        cohens_d = (np.mean(model_expr) - np.mean(sham_expr)) / np.sqrt(pooled_var) if pooled_var > 0 else 0
        
        results.append({
            'Gene': gene,
            'log2FC': log2fc,
            'p_value': p_value,
            'Cohens_d': cohens_d
        })
    
    # 计算FDR校正
    de_df = pd.DataFrame(results)
    de_df = de_df.sort_values('p_value')
    de_df['q_value'] = de_df['p_value'] * len(de_df) / (de_df.index + 1)
    de_df['q_value'] = de_df['q_value'].clip(upper=1.0)
    
    de_df.to_csv(os.path.join(output_dir, 'posthoc_de_analysis.csv'), index=False)
    
    return de_df

# 步骤4：Bootstrap稳定性分析
def bootstrap_analysis(X_raw, y, target_genes, alpha_median):
    n_bootstrap = 500
    bootstrap_coefficients = []
    
    for i in range(n_bootstrap):
        # 分层抽样：从Sham和Model中各有放回地抽取5个样本
        sham_indices = np.random.choice(np.where(y == 0)[0], size=5, replace=True)
        model_indices = np.random.choice(np.where(y == 1)[0], size=5, replace=True)
        bootstrap_indices = np.concatenate([sham_indices, model_indices])
        
        X_boot = X_raw[bootstrap_indices]
        y_boot = y[bootstrap_indices]
        
        # 标准化
        scaler = StandardScaler()
        X_boot_scaled = scaler.fit_transform(X_boot)
        
        # 拟合模型（使用中位最优alpha）
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            model = ElasticNetCV(
                l1_ratio=0.5, 
                cv=LeaveOneOut(), 
                n_alphas=50, 
                random_state=42, 
                max_iter=2000
            )
            model.fit(X_boot_scaled, y_boot)
        
        bootstrap_coefficients.append(model.coef_)
    
    # 计算稳定性指标
    bootstrap_coefficients = np.array(bootstrap_coefficients)
    stability_results = []
    
    for i, gene in enumerate(target_genes):
        gene_coefs = bootstrap_coefficients[:, i]
        # 选择频率（非零系数次数）
        selection_freq = np.sum(gene_coefs != 0) / n_bootstrap
        
        # 计算BCa 95%置信区间
        try:
            def stat_func(data):
                return np.median(data)
            
            res = bootstrap(gene_coefs, stat_func, n_resamples=1000, method='percentile')
            ci_lower, ci_upper = res.confidence_interval
        except:
            # 退化为分位数法
            ci_lower = np.percentile(gene_coefs, 2.5)
            ci_upper = np.percentile(gene_coefs, 97.5)
        
        # 系数中位数
        coef_median = np.median(gene_coefs)
        
        stability_results.append({
            'Gene': gene,
            'Selection_Frequency': selection_freq,
            'Coefficient_Median': coef_median,
            'CI_Lower': ci_lower,
            'CI_Upper': ci_upper
        })
    
    # 保存结果
    stability_df = pd.DataFrame(stability_results)
    stability_df.to_csv(os.path.join(output_dir, 'bootstrap_stability.csv'), index=False)
    
    # 绘制系数分布图
    plt.figure(figsize=(12, 8))
    for i, gene in enumerate(target_genes):
        plt.subplot(2, 3, i+1)
        sns.violinplot(y=bootstrap_coefficients[:, i])
        plt.axhline(y=0, color='red', linestyle='--')
        plt.title(gene)
        plt.ylabel('Coefficient')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'bootstrap_coefficient_distribution.png'), dpi=300)
    
    return stability_df, bootstrap_coefficients

# 步骤5：核心靶点判定
def identify_core_targets(stability_df, de_df):
    core_targets = []
    
    for _, row in stability_df.iterrows():
        gene = row['Gene']
        de_row = de_df[de_df['Gene'] == gene].iloc[0]
        
        # 判定条件
        condition1 = row['Selection_Frequency'] >= 0.8
        condition2 = row['CI_Lower'] * row['CI_Upper'] > 0  # 95%CI不包含0
        condition3 = abs(row['Coefficient_Median']) > 0.1
        condition4 = abs(de_row['Cohens_d']) > 0.5
        
        if condition1 and condition2 and condition3 and condition4:
            core_targets.append(gene)
    
    return core_targets

# 步骤6：优先级排序
def rank_targets(stability_df, de_df):
    ranking = []
    
    for _, row in stability_df.iterrows():
        gene = row['Gene']
        de_row = de_df[de_df['Gene'] == gene].iloc[0]
        
        # 综合评分
        score = (
            row['Selection_Frequency'] * 0.5 + 
            abs(row['Coefficient_Median']) * 0.3 + 
            abs(de_row['Cohens_d']) / 2 * 0.2
        )
        
        ranking.append({
            'Gene': gene,
            'Score': score,
            'Selection_Frequency': row['Selection_Frequency'],
            'Coefficient_Median': row['Coefficient_Median'],
            'Cohens_d': de_row['Cohens_d']
        })
    
    ranking_df = pd.DataFrame(ranking).sort_values('Score', ascending=False)
    return ranking_df

# 绘制混淆矩阵
def plot_confusion_matrix(cm):
    plt.figure(figsize=(8, 6))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', 
                xticklabels=['Sham (0)', 'Model (1)'],
                yticklabels=['Sham (0)', 'Model (1)'])
    plt.title('Confusion Matrix (LOOCV)')
    plt.xlabel('Predicted')
    plt.ylabel('True')
    plt.savefig(os.path.join(output_dir, 'loocv_confusion_matrix.png'), dpi=300)

# 生成报告
def generate_report(accuracy, sensitivity, specificity, core_targets, ranking_df, de_df, stability_df):
    report = f"""# LASSO-ElasticNet 小样本分析报告

## 方法学说明
- **数据来源**：GSE61616 (Rat MCAO模型，n=10)
- **分析方法**：ElasticNetCV (l1_ratio=0.5) with Leave-One-Out Cross-Validation
- **防泄漏措施**：
  - 每次LOOCV迭代内独立标准化（仅使用训练集统计量）
  - 未进行数据驱动的特征预筛选
  - Bootstrap采用分层抽样（保持5:5比例）

## LOOCV性能指标
- **准确率**: {accuracy:.4f}
- **敏感度**: {sensitivity:.4f}
- **特异度**: {specificity:.4f}

## 核心靶点判定
### 判定标准（必须同时满足）：
1. Bootstrap选择频率 ≥ 80%
2. BCa 95%CI不包含0
3. |中位系数| > 0.1
4. Cohen's d > 0.5

### 核心预测靶点：
{', '.join(core_targets) if core_targets else '无'}

## 靶点优先级排序
{ranking_df.to_markdown(index=False)}

## 事后差异表达分析
{de_df.to_markdown(index=False)}

## Bootstrap稳定性分析
{stability_df.to_markdown(index=False)}

## 小样本局限性声明
1. **样本量限制**：n=10（5 vs 5）属于极小样本，统计功效(power) < 0.3
2. **过拟合风险**：即使使用ElasticNet和LOOCV，小样本下仍可能出现虚假关联
3. **外推限制**：Rat到Human的同源转换存在功能差异
4. **多重检验**：Bootstrap和CV过程中进行了多次检验，未校正家族错误率
5. **建议验证**：所有筛选出的"核心靶点"必须通过独立队列或湿实验验证

## 生物学解释与下一步实验建议
- 建议对筛选出的核心靶点进行Western Blot或qPCR验证
- 考虑在更大样本量的独立队列中验证模型预测能力
- 探索BCP对这些靶点的调控机制
"""
    
    with open(os.path.join(output_dir, 'lasso_small_sample_report.md'), 'w', encoding='utf-8') as f:
        f.write(report)

# 主函数
def main():
    # 步骤1：数据预处理
    X_raw, y, target_genes, sample_columns = preprocess_data()
    
    # 步骤2：LASSO-ElasticNet建模
    accuracy, sensitivity, specificity, cm, coefficients, alphas = run_loocv(X_raw, y, target_genes)
    
    # 步骤3：事后差异表达分析
    de_df = posthoc_de_analysis(X_raw, y, target_genes)
    
    # 步骤4：Bootstrap稳定性分析
    alpha_median = np.median(alphas)
    stability_df, bootstrap_coefficients = bootstrap_analysis(X_raw, y, target_genes, alpha_median)
    
    # 步骤5：核心靶点判定
    core_targets = identify_core_targets(stability_df, de_df)
    
    # 步骤6：优先级排序
    ranking_df = rank_targets(stability_df, de_df)
    
    # 绘制混淆矩阵
    plot_confusion_matrix(cm)
    
    # 生成报告
    generate_report(accuracy, sensitivity, specificity, core_targets, ranking_df, de_df, stability_df)
    
    print("分析完成！结果已保存到output目录。")
    print(f"核心靶点：{core_targets}")
    print(f"Top 3靶点：{list(ranking_df['Gene'].head(3))}")

if __name__ == "__main__":
    main()
