# -*- coding: utf-8 -*-
"""
多数据集合并 + 批次效应校正
合并GSE61616(大鼠)和GSE97537(大鼠)
使用ComBat校正批次效应
输出: merged_expr_matrix.csv
"""

import os
import gzip
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler

# 路径配置
DOWNLOAD_DIR = r"C:\Users\Jy-Mentor-7\Downloads"
RESULTS_DIR = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\results"
MERGED_DIR = os.path.join(RESULTS_DIR, "merged_datasets")
os.makedirs(MERGED_DIR, exist_ok=True)


def parse_series_matrix(filepath):
    """解析GEO series matrix文件"""
    print(f"解析: {os.path.basename(filepath)}")
    
    with gzip.open(filepath, 'rt', encoding='utf-8') as f:
        lines = f.readlines()
    
    # 提取样本信息
    sample_titles = []
    sample_ids = []
    characteristics = []
    
    for line in lines:
        if line.startswith('!Sample_title'):
            sample_titles = [x.strip().strip('"') for x in line.strip().split('\t')[1:]]
        elif line.startswith('!Sample_geo_accession'):
            sample_ids = [x.strip().strip('"') for x in line.strip().split('\t')[1:]]
        elif line.startswith('!Sample_characteristics_ch1'):
            characteristics.append([x.strip().strip('"') for x in line.strip().split('\t')[1:]])
        elif line.startswith('!series_matrix_table_begin'):
            break
    
    # 找到数据矩阵开始位置
    data_start_idx = None
    for i, line in enumerate(lines):
        if line.startswith('!series_matrix_table_begin'):
            data_start_idx = i + 1
            break
    
    if data_start_idx is None:
        raise ValueError("未找到数据矩阵")
    
    # 读取数据矩阵
    data_lines = []
    for i in range(data_start_idx, len(lines)):
        if lines[i].startswith('!series_matrix_table_end'):
            break
        data_lines.append(lines[i])
    
    # 解析为DataFrame
    from io import StringIO
    data_str = ''.join(data_lines)
    df = pd.read_csv(StringIO(data_str), sep='\t', index_col=0)
    
    # 设置列名
    df.columns = sample_ids
    
    print(f"  样本数: {len(sample_ids)}")
    print(f"  探针数: {len(df)}")
    print(f"  样本: {sample_titles[:3]}...")
    
    return df, sample_titles, sample_ids


def load_gse61616():
    """加载GSE61616已处理的表达矩阵"""
    expr_file = os.path.join(RESULTS_DIR, "stage1_rma_degs", "expr_matrix.csv")
    
    if not os.path.exists(expr_file):
        raise FileNotFoundError(f"GSE61616表达矩阵不存在: {expr_file}")
    
    df = pd.read_csv(expr_file, index_col=0)
    print(f"GSE61616: {df.shape[0]} 探针 x {df.shape[1]} 样本")
    
    return df


def merge_datasets(expr1, expr2):
    """合并两个数据集(取交集探针)"""
    print("\n合并数据集...")
    
    # 取交集探针
    common_probes = set(expr1.index) & set(expr2.index)
    print(f"  GSE61616探针数: {len(expr1)}")
    print(f"  GSE97537探针数: {len(expr2)}")
    print(f"  交集探针数: {len(common_probes)}")
    
    if len(common_probes) < 100:
        raise ValueError(f"交集探针过少: {len(common_probes)}")
    
    # 提取交集
    merged = pd.concat([
        expr1.loc[sorted(common_probes)],
        expr2.loc[sorted(common_probes)]
    ], axis=1)
    
    print(f"  合并后: {merged.shape[0]} 探针 x {merged.shape[1]} 样本")
    
    return merged


def combat_batch_correction(merged_expr, batch_labels):
    """ComBat批次效应校正(简化版)"""
    print("\nComBat批次效应校正...")
    
    try:
        import sklearn
        from sklearn.linear_model import LinearRegression
        
        # 转置: 样本x基因
        expr_t = merged_expr.T.copy()
        
        # 标准化每个基因
        scaler = StandardScaler()
        expr_scaled = scaler.fit_transform(expr_t)
        expr_scaled = pd.DataFrame(expr_scaled, 
                                   index=expr_t.index, 
                                   columns=expr_t.columns)
        
        # 按批次校正
        unique_batches = np.unique(batch_labels)
        corrected = expr_scaled.copy()
        
        for gene in expr_scaled.columns:
            gene_expr = expr_scaled[gene].values
            
            # 计算批次均值
            batch_means = []
            for batch in unique_batches:
                mask = np.array(batch_labels) == batch
                batch_means.append(gene_expr[mask].mean())
            
            # 全局均值
            global_mean = gene_expr.mean()
            
            # 校正
            for i, batch in enumerate(unique_batches):
                mask = np.array(batch_labels) == batch
                corrected.loc[mask, gene] = gene_expr[mask] - batch_means[i] + global_mean
        
        print(f"  批次: {unique_batches}")
        print(f"  校正完成")
        
        return corrected.T  # 转回: 基因x样本
        
    except Exception as e:
        print(f"  ComBat校正失败: {e}")
        print(f"  使用标准化替代...")
        
        # 备选: 简单标准化
        scaler = StandardScaler()
        corrected = pd.DataFrame(
            scaler.fit_transform(merged_expr.T).T,
            index=merged_expr.index,
            columns=merged_expr.columns
        )
        
        return corrected


def create_sample_annotation_61616(sample_names):
    """创建GSE61616样本注释"""
    annotations = []
    
    for sid in sample_names:
        sid_lower = sid.lower()
        
        if 'sham' in sid_lower:
            group = 'Sham'
        elif 'model' in sid_lower:
            group = 'Model'
        elif 'xst' in sid_lower:
            group = 'Treatment'
        else:
            group = 'Unknown'
        
        annotations.append({
            'SampleID': sid,
            'Title': sid,
            'Group': group,
            'Dataset': 'GSE61616'
        })
    
    return pd.DataFrame(annotations)


def create_sample_annotation_97537(sample_titles, sample_ids):
    """创建GSE97537样本注释"""
    annotations = []
    
    for title, sid in zip(sample_titles, sample_ids):
        title_lower = title.lower()
        
        if 'sham' in title_lower or 'control' in title_lower:
            group = 'Sham'
        elif 'mcao' in title_lower or 'ischemia' in title_lower or 'operated' in title_lower:
            group = 'Model'
        else:
            group = 'Unknown'
        
        annotations.append({
            'SampleID': sid,
            'Title': title,
            'Group': group,
            'Dataset': 'GSE97537'
        })
    
    return pd.DataFrame(annotations)


def main():
    """主流程"""
    print("=" * 60)
    print("多数据集合并 + 批次效应校正")
    print("仅使用Sham和Model样本（排除药物治疗组）")
    print("=" * 60)
    
    # 1. 加载GSE61616
    print("\n[1/5] 加载GSE61616...")
    gse61616_expr = load_gse61616()
    
    # 2. 加载GSE97537
    print("\n[2/5] 加载GSE97537...")
    gse97537_file = os.path.join(DOWNLOAD_DIR, "GSE97537_series_matrix (2).txt.gz")
    
    if not os.path.exists(gse97537_file):
        raise FileNotFoundError(f"GSE97537文件不存在: {gse97537_file}")
    
    gse97537_expr, titles97537, ids97537 = parse_series_matrix(gse97537_file)
    
    # 3. 过滤: 仅保留Sham和Model样本（排除XST药物治疗组）
    print("\n[3/5] 过滤样本（仅保留Sham和Model）...")
    
    # GSE61616过滤
    annot_61616 = create_sample_annotation_61616(gse61616_expr.columns.tolist())
    sham_model_61616 = annot_61616[annot_61616['Group'].isin(['Sham', 'Model'])]
    gse61616_filtered = gse61616_expr.loc[:, sham_model_61616['SampleID']]
    print(f"  GSE61616: {gse61616_expr.shape[1]} -> {gse61616_filtered.shape[1]} 样本（排除XST）")
    
    # GSE97537过滤
    annot_97537 = create_sample_annotation_97537(titles97537, ids97537)
    sham_model_97537 = annot_97537[annot_97537['Group'].isin(['Sham', 'Model'])]
    gse97537_filtered = gse97537_expr.loc[:, sham_model_97537['SampleID']]
    print(f"  GSE97537: {gse97537_expr.shape[1]} -> {gse97537_filtered.shape[1]} 样本")
    
    # 4. 合并数据集
    print("\n[4/5] 合并数据集...")
    merged_expr = merge_datasets(gse61616_filtered, gse97537_filtered)
    
    # 创建批次标签
    batch_labels = ['GSE61616'] * gse61616_filtered.shape[1] + ['GSE97537'] * gse97537_filtered.shape[1]
    
    # 5. 批次效应校正
    print("\n[5/5] 批次效应校正...")
    corrected_expr = combat_batch_correction(merged_expr, batch_labels)
    
    # 6. 保存结果
    print("\n保存结果...")
    
    # 表达矩阵
    output_file = os.path.join(MERGED_DIR, "merged_expr_matrix.csv")
    corrected_expr.to_csv(output_file)
    print(f"  表达矩阵: {output_file}")
    
    # 样本注释（仅Sham和Model）
    all_annotations = pd.concat([sham_model_61616, sham_model_97537], ignore_index=True)
    annot_file = os.path.join(MERGED_DIR, "sample_annotations.csv")
    all_annotations.to_csv(annot_file, index=False)
    print(f"  样本注释: {annot_file}")
    
    # 统计信息
    print("\n" + "=" * 60)
    print("合并统计（仅Sham和Model）:")
    print(f"  GSE61616: {gse61616_filtered.shape[1]} 样本")
    print(f"  GSE97537: {gse97537_filtered.shape[1]} 样本")
    print(f"  合并后: {corrected_expr.shape[1]} 样本")
    print(f"  探针数: {corrected_expr.shape[0]}")
    print(f"  分组分布:")
    for group in all_annotations['Group'].unique():
        count = len(all_annotations[all_annotations['Group'] == group])
        dataset_breakdown = all_annotations[all_annotations['Group'] == group]['Dataset'].value_counts()
        dataset_str = ', '.join([f"{d}: {c}" for d, c in dataset_breakdown.items()])
        print(f"    {group}: {count} ({dataset_str})")
    print("=" * 60)


if __name__ == "__main__":
    main()
