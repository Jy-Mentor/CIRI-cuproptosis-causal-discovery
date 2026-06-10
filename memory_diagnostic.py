#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
内存诊断脚本 - 定位内存暴涨点
"""

import os
import sys
import time
import psutil
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import coo_matrix

# 内存监控
class MemoryMonitor:
    def __init__(self):
        self.initial_mem = self.get_mem()
        self.peak_mem = self.initial_mem
        self.checkpoints = []
        
    def get_mem(self):
        process = psutil.Process(os.getpid())
        return process.memory_info().rss / 1024 / 1024
    
    def checkpoint(self, label):
        current = self.get_mem()
        delta = current - self.initial_mem
        self.peak_mem = max(self.peak_mem, current)
        self.checkpoints.append({
            'label': label,
            'memory_mb': current,
            'delta_mb': delta,
            'timestamp': time.time()
        })
        print(f"[内存检查点] {label}: {current:.1f} MB (Δ{delta:+.1f} MB)")
        return current
    
    def report(self):
        print("\n" + "="*70)
        print("内存诊断报告")
        print("="*70)
        print(f"初始内存: {self.initial_mem:.1f} MB")
        print(f"峰值内存: {self.peak_mem:.1f} MB")
        print(f"最大增量: {self.peak_mem - self.initial_mem:.1f} MB")
        print("\n详细检查点:")
        for cp in self.checkpoints:
            print(f"  {cp['label']:40s}: {cp['memory_mb']:8.1f} MB (Δ{cp['delta_mb']:+7.1f} MB)")

monitor = MemoryMonitor()

# 设置路径
WORK_DIR = r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙"
DATA_DIR = r"C:\Users\Jy-Mentor-7\Desktop\虚拟敲除"

print("="*70)
print("内存诊断 - GSE174574 SCISSOR分析")
print("="*70)

# 检查点1: 初始状态
monitor.checkpoint("脚本启动")

# 读取一个样本的数据
print("\n读取单个样本数据...")
sample_matrix = os.path.join(DATA_DIR, "GSM5319987_sham1_matrix.mtx")
sample_features = os.path.join(DATA_DIR, "GSM5319987_sham1_genes.tsv")
sample_barcodes = os.path.join(DATA_DIR, "GSM5319987_sham1_barcodes.tsv")

# 读取基因和barcodes
with open(sample_features, 'r') as f:
    genes = [line.strip().split('\t')[1] if '\t' in line else line.strip() for line in f]

with open(sample_barcodes, 'r') as f:
    barcodes = [line.strip() for line in f]

monitor.checkpoint("读取基因/barcodes")

# 读取MTX矩阵
with open(sample_matrix, 'r') as f:
    lines = f.readlines()

data_start = 0
for i, line in enumerate(lines):
    if line.startswith('%%'):
        data_start = i + 1
    elif line.strip() and not line.startswith('%'):
        dims = line.strip().split()
        n_rows, n_cols, n_entries = int(dims[0]), int(dims[1]), int(dims[2])
        data_start = i + 1
        break

row_indices = []
col_indices = []
data_values = []

for line in lines[data_start:]:
    if line.strip():
        parts = line.strip().split()
        row_idx = int(parts[0]) - 1
        col_idx = int(parts[1]) - 1
        value = int(parts[2])
        row_indices.append(row_idx)
        col_indices.append(col_idx)
        data_values.append(value)

matrix = coo_matrix((data_values, (row_indices, col_indices)), 
                    shape=(n_rows, n_cols)).T.tocsr()

monitor.checkpoint("读取MTX矩阵 (sparse)")

# 创建AnnData
adata = sc.AnnData(X=matrix)
adata.var_names = genes
monitor.checkpoint("创建AnnData (sparse)")

# 测试1: .toarray() 转换
print("\n测试1: sparse -> dense 转换...")
try:
    X_dense = adata.X.toarray()
    monitor.checkpoint("转换为dense数组")
    print(f"  Dense数组大小: {X_dense.nbytes / 1024 / 1024:.1f} MB")
    del X_dense
except Exception as e:
    print(f"  错误: {e}")

# 测试2: Module Score计算
print("\n测试2: Module Score计算...")
marker_genes = ['Cx3cr1', 'Tmem119', 'P2ry12']  # Microglia markers
available_markers = [g for g in marker_genes if g in adata.var_names]

if len(available_markers) > 0:
    monitor.checkpoint("Module Score前")
    try:
        sc.tl.score_genes(adata, available_markers, score_name='test_score')
        monitor.checkpoint("Module Score后")
    except Exception as e:
        print(f"  错误: {e}")
        monitor.checkpoint("Module Score失败")

# 测试3: 子集化减少内存
print("\n测试3: 子集化策略...")
if adata.n_vars > 3000:
    # 策略A: 选择高变基因
    sc.pp.filter_genes(adata, min_cells=3)
    monitor.checkpoint("过滤低表达基因")
    
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    monitor.checkpoint("标准化后")
    
    sc.pp.highly_variable_genes(adata, n_top_genes=3000)
    adata_hvg = adata[:, adata.var['highly_variable']].copy()
    monitor.checkpoint("子集化到3000 HVG")
    
    print(f"  HVG子集: {adata_hvg.n_obs} × {adata_hvg.n_vars}")
    
    # 测试在HVG上计算Module Score
    try:
        sc.tl.score_genes(adata_hvg, available_markers, score_name='hvg_score')
        monitor.checkpoint("HVG上Module Score成功")
    except Exception as e:
        print(f"  HVG Module Score错误: {e}")

# 测试4: 全量数据模拟
print("\n测试4: 模拟全量数据规模...")
total_cells = 57224
total_genes = 27933
estimated_dense_mb = (total_cells * total_genes * 8) / 1024 / 1024
print(f"  全量dense数组预估: {estimated_dense_mb:.1f} MB ({estimated_dense_mb/1024:.1f} GB)")

# 生成诊断报告
monitor.report()

print("\n" + "="*70)
print("诊断结论")
print("="*70)
print("""
问题根源: Module Score计算时scanpy内部尝试创建全基因集dense数组
解决方案:
  1. 预过滤基因池至3000个HVG
  2. 在子集上计算Module Score
  3. 将score写回原始adata

关键代码修改:
  - 先运行sc.pp.highly_variable_genes(adata, n_top_genes=3000)
  - 创建adata_subset = adata[:, adata.var['highly_variable']].copy()
  - 在adata_subset上计算module scores
  - 将结果转移回adata
""")
