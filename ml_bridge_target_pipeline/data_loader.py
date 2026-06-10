"""
数据加载模块
=============
负责从文件系统加载特征矩阵、药物靶点、疾病基因等数据，
构建训练标签并返回统一的数据结构。
"""

import os
import pandas as pd
import numpy as np
from dataclasses import dataclass, field
from typing import List, Set, Tuple

from .config import PathConfig
from .utils import log


@dataclass
class LoadedData:
    """统一的数据加载结果"""
    feature_df: pd.DataFrame = field(default_factory=pd.DataFrame)
    feature_cols: List[str] = field(default_factory=list)
    all_genes: Set[str] = field(default_factory=set)
    drug_targets: Set[str] = field(default_factory=set)
    disease_genes: Set[str] = field(default_factory=set)
    X: np.ndarray = field(default_factory=lambda: np.array([]))
    y_dt: np.ndarray = field(default_factory=lambda: np.array([]))
    y_dg: np.ndarray = field(default_factory=lambda: np.array([]))
    gene_index: List[str] = field(default_factory=list)
    unknown_mask: np.ndarray = field(default_factory=lambda: np.array([]))


class DataLoader:
    """
    数据加载器

    负责:
      1. 加载增强特征矩阵
      2. 加载药物靶点基因列表
      3. 加载疾病基因列表
      4. 可选子图基因过滤
      5. 构建 DT/DG 标签
    """

    def __init__(self, path_config: PathConfig):
        self.paths = path_config

    def load_data(self) -> LoadedData:
        """
        完整的数据加载与标签构建流程

        Returns:
            LoadedData 对象，包含特征矩阵、标签和元数据
        """
        log("=" * 70)
        log("STEP 1: 数据加载与标签构建")
        log("=" * 70)

        # 1a. 增强特征矩阵
        feat_df = self._load_feature_matrix()
        all_genes = set(feat_df.index)
        feature_cols = list(feat_df.columns)
        log(f"  特征矩阵: {len(all_genes)} genes × {len(feature_cols)} features")

        # 1b. 药物靶点
        drug_targets = self._load_gene_set(
            self.paths.drug_targets_path, "药物靶点"
        ) & all_genes

        # 1c. 疾病基因
        disease_genes = self._load_gene_set(
            self.paths.disease_genes_path, "疾病基因"
        ) & all_genes

        # 1d. 可选的子图基因过滤
        feat_df, all_genes, drug_targets, disease_genes = self._apply_subgraph_filter(
            feat_df, all_genes, drug_targets, disease_genes
        )

        # 1e. 构建标签
        feat_df['is_drug_target'] = feat_df.index.isin(drug_targets).astype(int)
        feat_df['is_disease_gene'] = feat_df.index.isin(disease_genes).astype(int)

        self._print_label_stats(feat_df)

        # 1f. 转换为训练数据
        X = feat_df[feature_cols].values.astype(np.float64)
        y_dt = feat_df['is_drug_target'].values.astype(int)
        y_dg = feat_df['is_disease_gene'].values.astype(int)
        gene_index = feat_df.index.tolist()
        unknown_mask = (y_dt == 0) & (y_dg == 0)

        return LoadedData(
            feature_df=feat_df,
            feature_cols=feature_cols,
            all_genes=all_genes,
            drug_targets=drug_targets,
            disease_genes=disease_genes,
            X=X,
            y_dt=y_dt,
            y_dg=y_dg,
            gene_index=gene_index,
            unknown_mask=unknown_mask,
        )

    def _load_feature_matrix(self) -> pd.DataFrame:
        """加载并预处理增强特征矩阵"""
        log(f"加载特征矩阵: {self.paths.feature_path}")
        feat_df = pd.read_csv(self.paths.feature_path)
        feat_df['gene_symbol'] = feat_df['gene_symbol'].str.upper()
        feat_df = feat_df.drop_duplicates(subset='gene_symbol', keep='first')
        feat_df = feat_df.set_index('gene_symbol')

        n_missing = feat_df.isnull().sum().sum()
        if n_missing > 0:
            log(f"  填充 {n_missing} 个缺失值 (列均值)")
            feat_df = feat_df.fillna(feat_df.mean())
            feat_df = feat_df.fillna(0.0)

        return feat_df

    def _load_gene_set(self, path: str, label: str) -> Set[str]:
        """从文本文件加载基因集合"""
        log(f"加载{label}: {path}")
        with open(path, 'r') as f:
            genes_raw = set(line.strip().upper() for line in f if line.strip())
        return genes_raw

    def _apply_subgraph_filter(
        self,
        feat_df: pd.DataFrame,
        all_genes: Set[str],
        drug_targets: Set[str],
        disease_genes: Set[str],
    ) -> Tuple[pd.DataFrame, Set[str], Set[str], Set[str]]:
        """可选的子图基因过滤"""
        subgraph_path = self.paths.subgraph_genes_path
        if not os.path.exists(subgraph_path):
            return feat_df, all_genes, drug_targets, disease_genes

        log(f"加载子图基因: {subgraph_path}")
        with open(subgraph_path, 'r') as f:
            sg = set(line.strip().upper() for line in f if line.strip())
            sg.discard('GENE')

        kept = all_genes & sg
        if not kept:
            log(f"  子图基因无交集，跳过过滤")
            return feat_df, all_genes, drug_targets, disease_genes

        log(f"  子图过滤: {len(all_genes)} → {len(kept)}")
        feat_df = feat_df.loc[list(kept)]
        new_genes = set(feat_df.index)
        drug_targets = drug_targets & new_genes
        disease_genes = disease_genes & new_genes
        return feat_df, new_genes, drug_targets, disease_genes

    def _print_label_stats(self, feat_df: pd.DataFrame):
        """打印标签统计信息"""
        n_dt = feat_df['is_drug_target'].sum()
        n_dg = feat_df['is_disease_gene'].sum()
        n_both = ((feat_df['is_drug_target'] == 1) & (feat_df['is_disease_gene'] == 1)).sum()
        n_unknown = len(feat_df) - n_dt - n_dg + n_both
        n_total = len(feat_df)

        log(f"  标签统计: DT+={n_dt} ({100 * n_dt / n_total:.1f}%), "
            f"DG+={n_dg} ({100 * n_dg / n_total:.1f}%), "
            f"Both={n_both}, Unknown={n_unknown}")