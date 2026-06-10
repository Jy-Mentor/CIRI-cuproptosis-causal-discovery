"""
模型集成模块
=============
负责模型集成的计算，包括桥梁得分计算和RRF倒数排名融合。
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional

from .config import TrainingConfig
from .utils import log


class EnsembleIntegrator:
    """模型集成器 — 计算桥梁得分和RRF排名"""

    def __init__(self, training_config: TrainingConfig):
        self.rrf_k = training_config.rrf_k

    def compute_bridge_scores(
        self,
        drug_probas: Dict[str, np.ndarray],
        disease_probas: Dict[str, np.ndarray],
    ) -> Dict[str, np.ndarray]:
        """
        计算桥梁得分 = P_drug × P_disease

        Args:
            drug_probas: {model_name: drug_prob_array}
            disease_probas: {model_name: disease_prob_array}

        Returns:
            {model_name: bridge_score_array}
        """
        bridge = {}
        all_models = set(drug_probas.keys()) & set(disease_probas.keys())
        for model in all_models:
            bridge[model] = drug_probas[model] * disease_probas[model]
        return bridge

    def rrf_integrate(
        self,
        bridge_scores: Dict[str, np.ndarray],
        n_genes: int,
        unknown_mask: np.ndarray,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        RRF (Reciprocal Rank Fusion) 集成排名

        RRF_score(g) = Σ 1/(k + rank_i(g))

        Args:
            bridge_scores: {model_name: score_array}
            n_genes: 总基因数
            unknown_mask: 未知样本掩码

        Returns:
            (rrf_scores, rrf_ranks):
            - rrf_scores: RRF得分数组 (长度=n_genes)
            - rrf_ranks: 最终排名数组 (NaN表示非未知基因)
        """
        rrf_sum = np.zeros(n_genes, dtype=np.float64)

        for model, scores in bridge_scores.items():
            uk_scores = scores[unknown_mask]
            order = np.argsort(-uk_scores)  # 降序
            ranks = np.empty(len(uk_scores), dtype=np.float64)
            ranks[order] = np.arange(1, len(uk_scores) + 1)
            rrf_sum[unknown_mask] += 1.0 / (self.rrf_k + ranks)

        # 计算最终排名
        unknown_idx = np.where(unknown_mask)[0]
        rrf_sort_order = unknown_idx[np.argsort(-rrf_sum[unknown_mask])]
        rrf_ranks = np.full(n_genes, np.nan)
        for rank, idx in enumerate(rrf_sort_order):
            rrf_ranks[idx] = rank + 1

        return rrf_sum, rrf_ranks

    @staticmethod
    def extract_importance(model_name: str, clf, feature_cols: List[str]) -> Optional[pd.DataFrame]:
        """
        提取模型的特征重要性 (支持树模型和线性模型)

        Args:
            model_name: 模型名称
            clf: 训练好的分类器
            feature_cols: 特征列名列表

        Returns:
            特征重要性 DataFrame，或 None
        """
        try:
            if hasattr(clf, 'feature_importances_'):
                imp = clf.feature_importances_
            elif hasattr(clf, 'coef_'):
                imp = np.abs(clf.coef_).flatten()
            else:
                return None

            if len(imp) != len(feature_cols):
                return None

            df = pd.DataFrame({
                'model': model_name,
                'feature': feature_cols,
                'importance': imp,
            }).sort_values('importance', ascending=False).head(30)
            return df
        except Exception:
            return None

    @staticmethod
    def build_metrics_df(
        model_names: List[str],
        dt_metrics: Dict[str, Dict[str, float]],
        dg_metrics: Dict[str, Dict[str, float]],
    ) -> pd.DataFrame:
        """
        构建模型指标 DataFrame

        Args:
            model_names: 所有模型名称列表
            dt_metrics: DT任务指标
            dg_metrics: DG任务指标

        Returns:
            指标 DataFrame
        """
        metrics_rows = []
        for model in model_names:
            if model in dt_metrics:
                metrics_rows.append({
                    'model': model, 'task': 'DT',
                    'auroc': dt_metrics[model]['auroc'],
                    'auprc': dt_metrics[model]['auprc'],
                })
            if model in dg_metrics:
                metrics_rows.append({
                    'model': model, 'task': 'DG',
                    'auroc': dg_metrics[model]['auroc'],
                    'auprc': dg_metrics[model]['auprc'],
                })
        return pd.DataFrame(metrics_rows)

    @staticmethod
    def compare_with_gat(
        df_unknown: pd.DataFrame,
        gat_bridge_path: str,
    ) -> Dict:
        """
        与 GAT 桥梁基因结果进行对比分析

        Args:
            df_unknown: ML预测结果 DataFrame
            gat_bridge_path: GAT桥梁基因文件路径

        Returns:
            对比结果字典
        """
        import os
        from scipy.stats import spearmanr

        results = {
            'gat_genes': [],
            'overlap': [],
            'n_overlap': 0,
            'spearman_corr': None,
            'spearman_pval': None,
        }

        if not os.path.exists(gat_bridge_path):
            log(f"  GAT 桥梁基因文件未找到，跳过")
            return results

        try:
            gat_df = pd.read_csv(gat_bridge_path)
            gat_genes = list(gat_df['gene_symbol'].str.upper())
            results['gat_genes'] = gat_genes
            log(f"  GAT Top-20: {gat_genes}")

            overlap = []
            for g in gat_genes:
                if g in df_unknown['gene_symbol'].values:
                    rank_info = df_unknown[df_unknown['gene_symbol'] == g]
                    rrf_rank = rank_info['final_rank'].values[0]
                    overlap.append((g, int(rrf_rank)))

            results['overlap'] = overlap
            results['n_overlap'] = len(overlap)
            log(f"  GAT 基因在 ML 候选池中: {len(overlap)}/20")
            for g, r in overlap[:5]:
                log(f"    {g}: GAT Top-20, ML RRF rank={r}")

            # Spearman 相关性
            ml_ranks = {}
            for i, (_, row) in enumerate(df_unknown.iterrows()):
                ml_ranks[row['gene_symbol']] = i + 1

            gat_r, ml_r = [], []
            for rank_gat, g in enumerate(gat_genes):
                if g in ml_ranks:
                    gat_r.append(rank_gat + 1)
                    ml_r.append(ml_ranks[g])

            if len(gat_r) >= 5:
                corr, pval = spearmanr(gat_r, ml_r)
                results['spearman_corr'] = corr
                results['spearman_pval'] = pval
                log(f"  Spearman ρ = {corr:.4f} (p = {pval:.4f})")
            else:
                log(f"  共同基因不足5个，跳过 Spearman")

        except Exception as e:
            log(f"  [WARN] GAT 对比失败: {e}")

        return results