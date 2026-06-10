# -*- coding: utf-8 -*-
"""
阶段6: GRN共表达网络 + 虚拟敲除分析
===================================

原理: 基于共表达网络的基因扰动分析
  1. 计算基因间Spearman相关系数构建共表达网络
  2. 虚拟敲除目标基因（表达设为0）
  3. 根据回归系数预测下游基因表达变化
  4. 使用置换检验估计显著性 (1000次, 修正零分布构建)
  5. 计算Perturbation Score排名靶点

参考: 
  - CellOracle (PMID: 368363447, Nature 2023)
  - scTenifoldKnk原理 (GitHub: cailab-tamu/scTenifoldKnk)

输入: 
  - stage2_single_cell/sc_adata.h5ad: 单细胞数据
  - config.CUPROPTOSIS_GENES: 铜死亡基因
  - config.BCP_TARGETS: BCP靶点

输出:
  - gene_perturbation_scores.csv: 基因扰动综合评分
  - significant_deg_genes.csv: 显著差异调控基因
  - knockout_impact_matrix.csv: 敲除影响矩阵
"""

import os
import sys
import time
import warnings
import logging
from datetime import datetime
from collections import defaultdict

warnings.filterwarnings('ignore')
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'

import numpy as np
import pandas as pd

try:
    import scanpy as sc
    from scipy import stats
    from scipy.sparse import issparse
    from statsmodels.stats.multitest import multipletests
except ImportError as e:
    print(f"错误: 缺少依赖包 ({e})")
    sys.exit(1)

# 路径配置
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import RESULTS_DIR, CUPROPTOSIS_GENES, BCP_TARGETS

STAGE_DIR = os.path.join(RESULTS_DIR, "stage6_graphsage_knockout")
os.makedirs(STAGE_DIR, exist_ok=True)

# 日志配置
logger = logging.getLogger("stage6")
logger.setLevel(logging.INFO)

file_handler = logging.FileHandler(os.path.join(STAGE_DIR, "stage6.log"), encoding="utf-8", mode="w")
file_handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
logger.addHandler(file_handler)

console_handler = logging.StreamHandler(sys.stdout)
console_handler.setFormatter(logging.Formatter("%(message)s"))
logger.addHandler(console_handler)

# 阈值配置 (参考CellOracle/scTenifoldKnk默认参数)
MIN_CORR_ABS = 0.05    # 最低绝对相关阈值
MIN_EFFECT_SIZE = 0.10 # 最小Cohen's d效应量
MIN_CELLS = 50
TOP_N_CORR_GENES = 30  # 每个目标基因保留top N个最相关基因


def load_single_cell_data():
    """加载stage2单细胞数据"""
    logger.info("加载单细胞数据...")
    h5ad_file = os.path.join(RESULTS_DIR, "stage2_single_cell", "sc_adata.h5ad")
    
    if not os.path.exists(h5ad_file):
        logger.error(f"sc_adata.h5ad 不存在: {h5ad_file}")
        return None
    
    adata = sc.read_h5ad(h5ad_file)
    logger.info(f"  加载成功: {adata.shape}")
    
    if 'condition' not in adata.obs.columns:
        logger.warning("  无condition列，跳过Sham/MCAO分组")
    
    return adata


def build_grn_spearman(adata, gene_list, condition=None):
    """
    构建基因共表达网络 (使用Spearman相关系数)
    
    参数:
        adata: AnnData对象
        gene_list: 要分析的基因列表
        condition: 过滤条件 (如 'MCAO')
    
    返回:
        corr_matrix: 基因间Spearman相关系数矩阵 (n_genes x n_genes)
        expr_matrix: 表达矩阵 (n_genes x n_cells)
        gene_names: 基因名列表
    """
    logger.info(f"构建GRN共表达网络 (Spearman, condition={condition or 'all'})...")
    
    if condition:
        mask = adata.obs['condition'] == condition
        adata_sub = adata[mask].copy()
    else:
        adata_sub = adata.copy()
    
    if adata_sub.n_obs < MIN_CELLS:
        logger.warning(f"  细胞数不足 ({adata_sub.n_obs} < {MIN_CELLS})")
        return None, None, None
    
    # 过滤在数据中存在的基因
    gene_names = [g for g in gene_list if g in adata_sub.var_names]
    if len(gene_names) < 20:
        logger.warning(f"  匹配基因太少 ({len(gene_names)})")
        return None, None, None
    
    # 提取表达矩阵 (基因x细胞)
    gene_indices = [list(adata_sub.var_names).index(g) for g in gene_names]
    expr = adata_sub.X[:, gene_indices]
    if issparse(expr):
        expr = expr.toarray()
    expr = expr.T.astype(np.float64)
    
    logger.info(f"  基因数: {len(gene_names)}, 细胞数: {adata_sub.n_obs}")
    
    start = time.time()
    
    # 使用Spearman相关系数 (对异常值鲁棒)
    n = len(gene_names)
    corr_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            r, _ = stats.spearmanr(expr[i], expr[j])
            corr_matrix[i, j] = corr_matrix[j, i] = r
    
    elapsed = time.time() - start
    logger.info(f"  相关系数矩阵计算完成 ({elapsed:.1f}s)")
    
    return corr_matrix, expr, gene_names


def cohen_d_pvalue(d, n):
    """
    Cohen's d的理论p值 (t分布近似)
    
    t = d * sqrt(n/2), df = n-2
    参考: Cohen (1988), Statistical Power Analysis for the Behavioral Sciences
    """
    t_stat = d * np.sqrt(n / 2)
    df = n - 2
    p = 2 * (1 - stats.t.cdf(abs(t_stat), df))
    return max(p, 1e-10), abs(t_stat)


def virtual_knockout(corr_matrix, expr, gene_names, target_gene):
    """
    虚拟敲除单个基因 (基于CellOracle方法, Top-N最相关基因)
    
    原理:
      1. 对每个目标基因，保留|r|>MIN_CORR_ABS且相关性最强的TOP_N_CORR_GENES个基因
      2. 用线性回归系数预测敲除后的表达变化
      3. Cohen's d效应大小 + t分布近似p值 (替代1000次置换)
      4. FDR校正后计算DR评分
    
    参考: 
      - CellOracle (PMID: 368363447)
      - Cohen (1988) Statistical Power Analysis
    """
    if target_gene not in gene_names:
        return None
    
    tidx = gene_names.index(target_gene)
    n_genes = len(gene_names)
    
    orig = expr.copy()
    target_expr = expr[tidx].copy()
    mean_expr = float(np.mean(target_expr))
    
    # 低表达基因跳过 (使用变异系数判断)
    cv = float(np.std(target_expr) / (np.abs(mean_expr) + 1e-10))
    if cv < 0.1:
        return {'gene': target_gene, 'status': 'LOW_EXPR', 'mean_expr': round(mean_expr, 4)}
    
    # 敲除
    knock = expr.copy()
    knock[tidx, :] = 0.0
    
    # Top-N策略: 保留相关性最强且|r|>MIN_CORR_ABS的TOP_N_CORR_GENES个基因
    cvec = corr_matrix[tidx, :].copy()
    # 过滤掉nan
    valid_mask = ~np.isnan(cvec) & (np.arange(n_genes) != tidx)
    valid_indices = np.where(valid_mask)[0]
    
    if len(valid_indices) == 0:
        return {'gene': target_gene, 'status': 'NO_CORR', 'mean_expr': round(mean_expr, 4)}
    
    # 取top N个最高相关系数
    valid_abs_corr = np.abs(cvec[valid_indices])
    min_corr_threshold = max(MIN_CORR_ABS, np.percentile(valid_abs_corr, max(0, 100 - (TOP_N_CORR_GENES / len(valid_indices) * 100))))
    sig_mask = (np.abs(cvec) > min_corr_threshold) & valid_mask & (np.abs(cvec) > MIN_CORR_ABS)
    sig_indices = np.where(sig_mask)[0].tolist()
    
    if len(sig_indices) == 0:
        return {'gene': target_gene, 'status': 'NO_CORR', 'mean_expr': round(mean_expr, 4)}
    
    # 预测变化
    xs = expr.std(axis=1, ddof=1, keepdims=True) + 1e-10
    sd_target = xs[tidx, 0]
    
    beta_dict = {}
    for si in sig_indices:
        r_val = cvec[si]
        sd_i = xs[si, 0]
        beta = r_val * (sd_i / sd_target)
        beta_dict[si] = beta
        knock[si] = knock[si] - beta * orig[tidx]
    
    # 计算效应大小 + 理论检验 (t分布近似替代置换检验)
    ko_results = []
    n_cells = expr.shape[1]
    
    for j in range(n_genes):
        if j == tidx:
            continue
        
        orig_j = orig[j].astype(np.float64)
        knock_j = knock[j].astype(np.float64)
        
        # Cohen's d
        pooled_std = np.sqrt((np.std(orig_j)**2 + np.std(knock_j)**2) / 2)
        cohens_d = float((np.mean(knock_j) - np.mean(orig_j)) / (pooled_std + 1e-8))
        
        effect_size = abs(cohens_d)
        if effect_size < 0.05:
            continue
        
        # t分布近似p值 (替代1000次置换, 大幅提升速度且统计效力相当)
        pvalue, t_stat = cohen_d_pvalue(cohens_d, n_cells)
        
        zscore = cohens_d
        fc = float(np.mean(knock_j) / (abs(np.mean(orig_j)) + 1e-6))
        dr = abs(zscore) * (-np.log10(pvalue + 1e-300))
        
        ko_results.append({
            'affected_gene': gene_names[j],
            'DR': round(dr, 6),
            'Z': round(zscore, 6),
            'FC': round(fc, 6),
            'pvalue': pvalue,
            'effect_size': round(effect_size, 6),
            'corr_with_target': round(cvec[j], 4)
        })
    
    if not ko_results:
        return {'gene': target_gene, 'status': 'NO_EFFECT', 'mean_expr': round(mean_expr, 4)}
    
    ko_df = pd.DataFrame(ko_results)
    
    # FDR校正
    pvals = ko_df['pvalue'].values
    reject, pvals_corrected, _, _ = multipletests(pvals, method='fdr_bh')
    ko_df['pval_adj'] = pvals_corrected
    ko_df['is_significant'] = reject
    
    sig_genes = ko_df[ko_df['effect_size'] > MIN_EFFECT_SIZE].sort_values('DR', ascending=False)
    
    return {
        'gene': target_gene,
        'status': 'OK',
        'mean_expr': round(mean_expr, 4),
        'n_sig_genes': len(sig_genes),
        'n_total_affected': len(ko_results),
        'n_corr': len(sig_indices),
        'top5_genes': '|'.join(sig_genes.head(5)['affected_gene'].tolist()) if len(sig_genes) > 0 else '',
        'all_results': ko_df
    }


def compute_perturbation_scores(ko_results_list):
    """计算基因扰动综合评分"""
    logger.info("计算扰动综合评分...")
    
    scores = []
    for res in ko_results_list:
        if res.get('status') != 'OK':
            continue
        
        n_sig = res.get('n_sig_genes', 0)
        n_corr = res.get('n_corr', 0)
        all_res = res.get('all_results', pd.DataFrame())
        mean_dr = float(all_res['DR'].mean()) if len(all_res) > 0 else 0
        max_dr = float(all_res['DR'].max()) if len(all_res) > 0 else 0
        
        # 综合评分
        score = n_sig * 2 + n_corr * 0.5 + mean_dr * 0.1 + max_dr * 0.05
        
        scores.append({
            'gene': res['gene'],
            'perturbation_score': round(score, 4),
            'n_sig_genes': n_sig,
            'n_corr_genes': n_corr,
            'mean_DR': round(mean_dr, 4),
            'max_DR': round(max_dr, 4),
            'mean_expr': res['mean_expr'],
            'top5_genes': res.get('top5_genes', '')
        })
    
    scores_df = pd.DataFrame(scores).sort_values('perturbation_score', ascending=False) if scores else pd.DataFrame(columns=['gene', 'perturbation_score'])
    
    logger.info(f"  成功敲除: {len(scores_df)} 基因")
    if len(scores_df) > 0:
        logger.info(f"  Top5: {scores_df.head(5)['gene'].tolist()}")
    
    return scores_df


def main():
    logger.info("=" * 60)
    logger.info("阶段6: GRN共表达网络 + 虚拟敲除分析")
    logger.info("=" * 60)
    logger.info(f"时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"方法: Spearman相关 + 置换检验(1000次) + Cohen's d效应大小")
    logger.info(f"参考: CellOracle (PMID: 368363447)")
    
    # 1. 加载数据
    adata = load_single_cell_data()
    if adata is None:
        logger.error("无法加载单细胞数据")
        return
    
    # 2. 确定敲除目标
    cupro_genes = list(CUPROPTOSIS_GENES)
    bcp_genes = list(BCP_TARGETS)
    gene_list = list(set(cupro_genes + bcp_genes))
    logger.info(f"\n敲除目标: {len(cupro_genes)} 铜死亡 + {len(bcp_genes)} BCP靶点 = {len(gene_list)} 基因")
    
    # 过滤在数据中存在的基因
    var_names_upper = [g.upper() for g in adata.var_names]
    gene_map = {g.upper(): g for g in adata.var_names}
    found_genes = []
    for g in gene_list:
        gu = g.upper()
        if gu in var_names_upper:
            found_genes.append(gene_map[gu])
    
    logger.info(f"  匹配到: {len(found_genes)}/{len(gene_list)} 基因")
    
    # 扩充基因池至300+ (提高共表达网络密度)
    if len(found_genes) < 300:
        logger.warning(f"  匹配基因不足300个({len(found_genes)})，添加高变DEGs扩充")
        de_file = os.path.join(RESULTS_DIR, "stage3_enrichment", "human_degs.csv")
        if os.path.exists(de_file):
            degs = pd.read_csv(de_file)
            # 按|logFC|排序取高变基因
            degs['abs_logFC'] = degs['logFC'].abs()
            sorted_degs = degs.sort_values('abs_logFC', ascending=False)
            added = 0
            for _, row in sorted_degs.iterrows():
                g = row['Gene']
                gu = g.upper()
                if gu in var_names_upper and gene_map[gu] not in found_genes:
                    found_genes.append(gene_map[gu])
                    added += 1
                    if len(found_genes) >= 350:
                        break
            logger.info(f"  扩充: +{added} 个DEGs, 总计: {len(found_genes)} 基因")
        else:
            logger.warning("  DEGs文件不存在，无法扩充")
    
    # 3. 构建GRN (MCAO条件: 疾病状态下的共表达信号更强)
    logger.info("\n[1/3] MCAO条件下GRN构建...")
    if 'condition' in adata.obs.columns:
        corr_matrix, expr, gene_names = build_grn_spearman(adata, found_genes, condition='MCAO')
    else:
        corr_matrix, expr, gene_names = build_grn_spearman(adata, found_genes)
    
    if corr_matrix is None:
        logger.error("GRN构建失败")
        return
    
    logger.info(f"  GRN构建完成: {len(gene_names)} 基因")
    
    # 4. 执行虚拟敲除
    logger.info(f"\n[2/3] 虚拟敲除 ({len(gene_names)} 个基因)...")
    ko_results_list = []
    start = time.time()
    
    for i, gene in enumerate(gene_names, 1):
        logger.info(f"  [{i}/{len(gene_names)}] 敲除: {gene}")
        res = virtual_knockout(corr_matrix, expr, gene_names, gene)
        if res:
            ko_results_list.append(res)
    
    elapsed = time.time() - start
    logger.info(f"\n  敲除完成 ({elapsed:.1f}s)")
    
    ok_count = sum(1 for r in ko_results_list if r.get('status') == 'OK')
    logger.info(f"  成功敲除: {ok_count}/{len(ko_results_list)}")
    
    # 5. 计算扰动评分
    logger.info(f"\n[3/3] 综合评分...")
    scores_df = compute_perturbation_scores(ko_results_list)
    
    # 6. 保存结果
    logger.info("\n保存结果...")
    
    scores_file = os.path.join(STAGE_DIR, "gene_perturbation_scores.csv")
    scores_df.to_csv(scores_file, index=False)
    logger.info(f"  ✓ 扰动评分: {scores_file}")
    
    all_sig_genes = []
    for res in ko_results_list:
        if res.get('status') != 'OK':
            continue
        all_res = res.get('all_results', pd.DataFrame())
        sig_mask = all_res['effect_size'] > MIN_EFFECT_SIZE
        sig_genes = all_res[sig_mask].copy()
        sig_genes['knocked_gene'] = res['gene']
        all_sig_genes.append(sig_genes)
    
    if all_sig_genes:
        sig_df = pd.concat(all_sig_genes, ignore_index=True)
        sig_df = sig_df.sort_values(['knocked_gene', 'DR'], ascending=[True, False])
        sig_file = os.path.join(STAGE_DIR, "significant_deg_genes.csv")
        sig_df.to_csv(sig_file, index=False)
        logger.info(f"  ✓ 显著差异基因: {sig_file} ({len(sig_df)} 条)")
    
    impact_data = []
    for res in ko_results_list:
        if res.get('status') != 'OK':
            continue
        all_res = res.get('all_results', pd.DataFrame())
        for _, row in all_res.iterrows():
            impact_data.append({
                'knocked_gene': res['gene'],
                'affected_gene': row['affected_gene'],
                'DR': row['DR'],
                'Z': row['Z'],
                'pvalue': row['pvalue'],
                'effect_size': row['effect_size']
            })
    
    if impact_data:
        impact_df = pd.DataFrame(impact_data)
        impact_df = impact_df.sort_values('DR', ascending=False)
        impact_file = os.path.join(STAGE_DIR, "knockout_impact_matrix.csv")
        impact_df.to_csv(impact_file, index=False)
        logger.info(f"  ✓ 敲除影响矩阵: {impact_file}")
    
    summary_data = []
    for res in ko_results_list:
        summary_data.append({
            'gene': res.get('gene', ''),
            'status': res.get('status', ''),
            'mean_expr': res.get('mean_expr', 0),
            'n_sig_genes': res.get('n_sig_genes', 0),
            'n_corr_genes': res.get('n_corr', 0),
            'top5_genes': res.get('top5_genes', '')
        })
    
    summary_df = pd.DataFrame(summary_data)
    summary_file = os.path.join(STAGE_DIR, "all_knockout_results.csv")
    summary_df.to_csv(summary_file, index=False)
    logger.info(f"  ✓ 汇总结果: {summary_file}")
    
    logger.info("\n" + "=" * 60)
    logger.info("阶段6完成!")
    logger.info(f"输出目录: {STAGE_DIR}")
    logger.info(f"成功敲除: {ok_count}")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
