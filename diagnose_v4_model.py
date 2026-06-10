#!/usr/bin/env python3
"""
V4模型数据深度诊断脚本（无torch依赖版）
=========================================
目标：定位模型AUPRC=0.075（接近随机）的根本原因
诊断维度：
  1. 特征泄漏检测
  2. 特征方差分析
  3. 标签分布分析
  4. 铜死亡基因检查
  5. 训练日志分析
"""

import os
import sys
import json
import pickle
import logging
import warnings
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings('ignore')

# ---------------------------------------------------------------------------
# 配置
# ---------------------------------------------------------------------------
ROOT = Path(os.path.dirname(os.path.abspath(__file__)))
PROCESSED_DIR = ROOT / "processed"
DIAG_DIR = ROOT / "diagnostics_v4"
LOGS_DIR = ROOT / "logs_v4"
DIAG_DIR.mkdir(exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    handlers=[
        logging.FileHandler(DIAG_DIR / "diagnostic.log", encoding='utf-8'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

CUPROPTOSIS_GENES = {
    "FDX1", "LIAS", "LIPT1", "DLAT", "PDHB", "PDHX", "SLC31A1",
    "ATP7B", "ATOX1", "MTF1", "GLS", "CDKN2A", "NFKB1", "HIF1A", 
    "NFE2L2", "GPX4", "ATP7A", "MTOR"
}

# ---------------------------------------------------------------------------
# 1. 特征泄漏检测
# ---------------------------------------------------------------------------
def diagnose_feature_leakage(features_df: pd.DataFrame, labels_df: pd.DataFrame) -> Dict:
    """检测特征是否泄漏标签信息"""
    logger.info("\n" + "="*60)
    logger.info("【诊断1】特征泄漏检测")
    logger.info("="*60)
    
    feature_cols = [c for c in features_df.columns if c != "GeneSymbol"]
    labels = labels_df["Label"].values
    
    # 识别可疑特征名称
    leakage_keywords = ['dist_to', 'sim_', 'distance_to', 'similarity', 
                        'label', 'dist_']
    suspicious_cols = []
    for col in feature_cols:
        if any(kw in col.lower() for kw in leakage_keywords):
            suspicious_cols.append(col)
    
    logger.info(f"\n可疑泄漏特征 ({len(suspicious_cols)}个):")
    for col in suspicious_cols:
        logger.info(f"  ⚠️ {col}")
    
    # 计算每个特征与标签的相关性（仅对有标签节点）
    labeled_mask = labels >= 0
    y_binary = (labels[labeled_mask] == 1).astype(int)
    
    correlations = {}
    X_labeled = features_df.loc[labeled_mask, feature_cols].fillna(0).values
    
    for i, col in enumerate(feature_cols):
        if np.std(X_labeled[:, i]) > 1e-8 and len(X_labeled[:, i]) > 0:
            corr = np.corrcoef(X_labeled[:, i], y_binary)[0, 1]
            if not np.isnan(corr):
                correlations[col] = abs(corr)
    
    # 高相关特征（可能泄漏）
    high_corr = {k: v for k, v in sorted(correlations.items(), key=lambda x: x[1], reverse=True) if v > 0.15}
    
    logger.info(f"\n高相关特征 (|r|>0.15, {len(high_corr)}个):")
    for col, corr in list(high_corr.items())[:15]:
        marker = "⚠️ LEAKAGE" if col in suspicious_cols else ""
        logger.info(f"  {col}: |r|={corr:.4f} {marker}")
    
    # 详细对比可疑特征在不同标签组的分布
    if suspicious_cols:
        for col in suspicious_cols:
            if col not in feature_cols:
                continue
            
            col_idx = feature_cols.index(col)
            X_all = features_df[feature_cols].fillna(0).values
            
            mask_pos = (labels == 1)
            mask_neg = (labels == 0)
            mask_unknown = (labels == -1)
            
            val_pos = X_all[mask_pos, col_idx]
            val_neg = X_all[mask_neg, col_idx]
            val_unknown = X_all[mask_unknown, col_idx]
            
            logger.info(f"\n{col} 分布对比:")
            logger.info(f"  阳性(n={mask_pos.sum()}): mean={val_pos.mean():.6f}, std={val_pos.std():.6f}")
            logger.info(f"  阴性(n={mask_neg.sum()}): mean={val_neg.mean():.6f}, std={val_neg.std():.6f}")
            logger.info(f"  未知(n={mask_unknown.sum()}): mean={val_unknown.mean():.6f}, std={val_unknown.std():.6f}")
            
            # 统计检验
            if len(val_pos) > 1 and len(val_neg) > 1:
                try:
                    t_stat, p_val = stats.ttest_ind(val_pos, val_neg)
                    logger.info(f"  阳性vs阴性: t={t_stat:.4f}, p={p_val:.2e}")
                    if p_val < 0.05:
                        logger.info(f"  ⚠️ 两组存在显著差异 (p<0.05) - 可能是标签泄漏！")
                except:
                    pass
    
    return {
        "suspicious_cols": suspicious_cols,
        "high_corr_features": high_corr,
        "correlations": correlations
    }


# ---------------------------------------------------------------------------
# 2. 特征方差分析
# ---------------------------------------------------------------------------
def diagnose_feature_variance(features_df: pd.DataFrame) -> Dict:
    """分析特征的方差分布"""
    logger.info("\n" + "="*60)
    logger.info("【诊断2】特征方差分析")
    logger.info("="*60)
    
    feature_cols = [c for c in features_df.columns if c != "GeneSymbol"]
    X = features_df[feature_cols].fillna(0).values
    
    variances = np.var(X, axis=0)
    
    logger.info(f"\n方差统计:")
    logger.info(f"  最小方差: {variances.min():.8f}")
    logger.info(f"  最大方差: {variances.max():.6f}")
    logger.info(f"  中位数方差: {np.median(variances):.6f}")
    
    # 低方差特征
    low_var_threshold = 1e-6
    low_var_mask = variances < low_var_threshold
    low_var_features = [f for f, mask in zip(feature_cols, low_var_mask) if mask]
    
    if low_var_features:
        logger.warning(f"\n⚠️ 低方差特征 (var<{low_var_threshold}): {len(low_var_features)}个")
        for feat in low_var_features:
            logger.warning(f"  {feat}")
    
    # 特征分布类型
    logger.info(f"\n特征分布类型:")
    for i, col in enumerate(feature_cols):
        vals = X[:, i]
        unique_vals = len(np.unique(vals))
        unique_ratio = unique_vals / len(vals)
        
        if unique_ratio < 0.01:
            dist_type = "离散/常量"
        elif unique_ratio < 0.1:
            dist_type = "低多样性"
        else:
            dist_type = "连续"
        
        if col in ['degree', 'PageRank', 'ClusteringCoefficient', 'Betweenness', 'Closeness']:
            logger.info(f"  {col}: {dist_type} (unique={unique_vals}, ratio={unique_ratio:.4f})")
    
    return {
        "variances": dict(zip(feature_cols, variances)),
        "low_variance_features": low_var_features
    }


# ---------------------------------------------------------------------------
# 3. 标签分布分析
# ---------------------------------------------------------------------------
def diagnose_label_distribution(labels_df: pd.DataFrame, features_df: pd.DataFrame) -> Dict:
    """分析标签分布和类别不平衡问题"""
    logger.info("\n" + "="*60)
    logger.info("【诊断3】标签分布分析")
    logger.info("="*60)
    
    labels = labels_df["Label"].values
    
    total = len(labels)
    n_pos = (labels == 1).sum()
    n_neg = (labels == 0).sum()
    n_unknown = (labels == -1).sum()
    
    logger.info(f"\n总基因数: {total}")
    logger.info(f"阳性 (1): {n_pos} ({n_pos/total*100:.2f}%)")
    logger.info(f"阴性 (0): {n_neg} ({n_neg/total*100:.2f}%)")
    logger.info(f"未知 (-1): {n_unknown} ({n_unknown/total*100:.2f}%)")
    
    # 类别不平衡比
    if n_pos > 0:
        pos_neg_ratio = n_neg / n_pos
        logger.info(f"\n⚠️ 阴阳性比例: 1:{pos_neg_ratio:.1f}")
        if pos_neg_ratio > 10:
            logger.warning("严重类别不平衡！建议使用过采样或调整损失函数权重")
    
    # 阳性基因列表
    pos_genes = labels_df[labels_df["Label"] == 1]["GeneSymbol"].tolist()
    logger.info(f"\n阳性基因 (Top 30): {pos_genes[:30]}")
    
    # 阴性基因采样
    neg_genes = labels_df[labels_df["Label"] == 0]["GeneSymbol"].tolist()
    logger.info(f"阴性基因 (样本10): {neg_genes[:10]}")
    
    return {
        "total": total,
        "n_pos": n_pos,
        "n_neg": n_neg,
        "n_unknown": n_unknown,
        "pos_neg_ratio": n_neg / n_pos if n_pos > 0 else float('inf')
    }


# ---------------------------------------------------------------------------
# 4. 铜死亡基因专项诊断
# ---------------------------------------------------------------------------
def diagnose_cuproptosis_genes(features_df: pd.DataFrame, labels_df: pd.DataFrame) -> Dict:
    """诊断铜死亡基因的标签和特征情况"""
    logger.info("\n" + "="*60)
    logger.info("【诊断4】铜死亡基因专项诊断")
    logger.info("="*60)
    
    gene_to_label = dict(zip(labels_df["GeneSymbol"], labels_df["Label"]))
    
    results = []
    missing_genes = []
    
    for gene in sorted(CUPROPTOSIS_GENES):
        if gene in gene_to_label:
            label = gene_to_label[gene]
            label_name = "阳性" if label == 1 else ("阴性" if label == 0 else "未知")
            results.append({
                "gene": gene,
                "label": label,
                "label_name": label_name
            })
        else:
            missing_genes.append(gene)
    
    logger.info(f"\n铜死亡基因标签分布:")
    for result in results:
        logger.info(f"  {result['gene']}: {result['label_name']} (label={result['label']})")
    
    if missing_genes:
        logger.warning(f"\n⚠️ 缺失的铜死亡基因 ({len(missing_genes)}个): {missing_genes}")
    
    # 统计
    n_known = sum(1 for r in results if r['label'] != -1)
    n_unknown = sum(1 for r in results if r['label'] == -1)
    
    logger.info(f"\n铜死亡基因已知标签: {n_known}")
    logger.info(f"铜死亡基因未知标签: {n_unknown}")
    
    if n_unknown > 0:
        logger.info("⚠️ 部分铜死亡基因标签为-1（未知），需要检查是否应该加入训练集")
    
    return {
        "results": results,
        "missing_genes": missing_genes
    }


# ---------------------------------------------------------------------------
# 5. 训练日志分析
# ---------------------------------------------------------------------------
def diagnose_training_log(log_file: Path) -> Dict:
    """分析训练日志，诊断训练动态"""
    logger.info("\n" + "="*60)
    logger.info("【诊断5】训练日志分析")
    logger.info("="*60)
    
    if not log_file.exists():
        logger.warning(f"训练日志不存在: {log_file}")
        return {}
    
    log_df = pd.read_csv(log_file)
    
    logger.info(f"\n训练总Epoch: {len(log_df)}")
    logger.info(f"初始Train Loss: {log_df['train_loss'].iloc[0]:.6f}")
    logger.info(f"最终Train Loss: {log_df['train_loss'].iloc[-1]:.6f}")
    logger.info(f"最终Val AUPRC: {log_df['val_auprc'].iloc[-1]:.6f}")
    logger.info(f"最终Val F1: {log_df['val_f1'].iloc[-1]:.6f}")
    
    # 最佳性能
    best_auprc = log_df['val_auprc'].max()
    best_f1 = log_df['val_f1'].max()
    best_auprc_epoch = log_df['val_auprc'].idxmax() + 1
    best_f1_epoch = log_df['val_f1'].idxmax() + 1
    
    logger.info(f"\n最佳AUPRC: {best_auprc:.6f} (Epoch {best_auprc_epoch})")
    logger.info(f"最佳F1: {best_f1:.6f} (Epoch {best_f1_epoch})")
    
    # 随机基线对比
    random_auprc = 0.059  # 近似阳性比例
    logger.info(f"\n随机基线AUPRC: {random_auprc:.4f}")
    logger.info(f"最佳AUPRC/随机基线: {best_auprc/random_auprc:.2f}x")
    
    if best_auprc < random_auprc * 1.5:
        logger.warning("⚠️ 模型AUPRC仅略高于随机猜测，几乎未学到有效模式！")
    
    # AUPRC趋势
    first_10_auprc = log_df['val_auprc'].head(10).mean()
    last_10_auprc = log_df['val_auprc'].tail(10).mean()
    
    logger.info(f"\nAUPRC趋势:")
    logger.info(f"  前10个epoch平均: {first_10_auprc:.6f}")
    logger.info(f"  后10个epoch平均: {last_10_auprc:.6f}")
    
    if last_10_auprc < first_10_auprc * 0.9:
        logger.warning("⚠️ AUPRC在训练后期显著下降 - 严重过拟合！")
    
    # Loss趋势
    loss_values = log_df['train_loss'].values
    if loss_values[-1] > loss_values[0]:
        logger.warning("⚠️ Loss在训练结束时高于初始值 - 模型可能未正确优化！")
    
    # 早停分析
    patience_counter = 0
    early_stop_epoch = len(log_df)
    for i in range(1, len(log_df)):
        if log_df['val_auprc'].iloc[i] <= log_df['val_auprc'].iloc[:i].max():
            patience_counter += 1
        else:
            patience_counter = 0
        
        if patience_counter >= 50:
            early_stop_epoch = i + 1
            break
    
    logger.info(f"\n早停于Epoch: {early_stop_epoch}")
    logger.info(f"早停耐心: 50")
    
    return {
        "best_auprc": best_auprc,
        "best_f1": best_f1,
        "training_epochs": len(log_df),
        "early_stop_epoch": early_stop_epoch,
        "auprc_trend_declining": last_10_auprc < first_10_auprc * 0.9
    }


# ---------------------------------------------------------------------------
# 6. 选定特征分析
# ---------------------------------------------------------------------------
def diagnose_selected_features(features_df: pd.DataFrame, selected_features: List[str]) -> Dict:
    """分析选定特征的质量和合理性"""
    logger.info("\n" + "="*60)
    logger.info("【诊断6】选定特征分析")
    logger.info("="*60)
    
    logger.info(f"\n选定特征 ({len(selected_features)}个):")
    for i, feat in enumerate(selected_features, 1):
        logger.info(f"  {i}. {feat}")
    
    # 特征类型分类
    topo_features = []
    expr_features = []
    enrichment_features = []
    leakage_features = []
    
    for feat in selected_features:
        if any(kw in feat.lower() for kw in ['deg', 'centrality', 'page', 'cluster', 'coreness', 'harmonic']):
            topo_features.append(feat)
        elif any(kw in feat.lower() for kw in ['logfc', 'neg_log', 'p_', 'pval']):
            expr_features.append(feat)
        elif any(kw in feat.lower() for kw in ['go_', 'kegg_', 'domain', 'pathway']):
            enrichment_features.append(feat)
        elif any(kw in feat.lower() for kw in ['dist_', 'sim_']):
            leakage_features.append(feat)
    
    logger.info(f"\n特征类型分布:")
    logger.info(f"  拓扑特征: {len(topo_features)}个 - {topo_features}")
    logger.info(f"  表达特征: {len(expr_features)}个 - {expr_features}")
    logger.info(f"  富集特征: {len(enrichment_features)}个 - {enrichment_features}")
    logger.info(f"  ⚠️ 泄漏特征: {len(leakage_features)}个 - {leakage_features}")
    
    if leakage_features:
        logger.warning("\n⚠️ 发现基于标签计算的泄漏特征！这些特征会导致模型学习到虚假关联！")
        for feat in leakage_features:
            logger.warning(f"  {feat}")
    
    # 检查特征选择方法
    feature_cols = [c for c in features_df.columns if c != "GeneSymbol"]
    X = features_df[selected_features].fillna(0).values
    
    variances = np.var(X, axis=0)
    
    logger.info(f"\n选定特征方差:")
    for feat, var in sorted(zip(selected_features, variances), key=lambda x: x[1]):
        logger.info(f"  {feat}: {var:.6f}")
    
    return {
        "selected_features": selected_features,
        "topo_features": topo_features,
        "expr_features": expr_features,
        "enrichment_features": enrichment_features,
        "leakage_features": leakage_features,
        "variances": dict(zip(selected_features, variances))
    }


# ---------------------------------------------------------------------------
# 7. 综合诊断总结
# ---------------------------------------------------------------------------
def generate_summary_report(results: Dict) -> str:
    """生成诊断总结报告"""
    report = []
    report.append("="*60)
    report.append("V4模型诊断总结报告")
    report.append("="*60)
    report.append("")
    
    critical_issues = []
    warnings = []
    recommendations = []
    
    # 特征泄漏检查
    if 'leakage' in results and results['leakage'].get('suspicious_cols'):
        n_leakage = len(results['leakage']['suspicious_cols'])
        critical_issues.append(f"发现{n_leakage}个疑似标签泄漏特征")
        recommendations.append("移除dist_to_positive, go_sim_cupro等基于标签计算的特征")
    
    # 类别不平衡检查
    if 'labels' in results and results['labels'].get('pos_neg_ratio', 0) > 10:
        ratio = results['labels']['pos_neg_ratio']
        critical_issues.append(f"严重类别不平衡 (1:{ratio:.0f})")
        recommendations.append("使用Focal Loss (alpha>0.9)或阳性过采样")
    
    # 训练性能检查
    if 'training' in results:
        if results['training'].get('best_auprc', 0) < 0.1:
            critical_issues.append("模型AUPRC极低 (<0.1)，几乎未学到有效模式")
        
        if results['training'].get('auprc_trend_declining', False):
            warnings.append("AUPRC在训练后期下降，存在过拟合")
            recommendations.append("增加dropout、减少模型容量、使用早停")
    
    # 特征选择检查
    if 'features' in results and results['features'].get('leakage_features'):
        n = len(results['features']['leakage_features'])
        critical_issues.append(f"选定特征中包含{n}个泄漏特征")
    
    # 生成报告
    if critical_issues:
        report.append("🔴 严重问题:")
        for issue in critical_issues:
            report.append(f"  • {issue}")
        report.append("")
    
    if warnings:
        report.append("🟡 警告:")
        for warn in warnings:
            report.append(f"  • {warn}")
        report.append("")
    
    if recommendations:
        report.append("💡 修复建议:")
        for rec in recommendations:
            report.append(f"  • {rec}")
        report.append("")
    
    if not critical_issues and not warnings:
        report.append("✅ 未发现明显问题")
    
    return "\n".join(report)


# ---------------------------------------------------------------------------
# 主函数
# ---------------------------------------------------------------------------
def main():
    logger.info("="*60)
    logger.info("V4模型数据深度诊断开始")
    logger.info("="*60)
    
    results = {}
    
    # 加载数据
    logger.info("\n加载数据...")
    features = pd.read_csv(PROCESSED_DIR / "node_features.csv")
    labels = pd.read_csv(PROCESSED_DIR / "labels.csv")
    
    logger.info(f"基因数: {len(features)}")
    logger.info(f"特征数: {len(features.columns) - 1}")
    logger.info(f"标签分布: {labels['Label'].value_counts().sort_index().to_dict()}")
    
    # 加载选定特征
    try:
        with open(PROCESSED_DIR / "feature_cols_v4.pkl", 'rb') as f:
            selected_features = pickle.load(f)
        logger.info(f"选定特征数: {len(selected_features)}")
    except:
        selected_features = None
        logger.warning("无法加载选定特征")
    
    # 执行诊断
    # 1. 特征泄漏
    results['leakage'] = diagnose_feature_leakage(features, labels)
    
    # 2. 特征方差
    results['variance'] = diagnose_feature_variance(features)
    
    # 3. 标签分布
    results['labels'] = diagnose_label_distribution(labels, features)
    
    # 4. 铜死亡基因
    results['cuproptosis'] = diagnose_cuproptosis_genes(features, labels)
    
    # 5. 训练日志
    log_file = LOGS_DIR / "training_log_v4.csv"
    results['training'] = diagnose_training_log(log_file)
    
    # 6. 选定特征
    if selected_features:
        results['features'] = diagnose_selected_features(features, selected_features)
    
    # 保存结果
    with open(DIAG_DIR / "diagnostic_results.json", 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2, default=str, ensure_ascii=False)
    logger.info(f"\n诊断结果已保存: {DIAG_DIR / 'diagnostic_results.json'}")
    
    # 生成总结报告
    report = generate_summary_report(results)
    logger.info("\n" + report)
    
    # 保存报告
    with open(DIAG_DIR / "diagnostic_report.txt", 'w', encoding='utf-8') as f:
        f.write(report)
    logger.info(f"\n诊断报告已保存: {DIAG_DIR / 'diagnostic_report.txt'}")
    
    logger.info("\n诊断完成！")


if __name__ == "__main__":
    main()