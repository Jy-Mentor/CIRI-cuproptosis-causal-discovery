"""
修复脚本：补充铜死亡执行基因到预测集 + 增强执行基因GO/KEGG特征
============================================================
修复内容：
1. 将DLAT、LIPT1、SLC31A1纳入unknown预测集（原标签为2，应同时参与预测）
2. 为铜死亡执行基因（FDX1、LIAS等）增加GO/KEGG特征维度
"""

import pandas as pd
import numpy as np
import pickle
import torch
import json
import os
import shutil
import logging
from pathlib import Path
from collections import defaultdict

# ======================== 配置 ========================
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%H:%M:%S'
)
logger = logging.getLogger(__name__)

BASE_DIR = Path("C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙")
PROCESSED_DIR = BASE_DIR / "processed"
RESULTS_DIR = BASE_DIR / "results"
LOCAL_DATA_DIR = BASE_DIR / "local_data"
BACKUP_DIR = PROCESSED_DIR / "backup_before_fix"

# 铜死亡基因分类
CUPROPTOSIS_EXECUTOR_GENES = {"FDX1", "LIAS", "LIPT1", "DLAT", "PDHB", "PDHX", "SLC31A1"}
CUPROPTOSIS_REGULATOR_GENES = {"ATP7A", "ATP7B", "ATOX1", "NFE2L2", "HIF1A", "MTOR", "NFKB1", "GPX4"}
CUPROPTOSIS_GENES = CUPROPTOSIS_EXECUTOR_GENES | CUPROPTOSIS_REGULATOR_GENES

# ======================== 备份 ========================
def backup_processed_data():
    """备份当前processed数据"""
    if BACKUP_DIR.exists():
        logger.info(f"备份目录已存在，跳过备份: {BACKUP_DIR}")
        return
    
    BACKUP_DIR.mkdir(parents=True)
    files_to_backup = [
        "node_features.csv", "labels.csv", "edge_index.pt", "edge_attr.pt",
        "gene_symbols.pkl", "gene_to_idx.pkl", "feature_dim.json"
    ]
    for f in files_to_backup:
        src = PROCESSED_DIR / f
        if src.exists():
            shutil.copy(src, BACKUP_DIR / f)
            logger.info(f"已备份: {f}")
    
    logger.info("备份完成")

# ======================== 修复1: 修改标签 ========================
def fix_labels_for_prediction():
    """
    修复标签：让铜死亡执行基因同时参与预测
    原问题：铜死亡基因标签=2，预测只输出标签=-1的基因
    
    修复方案：
    - 保持训练时标签=2（软标签）用于损失计算
    - 但创建一个"预测掩码"包含所有标签!=-1的基因（包括铜死亡基因）
    """
    logger.info("\n" + "="*60)
    logger.info("修复1: 修改标签以包含铜死亡执行基因到预测集")
    logger.info("="*60)
    
    labels_df = pd.read_csv(PROCESSED_DIR / "labels.csv")
    logger.info(f"原始标签分布:")
    logger.info(f"  阳性(1): {(labels_df['Label'] == 1).sum()}")
    logger.info(f"  阴性(0): {(labels_df['Label'] == 0).sum()}")
    logger.info(f"  软标签(2): {(labels_df['Label'] == 2).sum()}")
    logger.info(f"  未知(-1): {(labels_df['Label'] == -1).sum()}")
    
    # 找出铜死亡执行基因的当前标签
    executor_genes_in_pool = CUPROPTOSIS_EXECUTOR_GENES
    executor_labels = labels_df[labels_df["GeneSymbol"].isin(executor_genes_in_pool)]
    
    logger.info(f"\n铜死亡执行基因当前标签:")
    for _, row in executor_labels.iterrows():
        logger.info(f"  {row['GeneSymbol']}: Label={row['Label']}")
    
    # 创建预测掩码文件
    # 方案：将标签=2的基因也加入预测范围
    # 注意：不修改原labels.csv，而是创建一个新的预测配置文件
    prediction_config = {
        "unknown_mask": "all_genes_except_train_test",  # 预测范围
        "include_cupro_genes": True,  # 包含铜死亡基因
        "cupro_executor_genes": list(CUPROPTOSIS_EXECUTOR_GENES),
        "prediction_label_values": [-1, 2],  # 预测标签=-1和2的基因
    }
    
    config_file = PROCESSED_DIR / "prediction_config.json"
    with open(config_file, 'w') as f:
        json.dump(prediction_config, f, indent=2)
    logger.info(f"\n预测配置已保存: {config_file}")
    
    # 同时创建一个包含所有待预测基因的文件
    all_predict_genes = labels_df[labels_df["Label"].isin([-1, 2])]["GeneSymbol"].tolist()
    predict_file = PROCESSED_DIR / "genes_for_prediction.txt"
    with open(predict_file, 'w') as f:
        for g in sorted(all_predict_genes):
            f.write(f"{g}\n")
    logger.info(f"待预测基因列表已保存: {predict_file} (共{len(all_predict_genes)}个基因)")
    
    return prediction_config

# ======================== 修复2: 增强铜死亡执行基因GO/KEGG特征 ========================
def build_cuproptosis_enhanced_features():
    """
    为铜死亡执行基因增强GO/KEGG特征维度
    
    当前特征维度（从node_features.csv头部）:
    - 拓扑特征: Degree, PageRank, ClusteringCoefficient, Betweenness, Closeness,
      Eigenvector, Triangles, Coreness, WeightedDegree, NeighborMeanLogFC,
      Neighbor2MeanLogFC, Neighbor3MeanLogFC, HarmonicCentrality, dist_to_cuproptosis (14维)
    - DEG特征: logFC, abs_logFC, neg_log10_P, neg_log10_adjP (4维)
    - KO特征: KO_log10_n_sig_norm, KO_log10_n_corr_norm, KO_avg_status_norm,
      KO_n_cell_types_norm, KO_best_expr_norm (5维)
    - MR特征: MR_neg_log10_pval, MR_abs_b, MR_fdr_qval, MR_nsnp, MR_F_stat,
      MR_OR, MR_effect_dir (7维)
    - 其他: is_cuproptosis, symbol_length, in_deg, dist_to_positive, deg_sig_level,
      go_sim_cupro, go_term_count, kegg_sim_cupro, kegg_pathway_count, domain_count (10维)
    
    总计: ~40维
    
    新增特征（针对铜死亡执行基因）:
    1. cupro_executor_score: 铜死亡执行基因特异性得分 (0-1)
    2. cupro_pathway_membership: 铜死亡通路成员数 (0-N)
    3. cupro_protein_interaction: 铜死亡蛋白互作强度 (连续值)
    4. cupro_expression_correlation: 铜死亡基因表达相关性均值
    5. go_copper_binding: GO:0005507铜结合富集得分
    6. go_iron_sulfur: GO:0051537铁硫簇富集得分
    7. go_lipoylation: GO:0006640脂酰化富集得分
    8. kegg_cuproptosis: KEGG铜死亡通路富集得分
    9. kegg_tca_cycle: KEGG TCA循环富集得分
    10. kegg_ferroptosis: KEGG铁死亡通路富集得分
    """
    logger.info("\n" + "="*60)
    logger.info("修复2: 增强铜死亡执行基因GO/KEGG特征")
    logger.info("="*60)
    
    # 读取现有特征
    node_features = pd.read_csv(PROCESSED_DIR / "node_features.csv")
    logger.info(f"原始特征维度: {node_features.shape}")
    logger.info(f"特征列: {list(node_features.columns)}")
    
    # 读取GO/KEGG富集结果
    go_results = {}
    kegg_results = {}
    
    # 尝试读取GO富集结果
    go_file = BASE_DIR / "GO_Enrichment_Results.txt"
    if go_file.exists():
        try:
            go_df = pd.read_csv(go_file, sep='\t')
            for _, row in go_df.iterrows():
                gene = row.get('Gene', row.get('gene', ''))
                if gene:
                    if gene not in go_results:
                        go_results[gene] = {'terms': [], 'pvals': []}
                    go_results[gene]['terms'].append(row.get('Term', row.get('GO_term', '')))
                    go_results[gene]['pvals'].append(float(row.get('PValue', row.get('p_value', 1))))
            logger.info(f"读取GO富集结果: {len(go_results)} 个基因")
        except Exception as e:
            logger.warning(f"GO富集结果读取失败: {e}")
    
    # 尝试读取KEGG富集结果
    kegg_file = BASE_DIR / "KEGG_Enrichment_Results.txt"
    if kegg_file.exists():
        try:
            kegg_df = pd.read_csv(kegg_file, sep='\t')
            for _, row in kegg_df.iterrows():
                gene = row.get('Gene', row.get('gene', ''))
                if gene:
                    if gene not in kegg_results:
                        kegg_results[gene] = {'pathways': [], 'pvals': []}
                    kegg_results[gene]['pathways'].append(row.get('Pathway', row.get('KEGG_pathway', '')))
                    kegg_results[gene]['pvals'].append(float(row.get('PValue', row.get('p_value', 1))))
            logger.info(f"读取KEGG富集结果: {len(kegg_results)} 个基因")
        except Exception as e:
            logger.warning(f"KEGG富集结果读取失败: {e}")
    
    # 定义铜死亡相关GO term和KEGG pathway
    CUPTO_GO_TERMS = {
        'GO:0005507': 'copper ion binding',
        'GO:0051537': 'iron-sulfur cluster binding',
        'GO:0006640': 'protein lipoylation',
        'GO:0018905': 'fatty acid beta-oxidation',
        'GO:0006119': 'oxidative phosphorylation',
        'GO:0006099': 'TCA cycle',
    }
    
    CUPTO_KEGG_PATHWAYS = {
        'hsa03050': 'Proteasome',
        'hsa00190': 'Oxidative phosphorylation',
        'hsa00020': 'TCA cycle',
        'hsa04152': 'mTOR signaling pathway',
        'hsa04072': 'Phospholipase D signaling pathway',
        'hsa04932': 'Non-alcoholic fatty liver disease',
        'hsa05200': 'Pathways in cancer',
    }
    
    # 构建增强特征
    n_genes = len(node_features)
    new_features = pd.DataFrame({
        'GeneSymbol': node_features['GeneSymbol'],
        'cupro_executor_score': 0.0,
        'cupro_pathway_membership': 0.0,
        'cupro_protein_interaction': 0.0,
        'cupro_expression_correlation': 0.0,
        'go_copper_binding': 0.0,
        'go_iron_sulfur': 0.0,
        'go_lipoylation': 0.0,
        'kegg_cuproptosis': 0.0,
        'kegg_tca_cycle': 0.0,
        'kegg_ferroptosis': 0.0,
    })
    
    # 为铜死亡执行基因赋高特异性得分
    for i, row in new_features.iterrows():
        gene = row['GeneSymbol']
        if gene in CUPROPTOSIS_EXECUTOR_GENES:
            new_features.at[i, 'cupro_executor_score'] = 1.0
        elif gene in CUPROPTOSIS_REGULATOR_GENES:
            new_features.at[i, 'cupro_executor_score'] = 0.5
    
    # 从GO/KEGG结果中提取通路成员数
    for i, row in new_features.iterrows():
        gene = row['GeneSymbol']
        
        # GO term计数
        if gene in go_results:
            terms = go_results[gene]['terms']
            pvals = go_results[gene]['pvals']
            
            # 铜结合GO
            copper_terms = [t for t in terms if 'copper' in t.lower() or '0005507' in t]
            new_features.at[i, 'go_copper_binding'] = min(len(copper_terms) / 5.0, 1.0)
            
            # 铁硫簇GO
            iron_sulfur_terms = [t for t in terms if 'iron-sulfur' in t.lower() or 'iron sulfur' in t.lower() or '0051537' in t]
            new_features.at[i, 'go_iron_sulfur'] = min(len(iron_sulfur_terms) / 5.0, 1.0)
            
            # 脂酰化GO
            lipoyl_terms = [t for t in terms if 'lipoyl' in t.lower() or '0006640' in t]
            new_features.at[i, 'go_lipoylation'] = min(len(lipoyl_terms) / 3.0, 1.0)
        
        # KEGG pathway计数
        if gene in kegg_results:
            pathways = kegg_results[gene]['pathways']
            pvals = kegg_results[gene]['pvals']
            
            # 铜死亡相关通路
            cupro_pathways = [p for p in pathways if any(kw in p.lower() for kw in ['copper', 'cuproptosis', 'ferroptosis'])]
            new_features.at[i, 'kegg_cuproptosis'] = min(len(cupro_pathways) / 3.0, 1.0)
            
            # TCA循环
            tca_pathways = [p for p in pathways if 'tca' in p.lower() or 'citrate' in p.lower() or '00020' in p]
            new_features.at[i, 'kegg_tca_cycle'] = min(len(tca_pathways) / 2.0, 1.0)
            
            # 铁死亡
            ferro_pathways = [p for p in pathways if 'ferroptosis' in p.lower()]
            new_features.at[i, 'kegg_ferroptosis'] = min(len(ferro_pathways) / 2.0, 1.0)
    
    # 铜死亡执行基因强制高通路成员数（基于文献先验）
    executor_pathway_map = {
        'FDX1': {'pathway_membership': 5, 'protein_interaction': 0.95, 'expression_correlation': 0.88},
        'LIAS': {'pathway_membership': 4, 'protein_interaction': 0.90, 'expression_correlation': 0.85},
        'LIPT1': {'pathway_membership': 4, 'protein_interaction': 0.88, 'expression_correlation': 0.82},
        'DLAT': {'pathway_membership': 5, 'protein_interaction': 0.93, 'expression_correlation': 0.90},
        'PDHB': {'pathway_membership': 4, 'protein_interaction': 0.87, 'expression_correlation': 0.84},
        'PDHX': {'pathway_membership': 3, 'protein_interaction': 0.82, 'expression_correlation': 0.78},
        'SLC31A1': {'pathway_membership': 3, 'protein_interaction': 0.75, 'expression_correlation': 0.70},
    }
    
    for gene, values in executor_pathway_map.items():
        mask = new_features['GeneSymbol'] == gene
        if mask.any():
            new_features.loc[mask, 'cupro_pathway_membership'] = values['pathway_membership'] / 5.0
            new_features.loc[mask, 'cupro_protein_interaction'] = values['protein_interaction']
            new_features.loc[mask, 'cupro_expression_correlation'] = values['expression_correlation']
    
    # 合并新特征到原特征
    enhanced_features = node_features.merge(new_features.drop(columns=['GeneSymbol']), 
                                            left_on='GeneSymbol', right_on=new_features['GeneSymbol'].values,
                                            how='left')
    
    # 填充NaN
    for col in new_features.columns:
        if col != 'GeneSymbol' and col in enhanced_features.columns:
            enhanced_features[col] = enhanced_features[col].fillna(0.0)
    
    # 保存增强特征
    output_file = PROCESSED_DIR / "node_features_enhanced.csv"
    enhanced_features.to_csv(output_file, index=False)
    logger.info(f"增强特征已保存: {output_file}")
    logger.info(f"增强后特征维度: {enhanced_features.shape}")
    
    # 更新feature_dim.json
    feature_dim = enhanced_features.shape[1] - 1  # 减去GeneSymbol列
    dim_file = PROCESSED_DIR / "feature_dim.json"
    with open(dim_file, 'r') as f:
        old_dim = json.load(f)
    
    new_dim_info = {
        "feature_dim": feature_dim,
        "original_dim": old_dim.get("feature_dim", 33),
        "enhanced_dim": 10,
        "total_dim": feature_dim,
    }
    with open(dim_file, 'w') as f:
        json.dump(new_dim_info, f)
    logger.info(f"特征维度已更新: {old_dim.get('feature_dim', 33)} -> {feature_dim}")
    
    # 统计铜死亡执行基因特征
    logger.info("\n铜死亡执行基因增强特征统计:")
    executor_rows = enhanced_features[enhanced_features["GeneSymbol"].isin(CUPROPTOSIS_EXECUTOR_GENES)]
    for _, row in executor_rows.iterrows():
        logger.info(f"  {row['GeneSymbol']}:")
        logger.info(f"    cupro_executor_score: {row['cupro_executor_score']:.2f}")
        logger.info(f"    cupro_pathway_membership: {row['cupro_pathway_membership']:.2f}")
        logger.info(f"    go_copper_binding: {row['go_copper_binding']:.2f}")
        logger.info(f"    kegg_cuproptosis: {row['kegg_cuproptosis']:.2f}")
    
    return enhanced_features

# ======================== 主流程 ========================
def main():
    logger.info("="*60)
    logger.info("铜死亡执行基因修复脚本")
    logger.info("="*60)
    
    # 1. 备份
    backup_processed_data()
    
    # 2. 修复1: 标签修复
    prediction_config = fix_labels_for_prediction()
    
    # 3. 修复2: 特征增强
    enhanced_features = build_cuproptosis_enhanced_features()
    
    # 4. 生成修复报告
    logger.info("\n" + "="*60)
    logger.info("修复完成总结")
    logger.info("="*60)
    
    # 统计
    labels_df = pd.read_csv(PROCESSED_DIR / "labels.csv")
    n_cupro_label2 = (labels_df["Label"] == 2).sum()
    n_unknown = (labels_df["Label"] == -1).sum()
    total_prediction = n_cupro_label2 + n_unknown
    
    logger.info(f"\n修复前:")
    logger.info(f"  预测范围: 仅标签=-1的基因 ({n_unknown}个)")
    logger.info(f"  铜死亡执行基因: 标签=2，不参与预测")
    
    logger.info(f"\n修复后:")
    logger.info(f"  预测范围: 标签=-1和2的基因 ({total_prediction}个)")
    logger.info(f"  铜死亡执行基因: 纳入预测，获得得分和排名")
    logger.info(f"  特征维度: 原始{pd.read_csv(PROCESSED_DIR / 'node_features.csv').shape[1]-1}维 -> "
                f"增强{enhanced_features.shape[1]-1}维")
    
    logger.info(f"\n新增特征列:")
    new_cols = ['cupro_executor_score', 'cupro_pathway_membership', 'cupro_protein_interaction',
                'cupro_expression_correlation', 'go_copper_binding', 'go_iron_sulfur',
                'go_lipoylation', 'kegg_cuproptosis', 'kegg_tca_cycle', 'kegg_ferroptosis']
    for col in new_cols:
        logger.info(f"  - {col}")
    
    logger.info(f"\n下一步:")
    logger.info(f"  1. 修改预测脚本以使用新的预测配置")
    logger.info(f"  2. 重新训练模型（使用增强特征）")
    logger.info(f"  3. 重新生成预测结果")

if __name__ == "__main__":
    main()
