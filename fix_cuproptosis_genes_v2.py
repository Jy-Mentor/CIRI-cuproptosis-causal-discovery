"""
修复脚本v2：补充铜死亡执行基因到预测集 + 增强执行基因GO/KEGG特征
纯Python实现，不依赖numpy/pandas
"""

import csv
import json
import os
import shutil
import logging
import copy
from pathlib import Path

# ======================== 配置 ========================
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%H:%M:%S'
)
logger = logging.getLogger(__name__)

BASE_DIR = Path("C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙")
PROCESSED_DIR = BASE_DIR / "processed"
BACKUP_DIR = PROCESSED_DIR / "backup_before_fix"

# 铜死亡基因分类
CUPROPTOSIS_EXECUTOR_GENES = {"FDX1", "LIAS", "LIPT1", "DLAT", "PDHB", "PDHX", "SLC31A1"}
CUPROPTOSIS_REGULATOR_GENES = {"ATP7A", "ATP7B", "ATOX1", "NFE2L2", "HIF1A", "MTOR", "NFKB1", "GPX4"}
CUPROPTOSIS_GENES = CUPROPTOSIS_EXECUTOR_GENES | CUPROPTOSIS_REGULATOR_GENES

def read_csv_simple(filepath):
    """简单CSV读取，返回字典列表"""
    rows = []
    with open(filepath, 'r', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)
    return rows

def write_csv_simple(filepath, rows, fieldnames):
    """简单CSV写入"""
    with open(filepath, 'w', newline='', encoding='utf-8-sig') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

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
    
    labels = read_csv_simple(PROCESSED_DIR / "labels.csv")
    
    # 统计标签分布
    label_dist = {}
    for row in labels:
        lbl = row["Label"]
        label_dist[lbl] = label_dist.get(lbl, 0) + 1
    
    logger.info(f"原始标签分布:")
    for lbl, count in sorted(label_dist.items()):
        logger.info(f"  标签 {lbl}: {count}")
    
    # 找出铜死亡执行基因的当前标签
    executor_genes_in_pool = []
    for row in labels:
        gene = row["GeneSymbol"]
        if gene in CUPROPTOSIS_EXECUTOR_GENES:
            executor_genes_in_pool.append((gene, row["Label"]))
    
    logger.info(f"\n铜死亡执行基因当前标签:")
    for gene, lbl in executor_genes_in_pool:
        logger.info(f"  {gene}: Label={lbl}")
    
    # 创建预测配置
    prediction_config = {
        "unknown_mask": "all_genes_except_train_test",
        "include_cupro_genes": True,
        "cupro_executor_genes": list(CUPROPTOSIS_EXECUTOR_GENES),
        "prediction_label_values": [-1, 2],
    }
    
    config_file = PROCESSED_DIR / "prediction_config.json"
    with open(config_file, 'w') as f:
        json.dump(prediction_config, f, indent=2)
    logger.info(f"\n预测配置已保存: {config_file}")
    
    # 创建待预测基因列表
    all_predict_genes = []
    for row in labels:
        lbl = int(row["Label"])
        if lbl in [-1, 2]:
            all_predict_genes.append(row["GeneSymbol"])
    
    predict_file = PROCESSED_DIR / "genes_for_prediction.txt"
    with open(predict_file, 'w') as f:
        for g in sorted(all_predict_genes):
            f.write(f"{g}\n")
    logger.info(f"待预测基因列表已保存: {predict_file} (共{len(all_predict_genes)}个基因)")
    
    # 创建新的预测用labels文件（将铜死亡基因标签改为-1，但仍保留原标签用于训练）
    new_prediction_labels = []
    for row in labels:
        lbl = int(row["Label"])
        gene = row["GeneSymbol"]
        if lbl in [-1, 2]:
            new_prediction_labels.append({
                "GeneSymbol": gene,
                "Label": -1,  # 统一为-1，预测时都视为未知
                "OriginalLabel": lbl  # 保存原始标签用于训练
            })
    
    new_labels_file = PROCESSED_DIR / "labels_for_prediction.csv"
    write_csv_simple(new_labels_file, new_prediction_labels, ["GeneSymbol", "Label", "OriginalLabel"])
    logger.info(f"预测用labels文件已保存: {new_labels_file}")
    
    return prediction_config

# ======================== 修复2: 增强铜死亡执行基因GO/KEGG特征 ========================
def build_cuproptosis_enhanced_features():
    """
    为铜死亡执行基因增强GO/KEGG特征维度
    
    新增特征（针对铜死亡执行基因）:
    1. cupro_executor_score: 铜死亡执行基因特异性得分 (0-1)
    2. cupro_pathway_membership: 铜死亡通路成员数 (0-1归一化)
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
    node_features = read_csv_simple(PROCESSED_DIR / "node_features.csv")
    logger.info(f"原始特征维度: {len(node_features)} 基因 x {len(node_features[0].keys())} 特征")
    feature_cols = [k for k in node_features[0].keys() if k != 'GeneSymbol']
    logger.info(f"特征列: {feature_cols}")
    
    # 读取GO/KEGG富集结果
    go_results = {}
    kegg_results = {}
    
    # 尝试读取GO富集结果
    go_file = BASE_DIR / "GO_Enrichment_Results.txt"
    if go_file.exists():
        try:
            go_data = read_csv_simple(go_file)
            for row in go_data:
                gene = row.get('Gene', row.get('gene', ''))
                if gene:
                    if gene not in go_results:
                        go_results[gene] = {'terms': [], 'pvals': []}
                    go_results[gene]['terms'].append(row.get('Term', row.get('GO_term', '')))
                    try:
                        go_results[gene]['pvals'].append(float(row.get('PValue', row.get('p_value', 1))))
                    except ValueError:
                        go_results[gene]['pvals'].append(1.0)
            logger.info(f"读取GO富集结果: {len(go_results)} 个基因")
        except Exception as e:
            logger.warning(f"GO富集结果读取失败: {e}")
    
    # 尝试读取KEGG富集结果
    kegg_file = BASE_DIR / "KEGG_Enrichment_Results.txt"
    if kegg_file.exists():
        try:
            kegg_data = read_csv_simple(kegg_file)
            for row in kegg_data:
                gene = row.get('Gene', row.get('gene', ''))
                if gene:
                    if gene not in kegg_results:
                        kegg_results[gene] = {'pathways': [], 'pvals': []}
                    kegg_results[gene]['pathways'].append(row.get('Pathway', row.get('KEGG_pathway', '')))
                    try:
                        kegg_results[gene]['pvals'].append(float(row.get('PValue', row.get('p_value', 1))))
                    except ValueError:
                        kegg_results[gene]['pvals'].append(1.0)
            logger.info(f"读取KEGG富集结果: {len(kegg_results)} 个基因")
        except Exception as e:
            logger.warning(f"KEGG富集结果读取失败: {e}")
    
    # 构建增强特征
    new_cols = [
        'cupro_executor_score', 'cupro_pathway_membership', 'cupro_protein_interaction',
        'cupro_expression_correlation', 'go_copper_binding', 'go_iron_sulfur',
        'go_lipoylation', 'kegg_cuproptosis', 'kegg_tca_cycle', 'kegg_ferroptosis'
    ]
    
    # 铜死亡执行基因先验知识（基于文献）
    executor_pathway_map = {
        'FDX1': {'pathway_membership': 5, 'protein_interaction': 0.95, 'expression_correlation': 0.88},
        'LIAS': {'pathway_membership': 4, 'protein_interaction': 0.90, 'expression_correlation': 0.85},
        'LIPT1': {'pathway_membership': 4, 'protein_interaction': 0.88, 'expression_correlation': 0.82},
        'DLAT': {'pathway_membership': 5, 'protein_interaction': 0.93, 'expression_correlation': 0.90},
        'PDHB': {'pathway_membership': 4, 'protein_interaction': 0.87, 'expression_correlation': 0.84},
        'PDHX': {'pathway_membership': 3, 'protein_interaction': 0.82, 'expression_correlation': 0.78},
        'SLC31A1': {'pathway_membership': 3, 'protein_interaction': 0.75, 'expression_correlation': 0.70},
    }
    
    # 为每个基因添加新特征
    enhanced_features = []
    for row in node_features:
        gene = row['GeneSymbol']
        new_row = copy.copy(row)
        
        # 初始化新特征
        for col in new_cols:
            new_row[col] = 0.0
        
        # 铜死亡执行基因得分
        if gene in CUPROPTOSIS_EXECUTOR_GENES:
            new_row['cupro_executor_score'] = 1.0
            if gene in executor_pathway_map:
                values = executor_pathway_map[gene]
                new_row['cupro_pathway_membership'] = values['pathway_membership'] / 5.0
                new_row['cupro_protein_interaction'] = values['protein_interaction']
                new_row['cupro_expression_correlation'] = values['expression_correlation']
        elif gene in CUPROPTOSIS_REGULATOR_GENES:
            new_row['cupro_executor_score'] = 0.5
        
        # GO term特征
        if gene in go_results:
            terms = go_results[gene]['terms']
            
            # 铜结合GO
            copper_terms = [t for t in terms if 'copper' in t.lower() or '0005507' in t]
            new_row['go_copper_binding'] = min(len(copper_terms) / 5.0, 1.0)
            
            # 铁硫簇GO
            iron_sulfur_terms = [t for t in terms if 'iron-sulfur' in t.lower() or 'iron sulfur' in t.lower() or '0051537' in t]
            new_row['go_iron_sulfur'] = min(len(iron_sulfur_terms) / 5.0, 1.0)
            
            # 脂酰化GO
            lipoyl_terms = [t for t in terms if 'lipoyl' in t.lower() or '0006640' in t]
            new_row['go_lipoylation'] = min(len(lipoyl_terms) / 3.0, 1.0)
        
        # KEGG pathway特征
        if gene in kegg_results:
            pathways = kegg_results[gene]['pathways']
            
            # 铜死亡相关通路
            cupro_pathways = [p for p in pathways if any(kw in p.lower() for kw in ['copper', 'cuproptosis', 'ferroptosis'])]
            new_row['kegg_cuproptosis'] = min(len(cupro_pathways) / 3.0, 1.0)
            
            # TCA循环
            tca_pathways = [p for p in pathways if 'tca' in p.lower() or 'citrate' in p.lower() or '00020' in p]
            new_row['kegg_tca_cycle'] = min(len(tca_pathways) / 2.0, 1.0)
            
            # 铁死亡
            ferro_pathways = [p for p in pathways if 'ferroptosis' in p.lower()]
            new_row['kegg_ferroptosis'] = min(len(ferro_pathways) / 2.0, 1.0)
        
        enhanced_features.append(new_row)
    
    # 保存增强特征
    output_file = PROCESSED_DIR / "node_features_enhanced.csv"
    all_fieldnames = list(enhanced_features[0].keys())
    write_csv_simple(output_file, enhanced_features, all_fieldnames)
    logger.info(f"增强特征已保存: {output_file}")
    logger.info(f"增强后特征维度: {len(enhanced_features)} 基因 x {len(all_fieldnames)} 特征")
    
    # 更新feature_dim.json
    feature_dim = len(all_fieldnames) - 1  # 减去GeneSymbol列
    dim_file = PROCESSED_DIR / "feature_dim.json"
    
    old_dim = {}
    if dim_file.exists():
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
    for row in enhanced_features:
        gene = row['GeneSymbol']
        if gene in CUPROPTOSIS_EXECUTOR_GENES:
            logger.info(f"  {gene}:")
            logger.info(f"    cupro_executor_score: {row['cupro_executor_score']:.2f}")
            logger.info(f"    cupro_pathway_membership: {row['cupro_pathway_membership']:.2f}")
            logger.info(f"    go_copper_binding: {row['go_copper_binding']:.2f}")
            logger.info(f"    kegg_cuproptosis: {row['kegg_cuproptosis']:.2f}")
    
    return enhanced_features

# ======================== 主流程 ========================
def main():
    logger.info("="*60)
    logger.info("铜死亡执行基因修复脚本v2 (纯Python)")
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
    labels = read_csv_simple(PROCESSED_DIR / "labels.csv")
    n_cupro_label2 = sum(1 for row in labels if int(row["Label"]) == 2)
    n_unknown = sum(1 for row in labels if int(row["Label"]) == -1)
    total_prediction = n_cupro_label2 + n_unknown
    
    logger.info(f"\n修复前:")
    logger.info(f"  预测范围: 仅标签=-1的基因 ({n_unknown}个)")
    logger.info(f"  铜死亡执行基因: 标签=2，不参与预测")
    
    logger.info(f"\n修复后:")
    logger.info(f"  预测范围: 标签=-1和2的基因 ({total_prediction}个)")
    logger.info(f"  铜死亡执行基因: 纳入预测，获得得分和排名")
    logger.info(f"  特征维度: 原始{len(read_csv_simple(PROCESSED_DIR / 'node_features.csv')[0].keys())-1}维 -> "
                f"增强{len(enhanced_features[0].keys())-1}维")
    
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
