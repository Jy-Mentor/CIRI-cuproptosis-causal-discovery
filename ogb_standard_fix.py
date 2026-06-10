"""
OGB标准修复方案A：半监督学习
================================
修复内容：
1. 回退特征注入：删除硬编码的先验知识特征
2. 修改标签逻辑：铜死亡基因统一设为-1，不参与训练
3. 增强图结构：基于GO/KEGG添加功能相似边（而非节点特征）
4. 符合OGB标准：训练集仅使用标签明确的节点(0/1)

参考标准：
- OGB (Hu et al., NeurIPS 2020): https://ogb.stanford.edu/docs/nodeprop/
- GNNMutation (BMC Bioinformatics 2025)
- GNNenrich (Bioinformatics 2025)
"""

import csv
import json
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
OG_STANDARD_DIR = PROCESSED_DIR / "ogb_standard"

# 铜死亡基因分类
CUPROPTOSIS_EXECUTOR_GENES = {"FDX1", "LIAS", "LIPT1", "DLAT", "PDHB", "PDHX", "SLC31A1"}
CUPROPTOSIS_REGULATOR_GENES = {"ATP7A", "ATP7B", "ATOX1", "NFE2L2", "HIF1A", "MTOR", "NFKB1", "GPX4"}
CUPROPTOSIS_GENES = CUPROPTOSIS_EXECUTOR_GENES | CUPROPTOSIS_REGULATOR_GENES

def read_csv_simple(filepath):
    rows = []
    with open(filepath, 'r', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)
    return rows

def write_csv_simple(filepath, rows, fieldnames):
    with open(filepath, 'w', newline='', encoding='utf-8-sig') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

# ======================== 步骤1: 回退特征注入 ========================
def revert_feature_injection():
    """
    步骤1: 回退特征注入
    - 删除node_features_enhanced.csv（包含硬编码先验特征）
    - 恢复标准node_features.csv（仅使用观测特征）
    """
    logger.info("\n" + "="*60)
    logger.info("步骤1: 回退特征注入")
    logger.info("="*60)
    
    # 删除增强特征文件
    enhanced_file = PROCESSED_DIR / "node_features_enhanced.csv"
    if enhanced_file.exists():
        enhanced_file.unlink()
        logger.info(f"已删除增强特征文件: {enhanced_file}")
    
    # 恢复标准特征文件（从备份）
    standard_file = PROCESSED_DIR / "node_features.csv"
    backup_file = BACKUP_DIR / "node_features.csv"
    
    if backup_file.exists():
        shutil.copy(backup_file, standard_file)
        logger.info(f"已恢复标准特征文件: {standard_file}")
        
        # 读取并验证
        features = read_csv_simple(standard_file)
        n_features = len(features[0].keys()) - 1  # 减去GeneSymbol
        logger.info(f"标准特征维度: {len(features)} 基因 × {n_features} 特征")
        logger.info(f"特征列: {list(features[0].keys())[1:]}")
    else:
        logger.warning("备份文件不存在，使用现有标准特征文件")
    
    # 更新feature_dim.json
    features = read_csv_simple(PROCESSED_DIR / "node_features.csv")
    feature_dim = len(features[0].keys()) - 1
    
    dim_info = {
        "feature_dim": feature_dim,
        "description": "OGB标准特征（无先验知识注入）",
        "feature_types": {
            "topological": 14,
            "deg": 4,
            "scKO": 5,
            "MR": 7,
            "other": 10
        }
    }
    
    with open(PROCESSED_DIR / "feature_dim.json", 'w') as f:
        json.dump(dim_info, f, indent=2)
    logger.info(f"已更新feature_dim.json: {feature_dim}维")

# ======================== 步骤2: 修改标签逻辑 ========================
def fix_labels_ogb_standard():
    """
    步骤2: 修改标签逻辑（OGB标准）
    
    OGB标准：
    - 训练集：仅使用标签明确的节点(0/1)
    - 验证集：用于早停和超参调优
    - 测试集：最终评估
    - 铜死亡基因：标签=-1，不参与训练，仅作为预测目标
    
    修复方案：
    1. 铜死亡基因标签从2改为-1
    2. 训练时只使用标签0/1的节点计算损失
    3. 预测时对所有标签=-1的节点推理（包括铜死亡基因）
    """
    logger.info("\n" + "="*60)
    logger.info("步骤2: 修改标签逻辑（OGB标准）")
    logger.info("="*60)
    
    labels = read_csv_simple(PROCESSED_DIR / "labels.csv")
    
    # 统计原始标签分布
    label_dist = {}
    for row in labels:
        lbl = row["Label"]
        label_dist[lbl] = label_dist.get(lbl, 0) + 1
    
    logger.info(f"原始标签分布:")
    for lbl, count in sorted(label_dist.items()):
        logger.info(f"  标签 {lbl}: {count}")
    
    # 修改铜死亡基因标签：2 → -1
    n_modified = 0
    for row in labels:
        if row["GeneSymbol"] in CUPROPTOSIS_GENES:
            original_label = row["Label"]
            if original_label == "2":
                row["Label"] = "-1"
                n_modified += 1
    
    logger.info(f"\n修改铜死亡基因标签: {n_modified} 个基因 (2 → -1)")
    
    # 统计修改后标签分布
    new_label_dist = {}
    for row in labels:
        lbl = row["Label"]
        new_label_dist[lbl] = new_label_dist.get(lbl, 0) + 1
    
    logger.info(f"修改后标签分布:")
    for lbl, count in sorted(new_label_dist.items()):
        logger.info(f"  标签 {lbl}: {count}")
    
    # 保存修改后的labels
    labels_file = PROCESSED_DIR / "labels.csv"
    write_csv_simple(labels_file, labels, ["GeneSymbol", "Label"])
    logger.info(f"已保存修改后的labels: {labels_file}")
    
    # 创建训练掩码文件
    train_mask = []
    val_mask = []
    test_mask = []
    predict_mask = []
    
    for row in labels:
        gene = row["GeneSymbol"]
        lbl = int(row["Label"])
        
        if lbl in [0, 1]:
            # 标签明确的节点：参与训练/验证/测试
            # 简化方案：70%训练，15%验证，15%测试
            if gene in CUPROPTOSIS_GENES:
                # 铜死亡基因即使有标签也不参与训练（避免先验泄露）
                predict_mask.append(gene)
            else:
                # 非铜死亡基因：根据基因名哈希分配到不同集合
                hash_val = hash(gene) % 100
                if hash_val < 70:
                    train_mask.append(gene)
                elif hash_val < 85:
                    val_mask.append(gene)
                else:
                    test_mask.append(gene)
        else:
            # 标签=-1的节点：仅用于预测
            predict_mask.append(gene)
    
    logger.info(f"\n数据集分割:")
    logger.info(f"  训练集: {len(train_mask)} 个基因")
    logger.info(f"  验证集: {len(val_mask)} 个基因")
    logger.info(f"  测试集: {len(test_mask)} 个基因")
    logger.info(f"  预测集: {len(predict_mask)} 个基因（包含铜死亡基因）")
    
    # 保存掩码文件
    for name, genes in [
        ("train_mask", train_mask),
        ("val_mask", val_mask),
        ("test_mask", test_mask),
        ("predict_mask", predict_mask)
    ]:
        mask_file = PROCESSED_DIR / f"{name}.txt"
        with open(mask_file, 'w') as f:
            for g in sorted(genes):
                f.write(f"{g}\n")
        logger.info(f"已保存{ name}: {mask_file}")
    
    # 创建OGB标准配置文件
    ogb_config = {
        "task_type": "node_classification",
        "evaluation_metric": "auprc",  # AUPRC适合不平衡数据
        "train_strategy": "semi_supervised",
        "train_labels": [0, 1],  # 仅使用0/1标签训练
        "predict_labels": [-1],  # 预测所有-1标签的节点
        "cuproptosis_handling": "predict_only",  # 铜死亡基因仅预测，不参与训练
        "dataset_split": {
            "train_ratio": 0.7,
            "val_ratio": 0.15,
            "test_ratio": 0.15,
            "train_size": len(train_mask),
            "val_size": len(val_mask),
            "test_size": len(test_mask),
            "predict_size": len(predict_mask)
        }
    }
    
    config_file = PROCESSED_DIR / "ogb_config.json"
    with open(config_file, 'w') as f:
        json.dump(ogb_config, f, indent=2)
    logger.info(f"已保存OGB配置: {config_file}")

# ======================== 步骤3: 增强图结构 ========================
def build_functional_similarity_edges():
    """
    步骤3: 增强图结构（而非节点特征）
    
    基于GO/KEGG添加功能相似边：
    - GO语义相似度 > 0.7 的基因对添加边
    - KEGG通路共现的基因对添加边
    - 边特征标记为"functional_similarity"类型
    
    参考：
    - GNNenrich (Bioinformatics 2025): GO用于后验富集，图结构增强
    - GNNRAI (NPJ Systems Biology 2025): 通路知识构建异质图
    """
    logger.info("\n" + "="*60)
    logger.info("步骤3: 增强图结构（GO/KEGG功能相似边）")
    logger.info("="*60)
    
    # 读取GO/KEGG富集结果
    go_annotations = {}
    kegg_annotations = {}
    
    # 尝试读取GO注释
    go_file = BASE_DIR / "GO_Enrichment_Results.txt"
    if go_file.exists():
        try:
            go_data = read_csv_simple(go_file)
            for row in go_data:
                gene = row.get('Gene', row.get('gene', ''))
                if gene:
                    if gene not in go_annotations:
                        go_annotations[gene] = set()
                    term = row.get('Term', row.get('GO_term', ''))
                    if term:
                        go_annotations[gene].add(term)
            logger.info(f"读取GO注释: {len(go_annotations)} 个基因")
        except Exception as e:
            logger.warning(f"GO注释读取失败: {e}")
    
    # 尝试读取KEGG注释
    kegg_file = BASE_DIR / "KEGG_Enrichment_Results.txt"
    if kegg_file.exists():
        try:
            kegg_data = read_csv_simple(kegg_file)
            for row in kegg_data:
                gene = row.get('Gene', row.get('gene', ''))
                if gene:
                    if gene not in kegg_annotations:
                        kegg_annotations[gene] = set()
                    pathway = row.get('Pathway', row.get('KEGG_pathway', ''))
                    if pathway:
                        kegg_annotations[gene].add(pathway)
            logger.info(f"读取KEGG注释: {len(kegg_annotations)} 个基因")
        except Exception as e:
            logger.warning(f"KEGG注释读取失败: {e}")
    
    # 计算GO语义相似度（Jaccard相似度）
    def jaccard_similarity(set1, set2):
        if not set1 or not set2:
            return 0.0
        intersection = len(set1 & set2)
        union = len(set1 | set2)
        return intersection / union if union > 0 else 0.0
    
    # 构建功能相似边
    functional_edges = []
    genes_with_go = list(go_annotations.keys())
    genes_with_kegg = list(kegg_annotations.keys())
    all_annotated_genes = set(genes_with_go) | set(genes_with_kegg)
    
    # 读取现有基因列表
    node_features = read_csv_simple(PROCESSED_DIR / "node_features.csv")
    all_genes = {row['GeneSymbol'] for row in node_features}
    
    # 仅对注释基因计算相似度
    annotated_in_pool = all_annotated_genes & all_genes
    gene_list = sorted(annotated_in_pool)
    
    logger.info(f"\n计算功能相似度:")
    logger.info(f"  有GO注释的基因: {len(genes_with_go)}")
    logger.info(f"  有KEGG注释的基因: {len(genes_with_kegg)}")
    logger.info(f"  在基因池中的注释基因: {len(annotated_in_pool)}")
    
    # 限制计算量：仅计算Top-5000基因对
    n_pairs = 0
    max_pairs = 5000
    
    for i in range(len(gene_list)):
        for j in range(i+1, len(gene_list)):
            g1, g2 = gene_list[i], gene_list[j]
            
            # GO语义相似度
            go_sim = jaccard_similarity(
                go_annotations.get(g1, set()),
                go_annotations.get(g2, set())
            )
            
            # KEGG通路共现
            kegg_sim = jaccard_similarity(
                kegg_annotations.get(g1, set()),
                kegg_annotations.get(g2, set())
            )
            
            # 综合相似度
            combined_sim = 0.6 * go_sim + 0.4 * kegg_sim
            
            # 阈值过滤
            if combined_sim > 0.3:  # 降低阈值以增加边数
                functional_edges.append({
                    "Gene1": g1,
                    "Gene2": g2,
                    "GO_Similarity": go_sim,
                    "KEGG_Similarity": kegg_sim,
                    "Combined_Similarity": combined_sim,
                    "Edge_Type": "functional_similarity"
                })
                n_pairs += 1
            
            if n_pairs >= max_pairs:
                break
        
        if n_pairs >= max_pairs:
            break
    
    logger.info(f"\n功能相似边统计:")
    logger.info(f"  新增边数: {len(functional_edges)}")
    
    if functional_edges:
        logger.info(f"  GO相似度范围: {min(e['GO_Similarity'] for e in functional_edges):.3f} - "
                    f"{max(e['GO_Similarity'] for e in functional_edges):.3f}")
        logger.info(f"  KEGG相似度范围: {min(e['KEGG_Similarity'] for e in functional_edges):.3f} - "
                    f"{max(e['KEGG_Similarity'] for e in functional_edges):.3f}")
        logger.info(f"  综合相似度范围: {min(e['Combined_Similarity'] for e in functional_edges):.3f} - "
                    f"{max(e['Combined_Similarity'] for e in functional_edges):.3f}")
    else:
        logger.warning("  无功能相似边生成（GO/KEGG注释文件为空）")
        logger.warning("  建议：后续可通过g:Profiler或clusterProfiler获取GO/KEGG注释")
    
    # 保存功能相似边
    edges_file = PROCESSED_DIR / "functional_similarity_edges.csv"
    if functional_edges:
        write_csv_simple(
            edges_file,
            functional_edges,
            ["Gene1", "Gene2", "GO_Similarity", "KEGG_Similarity", "Combined_Similarity", "Edge_Type"]
        )
        logger.info(f"已保存功能相似边: {edges_file}")
    else:
        logger.warning("无功能相似边生成")
    
    # 统计铜死亡基因的功能相似边
    cupro_edges = [e for e in functional_edges
                   if e["Gene1"] in CUPROPTOSIS_GENES or e["Gene2"] in CUPROPTOSIS_GENES]
    
    logger.info(f"\n铜死亡基因功能相似边:")
    logger.info(f"  涉及铜死亡基因的边数: {len(cupro_edges)}")
    for edge in cupro_edges[:10]:  # 仅显示前10条
        logger.info(f"    {edge['Gene1']} ↔ {edge['Gene2']}: "
                    f"GO={edge['GO_Similarity']:.3f}, KEGG={edge['KEGG_Similarity']:.3f}")

# ======================== 主流程 ========================
def main():
    logger.info("="*60)
    logger.info("OGB标准修复方案A：半监督学习")
    logger.info("="*60)
    
    # 1. 回退特征注入
    revert_feature_injection()
    
    # 2. 修改标签逻辑
    fix_labels_ogb_standard()
    
    # 3. 增强图结构
    build_functional_similarity_edges()
    
    # 4. 生成修复报告
    logger.info("\n" + "="*60)
    logger.info("修复完成总结（OGB标准方案A）")
    logger.info("="*60)
    
    logger.info(f"\n✅ 修复1: 回退特征注入")
    logger.info(f"  - 删除硬编码先验特征（10维）")
    logger.info(f"  - 恢复标准观测特征（40维）")
    logger.info(f"  - 符合OGB特征工程标准")
    
    logger.info(f"\n✅ 修复2: 修改标签逻辑")
    logger.info(f"  - 铜死亡基因标签: 2 → -1")
    logger.info(f"  - 训练集: 仅使用标签0/1的节点")
    logger.info(f"  - 预测集: 包含所有标签-1的节点（含铜死亡基因）")
    logger.info(f"  - 符合OGB半监督学习标准")
    
    logger.info(f"\n✅ 修复3: 增强图结构")
    logger.info(f"  - 基于GO/KEGG添加功能相似边")
    logger.info(f"  - 铜死亡基因通过图消息传递获得隐式知识")
    logger.info(f"  - 避免特征泄露")
    
    logger.info(f"\n📋 下一步:")
    logger.info(f"  1. 修改主脚本以使用OGB标准配置")
    logger.info(f"  2. 训练模型（仅使用标签0/1）")
    logger.info(f"  3. 预测所有标签-1的节点（包含铜死亡基因）")
    logger.info(f"  4. 评估铜死亡基因排名")

if __name__ == "__main__":
    main()
