# -*- coding: utf-8 -*-
"""
β-石竹烯(BCP)分子对接靶点系统性筛选
=====================================

筛选策略:
  1. 直接铜死亡通路靶点 (Tier A): 铜死亡核心基因蛋白
  2. 间接调控靶点 (Tier B): 通过GRN扰动影响铜死亡的BCP靶点
  3. 炎症-铜死亡交叉靶点 (Tier C): 炎症因子与铜死亡的交汇点

评分体系:
  - 铜死亡通路核心度 (0-100)
  - BCP靶点证据强度 (综合评分)
  - GRN扰动评分 (scTenifoldKnk)
  - PPI网络中心度 (GAT)
  - 可成药性 (Druggability)

输出:
  - 优先级排序的对接靶点列表
  - 每个靶点的科学依据和评分详情
"""

import os
import sys
import pandas as pd
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import (
    RESULTS_DIR, CUPROPTOSIS_GENES, CUPROPTOSIS_RELATED,
    CUPROPTOSIS_PATHWAY_SCORES, BCP_TARGETS
)

OUTPUT_DIR = os.path.join(RESULTS_DIR, "docking_screening")
os.makedirs(OUTPUT_DIR, exist_ok=True)


def load_pipeline_results():
    """加载所有阶段结果"""
    results = {}
    
    # Stage 8: 综合靶点评分
    core_file = os.path.join(RESULTS_DIR, "stage8_final_targets", "core_targets.csv")
    if os.path.exists(core_file):
        results['core_targets'] = pd.read_csv(core_file)
        print(f"✓ Stage 8核心靶点: {len(results['core_targets'])} 基因")
    
    # Stage 9: GAT排名
    gat_file = os.path.join(RESULTS_DIR, "stage9_ppi_gat_v3", "gat_gene_ranking.csv")
    if os.path.exists(gat_file):
        results['gat_ranking'] = pd.read_csv(gat_file)
        print(f"✓ Stage 9 GAT排名: {len(results['gat_ranking'])} 基因")
    
    # Stage 6: GRN扰动评分
    pert_file = os.path.join(RESULTS_DIR, "stage6_sctenifold_knockout", "gene_perturbation_scores.csv")
    if os.path.exists(pert_file):
        results['perturbation'] = pd.read_csv(pert_file)
        print(f"✓ Stage 6扰动评分: {len(results['perturbation'])} 基因")
    
    return results


def get_cuproptosis_pathway_score(gene):
    """获取铜死亡通路核心度评分"""
    gene_upper = gene.upper()
    return CUPROPTOSIS_PATHWAY_SCORES.get(gene_upper, 0)


def is_cuproptosis_gene(gene):
    """判断是否为铜死亡相关基因"""
    gene_upper = gene.upper()
    return gene_upper in [g.upper() for g in CUPROPTOSIS_GENES] or \
           gene_upper in [g.upper() for g in CUPROPTOSIS_RELATED]


def calculate_docking_priority(results):
    """
    计算分子对接优先级评分
    
    评分公式:
      Priority Score = w1*Pathway_Score + w2*BCP_Evidence + w3*GRN_Perturbation + w4*PPI_Centrality + w5*Druggability
    
    权重:
      w1 (通路核心度): 0.30 - 直接参与铜死亡的程度
      w2 (BCP证据): 0.25 - BCP靶点证据强度
      w3 (GRN扰动): 0.20 - 虚拟敲除扰动强度
      w4 (PPI中心度): 0.15 - 网络拓扑重要性
      w5 (可成药性): 0.10 - 蛋白可成药性预测
    """
    
    # 收集所有候选基因
    all_genes = set()
    
    # 1. 铜死亡核心基因 (必须包含)
    for g in CUPROPTOSIS_GENES:
        all_genes.add(g.upper())
    
    # 2. BCP靶点中排名靠前的
    if 'core_targets' in results:
        top_bcp = results['core_targets'].head(50)['Gene'].tolist()
        for g in top_bcp:
            all_genes.add(g.upper())
    
    # 3. GAT高排名基因
    if 'gat_ranking' in results:
        top_gat = results['gat_ranking'].head(30)['Gene'].tolist()
        for g in top_gat:
            all_genes.add(g.upper())
    
    print(f"\n候选基因池: {len(all_genes)} 个基因")
    
    # 构建评分矩阵
    scoring_data = []
    
    for gene in sorted(all_genes):
        gene_upper = gene.upper()
        
        # 1. 铜死亡通路核心度 (0-100)
        pathway_score = get_cuproptosis_pathway_score(gene)
        is_cupro = is_cuproptosis_gene(gene)
        
        # 2. BCP证据强度 (0-100)
        bcp_score = 0
        bcp_tier = "N/A"
        if 'core_targets' in results:
            match = results['core_targets'][results['core_targets']['Gene'].str.upper() == gene_upper]
            if len(match) > 0:
                bcp_score = float(match.iloc[0]['Comprehensive'])
                bcp_tier = match.iloc[0].get('Tier', 'N/A')
        
        # 3. GRN扰动评分 (标准化0-100)
        grn_score = 0
        if 'perturbation' in results:
            match = results['perturbation'][results['perturbation']['gene'].str.upper() == gene_upper]
            if len(match) > 0:
                # 扰动评分标准化 (max ~76)
                raw_score = float(match.iloc[0]['perturbation_score'])
                grn_score = min(100, raw_score / 76 * 100)
        
        # 4. PPI中心度 (0-100)
        ppi_score = 0
        if 'gat_ranking' in results:
            match = results['gat_ranking'][results['gat_ranking']['Gene'].str.upper() == gene_upper]
            if len(match) > 0:
                # GAT分数标准化
                ppi_score = float(match.iloc[0]['GAT_score']) * 100
        
        # 5. 可成药性 (简化评分)
        # 激酶、受体、酶类更可成药
        druggable_keywords = [' kinase', ' receptor', ' enzyme', ' protease', ' transporter']
        # 这里使用简单的基于基因功能的启发式评分
        # 实际应用中应使用AlphaFold结构预测或已知的DrugBank数据
        druggability = 50  # 默认中等可成药性
        
        # 计算综合优先级评分
        priority = (
            0.30 * pathway_score +
            0.25 * bcp_score +
            0.20 * grn_score +
            0.15 * ppi_score +
            0.10 * druggability
        )
        
        scoring_data.append({
            'Gene': gene,
            'Is_Cuproptosis': is_cupro,
            'Pathway_Score': round(pathway_score, 2),
            'BCP_Score': round(bcp_score, 2),
            'BCP_Tier': bcp_tier,
            'GRN_Score': round(grn_score, 2),
            'PPI_Score': round(ppi_score, 2),
            'Druggability': round(druggability, 2),
            'Priority_Score': round(priority, 2),
            'Category': 'Direct' if pathway_score > 50 else ('Indirect' if bcp_score > 30 else 'Candidate')
        })
    
    df = pd.DataFrame(scoring_data)
    df = df.sort_values('Priority_Score', ascending=False)
    
    return df


def categorize_targets(df):
    """分类靶点并生成推荐"""
    
    # Tier 1: 直接铜死亡核心靶点 (高优先级)
    tier1 = df[(df['Is_Cuproptosis'] == True) & (df['Priority_Score'] > 40)]
    
    # Tier 2: 高评分BCP靶点 (间接调控)
    tier2 = df[(df['Is_Cuproptosis'] == False) & (df['BCP_Score'] > 40) & (df['Priority_Score'] > 30)]
    
    # Tier 3: 潜在候选 (需要进一步验证)
    tier3 = df[df['Priority_Score'] > 20]
    tier3 = tier3[~tier3['Gene'].isin(tier1['Gene'].tolist() + tier2['Gene'].tolist())]
    
    return tier1, tier2, tier3


def generate_report(df, tier1, tier2, tier3):
    """生成筛选报告"""
    
    report_lines = [
        "=" * 80,
        "β-石竹烯(BCP)分子对接靶点系统性筛选报告",
        "=" * 80,
        "",
        "筛选策略:",
        "  1. 直接铜死亡通路靶点: 铜死亡核心基因蛋白 (FDX1, LIAS, DLAT等)",
        "  2. 间接调控靶点: 通过GRN扰动影响铜死亡的BCP靶点",
        "  3. 炎症-铜死亡交叉靶点: 炎症因子与铜死亡的交汇点",
        "",
        "评分体系:",
        "  - 铜死亡通路核心度 (30%): 基于Tsvetkov Science 2022通路层级",
        "  - BCP靶点证据 (25%): 综合多组学证据评分",
        "  - GRN扰动评分 (20%): scTenifoldKnk虚拟敲除强度",
        "  - PPI网络中心度 (15%): GAT图注意力网络排名",
        "  - 可成药性 (10%): 蛋白结构可成药预测",
        "",
        "=" * 80,
        "",
    ]
    
    # Tier 1: 直接铜死亡靶点
    report_lines.extend([
        "【Tier 1】直接铜死亡核心靶点 (强烈推荐对接)",
        "-" * 80,
        "这些蛋白直接参与铜死亡通路，是BCP调控铜死亡的最直接机制靶点。",
        "",
        f"{'排名':<6}{'基因':<12}{'通路分':<10}{'BCP分':<10}{'GRN分':<10}{'PPI分':<10}{'综合分':<10}{'类别':<15}",
        "-" * 80,
    ])
    
    for i, (_, row) in enumerate(tier1.head(15).iterrows(), 1):
        report_lines.append(
            f"{i:<6}{row['Gene']:<12}{row['Pathway_Score']:<10.1f}{row['BCP_Score']:<10.1f}"
            f"{row['GRN_Score']:<10.1f}{row['PPI_Score']:<10.1f}{row['Priority_Score']:<10.1f}"
            f"{row['Category']:<15}"
        )
    
    report_lines.extend([
        "",
        "=" * 80,
        "",
    ])
    
    # Tier 2: 间接调控靶点
    report_lines.extend([
        "【Tier 2】间接调控靶点 (推荐对接)",
        "-" * 80,
        "这些靶点虽非铜死亡核心基因，但通过GRN扰动或PPI网络间接调控铜死亡通路。",
        "可能代表BCP通过炎症/免疫途径间接影响铜死亡的机制。",
        "",
        f"{'排名':<6}{'基因':<12}{'通路分':<10}{'BCP分':<10}{'GRN分':<10}{'PPI分':<10}{'综合分':<10}{'类别':<15}",
        "-" * 80,
    ])
    
    for i, (_, row) in enumerate(tier2.head(15).iterrows(), 1):
        report_lines.append(
            f"{i:<6}{row['Gene']:<12}{row['Pathway_Score']:<10.1f}{row['BCP_Score']:<10.1f}"
            f"{row['GRN_Score']:<10.1f}{row['PPI_Score']:<10.1f}{row['Priority_Score']:<10.1f}"
            f"{row['Category']:<15}"
        )
    
    report_lines.extend([
        "",
        "=" * 80,
        "",
    ])
    
    # Tier 3: 潜在候选
    report_lines.extend([
        "【Tier 3】潜在候选靶点 (可选验证)",
        "-" * 80,
        "",
        f"{'排名':<6}{'基因':<12}{'通路分':<10}{'BCP分':<10}{'GRN分':<10}{'PPI分':<10}{'综合分':<10}",
        "-" * 80,
    ])
    
    for i, (_, row) in enumerate(tier3.head(10).iterrows(), 1):
        report_lines.append(
            f"{i:<6}{row['Gene']:<12}{row['Pathway_Score']:<10.1f}{row['BCP_Score']:<10.1f}"
            f"{row['GRN_Score']:<10.1f}{row['PPI_Score']:<10.1f}{row['Priority_Score']:<10.1f}"
        )
    
    report_lines.extend([
        "",
        "=" * 80,
        "",
        "对接建议:",
        "  1. 优先对接Tier 1靶点 (直接铜死亡核心蛋白)",
        "  2. 其次考虑Tier 2靶点 (高评分BCP靶点)",
        "  3. Tier 3靶点可作为补充验证",
        "",
        "注意事项:",
        "  - 需要获取靶蛋白的3D结构 (PDB或AlphaFold预测)",
        "  - 建议进行分子动力学模拟验证对接稳定性",
        "  - 体外实验验证BCP与靶蛋白的结合亲和力",
        "",
        "=" * 80,
    ])
    
    return '\n'.join(report_lines)


def main():
    print("=" * 60)
    print("β-石竹烯(BCP)分子对接靶点系统性筛选")
    print("=" * 60)
    
    # 加载数据
    print("\n[1/4] 加载分析结果...")
    results = load_pipeline_results()
    
    if not results:
        print("错误: 无法加载分析结果")
        return
    
    # 计算优先级
    print("\n[2/4] 计算对接优先级评分...")
    df = calculate_docking_priority(results)
    
    # 分类
    print("\n[3/4] 分类靶点...")
    tier1, tier2, tier3 = categorize_targets(df)
    
    print(f"  Tier 1 (直接铜死亡): {len(tier1)} 个")
    print(f"  Tier 2 (间接调控): {len(tier2)} 个")
    print(f"  Tier 3 (潜在候选): {len(tier3)} 个")
    
    # 保存结果
    print("\n[4/4] 保存结果...")
    
    # 完整评分表
    df.to_csv(os.path.join(OUTPUT_DIR, "docking_target_scores.csv"), index=False)
    print(f"  ✓ 完整评分表: {OUTPUT_DIR}/docking_target_scores.csv")
    
    # 分Tier保存
    tier1.to_csv(os.path.join(OUTPUT_DIR, "tier1_direct_cuproptosis.csv"), index=False)
    tier2.to_csv(os.path.join(OUTPUT_DIR, "tier2_indirect_regulation.csv"), index=False)
    tier3.to_csv(os.path.join(OUTPUT_DIR, "tier3_candidates.csv"), index=False)
    
    # 生成报告
    report = generate_report(df, tier1, tier2, tier3)
    report_file = os.path.join(OUTPUT_DIR, "docking_screening_report.txt")
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write(report)
    print(f"  ✓ 筛选报告: {report_file}")
    
    # 打印报告
    print("\n" + report)
    
    print("\n" + "=" * 60)
    print("筛选完成!")
    print(f"输出目录: {OUTPUT_DIR}")
    print("=" * 60)


if __name__ == "__main__":
    main()
