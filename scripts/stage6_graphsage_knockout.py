# -*- coding: utf-8 -*-
"""
阶段6: scTenifoldKnk 虚拟敲除 (基于单细胞GRN)
方法: scTenifoldKnk (Osorio et al., 2022, Nature Communications)
流程:
  1. 从stage2加载单细胞count矩阵
  2. 对每个铜死亡基因调用scTenifoldKnk R包进行虚拟敲除
  3. scTenifoldKnk内部流程:
     a. PC回归构建GRN (基因调控网络)
     b. 移除目标基因 (虚拟敲除)
     c. 张量分解 (CP分解) 比较对照vs敲除
     d. 输出差异调控基因排名
  4. 综合排名分析 → 基因扰动评分
"""

import os
import sys
import subprocess
import logging
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import (
    BASE_DIR, RESULTS_DIR, DATA_DIR,
    CUPROPTOSIS_GENES, CUPROPTOSIS_RELATED,
    FIG_FORMAT, FIG_DPI
)
from scripts.utils import setup_logger, ensure_dir

STAGE_DIR = os.path.join(RESULTS_DIR, "stage6_graphsage_knockout")
ensure_dir(STAGE_DIR)

logger = setup_logger("stage6", os.path.join(STAGE_DIR, "stage6.log"))

R_EXE = r"C:\R\R-4.5.2\bin\Rscript.exe"
R_SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "stage6_scTenifoldKnk.R")


def check_stage2_outputs():
    """检查stage2输出文件是否存在"""
    logger.info("检查stage2输出...")

    # 优先检查mtx格式，回退到CSV格式
    count_file_mtx = os.path.join(RESULTS_DIR, "stage2_single_cell", "sc_count_matrix.mtx")
    count_file_csv = os.path.join(RESULTS_DIR, "stage2_single_cell", "sc_count_matrix.csv")
    annot_file = os.path.join(RESULTS_DIR, "stage2_single_cell", "sc_cell_annotations.csv")

    if os.path.exists(count_file_mtx):
        count_file = count_file_mtx
        logger.info(f"  使用mtx格式count矩阵: {count_file}")
    elif os.path.exists(count_file_csv):
        count_file = count_file_csv
        logger.info(f"  使用CSV格式count矩阵: {count_file}")
    else:
        logger.error(f"count矩阵不存在 (mtx或CSV)")
        logger.error("请先运行stage2_single_cell.py")
        return None, None

    if not os.path.exists(annot_file):
        logger.warning(f"细胞注释不存在: {annot_file}")
        annot_file = None

    logger.info(f"  细胞注释: {annot_file}")

    return count_file, annot_file


def run_sctenifoldknk_r(count_file, annot_file):
    """调用R脚本执行scTenifoldKnk虚拟敲除"""
    logger.info("=" * 50)
    logger.info("调用scTenifoldKnk R脚本...")
    logger.info("=" * 50)

    if not os.path.exists(R_SCRIPT):
        logger.error(f"R脚本不存在: {R_SCRIPT}")
        return False

    if not os.path.exists(R_EXE):
        logger.error(f"R可执行文件不存在: {R_EXE}")
        return False

    # 准备参数
    target_file = ""
    n_cores = "2"

    cmd = [
        R_EXE, R_SCRIPT,
        count_file,
        annot_file if annot_file else "",
        STAGE_DIR,
        target_file,
        n_cores
    ]

    logger.info(f"执行命令: {' '.join(cmd)}")
    logger.info("scTenifoldKnk正在运行，每个基因敲除需要5-15分钟...")
    logger.info(f"总敲除基因数: 约30-40个铜死亡相关基因")
    logger.info("预计总耗时: 2-6小时")

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=86400,
            encoding='utf-8',
            errors='replace'
        )

        if result.stdout:
            for line in result.stdout.strip().split('\n'):
                logger.info(f"  R: {line}")

        if result.returncode == 0:
            logger.info("scTenifoldKnk完成!")
            return True
        else:
            logger.error(f"scTenifoldKnk失败 (code={result.returncode})")
            if result.stderr:
                for line in result.stderr.strip().split('\n')[-20:]:
                    logger.error(f"  R ERR: {line}")
            return False

    except subprocess.TimeoutExpired:
        logger.error("scTenifoldKnk超时 (>24小时)")
        return False
    except Exception as e:
        logger.error(f"scTenifoldKnk执行失败: {e}")
        return False


def load_results():
    """加载scTenifoldKnk结果"""
    logger.info("加载scTenifoldKnk结果...")

    results = {}

    # 基因扰动评分
    score_file = os.path.join(STAGE_DIR, "gene_perturbation_scores.csv")
    if os.path.exists(score_file):
        results['scores'] = pd.read_csv(score_file)
        logger.info(f"  扰动评分: {len(results['scores'])} 个基因")
    else:
        logger.warning("  扰动评分文件不存在")

    # 所有敲除结果
    all_ko_file = os.path.join(STAGE_DIR, "all_knockout_results.csv")
    if os.path.exists(all_ko_file):
        results['all_ko'] = pd.read_csv(all_ko_file)
        logger.info(f"  所有敲除结果: {len(results['all_ko'])} 行")
    else:
        logger.warning("  所有敲除结果文件不存在")

    # 显著差异调控基因
    sig_file = os.path.join(STAGE_DIR, "significant_deg_genes.csv")
    if os.path.exists(sig_file):
        results['sig'] = pd.read_csv(sig_file)
        logger.info(f"  显著差异调控基因: {len(results['sig'])}")
    else:
        logger.warning("  显著差异调控基因文件不存在")

    # 敲除影响矩阵
    impact_file = os.path.join(STAGE_DIR, "knockout_impact_matrix.csv")
    if os.path.exists(impact_file):
        results['impact'] = pd.read_csv(impact_file, index_col=0)
        logger.info(f"  敲除影响矩阵: {results['impact'].shape}")
    else:
        logger.warning("  敲除影响矩阵文件不存在")

    return results


def plot_perturbation_scores(results):
    """绘制基因扰动评分图"""
    if 'scores' not in results:
        return

    logger.info("绘制基因扰动评分图...")

    scores = results['scores'].head(30)

    fig, ax = plt.subplots(figsize=(12, 10))

    colors = []
    cupro_core = set(CUPROPTOSIS_GENES)
    cupro_related = set(CUPROPTOSIS_RELATED)

    for g in scores['gene']:
        g_upper = g.upper()
        if g_upper in cupro_core:
            colors.append('#e74c3c')
        elif g_upper in cupro_related:
            colors.append('#e67e22')
        else:
            colors.append('#3498db')

    bars = ax.barh(range(len(scores)), scores['perturbation_score'], color=colors)
    ax.set_yticks(range(len(scores)))
    ax.set_yticklabels(scores['gene'])
    ax.invert_yaxis()
    ax.set_xlabel('Perturbation Score', fontsize=12)
    ax.set_title('scTenifoldKnk: Gene Perturbation Scores\n(GRN-based Virtual Knockout)',
                 fontsize=14, fontweight='bold')

    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#e74c3c', label='Cuproptosis Core'),
        Patch(facecolor='#e67e22', label='Cuproptosis Related'),
        Patch(facecolor='#3498db', label='Other Genes'),
    ]
    ax.legend(handles=legend_elements, loc='lower right')

    plt.tight_layout()
    fig.savefig(os.path.join(STAGE_DIR, f"perturbation_scores.{FIG_FORMAT}"),
                dpi=FIG_DPI, bbox_inches='tight')
    plt.close()
    logger.info("  扰动评分图已保存")


def plot_knockout_impact_heatmap(results):
    """绘制敲除影响热图"""
    if 'impact' not in results:
        return

    logger.info("绘制敲除影响热图...")

    impact = results['impact']

    # 限制大小：取影响最大的基因和敲除
    gene_impact = impact.sum(axis=1).sort_values(ascending=False)
    ko_impact = impact.sum(axis=0).sort_values(ascending=False)

    top_genes = gene_impact.head(20).index
    top_kos = ko_impact.head(20).index

    sub_impact = impact.loc[top_genes, top_kos]

    fig, ax = plt.subplots(figsize=(14, 10))
    sns.heatmap(sub_impact, cmap='Reds', ax=ax, cbar_kws={'label': 'Number of KO effects'})
    ax.set_xlabel('Knocked Out Gene', fontsize=12)
    ax.set_ylabel('Affected Gene', fontsize=12)
    ax.set_title('scTenifoldKnk: Knockout Impact Matrix\n(GRN-based)', fontsize=14, fontweight='bold')
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)

    plt.tight_layout()
    fig.savefig(os.path.join(STAGE_DIR, f"knockout_impact_heatmap.{FIG_FORMAT}"),
                dpi=FIG_DPI, bbox_inches='tight')
    plt.close()
    logger.info("  敲除影响热图已保存")


def plot_volcano(results):
    """绘制差异调控火山图"""
    if 'all_ko' not in results:
        return

    logger.info("绘制差异调控火山图...")

    df = results['all_ko']
    df['logFC'] = np.log2(df['FC'].clip(lower=1e-10))
    df['neg_log10_p'] = -np.log10(df['p.value'].clip(lower=1e-300))

    fig, ax = plt.subplots(figsize=(10, 8))

    colors = []
    for _, row in df.iterrows():
        if row['p.value'] < 0.05 and abs(row['logFC']) > 0.5:
            colors.append('#e74c3c')
        elif row['p.value'] < 0.05:
            colors.append('#3498db')
        else:
            colors.append('#95a5a6')

    ax.scatter(df['logFC'], df['neg_log10_p'], c=colors, alpha=0.5, s=5)
    ax.axhline(-np.log10(0.05), color='red', linestyle='--', alpha=0.5, label='p=0.05')
    ax.axvline(0.5, color='gray', linestyle='--', alpha=0.3)
    ax.axvline(-0.5, color='gray', linestyle='--', alpha=0.3)

    ax.set_xlabel('log2(Fold Change)', fontsize=12)
    ax.set_ylabel('-log10(p-value)', fontsize=12)
    ax.set_title('scTenifoldKnk: Differential Regulation Volcano Plot\n(All Knockouts Combined)',
                 fontsize=14, fontweight='bold')
    ax.legend()

    plt.tight_layout()
    fig.savefig(os.path.join(STAGE_DIR, f"deregulation_volcano.{FIG_FORMAT}"),
                dpi=FIG_DPI, bbox_inches='tight')
    plt.close()
    logger.info("  火山图已保存")


def save_summary(results):
    """保存结果摘要"""
    logger.info("保存结果摘要...")

    summary_lines = [
        "# ============================================================",
        "# scTenifoldKnk 虚拟敲除结果摘要",
        "# 方法: scTenifoldKnk (Osorio et al., 2022, Nature Comms)",
        "# 原理: 单细胞GRN推断 + 虚拟敲除 + 张量分解",
        "# ============================================================",
        "",
    ]

    if 'scores' in results:
        scores = results['scores']
        summary_lines.append(f"## 基因扰动评分 (Top 30)")
        summary_lines.append("")
        summary_lines.append(f"{'Rank':<6}{'Gene':<20}{'Score':<12}{'n_KO':<8}{'|logFC|':<10}")
        summary_lines.append("-" * 56)
        for i, (_, row) in enumerate(scores.head(30).iterrows()):
            summary_lines.append(
                f"{i+1:<6}{row['gene']:<20}{row['perturbation_score']:<12.2f}"
                f"{int(row['n_knockouts']):<8}{abs(row['logFC']):<10.3f}"
            )

    if 'sig' in results:
        summary_lines.append("")
        summary_lines.append(f"## 显著差异调控基因: {len(results['sig'])}")
        summary_lines.append("")

    if 'impact' in results:
        summary_lines.append("")
        summary_lines.append(f"## 敲除影响矩阵: {results['impact'].shape[0]} 基因 x {results['impact'].shape[1]} 敲除")
        summary_lines.append("")

    summary_lines.extend([
        "",
        "## 铜死亡核心基因敲除效应",
        "",
    ])

    if 'scores' in results:
        for g in CUPROPTOSIS_GENES:
            g_lower = g.lower()
            match = results['scores'][results['scores']['gene'].str.lower() == g_lower]
            if len(match) > 0:
                row = match.iloc[0]
                summary_lines.append(
                    f"  {g}: score={row['perturbation_score']:.2f}, "
                    f"n_KO={int(row['n_knockouts'])}, |logFC|={abs(row['logFC']):.3f}"
                )
            else:
                summary_lines.append(f"  {g}: 未在结果中")

    with open(os.path.join(STAGE_DIR, "knockout_summary.txt"), 'w', encoding='utf-8') as f:
        f.write('\n'.join(summary_lines))

    logger.info("  摘要已保存")


def main():
    logger.info("=" * 60)
    logger.info("阶段6: scTenifoldKnk 虚拟敲除 (基于单细胞GRN)")
    logger.info("=" * 60)

    # ---- 1. 检查stage2输出 ----
    logger.info("[1/4] 检查stage2输出...")
    count_file, annot_file = check_stage2_outputs()

    if count_file is None:
        logger.error("stage2输出不完整，终止")
        return None

    # ---- 2. 运行scTenifoldKnk ----
    logger.info("[2/4] 运行scTenifoldKnk虚拟敲除...")
    logger.info("  这将调用R脚本，对每个铜死亡基因进行GRN构建+虚拟敲除")
    logger.info("  每个基因约5-15分钟，总耗时取决于基因数量")

    success = run_sctenifoldknk_r(count_file, annot_file)

    if not success:
        logger.warning("scTenifoldKnk未完全成功，尝试加载部分结果...")

    # ---- 3. 加载结果 ----
    logger.info("[3/4] 加载并可视化结果...")
    results = load_results()

    if not results:
        logger.error("无可用结果，终止")
        return None

    # 可视化
    plot_perturbation_scores(results)
    plot_knockout_impact_heatmap(results)
    plot_volcano(results)
    save_summary(results)

    # ---- 4. 输出摘要 ----
    logger.info("\n" + "=" * 60)
    logger.info("阶段6完成! 摘要:")
    logger.info(f"  方法: scTenifoldKnk (单细胞GRN虚拟敲除)")

    if 'scores' in results:
        logger.info(f"  扰动评分基因: {len(results['scores'])}")
        logger.info(f"  Top 5:")
        for i, (_, row) in enumerate(results['scores'].head(5).iterrows()):
            logger.info(f"    {i+1}. {row['gene']}: score={row['perturbation_score']:.2f}")

    if 'sig' in results:
        logger.info(f"  显著差异调控基因: {len(results['sig'])}")

    if 'impact' in results:
        logger.info(f"  敲除影响矩阵: {results['impact'].shape}")

    logger.info("=" * 60)

    return results


if __name__ == "__main__":
    main()