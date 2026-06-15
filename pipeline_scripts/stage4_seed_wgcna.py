# -*- coding: utf-8 -*-
"""
阶段4: 三层种子池构建 + WGCNA共表达网络
输入: 阶段1表达矩阵 + 阶段2单细胞marker + 阶段3人类DEGs
输出: 种子池基因列表 + WGCNA模块
"""

import os
import sys
import json
import logging
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

# 路径配置 - 使用统一导入接口
from pipeline_scripts import (
    BASE_DIR, RESULTS_DIR, DATA_DIR,
    CUPROPTOSIS_GENES, CUPROPTOSIS_RELATED,
    BCP_TARGETS, FIG_FORMAT, FIG_DPI,
)
from config import RAT_MOUSE_HUMAN_MAP
from scripts.utils import setup_logger, ensure_dir

STAGE_DIR = os.path.join(RESULTS_DIR, "stage4_seed_wgcna")
ensure_dir(STAGE_DIR)

logger = setup_logger("stage4", os.path.join(STAGE_DIR, "stage4.log"))


def load_rat_human_mapping(map_file):
    """加载大鼠-人类基因映射库（双向）"""
    logger.info(f"加载映射库: {map_file}")
    
    rat_to_human = {}
    human_to_rat = {}
    
    with open(map_file, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                rat_symbol = parts[0].strip().upper()
                human_symbol = parts[3].strip().upper()
                if rat_symbol and human_symbol:
                    rat_to_human[rat_symbol] = human_symbol
                    if human_symbol not in human_to_rat:
                        human_to_rat[human_symbol] = []
                    human_to_rat[human_symbol].append(rat_symbol)
    
    logger.info(f"  大鼠→人类: {len(rat_to_human)}, 人类→大鼠: {len(human_to_rat)}")
    return rat_to_human, human_to_rat


def build_seed_pool(human_degs, sc_markers, rat_to_human, human_to_rat):
    """三层种子池构建"""
    logger.info("=" * 50)
    logger.info("三层种子池构建")
    logger.info("=" * 50)
    
    # ---- 第一层: Bulk RNA-seq DEGs ----
    layer1 = set(human_degs['Gene'].tolist())
    logger.info(f"第一层 (Bulk DEGs): {len(layer1)} 个基因")
    
    # ---- 第二层: 铜死亡 + 单细胞marker ----
    cupro_core = set(CUPROPTOSIS_GENES)
    cupro_related = set(CUPROPTOSIS_RELATED)
    
    sc_mouse_genes = set()
    if sc_markers is not None and len(sc_markers) > 0:
        sc_mouse_genes = set(sc_markers)
    
    layer2 = cupro_core | cupro_related | sc_mouse_genes
    logger.info(f"第二层 (铜死亡+单细胞): {len(layer2)} 个基因")
    logger.info(f"  铜死亡核心: {len(cupro_core)}")
    logger.info(f"  铜死亡相关: {len(cupro_related)}")
    logger.info(f"  单细胞marker: {len(sc_mouse_genes)}")
    
    # ---- 第三层: BCP靶点 ----
    layer3 = set(BCP_TARGETS)
    logger.info(f"第三层 (BCP靶点): {len(layer3)} 个基因")
    
    # ---- 三层并集 ----
    seed_pool = layer1 | layer2 | layer3
    logger.info(f"种子池 (并集): {len(seed_pool)} 个基因")
    
    # ---- 交集分析 ----
    inter_12 = layer1 & layer2
    inter_13 = layer1 & layer3
    inter_23 = layer2 & layer3
    inter_123 = layer1 & layer2 & layer3
    
    logger.info(f"  第一层∩第二层: {len(inter_12)}")
    logger.info(f"  第一层∩第三层: {len(inter_13)}")
    logger.info(f"  第二层∩第三层: {len(inter_23)}")
    logger.info(f"  三层交集: {len(inter_123)}")
    
    if inter_123:
        logger.info(f"  三层交集基因: {sorted(inter_123)}")
    
    # ---- 保存种子池 ----
    seed_info = {
        'layer1_count': len(layer1),
        'layer2_count': len(layer2),
        'layer3_count': len(layer3),
        'seed_pool_count': len(seed_pool),
        'inter_12': sorted(inter_12),
        'inter_13': sorted(inter_13),
        'inter_23': sorted(inter_23),
        'inter_123': sorted(inter_123),
        'seed_pool': sorted(seed_pool)
    }
    
    with open(os.path.join(STAGE_DIR, "seed_pool_info.json"), 'w', encoding='utf-8') as f:
        json.dump(seed_info, f, indent=2, ensure_ascii=False)
    
    with open(os.path.join(STAGE_DIR, "seed_pool_genes.txt"), 'w', encoding='utf-8') as f:
        f.write("# 种子池基因 (人类符号)\n")
        for g in sorted(seed_pool):
            f.write(f"{g}\n")
    
    # ---- 绘制Venn图 ----
    plot_seed_venn(layer1, layer2, layer3)
    
    return seed_pool, seed_info


def plot_seed_venn(layer1, layer2, layer3):
    """绘制三层种子池Venn图"""
    try:
        from matplotlib_venn import venn3
        
        fig, ax = plt.subplots(figsize=(8, 8))
        v = venn3([layer1, layer2, layer3],
                  ('Bulk DEGs', 'Cuproptosis\n+ scMarkers', 'BCP Targets'))
        
        if v.get_label_by_id('100'):
            v.get_label_by_id('100').set_text(str(len(layer1 - layer2 - layer3)))
        if v.get_label_by_id('010'):
            v.get_label_by_id('010').set_text(str(len(layer2 - layer1 - layer3)))
        if v.get_label_by_id('001'):
            v.get_label_by_id('001').set_text(str(len(layer3 - layer1 - layer2)))
        if v.get_label_by_id('110'):
            v.get_label_by_id('110').set_text(str(len((layer1 & layer2) - layer3)))
        if v.get_label_by_id('101'):
            v.get_label_by_id('101').set_text(str(len((layer1 & layer3) - layer2)))
        if v.get_label_by_id('011'):
            v.get_label_by_id('011').set_text(str(len((layer2 & layer3) - layer1)))
        if v.get_label_by_id('111'):
            v.get_label_by_id('111').set_text(str(len(layer1 & layer2 & layer3)))
        
        ax.set_title('Seed Gene Pool: Three-Layer Union', fontsize=14, fontweight='bold')
        
        plt.tight_layout()
        fig.savefig(os.path.join(STAGE_DIR, f"seed_pool_venn.{FIG_FORMAT}"),
                    dpi=FIG_DPI, bbox_inches='tight')
        plt.close()
        logger.info("  Venn图已保存")
        
    except ImportError:
        logger.warning("  matplotlib_venn未安装，跳过Venn图")
    except Exception as e:
        logger.warning(f"  Venn图绘制失败: {e}")


def map_seed_to_rat(seed_pool, human_to_rat):
    """将种子池人类基因映射回大鼠基因"""
    logger.info("映射种子池到 rat 基因...")
    
    rat_genes = set()
    unmapped = []
    
    for hg in seed_pool:
        if hg in human_to_rat:
            rat_genes.update(human_to_rat[hg])
        else:
            unmapped.append(hg)
    
    logger.info(f"  种子池人类基因: {len(seed_pool)}")
    logger.info(f"  映射到大鼠基因: {len(rat_genes)}")
    logger.info(f"  未映射: {len(unmapped)}")
    
    if unmapped:
        with open(os.path.join(STAGE_DIR, "unmapped_seed_genes.txt"), 'w', encoding='utf-8') as f:
            for g in sorted(unmapped):
                f.write(f"{g}\n")
    
    return rat_genes


def prepare_wgcna_input(expr_file, degs_file, rat_genes):
    """准备WGCNA输入：提取种子池基因的表达子矩阵"""
    logger.info("准备WGCNA输入数据...")
    
    expr = pd.read_csv(expr_file, index_col=0)
    logger.info(f"  表达矩阵: {expr.shape[0]} 探针 x {expr.shape[1]} 样本")
    
    # 加载样本注释（如果是合并数据集）
    annot_file = os.path.join(os.path.dirname(expr_file), "sample_annotations.csv")
    sample_groups = {}
    if os.path.exists(annot_file):
        annot = pd.read_csv(annot_file)
        for _, row in annot.iterrows():
            sample_groups[row['SampleID']] = row['Group']
        logger.info(f"  加载样本注释: {len(sample_groups)} 样本")
    
    # 仅使用Sham和Model样本（排除Treatment/XST）
    if sample_groups:
        sham_model_samples = [sid for sid, grp in sample_groups.items() 
                              if grp in ['Sham', 'Model'] and sid in expr.columns]
        if len(sham_model_samples) > 0:
            expr = expr[sham_model_samples]
            logger.info(f"  筛选后: {expr.shape[0]} 探针 x {expr.shape[1]} 样本 (Sham+Model)")
    
    degs = pd.read_csv(degs_file)
    degs['GeneSymbol'] = degs['GeneSymbol'].fillna('').astype(str).str.upper()
    degs = degs[degs['GeneSymbol'] != ''].copy()
    
    seed_probes = degs[degs['GeneSymbol'].isin(rat_genes)]
    
    if len(seed_probes) == 0:
        logger.warning("  种子池基因在表达矩阵中无匹配探针!")
        return None
    
    probe_ids = seed_probes['ProbeID'].unique()
    probe_ids_in_expr = [p for p in probe_ids if p in expr.index]
    
    logger.info(f"  种子池探针: {len(probe_ids)}, 在表达矩阵中: {len(probe_ids_in_expr)}")
    
    if len(probe_ids_in_expr) < 10:
        logger.warning(f"  探针数不足 ({len(probe_ids_in_expr)}), 无法运行WGCNA")
        return None
    
    wgcna_expr = expr.loc[probe_ids_in_expr]
    
    wgcna_file = os.path.join(STAGE_DIR, "wgcna_input_expr.csv")
    wgcna_expr.to_csv(wgcna_file)
    logger.info(f"  WGCNA输入已保存: {wgcna_file} ({wgcna_expr.shape[0]} 探针 x {wgcna_expr.shape[1]} 样本)")
    
    probe_gene_map = seed_probes[['ProbeID', 'GeneSymbol']].drop_duplicates('ProbeID')
    probe_gene_map.to_csv(os.path.join(STAGE_DIR, "wgcna_probe_gene_map.csv"), index=False)
    
    return wgcna_expr


def run_wgcna_r(wgcna_expr_file, n_samples=15):
    """调用R脚本运行WGCNA"""
    logger.info("运行WGCNA (R脚本)...")
    
    r_script = os.path.join(os.path.dirname(os.path.abspath(__file__)), "stage4_wgcna.R")
    
    if not os.path.exists(r_script):
        logger.warning(f"  WGCNA R脚本不存在: {r_script}, 跳过WGCNA")
        return None
    
    import subprocess
    
    r_exe = r"C:\R\R-4.5.2\bin\Rscript.exe"
    
    cmd = [
        r_exe, r_script,
        "--expr_file", wgcna_expr_file,
        "--out_dir", STAGE_DIR,
        "--n_samples", str(n_samples)
    ]
    
    logger.info(f"  执行: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600,
                                encoding='utf-8', errors='replace')
        
        if result.returncode == 0:
            logger.info("  WGCNA完成")
            if result.stdout:
                for line in result.stdout.strip().split('\n')[-10:]:
                    logger.info(f"    R: {line}")
        else:
            logger.error(f"  WGCNA失败 (code={result.returncode})")
            if result.stderr:
                logger.error(f"    {result.stderr[:500]}")
        
        return result.returncode == 0
        
    except subprocess.TimeoutExpired:
        logger.error("  WGCNA超时")
        return False
    except Exception as e:
        logger.error(f"  WGCNA执行失败: {e}")
        return False


def main():
    logger.info("=" * 60)
    logger.info("阶段4: 三层种子池构建 + WGCNA")
    logger.info("=" * 60)
    
    # ---- 1. 加载映射库 ----
    logger.info("[1/5] 加载基因映射库...")
    rat_to_human, human_to_rat = load_rat_human_mapping(RAT_MOUSE_HUMAN_MAP)
    
    # ---- 2. 加载人类DEGs ----
    logger.info("[2/5] 加载人类DEGs...")
    human_degs_file = os.path.join(RESULTS_DIR, "stage3_enrichment", "human_degs.csv")
    human_degs = pd.read_csv(human_degs_file)
    logger.info(f"  人类DEGs: {len(human_degs)}")
    
    # ---- 3. 加载单细胞marker ----
    logger.info("[3/5] 加载单细胞marker基因...")
    sc_marker_file = os.path.join(RESULTS_DIR, "stage2_single_cell", "sc_marker_genes.txt")
    sc_markers = None
    if os.path.exists(sc_marker_file):
        with open(sc_marker_file, 'r', encoding='utf-8') as f:
            sc_markers = [line.strip() for line in f if line.strip() and not line.startswith('#')]
        logger.info(f"  单细胞marker: {len(sc_markers)}")
    else:
        logger.warning(f"  单细胞marker文件不存在: {sc_marker_file}")
    
    # ---- 4. 构建种子池 ----
    logger.info("[4/5] 构建三层种子池...")
    seed_pool, seed_info = build_seed_pool(human_degs, sc_markers, rat_to_human, human_to_rat)
    
    # ---- 5. WGCNA ----
    logger.info("[5/5] WGCNA共表达网络...")
    rat_genes = map_seed_to_rat(seed_pool, human_to_rat)
    
    # 优先使用合并数据集
    merged_expr_file = os.path.join(RESULTS_DIR, "merged_datasets", "merged_expr_matrix.csv")
    degs_file = os.path.join(RESULTS_DIR, "stage1_rma_degs", "limma_degs.csv")
    
    if os.path.exists(merged_expr_file):
        logger.info(f"  使用合并数据集: {merged_expr_file}")
        expr_file = merged_expr_file
    else:
        expr_file = os.path.join(RESULTS_DIR, "stage1_rma_degs", "expr_matrix.csv")
        logger.info(f"  使用GSE61616: {expr_file}")
    
    wgcna_expr = prepare_wgcna_input(expr_file, degs_file, rat_genes)
    
    if wgcna_expr is not None:
        n_samples = wgcna_expr.shape[1]
        logger.info(f"  WGCNA样本数: {n_samples}")
        
        if n_samples < 15:
            logger.warning(f"  样本数不足({n_samples}), WGCNA可能过拟合")
        
        wgcna_success = run_wgcna_r(
            os.path.join(STAGE_DIR, "wgcna_input_expr.csv"),
            n_samples=n_samples
        )
    else:
        wgcna_success = False
    
    # ---- 输出摘要 ----
    logger.info("\n" + "=" * 60)
    logger.info("阶段4完成! 摘要:")
    logger.info(f"  种子池大小: {len(seed_pool)}")
    logger.info(f"  三层交集: {seed_info['inter_123']}")
    logger.info(f"  WGCNA: {'成功' if wgcna_success else '跳过/失败'}")
    logger.info("=" * 60)
    
    return seed_pool, seed_info


if __name__ == "__main__":
    main()