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

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import (
    BASE_DIR, RESULTS_DIR, DATA_DIR,
    RAT_MOUSE_HUMAN_MAP, CUPROPTOSIS_GENES, CUPROPTOSIS_RELATED,
    BCP_TARGETS, FIG_FORMAT, FIG_DPI
)
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
    
    # ---- v9: 移除第三层(BCP靶点)，仅保留两层 ----
    logger.info("第三层 (BCP靶点): 已移除 (v9无偏设计)")
    
    # ---- 两层并集 ----
    seed_pool = layer1 | layer2
    logger.info(f"种子池 (并集): {len(seed_pool)} 个基因")
    
    # ---- 交集分析 ----
    inter_12 = layer1 & layer2
    
    logger.info(f"  第一层∩第二层: {len(inter_12)}")
    
    if inter_12:
        logger.info(f"  两层交集基因: {sorted(inter_12)[:20]}...")
    
    # ---- 保存种子池 ----
    seed_info = {
        'layer1_count': len(layer1),
        'layer2_count': len(layer2),
        'layer3_count': 0,
        'seed_pool_count': len(seed_pool),
        'inter_12': sorted(inter_12),
        'inter_13': [],
        'inter_23': [],
        'inter_123': [],
        'seed_pool': sorted(seed_pool)
    }
    
    with open(os.path.join(STAGE_DIR, "seed_pool_info.json"), 'w', encoding='utf-8') as f:
        json.dump(seed_info, f, indent=2, ensure_ascii=False)
    
    with open(os.path.join(STAGE_DIR, "seed_pool_genes.txt"), 'w', encoding='utf-8') as f:
        f.write("# 种子池基因 (人类符号)\n")
        for g in sorted(seed_pool):
            f.write(f"{g}\n")
    
    # ---- Venn图 (v9: 两层) ----
    plot_seed_venn(layer1, layer2)
    
    return seed_pool, seed_info


def fetch_string_ppi_for_cupro(cupro_genes, output_dir, min_score=400):
    """Fetch a lightweight PPI from STRING for cuproptosis genes only
    
    FIX:[v5][break Stage4→Stage5 circular dependency]
    Downloads only cuproptosis gene neighbors from STRING, creating a minimal PPI
    file for seed pool filtering BEFORE full Stage5 runs.
    """
    import urllib.request
    import urllib.parse
    
    protein_names = list(set(cupro_genes))
    logger.info(f"  正在从STRING下载铜死亡基因PPI ({len(protein_names)} 基因)...")
    
    url = "https://string-db.org/api/tsv/network"
    body = "\n".join(protein_names)
    
    params = {
        "identifiers": body,
        "species": "9606",  # Human
        "limit": 500,       # Limit neighbors to keep it lightweight
        "network_type": "full",
        "required_score": str(min_score),
    }
    
    # Use a pre-built TSV format: preferredName_A\tpreferredName_B\tscore
    # For simplicity, use the string API with POST
    import http.client
    conn = http.client.HTTPSConnection("string-db.org")
    headers = {"Content-Type": "application/x-www-form-urlencoded"}
    body_encoded = urllib.parse.urlencode({
        "identifiers": "\r\n".join(protein_names),
        "species": "9606",
        "limit": 500,
        "required_score": str(min_score),
        "caller_identity": "bcp_cupro_ciri"
    })
    
    try:
        conn.request("POST", "/api/tsv/network", body_encoded, headers)
        response = conn.getresponse()
        if response.status == 200:
            data = response.read().decode("utf-8")
            lines = data.strip().split('\n')
            if len(lines) > 1:
                # Filter to only keep columns we need
                header = lines[0].split('\t')
                col_a = header.index('preferredName_A') if 'preferredName_A' in header else None
                col_b = header.index('preferredName_B') if 'preferredName_B' in header else None
                col_s = header.index('score') if 'score' in header else None
                
                if col_a and col_b:
                    import csv
                    os.makedirs(output_dir, exist_ok=True)
                    with open(os.path.join(output_dir, "string_ppi.tsv"), 'w', newline='', encoding='utf-8') as f:
                        writer = csv.writer(f, delimiter='\t')
                        writer.writerow(['preferredName_A', 'preferredName_B', 'score'])
                        for line in lines[1:]:
                            fields = line.split('\t')
                            a = fields[col_a] if col_a < len(fields) else ''
                            b = fields[col_b] if col_b < len(fields) else ''
                            s = fields[col_s] if col_s and col_s < len(fields) else '400'
                            if a and b:
                                writer.writerow([a, b, s])
                    logger.info(f"  STRING PPI下载成功")
                    return True
            else:
                logger.warning("  STRING返回空结果")
        else:
            logger.warning(f"  STRING请求失败 (status={response.status})")
    except Exception as e:
        logger.warning(f"  STRING PPI下载失败: {e}")
    finally:
        conn.close()
    
    return False


def filter_seed_by_cupro_ppi(seed_pool, ppi_file, cupro_core_genes):
    """仅保留与铜死亡核心基因有直接PPI关系的种子池基因
    
    FIX:[v4][从源头锁定靶向性，剔除SPP1等无直接连接的泛炎症基因]
    FIX:[v5][如果PPI文件不存在，尝试从STRING实时下载铜死亡基因PPI]
    FIX:[v6][优先使用本地预置PPI缓存，消除网络依赖]
    """
    cupro_set = {g.upper() for g in cupro_core_genes}
    
    if not os.path.exists(ppi_file):
        # FIX:[v6][first try local cache, then STRING API]
        cache_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "cuproptosis_copper_ppi_cache.tsv")
        output_dir = os.path.dirname(ppi_file)
        
        if os.path.exists(cache_file):
            logger.info(f"  使用本地预置PPI缓存: {cache_file}")
            import shutil
            os.makedirs(output_dir, exist_ok=True)
            shutil.copy(cache_file, ppi_file)
            logger.info("  本地PPI缓存加载成功")
        else:
            logger.warning("PPI缓存不存在，尝试从STRING下载铜死亡基因PPI...")
            if not fetch_string_ppi_for_cupro(cupro_core_genes, output_dir):
                logger.warning("STRING下载失败，跳过PPI紧邻过滤，保留全部种子池")
                return seed_pool
            logger.info("STRING PPI下载成功")
    
    ppi_df = pd.read_csv(ppi_file, sep='\t')
    
    # Extract direct neighbors of cuproptosis core genes
    mask_a = ppi_df["preferredName_A"].str.upper().isin(cupro_set)
    mask_b = ppi_df["preferredName_B"].str.upper().isin(cupro_set)
    partners = set(ppi_df.loc[mask_a, "preferredName_B"].str.upper()) | \
               set(ppi_df.loc[mask_b, "preferredName_A"].str.upper())
    
    # 保证铜死亡基因本身也在内
    filtered = (seed_pool & partners) | cupro_set
    logger.info(f"PPI紧邻过滤: {len(seed_pool)} -> {len(filtered)} 基因（仅保留铜死亡直接邻居）")
    logger.info(f"  过滤掉的基因数: {len(seed_pool) - len(filtered)}")
    
    # 统计过滤掉的基因
    removed = seed_pool - filtered
    if removed:
        logger.info(f"  被过滤的示例基因(前20): {sorted(removed)[:20]}")
    
    return filtered


def plot_seed_venn(layer1, layer2):
    """v9: 绘制两层种子池Venn图 (移除BCP靶点层)"""
    try:
        from matplotlib_venn import venn2
        
        fig, ax = plt.subplots(figsize=(8, 8))
        v = venn2([layer1, layer2],
                  ('Bulk DEGs', 'Cuproptosis\n+ scMarkers'))
        
        if v.get_label_by_id('10'):
            v.get_label_by_id('10').set_text(str(len(layer1 - layer2)))
        if v.get_label_by_id('01'):
            v.get_label_by_id('01').set_text(str(len(layer2 - layer1)))
        if v.get_label_by_id('11'):
            v.get_label_by_id('11').set_text(str(len(layer1 & layer2)))
        
        ax.set_title('Seed Gene Pool: Two-Layer Union (v9 Unbiased)', fontsize=14, fontweight='bold')
        
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
    
    # ---- v8: 无偏管线 - 保留全种子池，不进行铜死亡过滤 ----
    seed_info['seed_pool_count'] = len(seed_pool)
    seed_info['seed_pool'] = sorted(seed_pool)
    
    with open(os.path.join(STAGE_DIR, "seed_pool_info.json"), 'w', encoding='utf-8') as f:
        json.dump(seed_info, f, indent=2, ensure_ascii=False)
    
    with open(os.path.join(STAGE_DIR, "seed_pool_genes.txt"), 'w', encoding='utf-8') as f:
        f.write("# 种子池基因 (人类符号), v8无偏管线\n")
        for g in sorted(seed_pool):
            f.write(f"{g}\n")
    
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