# -*- coding: utf-8 -*-
"""
阶段5: STRING PPI网络 + MCODE模块识别
输入: 种子池基因 (人类符号)
输出: PPI网络 + MCODE模块 + 模块基因
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
import networkx as nx
import requests
from io import StringIO

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import (
    BASE_DIR, RESULTS_DIR, DATA_DIR,
    FIG_FORMAT, FIG_DPI
)
from scripts.utils import setup_logger, ensure_dir

STAGE_DIR = os.path.join(RESULTS_DIR, "stage5_ppi_mcode")
ensure_dir(STAGE_DIR)

logger = setup_logger("stage5", os.path.join(STAGE_DIR, "stage5.log"))

STRING_API_URL = "https://string-db.org/api"


def load_seed_pool():
    """加载种子池基因"""
    seed_file = os.path.join(RESULTS_DIR, "stage4_seed_wgcna", "seed_pool_genes.txt")
    
    if not os.path.exists(seed_file):
        logger.error(f"种子池文件不存在: {seed_file}")
        return []
    
    genes = []
    with open(seed_file, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                genes.append(line)
    
    logger.info(f"加载种子池: {len(genes)} 个基因")
    return genes


def query_string_network(genes, species=9606, score_threshold=400):
    """通过STRING API获取PPI网络"""
    logger.info(f"查询STRING PPI网络 (物种={species}, 阈值={score_threshold})...")
    
    # STRING API对URL长度有限制，分批查询
    batch_size = 100
    all_edges = []
    
    for i in range(0, len(genes), batch_size):
        batch = genes[i:i+batch_size]
        
        params = {
            'identifiers': '\r'.join(batch),
            'species': species,
            'required_score': score_threshold,
            'limit': 1,
            'caller_identity': 'BCP_Cuproptosis_CIRI'
        }
        
        try:
            resp = requests.post(f"{STRING_API_URL}/tsv/network", data=params, timeout=120)
            
            if resp.status_code != 200:
                logger.warning(f"  批次 {i//batch_size+1} 返回 {resp.status_code}")
                continue
            
            batch_df = pd.read_csv(StringIO(resp.text), sep='\t')
            all_edges.append(batch_df)
            
            logger.info(f"  批次 {i//batch_size+1}/{(len(genes)+batch_size-1)//batch_size}: {len(batch_df)} 边")
            
        except requests.exceptions.Timeout:
            logger.warning(f"  批次 {i//batch_size+1} 超时")
            continue
        except Exception as e:
            logger.warning(f"  批次 {i//batch_size+1} 错误: {e}")
            continue
    
    if len(all_edges) == 0:
        logger.warning("  所有批次均失败")
        return None
    
    try:
        df = pd.concat(all_edges, ignore_index=True)
    except ValueError as e:
        logger.error(f"  合并边列表失败: {e}")
        return None
    
    # 去重
    df = df.drop_duplicates(subset=['preferredName_A', 'preferredName_B'])
    
    logger.info(f"  总PPI边数: {len(df)}")
    
    if len(df) == 0:
        logger.warning("  PPI网络为空!")
        return None
    
    # 统计
    nodes = set(df['preferredName_A']) | set(df['preferredName_B'])
    logger.info(f"  PPI节点数: {len(nodes)}")
    
    return df


def build_networkx_graph(ppi_df):
    """构建NetworkX图"""
    logger.info("构建NetworkX图...")
    
    G = nx.Graph()
    
    for _, row in ppi_df.iterrows():
        node_a = row['preferredName_A']
        node_b = row['preferredName_B']
        score = float(row.get('score', 0)) / 1000.0
        
        G.add_edge(node_a, node_b, weight=score)
    
    logger.info(f"  图: {G.number_of_nodes()} 节点, {G.number_of_edges()} 边")
    
    return G


def run_mcode(G):
    """运行MCODE模块识别"""
    logger.info("运行MCODE模块识别...")
    
    try:
        import community as community_louvain
        
        # 使用Louvain社区检测作为MCODE替代
        partition = community_louvain.best_partition(G, weight='weight')
        
        modules = {}
        for node, mod_id in partition.items():
            if mod_id not in modules:
                modules[mod_id] = []
            modules[mod_id].append(node)
        
        # 按模块大小排序
        sorted_modules = sorted(modules.items(), key=lambda x: len(x[1]), reverse=True)
        
        logger.info(f"  识别 {len(sorted_modules)} 个模块:")
        for mod_id, mod_genes in sorted_modules[:10]:
            logger.info(f"    模块 {mod_id}: {len(mod_genes)} 基因")
        
        return sorted_modules
        
    except ImportError:
        logger.warning("  python-louvain未安装，使用连通分量")
        
        components = list(nx.connected_components(G))
        sorted_components = sorted(components, key=len, reverse=True)
        
        modules = [(i, list(comp)) for i, comp in enumerate(sorted_components)]
        
        logger.info(f"  识别 {len(modules)} 个连通分量:")
        for mod_id, mod_genes in modules[:10]:
            logger.info(f"    分量 {mod_id}: {len(mod_genes)} 基因")
        
        return modules


def identify_hub_genes(G, top_n=30):
    """识别PPI网络hub基因 (Degree + Betweenness)"""
    logger.info(f"识别hub基因 (Top {top_n})...")
    
    # Degree
    degrees = dict(G.degree())
    
    # Betweenness
    try:
        betweenness = nx.betweenness_centrality(G, weight='weight')
    except Exception:
        betweenness = {n: 0 for n in G.nodes()}
    
    # 综合评分
    deg_max = max(degrees.values()) if degrees else 1
    bet_max = max(betweenness.values()) if betweenness else 1
    
    scores = {}
    for node in G.nodes():
        deg_norm = degrees.get(node, 0) / deg_max if deg_max > 0 else 0
        bet_norm = betweenness.get(node, 0) / bet_max if bet_max > 0 else 0
        scores[node] = 0.5 * deg_norm + 0.5 * bet_norm
    
    sorted_hubs = sorted(scores.items(), key=lambda x: x[1], reverse=True)
    top_hubs = sorted_hubs[:top_n]
    
    for i, (gene, score) in enumerate(top_hubs[:10]):
        logger.info(f"  {i+1}. {gene}: score={score:.4f}, deg={degrees.get(gene, 0)}")
    
    return top_hubs


def plot_ppi_network(G, modules, hub_genes, top_n_modules=5):
    """绘制PPI网络图"""
    logger.info("绘制PPI网络图...")
    
    try:
        # 取前N个模块
        top_modules = modules[:top_n_modules]
        
        # 模块颜色映射
        module_colors = plt.cm.tab10(np.linspace(0, 1, len(top_modules)))
        node_color_map = {}
        for i, (mod_id, mod_genes) in enumerate(top_modules):
            for g in mod_genes:
                node_color_map[g] = module_colors[i]
        
        # 子图
        module_nodes = set()
        for _, mod_genes in top_modules:
            module_nodes.update(mod_genes)
        
        subG = G.subgraph(module_nodes)
        
        if subG.number_of_nodes() == 0:
            logger.warning("  子图为空，跳过绘图")
            return
        
        fig, ax = plt.subplots(figsize=(14, 12))
        
        pos = nx.spring_layout(subG, k=2, iterations=50, seed=42)
        
        node_colors = [node_color_map.get(n, 'gray') for n in subG.nodes()]
        node_sizes = [300 + 50 * subG.degree(n) for n in subG.nodes()]
        
        nx.draw_networkx_edges(subG, pos, alpha=0.2, edge_color='gray', ax=ax)
        nx.draw_networkx_nodes(subG, pos, node_color=node_colors, node_size=node_sizes,
                               alpha=0.8, ax=ax)
        
        # 标注hub基因
        hub_set = set(h[0] for h in hub_genes[:15])
        hub_labels = {n: n for n in subG.nodes() if n in hub_set}
        nx.draw_networkx_labels(subG, pos, labels=hub_labels, font_size=6,
                                font_weight='bold', ax=ax)
        
        ax.set_title(f'PPI Network: {subG.number_of_nodes()} nodes, {subG.number_of_edges()} edges',
                     fontsize=14, fontweight='bold')
        ax.axis('off')
        
        plt.tight_layout()
        fig.savefig(os.path.join(STAGE_DIR, f"ppi_network.{FIG_FORMAT}"),
                    dpi=FIG_DPI, bbox_inches='tight')
        plt.close()
        logger.info("  PPI网络图已保存")
        
    except Exception as e:
        logger.warning(f"  PPI绘图失败: {e}")


def main():
    logger.info("=" * 60)
    logger.info("阶段5: STRING PPI + MCODE模块识别")
    logger.info("=" * 60)
    
    # ---- 1. 加载种子池 ----
    logger.info("[1/5] 加载种子池...")
    seed_genes = load_seed_pool()
    
    if len(seed_genes) == 0:
        logger.error("种子池为空，终止")
        return None
    
    # ---- 2. STRING PPI ----
    logger.info("[2/5] 查询STRING PPI网络...")
    
    ppi_df = query_string_network(seed_genes, species=9606, score_threshold=400)
    
    if ppi_df is None or len(ppi_df) == 0:
        logger.warning("PPI网络为空，使用本地备份方案")
        # 保存空结果
        with open(os.path.join(STAGE_DIR, "ppi_status.txt"), 'w') as f:
            f.write("PPI_FAILED")
        return None
    
    ppi_df.to_csv(os.path.join(STAGE_DIR, "string_ppi.tsv"), sep='\t', index=False)
    
    # ---- 3. 构建NetworkX图 ----
    logger.info("[3/5] 构建NetworkX图...")
    G = build_networkx_graph(ppi_df)
    
    # ---- 4. MCODE模块识别 ----
    logger.info("[4/5] MCODE模块识别...")
    modules = run_mcode(G)
    
    # 保存模块
    module_data = {}
    for mod_id, mod_genes in modules:
        module_data[str(mod_id)] = {
            'size': len(mod_genes),
            'genes': mod_genes
        }
    
    with open(os.path.join(STAGE_DIR, "mcode_modules.json"), 'w', encoding='utf-8') as f:
        json.dump(module_data, f, indent=2, ensure_ascii=False)
    
    # ---- 5. Hub基因识别 ----
    logger.info("[5/5] Hub基因识别...")
    hub_genes = identify_hub_genes(G, top_n=30)
    
    with open(os.path.join(STAGE_DIR, "ppi_hub_genes.txt"), 'w', encoding='utf-8') as f:
        f.write("# PPI Hub基因 (Degree+Betweenness综合评分)\n")
        for gene, score in hub_genes:
            f.write(f"{gene}\t{score:.4f}\n")
    
    # ---- 绘图 ----
    plot_ppi_network(G, modules, hub_genes)
    
    # ---- 输出摘要 ----
    logger.info("\n" + "=" * 60)
    logger.info("阶段5完成! 摘要:")
    logger.info(f"  PPI节点: {G.number_of_nodes()}")
    logger.info(f"  PPI边: {G.number_of_edges()}")
    logger.info(f"  MCODE模块: {len(modules)}")
    logger.info(f"  Hub基因: {len(hub_genes)}")
    logger.info("=" * 60)
    
    return G, modules, hub_genes


if __name__ == "__main__":
    main()