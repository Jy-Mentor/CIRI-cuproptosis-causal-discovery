# -*- coding: utf-8 -*-
import os
from collections import defaultdict

SEED_FILE = r"C:\Users\Jy-Mentor-7\Desktop\种子基因.txt"
PPI_FILE = r"C:\Users\Jy-Mentor-7\Desktop\9606蛋白质\9606_human_ppi_symbol.txt"

OUTPUT_DIR = os.path.dirname(SEED_FILE)
OUTPUT_PPI = os.path.join(OUTPUT_DIR, "ppi_subgraph.csv")
OUTPUT_GENES = os.path.join(OUTPUT_DIR, "subgraph_genes.txt")

MIN_SCORE = 700
MAX_DISTANCE = 2
MIN_DEGREE = 0

print("=" * 60)
print("PPI 子图提取 (≤2步, score≥700)")
print("=" * 60)

print(f"\n[1/5] 读取种子基因: {SEED_FILE}")
with open(SEED_FILE, "r", encoding="utf-8") as f:
    raw_seeds = [line.strip().upper() for line in f if line.strip()]
raw_seeds = list(dict.fromkeys(raw_seeds))
print(f"  -> 原始种子基因数: {len(raw_seeds)}")

print(f"\n[2/5] 构建 PPI 邻接表: {PPI_FILE}")
adjacency = defaultdict(set)
edge_count_total = 0
self_loop_count = 0

with open(PPI_FILE, "r", encoding="utf-8") as f:
    header = f.readline().strip()
    for line in f:
        line = line.strip()
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        gene_a = parts[0].strip().upper()
        gene_b = parts[1].strip().upper()
        score_str = parts[2].strip()
        try:
            score = float(score_str)
            if score < MIN_SCORE:
                continue
        except ValueError:
            continue
        if gene_a == gene_b:
            self_loop_count += 1
            continue
        adjacency[gene_a].add(gene_b)
        adjacency[gene_b].add(gene_a)
        edge_count_total += 1

print(f"  -> 总边数 (combined_score ≥ {MIN_SCORE}, 去自环): {edge_count_total}")
print(f"  -> 排除的自环数: {self_loop_count}")
print(f"  -> 唯一基因数: {len(adjacency)}")

print(f"\n[3/5] BFS 距离 ≤ {MAX_DISTANCE} 步扩展")

seeds_in_network = set()
for s in raw_seeds:
    if s in adjacency:
        seeds_in_network.add(s)

missing_seeds = [s for s in raw_seeds if s not in adjacency]
print(f"  -> 网络中存在的种子基因: {len(seeds_in_network)}")
print(f"  -> 未在网络中找到的种子基因: {len(missing_seeds)}")
if missing_seeds:
    print(f"     缺失基因 (前20个): {missing_seeds[:20]}")

visited = set(seeds_in_network)
frontier = set(seeds_in_network)

for step in range(1, MAX_DISTANCE + 1):
    next_frontier = set()
    for gene in frontier:
        for neighbor in adjacency[gene]:
            if neighbor not in visited:
                visited.add(neighbor)
                next_frontier.add(neighbor)
    print(f"  -> 距离 {step}: 新增 {len(next_frontier)} 个基因")
    frontier = next_frontier

subgraph_nodes = visited
print(f"  -> 子图总基因数 (度过滤前): {len(subgraph_nodes)}")

print(f"\n[4/6] 度过滤 (仅保留度 ≥ {MIN_DEGREE} 的基因)")

subgraph_degree = {}
for gene in subgraph_nodes:
    subgraph_degree[gene] = len(adjacency[gene] & subgraph_nodes)

pruned_count = 0
while True:
    low_degree = [g for g in subgraph_nodes if subgraph_degree[g] < MIN_DEGREE]
    if not low_degree:
        break
    for g in low_degree:
        subgraph_nodes.discard(g)
        for neighbor in (adjacency[g] & subgraph_nodes):
            subgraph_degree[neighbor] -= 1
        del subgraph_degree[g]
    pruned_count += len(low_degree)

print(f"  -> 移除的低度基因数: {pruned_count}")
print(f"  -> 子图总基因数 (度过滤后): {len(subgraph_nodes)}")

print(f"\n[5/6] 过滤边 (两端均在子图内)")

subgraph_edges = []
with open(PPI_FILE, "r", encoding="utf-8") as f:
    f.readline()
    for line in f:
        line = line.strip()
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        gene_a = parts[0].strip().upper()
        gene_b = parts[1].strip().upper()
        score_str = parts[2].strip()
        try:
            score = float(score_str)
            if score < MIN_SCORE:
                continue
        except ValueError:
            continue
        if gene_a == gene_b:
            continue
        if gene_a in subgraph_nodes and gene_b in subgraph_nodes:
            subgraph_edges.append((gene_a, gene_b, score_str))

print(f"  -> 子图边数: {len(subgraph_edges)}")

print(f"\n[6/6] 写入输出文件")

with open(OUTPUT_PPI, "w", encoding="utf-8") as f:
    f.write("gene_a\tgene_b\tcombined_score\n")
    for gene_a, gene_b, score in subgraph_edges:
        f.write(f"{gene_a}\t{gene_b}\t{score}\n")
print(f"  -> PPI 子图: {OUTPUT_PPI}")

sorted_genes = sorted(subgraph_nodes)
with open(OUTPUT_GENES, "w", encoding="utf-8") as f:
    f.write("gene\n")
    for gene in sorted_genes:
        f.write(f"{gene}\n")
print(f"  -> 基因列表: {OUTPUT_GENES}")

print("\n" + "=" * 60)
print("摘要")
print("=" * 60)
print(f"  原始 PPI 边数:      {edge_count_total}")
print(f"  原始 PPI 基因数:     {len(adjacency)}")
print(f"  种子基因数:          {len(raw_seeds)}")
print(f"  网络中种子基因数:    {len(seeds_in_network)}")
print(f"  子图基因数 (≤2步):   {len(subgraph_nodes)}")
print(f"  子图边数 (≤2步):     {len(subgraph_edges)}")
print("=" * 60)
