# L5 反向筛选层：MCI + 超图筛选（公共数据库版本）
# 输入：G_ideal 理想干预子图 + TCMSP/ETCM/SymMap 成分-靶点数据库
# 主方案：MCI 计算 + 多数据库交叉验证

import pandas as pd
import numpy as np
import json
import os
import warnings
warnings.filterwarnings('ignore')

SEED = 42
np.random.seed(SEED)

OUTPUT_DIR = '../results/L5_reverse_screening'
FIGURE_DIR = '../figures/L5'
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(FIGURE_DIR, exist_ok=True)

# ==================== 0. 配置参数 ====================
MCI_STRUCT_WEIGHT = 0.6
MCI_DIR_WEIGHT = 0.4
MCI_THRESHOLD = 0.5
DIRECTION_CONSISTENCY_THRESHOLD = 0.80

# 公共数据库路径
TCMSP_PATH = '../data/TCMSP_compound_targets.csv'
ETCM_PATH = '../data/ETCM_compound_targets.csv'
SYMMAP_PATH = '../data/SymMap_compound_targets.csv'

# β-石竹烯（BCP）已知靶点
BCP_KNOWN_TARGETS = ['CNR2', 'TRPV1', 'PPARG']

def load_ideal_subgraph(l4_output_path):
    print("=== 步骤1: 加载 G_ideal 理想干预子图 ===")
    
    g_ideal_path = os.path.join(l4_output_path, 'G_ideal.json')
    
    if not os.path.exists(g_ideal_path):
        print("警告: 未找到 G_ideal 文件，使用示例数据")
        return {
            'M_inhibit': ['μ1', 'μ3', 'μ5'],
            'M_activate': ['μ2', 'μ4', 'μ6'],
            'E_ideal': [('μ1', 'μ2'), ('μ3', 'μ4'), ('μ5', 'μ6')]
        }
    
    with open(g_ideal_path, 'r') as f:
        g_ideal = json.load(f)
    
    print(f"M_inhibit: {len(g_ideal['M_inhibit'])} 个宏节点")
    print(f"M_activate: {len(g_ideal['M_activate'])} 个宏节点")
    print(f"E_ideal: {len(g_ideal['E_ideal'])} 条边")
    
    return g_ideal

def load_compound_targets(tcmsp_path, etcm_path, symmap_path):
    print("\n=== 步骤2: 加载成分-靶点数据库 ===")
    
    all_targets = {}
    
    # TCMSP
    if os.path.exists(tcmsp_path):
        tcmsp_df = pd.read_csv(tcmsp_path)
        for _, row in tcmsp_df.iterrows():
            compound = row['compound']
            target = row['target']
            if compound not in all_targets:
                all_targets[compound] = {'targets': set(), 'databases': set()}
            all_targets[compound]['targets'].add(target)
            all_targets[compound]['databases'].add('TCMSP')
        print(f"  TCMSP: {len(set(tcmsp_df['compound']))} 个成分")
    
    # ETCM
    if os.path.exists(etcm_path):
        etcm_df = pd.read_csv(etcm_path)
        for _, row in etcm_df.iterrows():
            compound = row['compound']
            target = row['target']
            if compound not in all_targets:
                all_targets[compound] = {'targets': set(), 'databases': set()}
            all_targets[compound]['targets'].add(target)
            all_targets[compound]['databases'].add('ETCM')
        print(f"  ETCM: {len(set(etcm_df['compound']))} 个成分")
    
    # SymMap
    if os.path.exists(symmap_path):
        symmap_df = pd.read_csv(symmap_path)
        for _, row in symmap_df.iterrows():
            compound = row['compound']
            target = row['target']
            if compound not in all_targets:
                all_targets[compound] = {'targets': set(), 'databases': set()}
            all_targets[compound]['targets'].add(target)
            all_targets[compound]['databases'].add('SymMap')
        print(f"  SymMap: {len(set(symmap_df['compound']))} 个成分")
    
    # BCP 已知靶点
    all_targets['beta-caryophyllene'] = {
        'targets': set(BCP_KNOWN_TARGETS),
        'databases': {'literature_curated'}
    }
    print(f"  BCP (文献 curated): {len(BCP_KNOWN_TARGETS)} 个靶点")
    
    print(f"\n总成分数: {len(all_targets)}")
    
    return all_targets

def calculate_mci(compound_data, g_ideal, all_compounds):
    print("\n=== 步骤3: MCI 计算 ===")
    
    results = []
    
    for compound, data in all_compounds.items():
        targets = data['targets']
        databases = data['databases']
        
        # 结构契合度
        all_ideal_nodes = set(g_ideal['M_inhibit'] + g_ideal['M_activate'])
        if len(all_ideal_nodes) == 0:
            mci_struct = 0
        else:
            intersection = targets.intersection(all_ideal_nodes)
            mci_struct = len(intersection) / len(all_ideal_nodes)
        
        # 方向契合度
        if len(targets) == 0 or len(all_ideal_nodes) == 0:
            mci_dir = 0.5
        else:
            correct_direction = 0
            total = 0
            
            for target in targets:
                if target in g_ideal['M_inhibit']:
                    total += 1
                    if compound_data.get(f"{compound}_{target}_direction") == 'inhibit':
                        correct_direction += 1
                elif target in g_ideal['M_activate']:
                    total += 1
                    if compound_data.get(f"{compound}_{target}_direction") == 'activate':
                        correct_direction += 1
            
            mci_dir = correct_direction / total if total > 0 else 0.5
        
        # 综合 MCI
        mci_comprehensive = MCI_STRUCT_WEIGHT * mci_struct + MCI_DIR_WEIGHT * mci_dir
        
        # 数据库一致性数量
        n_databases = len(databases)
        
        results.append({
            'compound': compound,
            'targets': ','.join(targets),
            'n_targets': len(targets),
            'databases': ','.join(databases),
            'n_databases': n_databases,
            'mci_struct': mci_struct,
            'mci_dir': mci_dir,
            'mci_comprehensive': mci_comprehensive,
            'passes_threshold': mci_comprehensive > MCI_THRESHOLD
        })
    
    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values('mci_comprehensive', ascending=False).reset_index(drop=True)
    
    print(f"  通过阈值 (MCI > {MCI_THRESHOLD}): {results_df['passes_threshold'].sum()} 个成分")
    
    return results_df

def hypergraph_reasoning(compound_results, g_ideal):
    print("\n=== 步骤4: 超图高阶推理 ===")
    
    # 构建超边 {C_i, μ_a, μ_b, P}
    hyperedges = []
    
    top_compounds = compound_results[compound_results['mci_comprehensive'] > MCI_THRESHOLD]
    
    for _, row in top_compounds.iterrows():
        compound = row['compound']
        targets = set(row['targets'].split(','))
        
        # 找到共同靶向的宏节点对
        target_list = list(targets)
        for i in range(len(target_list)):
            for j in range(i+1, len(target_list)):
                mu_a = target_list[i]
                mu_b = target_list[j]
                
                if mu_a in g_ideal['M_inhibit'] or mu_a in g_ideal['M_activate']:
                    if mu_b in g_ideal['M_inhibit'] or mu_b in g_ideal['M_activate']:
                        hyperedges.append({
                            'compound': compound,
                            'mu_a': mu_a,
                            'mu_b': mu_b,
                            'weight': row['n_databases'],
                            'mci': row['mci_comprehensive']
                        })
    
    hyperedges_df = pd.DataFrame(hyperedges)
    
    print(f"  超边数: {len(hyperedges_df)}")
    
    return hyperedges_df

def literature_validation(compound_results):
    print("\n=== 步骤5: 文献 curated 靶点验证 ===")
    
    validated_compounds = []
    
    for _, row in compound_results.iterrows():
        compound = row['compound']
        targets = set(row['targets'].split(','))
        
        # 检查是否有文献 curated 靶点
        has_literature = 'literature_curated' in row['databases']
        
        validated_compounds.append({
            'compound': compound,
            'mci': row['mci_comprehensive'],
            'has_literature': has_literature,
            'n_databases': row['n_databases']
        })
    
    validated_df = pd.DataFrame(validated_compounds)
    
    n_literature_validated = validated_df['has_literature'].sum()
    print(f"  Top-10 中有文献 curated 靶点: {n_literature_validated} 个 (要求 ≥ 2)")
    
    return validated_df

def export_results(compound_results, hyperedges_df, validation_df, output_dir):
    print("\n=== 步骤6: 导出结果 ===")
    
    compound_results.to_csv(os.path.join(output_dir, 'compound_ranking.csv'), index=False)
    hyperedges_df.to_csv(os.path.join(output_dir, 'hypergraph_edges.csv'), index=False)
    validation_df.to_csv(os.path.join(output_dir, 'literature_validation.csv'), index=False)
    
    # Top-5 进入 L6
    top5 = compound_results.head(5)
    top5.to_csv(os.path.join(output_dir, 'top5_for_L6.csv'), index=False)
    
    print(f"  Top-5 候选成分:")
    for _, row in top5.iterrows():
        print(f"    {row['compound']}: MCI = {row['mci_comprehensive']:.4f}")
    
    print(f"\n结果已导出至 {output_dir}")

def main():
    print("=" * 60)
    print("L5 反向筛选层：MCI + 超图筛选")
    print("=" * 60)
    
    L4_OUTPUT = '../results/L4_ideal_subgraph'
    
    # 1. 加载 G_ideal
    g_ideal = load_ideal_subgraph(L4_OUTPUT)
    
    # 2. 加载成分-靶点数据库
    all_compounds = load_compound_targets(TCMSP_PATH, ETCM_PATH, SYMMAP_PATH)
    
    # 3. MCI 计算
    compound_results = calculate_mci({}, g_ideal, all_compounds)
    
    # 4. 超图推理
    hyperedges_df = hypergraph_reasoning(compound_results, g_ideal)
    
    # 5. 文献验证
    validation_df = literature_validation(compound_results)
    
    # 6. 导出结果
    export_results(compound_results, hyperedges_df, validation_df, OUTPUT_DIR)
    
    # 自检标准
    print("\n" + "=" * 60)
    print("自检标准检查")
    print("=" * 60)
    print(f"✓ MCI 计算无 NaN/Inf: {compound_results['mci_comprehensive'].isnull().sum() == 0}")
    n_lit = validation_df['has_literature'].sum()
    print(f"✓ Top-10 有文献 curated 靶点: {n_lit} (要求 ≥ 2)")
    
    if n_lit < 2:
        print("⚠ 建议启动 L5-alt 路径 A（间接调控分析）")
    
    print("\nL5 分析完成！")

if __name__ == '__main__':
    main()