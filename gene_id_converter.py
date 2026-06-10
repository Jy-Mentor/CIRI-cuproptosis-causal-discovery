#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
基因ID转换工具
支持大鼠、小鼠和人类基因之间的相互转换
"""

import os
import sys
import argparse
import pandas as pd
from typing import Dict, List, Tuple

def install_package(package: str) -> None:
    """
    智能安装包，只安装未安装的包
    """
    try:
        __import__(package)
    except ImportError:
        print(f"正在安装 {package}...")
        import subprocess
        subprocess.check_call([sys.executable, "-m", "pip", "install", package])

# 智能安装必要的包
install_package("pandas")

def parse_args() -> argparse.Namespace:
    """
    解析命令行参数
    """
    parser = argparse.ArgumentParser(description="基因ID转换工具")
    parser.add_argument("--input", required=True, help="输入基因列表文件路径")
    parser.add_argument("--output", required=True, help="输出转换结果文件路径")
    parser.add_argument("--from-species", required=True, choices=["rat", "mouse", "human"], help="源物种")
    parser.add_argument("--to-species", required=True, choices=["rat", "mouse", "human"], help="目标物种")
    parser.add_argument("--mapping-file", default="C:\\Users\\Jy-Mentor-7\\Desktop\\大创\\大鼠 小鼠 人类映射库.txt", help="映射库文件路径")
    parser.add_argument("--unmapped-output", default="unmapped_genes.txt", help="未映射基因输出文件路径")
    return parser.parse_args()

def load_mapping_file(mapping_file: str) -> Dict[str, Dict[str, str]]:
    """
    加载映射库文件，构建物种间的映射关系
    """
    mapping_dict = {}
    try:
        df = pd.read_csv(mapping_file, sep="\t", header=0)
        # 确保列名正确
        columns = df.columns.tolist()
        print(f"映射库文件列名: {columns}")
        
        # 构建映射关系
        for _, row in df.iterrows():
            # 提取各个物种的基因ID
            rat_gene = str(row.get("大鼠基因", "")) if pd.notna(row.get("大鼠基因")) else ""
            mouse_gene = str(row.get("小鼠基因", "")) if pd.notna(row.get("小鼠基因")) else ""
            human_gene = str(row.get("人类基因", "")) if pd.notna(row.get("人类基因")) else ""
            
            # 构建双向映射
            if rat_gene:
                if rat_gene not in mapping_dict:
                    mapping_dict[rat_gene] = {}
                if mouse_gene:
                    mapping_dict[rat_gene]["mouse"] = mouse_gene
                if human_gene:
                    mapping_dict[rat_gene]["human"] = human_gene
            
            if mouse_gene:
                if mouse_gene not in mapping_dict:
                    mapping_dict[mouse_gene] = {}
                if rat_gene:
                    mapping_dict[mouse_gene]["rat"] = rat_gene
                if human_gene:
                    mapping_dict[mouse_gene]["human"] = human_gene
            
            if human_gene:
                if human_gene not in mapping_dict:
                    mapping_dict[human_gene] = {}
                if rat_gene:
                    mapping_dict[human_gene]["rat"] = rat_gene
                if mouse_gene:
                    mapping_dict[human_gene]["mouse"] = mouse_gene
        
        print(f"成功加载映射库，共 {len(mapping_dict)} 个基因条目")
        return mapping_dict
    except Exception as e:
        print(f"加载映射库失败: {e}")
        sys.exit(1)

def load_input_genes(input_file: str) -> List[str]:
    """
    加载输入基因列表
    """
    genes = []
    try:
        with open(input_file, 'r', encoding='utf-8') as f:
            for line in f:
                gene = line.strip()
                if gene:
                    genes.append(gene)
        print(f"成功加载输入基因列表，共 {len(genes)} 个基因")
        return genes
    except Exception as e:
        print(f"加载输入文件失败: {e}")
        sys.exit(1)

def convert_genes(genes: List[str], from_species: str, to_species: str, mapping_dict: Dict[str, Dict[str, str]]) -> Tuple[List[Tuple[str, str]], List[str]]:
    """
    转换基因ID
    返回 (转换结果, 未映射基因)
    """
    results = []
    unmapped = []
    
    for gene in genes:
        if gene in mapping_dict and to_species in mapping_dict[gene]:
            converted_gene = mapping_dict[gene][to_species]
            results.append((gene, converted_gene))
        else:
            unmapped.append(gene)
    
    print(f"转换完成: {len(results)} 个基因成功转换, {len(unmapped)} 个基因未映射")
    return results, unmapped

def save_results(results: List[Tuple[str, str]], output_file: str) -> None:
    """
    保存转换结果
    """
    try:
        df = pd.DataFrame(results, columns=["原始基因", "转换基因"])
        df.to_csv(output_file, sep="\t", index=False, encoding='utf-8')
        print(f"转换结果已保存到: {output_file}")
    except Exception as e:
        print(f"保存结果失败: {e}")
        sys.exit(1)

def save_unmapped(unmapped: List[str], unmapped_output: str) -> None:
    """
    保存未映射基因
    """
    try:
        with open(unmapped_output, 'w', encoding='utf-8') as f:
            for gene in unmapped:
                f.write(f"{gene}\n")
        print(f"未映射基因已保存到: {unmapped_output}")
    except Exception as e:
        print(f"保存未映射基因失败: {e}")
        sys.exit(1)

def main():
    """
    主函数
    """
    # 解析命令行参数
    args = parse_args()
    
    # 加载映射库
    mapping_dict = load_mapping_file(args.mapping_file)
    
    # 加载输入基因列表
    genes = load_input_genes(args.input)
    
    # 转换基因
    results, unmapped = convert_genes(genes, args.from_species, args.to_species, mapping_dict)
    
    # 保存结果
    save_results(results, args.output)
    
    # 保存未映射基因
    if unmapped:
        save_unmapped(unmapped, args.unmapped_output)
    else:
        print("所有基因都成功映射")

if __name__ == "__main__":
    main()
