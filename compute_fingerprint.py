# -*- coding: utf-8 -*-
"""
计算 β-石竹烯 (Beta-caryophyllene) 的 Morgan 指纹 (ECFP4)
并保存为 CSV 文件
"""

import csv
import sys


def main():
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        print("RDKit 未安装。请使用以下命令安装：")
        print("  conda install -c conda-forge rdkit")
        print("  或")
        print("  pip install rdkit-pypi")
        sys.exit(1)

    # β-石竹烯的 SMILES
    smiles = "CC1=CCCC(=C)C2CC(C2(C)C)CC1"
    drug_name = "β-石竹烯 (Beta-caryophyllene)"

    print(f"正在处理药物: {drug_name}")
    print(f"SMILES: {smiles}")

    # 从 SMILES 构建分子对象
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("错误: 无法从 SMILES 构建分子，请检查 SMILES 是否正确。")
        sys.exit(1)

    print("分子构建成功。")

    # 计算 Morgan 指纹 (ECFP4)
    # radius=2 对应 ECFP4（直径=4），nBits=1024 生成 1024 位指纹向量
    radius = 2
    n_bits = 1024
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)

    # 转换为位列表 (0/1)
    fp_list = list(fp)
    print(f"指纹计算完成。位数: {len(fp_list)}")

    # 保存为 CSV 文件
    output_file = "drug_fingerprint.csv"
    columns = [f"fp_{i}" for i in range(n_bits)]

    try:
        with open(output_file, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(columns)
            writer.writerow(fp_list)
        print(f"\n文件保存成功: {output_file}")
    except IOError as e:
        print(f"错误: 无法写入文件 {output_file}: {e}")
        sys.exit(1)

    # 打印前 10 个指纹位作为检查
    print(f"\n前 10 个指纹位 (fp_0 ~ fp_9):")
    print("-" * 40)
    for i in range(10):
        print(f"  fp_{i:4d} = {fp_list[i]}")
    print("-" * 40)
    print(f"共 {len(fp_list)} 位，其中值为 1 的位数: {sum(fp_list)}")
    print("脚本执行完毕。")


if __name__ == "__main__":
    main()