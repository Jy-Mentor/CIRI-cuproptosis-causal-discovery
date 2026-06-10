import pandas as pd
from pathlib import Path

PROCESSED_DIR = Path("processed")

# 读取当前标签
labels = pd.read_csv(PROCESSED_DIR / "labels.csv")

# 转录因子/结构蛋白/血浆蛋白黑名单（无成药性）
BLACKLIST = {
    "FGA", "SNAI1", "ZEB2", "EGR1", "RELA",
    "KEAP1", "DDIT3", "CAT", "SOD1", "CP",
    "BAX", "TERT", "IGF1R", "NFE2L2",
    # 边缘案例：STAT3/NR1H3/F3/GSR — 暂时保留，由用户决定是否加入
}

print("=" * 60)
print("阳性标签黑名单剔除补丁")
print("=" * 60)

# 剔除前统计
n_pos_before = (labels["Label"] == 1).sum()
print(f"\n剔除前阳性标签: {n_pos_before} 个")

# 从阳性中剔除黑名单
demoted_mask = (labels["Label"] == 1) & (labels["GeneSymbol"].isin(BLACKLIST))
demoted_genes = labels.loc[demoted_mask, "GeneSymbol"].tolist()
labels.loc[demoted_mask, "Label"] = -1

# 保存为最终标签
labels.to_csv(PROCESSED_DIR / "labels_final.csv", index=False)

# 最终统计
n_pos = (labels["Label"] == 1).sum()
n_neg = (labels["Label"] == 0).sum()
n_unk = (labels["Label"] == -1).sum()

print(f"\n最终标签分布:")
print(f"  阳性(1): {n_pos}")
print(f"  阴性(0): {n_neg}")
print(f"  未知(-1): {n_unk}")
print(f"  阳性/阴性比: 1:{n_neg/max(n_pos,1):.1f}")

print(f"\n被降级的基因 ({len(demoted_genes)} 个):")
print(f"  {sorted(demoted_genes)}")

# 列出最终阳性基因
final_positive = labels[labels["Label"] == 1]["GeneSymbol"].tolist()
print(f"\n最终阳性基因列表 ({len(final_positive)} 个):")
print(f"  {sorted(final_positive)}")

print("\n" + "=" * 60)
print("labels_final.csv 已生成。请人工确认阳性列表后")
print("将其重命名为 labels.csv 或让 GAT 模块读取 labels_final.csv")
print("=" * 60)
