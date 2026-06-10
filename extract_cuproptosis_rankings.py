"""输出铜死亡基因的完整排名和得分"""
import pandas as pd
from pathlib import Path

RESULTS_DIR = Path("results")
CUPROPTOSIS_GENES = {
    # 铜死亡执行基因
    "FDX1", "LIAS", "LIPT1", "DLAT", "PDHB", "PDHX", "SLC31A1",
    # 铜死亡调控基因
    "ATP7B", "ATOX1", "MTF1", "GLS", "CDKN2A",
    # 铜离子结合相关
    "ATP7A", "COMMD1", "SCO1", "SCO2", "COX17", "COX11"
}

def main():
    # 读取预测结果
    pred_file = RESULTS_DIR / "all_unknown_predictions.csv"
    if not pred_file.exists():
        print(f"✗ 预测结果文件不存在: {pred_file}")
        return
    
    df = pd.read_csv(pred_file)
    
    # 筛选铜死亡基因
    cupro_df = df[df["GeneSymbol"].isin(CUPROPTOSIS_GENES)].copy()
    cupro_df = cupro_df.sort_values("Rank")
    
    print("="*80)
    print("铜死亡基因完整排名")
    print("="*80)
    print(f"\n总预测基因数: {len(df)}")
    print(f"铜死亡基因数: {len(cupro_df)}/{len(CUPROPTOSIS_GENES)}")
    
    print(f"\n{'排名':<8} {'基因名':<12} {'P_target':<12} {'铜死亡':<8} {'距离铜死亡':<12} {'分类'}")
    print("-"*80)
    
    for _, row in cupro_df.iterrows():
        gene = row["GeneSymbol"]
        rank = int(row["Rank"])
        p_target = row["P_target"]
        is_cupro = "✓" if row["is_cuproptosis"] == 1 else "✗"
        dist = row["dist_to_cuproptosis"]
        
        # 分类
        if gene in {"FDX1", "LIAS", "LIPT1", "DLAT", "PDHB", "PDHX", "SLC31A1"}:
            category = "执行基因"
        elif gene in {"ATP7B", "ATOX1", "MTF1", "GLS", "CDKN2A", "ATP7A", "COMMD1"}:
            category = "调控基因"
        elif gene in {"SCO1", "SCO2", "COX17", "COX11"}:
            category = "铜代谢"
        else:
            category = "其他"
        
        print(f"{rank:<8} {gene:<12} {p_target:<12.4f} {is_cupro:<8} {dist:<12.2f} {category}")
    
    # 统计
    print(f"\n{'='*80}")
    print("统计分析")
    print("="*80)
    
    if len(cupro_df) > 0:
        print(f"  平均排名: {cupro_df['Rank'].mean():.1f}")
        print(f"  中位数排名: {cupro_df['Rank'].median():.0f}")
        print(f"  最佳排名: {cupro_df['Rank'].min():.0f} ({cupro_df.loc[cupro_df['Rank'].idxmin(), 'GeneSymbol']})")
        print(f"  最差排名: {cupro_df['Rank'].max():.0f} ({cupro_df.loc[cupro_df['Rank'].idxmax(), 'GeneSymbol']})")
        print(f"  平均P_target: {cupro_df['P_target'].mean():.4f}")
        print(f"  P_target范围: {cupro_df['P_target'].min():.4f} - {cupro_df['P_target'].max():.4f}")
        
        # 按分类统计
        print(f"\n  按分类统计:")
        for category in ["执行基因", "调控基因", "铜代谢", "其他"]:
            cat_df = cupro_df[cupro_df.apply(lambda r: category in r.get("分类", ""), axis=1)]
            if len(cat_df) > 0:
                print(f"    {category}: {len(cat_df)}个, 平均排名={cat_df['Rank'].mean():.0f}, 平均P={cat_df['P_target'].mean():.4f}")
    
    # 保存结果
    output_file = RESULTS_DIR / "cuproptosis_gene_rankings.csv"
    cupro_df.to_csv(output_file, index=False)
    print(f"\n铜死亡基因排名已保存至: {output_file}")

if __name__ == "__main__":
    main()
