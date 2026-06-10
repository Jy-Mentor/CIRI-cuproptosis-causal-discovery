"""输出完整的预测结果统计"""
import pandas as pd
from pathlib import Path

RESULTS_DIR = Path("results")
CUPROPTOSIS_GENES = {
    "FDX1", "LIAS", "LIPT1", "DLAT", "PDHB", "PDHX", "SLC31A1",
    "ATP7B", "ATOX1", "MTF1", "GLS", "CDKN2A",
    "ATP7A", "COMMD1", "SCO1", "SCO2", "COX17", "COX11"
}

def main():
    # 读取预测结果
    pred_file = RESULTS_DIR / "all_unknown_predictions.csv"
    if not pred_file.exists():
        print(f"✗ 预测结果文件不存在: {pred_file}")
        return
    
    df = pd.read_csv(pred_file)
    
    print("="*80)
    print("石竹烯-CIRI 靶点预测 - 完整结果统计")
    print("="*80)
    
    # 基本统计
    print(f"\n【数据集概览】")
    print(f"  总预测基因数: {len(df):,}")
    print(f"  铜死亡基因数: {df['is_cuproptosis'].sum():.0f}")
    print(f"  P_target范围: {df['P_target'].min():.4f} - {df['P_target'].max():.4f}")
    print(f"  P_target均值: {df['P_target'].mean():.4f}")
    print(f"  P_target中位数: {df['P_target'].median():.4f}")
    
    # 分数分布
    print(f"\n【分数分布】")
    thresholds = [0.9, 0.85, 0.8, 0.75, 0.7, 0.6, 0.5]
    for thresh in thresholds:
        count = (df["P_target"] >= thresh).sum()
        pct = count / len(df) * 100
        print(f"  P >= {thresh:.2f}: {count:>6,} 个 ({pct:>5.2f}%)")
    
    # Top 50
    print(f"\n【Top 50 候选靶点】")
    top50 = df.head(50)
    print(f"  {'排名':<6} {'基因名':<12} {'P_target':<10} {'铜死亡':<6}")
    print(f"  {'-'*40}")
    for _, row in top50.iterrows():
        flag = "✓" if row["is_cuproptosis"] == 1 else ""
        print(f"  {int(row['Rank']):<6} {row['GeneSymbol']:<12} {row['P_target']:<10.4f} {flag}")
    
    # 铜死亡基因统计
    cupro_df = df[df["GeneSymbol"].isin(CUPROPTOSIS_GENES)]
    print(f"\n【铜死亡基因统计】")
    print(f"  总数: {len(cupro_df)}/{len(CUPROPTOSIS_GENES)}")
    if len(cupro_df) > 0:
        print(f"  平均排名: {cupro_df['Rank'].mean():.0f}")
        print(f"  最佳排名: {int(cupro_df['Rank'].min())} ({cupro_df.loc[cupro_df['Rank'].idxmin(), 'GeneSymbol']})")
        print(f"  平均P_target: {cupro_df['P_target'].mean():.4f}")
        
        # 按分数区间统计铜死亡基因
        print(f"\n  铜死亡基因分数分布:")
        for thresh in [0.82, 0.80, 0.78, 0.75]:
            count = (cupro_df["P_target"] >= thresh).sum()
            print(f"    P >= {thresh:.2f}: {count} 个")
    
    # 保存完整统计
    output_file = RESULTS_DIR / "full_prediction_summary.csv"
    df.to_csv(output_file, index=False)
    print(f"\n完整预测结果已保存至: {output_file}")

if __name__ == "__main__":
    main()
