"""诊断模型性能弱的原因"""
import pandas as pd
from pathlib import Path

LOG_FILE = Path("logs/training_log.csv")

def main():
    if not LOG_FILE.exists():
        print("训练日志不存在")
        return
    
    df = pd.read_csv(LOG_FILE)
    
    print("="*60)
    print("模型性能诊断")
    print("="*60)
    
    # 找到最佳epoch
    best_idx = df["val_pr_auc"].idxmax()
    best_epoch = int(df.loc[best_idx, "epoch"])
    
    print(f"\n最佳验证集指标 (Epoch {best_epoch}):")
    print(f"  Val_PR_AUC (AUPRC): {df.loc[best_idx, 'val_pr_auc']:.4f}")
    print(f"  Val_F1: {df.loc[best_idx, 'val_f1']:.4f}")
    print(f"  Train_Loss: {df.loc[best_idx, 'train_loss']:.4f}")
    
    # 训练动态
    print(f"\n训练动态分析:")
    print(f"  Epoch 1: AUPRC={df.loc[0, 'val_pr_auc']:.4f}, F1={df.loc[0, 'val_f1']:.4f}")
    print(f"  Epoch 10: AUPRC={df.loc[9, 'val_pr_auc']:.4f}, F1={df.loc[9, 'val_f1']:.4f}")
    print(f"  Best Epoch {best_epoch}: AUPRC={df.loc[best_idx, 'val_pr_auc']:.4f}, F1={df.loc[best_idx, 'val_f1']:.4f}")
    
    # 检查训练动态
    train_loss_trend = df["train_loss"].iloc[:best_epoch+1]
    print(f"  训练Loss趋势: {train_loss_trend.iloc[0]:.4f} → {train_loss_trend.iloc[-1]:.4f} ({'下降' if train_loss_trend.iloc[-1] < train_loss_trend.iloc[0] else '上升'})")
    
    # 类别不平衡
    print(f"\n类别不平衡分析:")
    n_pos, n_neg = 175, 2784
    print(f"  阳性样本: {n_pos}")
    print(f"  阴性样本: {n_neg}")
    print(f"  阳性比例: {n_pos/(n_pos+n_neg)*100:.2f}%")
    print(f"  类别比: 1:{n_neg/n_pos:.0f}")
    
    # Baseline
    print(f"\n理论Baseline:")
    print(f"  随机猜测AUPRC: {n_pos/(n_pos+n_neg):.4f}")
    print(f"  全预测阴性F1: 0.0000")
    print(f"  全预测阳性Recall=1.0, Precision={n_pos/(n_pos+n_neg):.4f}, F1={2*n_pos/(2*n_pos+n_neg):.4f}")
    
    # 诊断
    model_auprc = df.loc[best_idx, 'val_pr_auc']
    model_f1 = df.loc[best_idx, 'val_f1']
    
    print(f"\n⚠ 问题诊断:")
    print(f"  1. AUPRC={model_auprc:.4f}, F1={model_f1:.4f}")
    print(f"  2. AUPRC仅比随机baseline({n_pos/(n_pos+n_neg):.4f})高{(model_auprc - n_pos/(n_pos+n_neg))/n_pos/(n_pos+n_neg)*100:.1f}倍")
    print(f"  3. 模型判别能力弱的原因：")
    print(f"     - 标签噪音：阳性标签来自计算预测，缺乏实验验证")
    print(f"     - 类别极度不平衡：1:16，正类仅175个样本")
    print(f"     - 图结构噪音：STRING PPI包含大量间接相互作用")
    print(f"     - 特征工程：40维特征可能包含噪音")
    print(f"     - 正类权重不足：2.0对1:16的比不够")
    print(f"\n  4. 优化建议：")
    print(f"     - 增加正类权重至10-16")
    print(f"     - 使用Focal Loss")
    print(f"     - 特征选择：去除低方差/高相关特征")
    print(f"     - 图数据增强：BAT")
    print(f"     - 集成学习：多模型投票")

if __name__ == "__main__":
    main()
