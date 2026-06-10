# -*- coding: utf-8 -*-
"""
汇总最新正确结果到Excel (使用CIRI和BCP靶点文件夹的数据)
"""
import os
import sys
import pandas as pd
from datetime import datetime
import shutil

# 正确的源数据路径
SOURCE_RESULTS = r"c:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\CIRI和BCP靶点\results"
SOURCE_SCRIPTS = r"c:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\CIRI和BCP靶点\scripts_v9"

# 目标路径
TARGET_DIR = r"c:\Users\Jy-Mentor-7\Desktop\BCP∩CIRI"
OUTPUT_FILE = os.path.join(TARGET_DIR, "CIRI_BCP_多组学分析结果汇总.xlsx")

def load_csv_safe(filepath):
    """安全加载CSV文件"""
    try:
        if os.path.exists(filepath):
            return pd.read_csv(filepath)
    except Exception as e:
        print(f"  警告: 无法加载 {filepath}: {e}")
    return None

def create_summary_excel():
    """创建汇总Excel文件 - 只包含有用的核心结果"""
    print("=" * 60)
    print("创建多组学分析结果汇总Excel (正确版本)")
    print("=" * 60)
    print(f"数据源: {SOURCE_RESULTS}")
    print(f"输出: {OUTPUT_FILE}")
    print()
    
    os.makedirs(TARGET_DIR, exist_ok=True)
    
    with pd.ExcelWriter(OUTPUT_FILE, engine='openpyxl') as writer:
        
        # 1. 分析概览
        print("[0] 创建: 分析概览")
        overview_data = {
            '分析项目': [
                '分析日期',
                '管线版本',
                '总评估基因',
                'Tier1靶点',
                'GRN扰动基因',
                'ML特征基因',
                'PPI网络节点',
                'GAT验证基因',
                '差异表达基因',
                'WGCNA模块基因',
                '铜死亡核心基因在Tier1',
                'ML模型AUC',
                'GAT R²',
                '框架类型',
            ],
            '结果': [
                datetime.now().strftime('%Y-%m-%d'),
                'v8.0 无偏管线',
                '9,853',
                '40',
                '105',
                '50',
                '1,998',
                '1,998',
                '1,836',
                '2,029',
                '0 (数据驱动)',
                '0.85',
                '0.7778',
                '无先验权重 (BCP_prior=0, Cuproptosis_prior=0)',
            ]
        }
        overview_df = pd.DataFrame(overview_data)
        overview_df.to_excel(writer, sheet_name='0.分析概览', index=False)
        print("  ✓ 分析概览")
        
        # 2. Tier1核心靶点 (最重要的结果)
        print("\n[1] 处理: Tier1核心靶点")
        df = load_csv_safe(os.path.join(SOURCE_RESULTS, "stage8_final_targets", "tier1_targets.csv"))
        if df is not None:
            # 检查并删除冗余列
            redundant_cols = ['BCP_prior', 'Cuproptosis_prior', 'Pathway_coreness']
            df_clean = df.drop(columns=[c for c in redundant_cols if c in df.columns], errors='ignore')
            df_clean.to_excel(writer, sheet_name='1.Tier1核心靶点', index=False)
            print(f"  ✓ Tier1靶点: {len(df_clean)}个基因")
            print(f"     列: {list(df_clean.columns)}")
        
        # 3. 全部靶点排名 (Top 200)
        print("\n[2] 处理: 全部靶点排名")
        df = load_csv_safe(os.path.join(SOURCE_RESULTS, "stage8_final_targets", "core_targets.csv"))
        if df is not None:
            redundant_cols = ['BCP_prior', 'Cuproptosis_prior', 'Pathway_coreness']
            df_clean = df.drop(columns=[c for c in redundant_cols if c in df.columns], errors='ignore')
            df_clean.head(200).to_excel(writer, sheet_name='2.全部靶点排名Top200', index=False)
            print(f"  ✓ 全部靶点: Top 200 / {len(df_clean)} 总基因")
        
        # 4. GRN扰动评分
        print("\n[3] 处理: GRN扰动评分")
        df = load_csv_safe(os.path.join(SOURCE_RESULTS, "stage6_graphsage_knockout", "gene_perturbation_scores.csv"))
        if df is not None:
            df.to_excel(writer, sheet_name='3.GRN扰动评分', index=False)
            print(f"  ✓ GRN扰动: {len(df)}个基因")
        
        # 5. ML重要性
        print("\n[4] 处理: ML重要性")
        df = load_csv_safe(os.path.join(SOURCE_RESULTS, "stage7_ml_shap", "gene_shap_importance.csv"))
        if df is not None:
            df.to_excel(writer, sheet_name='4.ML重要性', index=False)
            print(f"  ✓ ML重要性: {len(df)}个基因")
        
        # 6. PPI网络度排名
        print("\n[5] 处理: PPI网络度排名")
        df = load_csv_safe(os.path.join(SOURCE_RESULTS, "stage5_string_ppi", "node_degree_ranking.csv"))
        if df is not None:
            df.to_excel(writer, sheet_name='5.PPI度排名', index=False)
            print(f"  ✓ PPI度排名: {len(df)}个基因")
        
        # 7. GAT排名
        print("\n[6] 处理: GAT排名")
        df = load_csv_safe(os.path.join(SOURCE_RESULTS, "stage9_ppi_gat", "gat_gene_ranking.csv"))
        if df is not None:
            df.to_excel(writer, sheet_name='6.GAT排名', index=False)
            print(f"  ✓ GAT排名: {len(df)}个基因")
        
        # 8. GAT铜死亡验证
        print("\n[7] 处理: GAT铜死亡验证")
        df = load_csv_safe(os.path.join(SOURCE_RESULTS, "stage9_ppi_gat", "cuproptosis_validation.csv"))
        if df is not None:
            df.to_excel(writer, sheet_name='7.GAT铜死亡验证', index=False)
            print(f"  ✓ 铜死亡验证: {len(df)}个基因")
        
        # 9. 差异表达基因 (只取显著的前500)
        print("\n[8] 处理: 差异表达基因")
        df = load_csv_safe(os.path.join(SOURCE_RESULTS, "stage1_rma_degs", "limma_degs.csv"))
        if df is not None:
            # 按adjPVal排序取前500
            if 'adjPVal' in df.columns:
                df_sig = df.nsmallest(500, 'adjPVal')
            else:
                df_sig = df.head(500)
            df_sig.to_excel(writer, sheet_name='8.差异表达基因Top500', index=False)
            print(f"  ✓ DEGs: Top 500 / {len(df)} 总基因")
        
        # 10. WGCNA模块
        print("\n[9] 处理: WGCNA模块")
        df = load_csv_safe(os.path.join(SOURCE_RESULTS, "stage4_seed_wgcna", "wgcna_modules.csv"))
        if df is not None:
            df.to_excel(writer, sheet_name='9.WGCNA模块', index=False)
            print(f"  ✓ WGCNA: {len(df)}个基因")
        
        # 11. ML模型性能
        print("\n[10] 处理: ML模型性能")
        df = load_csv_safe(os.path.join(SOURCE_RESULTS, "stage7_ml_shap", "ml_model_performance.csv"))
        if df is not None:
            df.to_excel(writer, sheet_name='10.ML模型性能', index=False)
            print(f"  ✓ ML性能: {len(df)}个模型")
        
        # 12. GAT训练性能
        print("\n[11] 处理: GAT训练性能")
        df = load_csv_safe(os.path.join(SOURCE_RESULTS, "stage9_ppi_gat", "gat_training_loss.csv"))
        if df is not None:
            df.to_excel(writer, sheet_name='11.GAT训练曲线', index=False)
            print(f"  ✓ GAT训练: {len(df)}个epoch")
    
    print("\n" + "=" * 60)
    print(f"✓ Excel文件已保存: {OUTPUT_FILE}")
    print("=" * 60)

def copy_scripts():
    """复制脚本到目标文件夹"""
    print("\n" + "=" * 60)
    print("复制脚本到目标文件夹")
    print("=" * 60)
    
    target_scripts = os.path.join(TARGET_DIR, "scripts")
    os.makedirs(target_scripts, exist_ok=True)
    
    script_count = 0
    for filename in os.listdir(SOURCE_SCRIPTS):
        if filename.endswith('.py') or filename.endswith('.R') or filename == 'config.py':
            src = os.path.join(SOURCE_SCRIPTS, filename)
            dst = os.path.join(target_scripts, filename)
            shutil.copy2(src, dst)
            script_count += 1
            print(f"  ✓ {filename}")
    
    # 创建脚本清单
    readme = f"""# 脚本清单

生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
管线版本: v8.0 无偏多组学整合

## 核心Stage脚本 (Stage 0-9)

| 脚本 | 功能 |
|-----|------|
| stage0_data_validation.py | 数据验证层 |
| stage1_rma_degs.R | Bulk RNA-seq差异分析 |
| stage2_single_cell.py | 单细胞数据处理 |
| stage3_enrichment.py | 功能富集分析 |
| stage4_seed_wgcna.py | 种子池构建+WGCNA |
| stage5_string_ppi_local.py | PPI网络构建(本地STRING) |
| stage6_graphsage_knockout.py | GRN虚拟敲除 |
| stage7_ml_shap.py | 机器学习+SHAP |
| stage8_final_targets.py | 多组学融合排名 (无偏) |
| stage9_ppi_gat.py | GAT图神经网络验证 |

## 铜死亡验证模块 (M1-M6)

| 脚本 | 功能 |
|-----|------|
| cuproptosis_gsva.py | 铜死亡通路GSVA |
| cuproptosis_gsea.py | 铜死亡GSEA |
| cuproptosis_wgcna.py | 铜死亡WGCNA分析 |
| cuproptosis_singlecell.py | 单细胞铜死亡 |
| cuproptosis_ppi_neighbors.py | PPI邻居分析 |
| cuproptosis_immunology.py | 免疫浸润分析 |

## 配置文件

- config.py - 主配置
- utils.py - 工具函数

## 总脚本数: {script_count}个

## 关键特性

✅ 无偏管线: BCP_prior=0, Cuproptosis_prior=0  
✅ 三维度自由竞争: GRN(45.9%) + ML(13.1%) + PPI(40.9%)  
✅ 数据驱动: 不预设铜死亡优先  
✅ 独立验证: GAT R²=0.7778
"""
    
    readme_path = os.path.join(target_scripts, "脚本清单.md")
    with open(readme_path, 'w', encoding='utf-8') as f:
        f.write(readme)
    
    print(f"\n  ✓ 脚本清单: 脚本清单.md")
    print(f"\n完成! 共复制 {script_count} 个脚本")

if __name__ == "__main__":
    create_summary_excel()
    copy_scripts()
    
    print("\n" + "=" * 60)
    print("所有任务完成!")
    print(f"目标文件夹: {TARGET_DIR}")
    print("=" * 60)
