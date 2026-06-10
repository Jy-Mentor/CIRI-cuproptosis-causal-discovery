"""
GSE23160与GSE174574数据整合分析
用于L2c NeuralODE模块的时间校准
"""

import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import json
import warnings
warnings.filterwarnings('ignore')

plt.rcParams['font.size'] = 12
plt.rcParams['axes.linewidth'] = 1.5

BULK_DIR = Path(r'D:\反向网络药理学\L1 数据集\bulk\GSE23160(主验证集时序差异分析，2h,8h,24h)\GSE23160_limma_results')
SCRNA_DIR = Path(r'C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\GSE174574_analysis')

CUPROPTOSIS_CORE = ['Fdx1', 'Lias', 'Lipt1', 'Dlat', 'Dld', 'Pdha1', 'Pdhb', 'Mtf1', 'Gls', 'Cdkn2a']
CUPROPTOSIS_EXT = ['Sirt7', 'Atp7b', 'Slc31a1', 'Cox17', 'Atox1', 'Ccs']
CUPROPTOSIS_ALL = CUPROPTOSIS_CORE + CUPROPTOSIS_EXT

def load_bulk_data():
    """加载Bulk时序数据"""
    print("="*60)
    print("加载GSE23160 Bulk时序数据")
    print("="*60)
    
    copper_deg = pd.read_csv(BULK_DIR / 'copper_death_genes_DEG.txt', sep='\t')
    print(f"\ncopper_death_genes_DEG.txt:")
    print(f"  行数: {len(copper_deg)}")
    print(f"  列: {copper_deg.columns.tolist()}")
    
    sig_2h = pd.read_csv(BULK_DIR / 'sig_DEGs_2h.txt', sep='\t')
    sig_8h = pd.read_csv(BULK_DIR / 'sig_DEGs_8h.txt', sep='\t')
    sig_24h = pd.read_csv(BULK_DIR / 'sig_DEGs_24h.txt', sep='\t')
    
    print(f"\n显著差异基因:")
    print(f"  2h: {len(sig_2h)} genes")
    print(f"  8h: {len(sig_8h)} genes")
    print(f"  24h: {len(sig_24h)} genes")
    
    return copper_deg, {'2h': sig_2h, '8h': sig_8h, '24h': sig_24h}

def analyze_cuproptosis_temporal(copper_deg, sig_degs):
    """分析铜死亡基因的时序表达模式"""
    print("\n" + "="*60)
    print("铜死亡基因时序分析")
    print("="*60)
    
    copper_deg['Gene_Symbol_upper'] = copper_deg['Gene_Symbol'].str.upper()
    
    results = {}
    for gene in CUPROPTOSIS_ALL:
        gene_upper = gene.upper()
        gene_data = copper_deg[copper_deg['Gene_Symbol_upper'] == gene_upper]
        
        if len(gene_data) > 0:
            time_points = {}
            for tp in ['2h', '8h', '24h']:
                tp_data = gene_data[gene_data['time_point'] == tp]
                if len(tp_data) > 0:
                    time_points[tp] = {
                        'logFC': float(tp_data['logFC'].iloc[0]),
                        'pvalue': float(tp_data['P.Value'].iloc[0]),
                        'adj_pval': float(tp_data['adj.P.Val'].iloc[0]),
                        'significance': str(tp_data['significance'].iloc[0])
                    }
            
            if time_points:
                results[gene] = time_points
    
    print(f"\n铜死亡基因在Bulk数据中检出: {len(results)}/{len(CUPROPTOSIS_ALL)}")
    for gene, tp_data in results.items():
        print(f"\n  {gene}:")
        for tp, vals in tp_data.items():
            print(f"    {tp}: logFC={vals['logFC']:.4f}, adj.P={vals['adj_pval']:.6f}, {vals['significance']}")
    
    return results

def create_temporal_figures(copper_results, sig_degs):
    """生成时序可视化图"""
    print("\n" + "="*60)
    print("生成时序可视化图")
    print("="*60)
    
    fig_dir = SCRNA_DIR / 'figures'
    fig_dir.mkdir(parents=True, exist_ok=True)
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    time_points = ['2h', '8h', '24h']
    time_values = [2, 8, 24]
    
    for idx, (ax, gene) in enumerate(zip(axes.flat[:4], ['Fdx1', 'Lias', 'Dlat', 'Cdkn2a'])):
        if gene in copper_results:
            logfcs = []
            pvals = []
            tps = []
            for tp in time_points:
                if tp in copper_results[gene]:
                    logfcs.append(copper_results[gene][tp]['logFC'])
                    pvals.append(copper_results[gene][tp]['adj_pval'])
                    tps.append(tp)
            
            if logfcs:
                ax.plot(tps, logfcs, marker='o', linewidth=2, markersize=10, color='#2980b9')
                for i, (tp, lfc) in enumerate(zip(tps, logfcs)):
                    ax.annotate(f'{lfc:.2f}', (tp, lfc), textcoords="offset points",
                              xytext=(0, 10), ha='center', fontsize=10)
                
                ax.axhline(y=0, color='red', linestyle='--', alpha=0.5)
                ax.set_title(f'{gene} (Mouse)', fontweight='bold')
                ax.set_xlabel('Time Point')
                ax.set_ylabel('log2 Fold Change')
                ax.grid(True, alpha=0.3)
    
    for ax in axes.flat[4:]:
        ax.axis('off')
    
    plt.suptitle('GSE23160: Cuproptosis Genes Temporal Expression', fontsize=16, fontweight='bold')
    plt.tight_layout()
    fig_path = fig_dir / 'Figure_2_Cuproptosis_Temporal.png'
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  保存: Figure_2_Cuproptosis_Temporal.png")
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    for tp_name, tp_data in sig_degs.items():
        up = len(tp_data[tp_data['significance'] == 'Upregulated'])
        down = len(tp_data[tp_data['significance'] == 'Downregulated'])
        ax.bar([tp_name], [up], bottom=0, color='#e74c3c', label='Up' if tp_name == '2h' else "")
        ax.bar([tp_name], [down], bottom=up, color='#3498db', label='Down' if tp_name == '2h' else "")
    
    ax.legend(['Upregulated', 'Downregulated'])
    ax.set_title('GSE23160: DEG Counts by Time Point', fontweight='bold')
    ax.set_ylabel('Gene Count')
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    fig_path = fig_dir / 'Figure_3_DEG_Counts.png'
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  保存: Figure_3_DEG_Counts.png")
    
    return fig_path

def prepare_l2c_temporal_calibration(copper_results, sig_degs):
    """准备L2c时间校准参数"""
    print("\n" + "="*60)
    print("准备L2c时间校准参数")
    print("="*60)
    
    l2c_dir = SCRNA_DIR / 'l2c_interface'
    l2c_dir.mkdir(parents=True, exist_ok=True)
    
    bulk_temporal_data = {}
    for tp_name, tp_data in sig_degs.items():
        bulk_temporal_data[tp_name] = {
            'total_degs': int(len(tp_data)),
            'upregulated': int(len(tp_data[tp_data['significance'] == 'Upregulated'])),
            'downregulated': int(len(tp_data[tp_data['significance'] == 'Downregulated'])),
            'top_genes': tp_data.head(20)['Gene_Symbol'].tolist()
        }
    
    cupro_temporal_patterns = {}
    for gene, tp_data in copper_results.items():
        cupro_temporal_patterns[gene] = {
            'logFC_trajectory': {tp: vals['logFC'] for tp, vals in tp_data.items()},
            'significance': {tp: vals['adj_pval'] for tp, vals in tp_data.items()}
        }
    
    l2c_temporal_config = {
        'bulk_dataset': 'GSE23160',
        'platform': 'Illumina RatRef-12 v1.0',
        'species': 'Rat',
        'time_points': ['2h', '8h', '24h'],
        'temporal_deg_summary': bulk_temporal_data,
        'cuproptosis_temporal_patterns': cupro_temporal_patterns,
        'l2c_integration': {
            'scRNA_cell_types': ['Homeostatic', 'M2', 'M1', 'DAM'],
            'bulk_time_to_cell_mapping': {
                '2h': 'M2',
                '8h': 'M1',
                '24h': 'DAM'
            },
            'calibration_method': 'Use GSE23160 temporal logFC as L2c velocity field initial condition',
            'cuproptosis_genes_for_l2c': list(copper_results.keys())
        },
        'notes': [
            'Pseudotime from scRNA represents state ordering only',
            'Bulk time points (2h, 8h, 24h) are actual reperfusion times',
            'L2c should use Bulk temporal data for mathematical initialization'
        ]
    }
    
    with open(l2c_dir / 'l2c_temporal_config.json', 'w', encoding='utf-8') as f:
        json.dump(l2c_temporal_config, f, indent=2, ensure_ascii=False)
    print(f"  保存: l2c_temporal_config.json")
    
    return l2c_temporal_config

def generate_integration_report(copper_results, sig_degs, l2c_config):
    """生成整合分析报告"""
    print("\n" + "="*60)
    print("生成整合分析报告")
    print("="*60)
    
    report_path = SCRNA_DIR / 'GSE174574_GSE23160_Integration_Report.md'
    
    report = f"""# GSE174574 + GSE23160 整合分析报告

## 📊 数据集概览

| 数据集 | 类型 | 平台 | 物种 | 样本/时间 |
|--------|------|------|------|-----------|
| **GSE174574** | scRNA-seq | GPL21103 (10x) | Mouse | 6 (3 sham + 3 MCAO) |
| **GSE23160** | Bulk RNA-seq | Illumina RatRef-12 | Rat | 3 time points (2h, 8h, 24h) |

## 🔬 铜死亡基因时序表达

### GSE23160 Bulk数据中的铜死亡基因

| 基因 | 2h logFC | 8h logFC | 24h logFC | 模式 |
|------|----------|----------|-----------|------|
"""
    
    for gene, tp_data in copper_results.items():
        logfcs = []
        for tp in ['2h', '8h', '24h']:
            if tp in tp_data:
                logfcs.append(f"{tp_data[tp]['logFC']:.3f}")
            else:
                logfcs.append("N/A")
        report += f"| {gene} | {logfcs[0]} | {logfcs[1]} | {logfcs[2]} | {'→'.join(logfcs)} |\n"
    
    report += f"""
## 📈 显著差异基因统计

| 时间点 | 总DEGs | 上调 | 下调 |
|--------|--------|------|------|
"""
    
    for tp_name, tp_data in sig_degs.items():
        up = len(tp_data[tp_data['significance'] == 'Upregulated'])
        down = len(tp_data[tp_data['significance'] == 'Downregulated'])
        report += f"| {tp_name} | {len(tp_data)} | {up} | {down} |\n"
    
    report += f"""
## 🔄 L2c时间校准策略

### 映射关系

| Bulk时间 | 对应细胞状态 | 生物学意义 |
|----------|-------------|-----------|
| 2h | M2 (抗炎) | 早期炎症反应 |
| 8h | M1 (促炎) | 炎症峰值 |
| 24h | DAM (疾病相关) | 慢性损伤 |

### L2c NeuralODE配置

- **时间校准源**: GSE23160 Bulk时序数据
- **校准方法**: 使用Bulk temporal logFC作为L2c速度场初始条件
- **细胞状态路径**: Homeostatic → M2 → M1 → DAM
- **伪时间用途**: 仅用于可视化，不作为时间变量

## 📁 输出文件

| 文件 | 说明 |
|------|------|
| `l2c_interface/l2c_temporal_config.json` | L2c时间校准配置 |
| `l2c_interface/l2c_config.json` | L2c基础配置 |
| `figures/Figure_2_Cuproptosis_Temporal.png` | 铜死亡基因时序表达 |
| `figures/Figure_3_DEG_Counts.png` | 差异基因统计 |

## ⚠️ 重要说明

1. **GSE23160使用Rat模型**，GSE174574使用Mouse模型
2. 基因名需要物种间映射 (Rat ↔ Mouse)
3. 伪时间仅表示状态排序，不代表实际再灌注时间
4. L2c应使用Bulk数据的时间信息，而非scRNA的伪时间

---

**报告生成日期**: 2026-05-22
"""
    
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write(report)
    print(f"  保存: GSE174574_GSE23160_Integration_Report.md")
    
    return report_path

if __name__ == '__main__':
    print("="*60)
    print("GSE23160 + GSE174574 整合分析")
    print("="*60)
    
    copper_deg, sig_degs = load_bulk_data()
    
    copper_results = analyze_cuproptosis_temporal(copper_deg, sig_degs)
    
    fig_path = create_temporal_figures(copper_results, sig_degs)
    
    l2c_config = prepare_l2c_temporal_calibration(copper_results, sig_degs)
    
    report_path = generate_integration_report(copper_results, sig_degs, l2c_config)
    
    print("\n" + "="*60)
    print("✅ 整合分析完成!")
    print("="*60)
    print(f"\n输出文件:")
    print(f"  - {SCRNA_DIR / 'l2c_interface' / 'l2c_temporal_config.json'}")
    print(f"  - {SCRNA_DIR / 'figures' / 'Figure_2_Cuproptosis_Temporal.png'}")
    print(f"  - {SCRNA_DIR / 'figures' / 'Figure_3_DEG_Counts.png'}")
    print(f"  - {SCRNA_DIR / 'GSE174574_GSE23160_Integration_Report.md'}")
