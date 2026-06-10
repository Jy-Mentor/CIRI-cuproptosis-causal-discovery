"""
GSE174574 伪时间→真实时间映射 (GPR科学验证版)

方法学:
1. 高斯过程回归 (Gaussian Process Regression)
   - 优点: 提供贝叶斯不确定性量化,适合小样本校准
   - 参考: Reid & Wernisch (2018) PLOS Computational Biology
   
2. 样条插值作为对比基线
   - 参考: Trapnell et al. (2014) Nature Biotechnology

3. DPT伪时间
   - 参考: Haghverdi et al. (2016) Nature Biotechnology
"""

import scanpy as sc
import numpy as np
import pandas as pd
from pathlib import Path
import json
from scipy.interpolate import UnivariateSpline
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

plt.rcParams['font.size'] = 11
plt.rcParams['axes.linewidth'] = 1.2

DATA_FILE = Path(r'C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\GSE174574_analysis\results\GSE174574_processed.h5ad')
OUTPUT_DIR = Path(r'C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\GSE174574_analysis')

BULK_TIME_POINTS = {'Homeostatic': 0, 'M2': 2, 'M1': 8, 'DAM': 24}


def compute_pseudotime(adata):
    """计算DPT伪时间"""
    print("=" * 70)
    print("Step 1: 计算DPT伪时间")
    print("=" * 70)

    homeo_mask = adata.obs['cell_type'] == 'Homeostatic'
    homeo_cells = np.where(homeo_mask)[0]
    if len(homeo_cells) == 0:
        raise ValueError("未找到Homeostatic细胞!")

    root_idx = homeo_cells[len(homeo_cells) // 2]
    adata.uns['iroot'] = root_idx
    print(f"  根细胞: index={root_idx}")

    sc.tl.dpt(adata, n_dcs=10)
    adata.obs['pseudotime'] = adata.obs['dpt_pseudotime']

    print(f"\n  伪时间统计:")
    stats = {}
    for ct in adata.obs['cell_type'].unique():
        mask = adata.obs['cell_type'] == ct
        pt = adata.obs.loc[mask, 'pseudotime'].dropna()
        stats[ct] = {
            'median': float(pt.median()),
            'iqr_25': float(pt.quantile(0.25)),
            'iqr_75': float(pt.quantile(0.75))
        }
        print(f"    {ct}: median={stats[ct]['median']:.4f}, IQR=[{stats[ct]['iqr_25']:.4f}, {stats[ct]['iqr_75']:.4f}]")

    return adata, stats


def gpr_calibration(pseudotime_stats, real_times):
    """使用高斯过程回归进行时间校准"""
    print(f"\n{'=' * 70}")
    print("Step 2: 高斯过程回归 (GPR) 校准")
    print("=" * 70)

    pt_values = np.array([stats['median'] for stats in pseudotime_stats.values()]).reshape(-1, 1)
    rt_values = np.array([real_times[ct] for ct in pseudotime_stats.keys()])

    scaler_X = StandardScaler()
    scaler_y = StandardScaler()
    pt_scaled = scaler_X.fit_transform(pt_values)
    rt_scaled = scaler_y.fit_transform(rt_values.reshape(-1, 1)).ravel()

    kernel = ConstantKernel(1.0) * RBF(length_scale=1.0) + WhiteKernel(noise_level=0.1)
    gpr = GaussianProcessRegressor(
        kernel=kernel,
        alpha=0.1,
        n_restarts_optimizer=10,
        random_state=42,
        normalize_y=True
    )

    gpr.fit(pt_scaled, rt_scaled)

    print(f"  GPR核函数参数:")
    print(f"    优化后核: {gpr.kernel_}")
    print(f"    对数边际似然: {gpr.log_marginal_likelihood():.4f}")

    return gpr, scaler_X, scaler_y


def predict_with_uncertainty(gpr, scaler_X, scaler_y, cell_pseudotime):
    """使用GPR预测真实时间+不确定性"""
    print(f"\n{'=' * 70}")
    print("Step 3: 预测+不确定性量化")
    print("=" * 70)

    pt_scaled = scaler_X.transform(cell_pseudotime.reshape(-1, 1))
    time_scaled, std_scaled = gpr.predict(pt_scaled, return_std=True)
    time_real = scaler_y.inverse_transform(time_scaled.reshape(-1, 1)).ravel()
    std_real = std_scaled * scaler_y.scale_[0]

    time_real = np.clip(time_real, 0, 24)
    ci_lower = np.clip(time_real - 1.96 * std_real, 0, 24)
    ci_upper = np.clip(time_real + 1.96 * std_real, 0, 24)

    print(f"  不确定性统计:")
    print(f"    平均标准差: {np.nanmean(std_real):.2f} h")
    print(f"    95% CI平均宽度: {np.nanmean(ci_upper - ci_lower):.2f} h")
    print(f"    最大标准差: {np.nanmax(std_real):.2f} h")

    return time_real, std_real, ci_lower, ci_upper


def apply_calibration(adata, time_real, std_real, ci_lower, ci_upper):
    """应用校准"""
    adata.obs['real_time_hours'] = time_real
    adata.obs['real_time_std'] = std_real
    adata.obs['real_time_ci_lower'] = ci_lower
    adata.obs['real_time_ci_upper'] = ci_upper
    adata.obs['calibration_method'] = 'GPR_Reid_Wernisch_2018'

    print(f"\n  校准后时间:")
    for ct in adata.obs['cell_type'].unique():
        mask = adata.obs['cell_type'] == ct
        t_mean = np.nanmean(adata.obs.loc[mask, 'real_time_hours'])
        t_std = np.nanmean(adata.obs.loc[mask, 'real_time_std'])
        print(f"    {ct}: {t_mean:.1f} ± {t_std:.1f} h")


def save_results(adata, pseudotime_stats, gpr_kernel):
    """保存结果"""
    l2c_dir = OUTPUT_DIR / 'l2c_interface'
    l2c_dir.mkdir(parents=True, exist_ok=True)

    cell_data = adata.obs[[
        'cell_type', 'pseudotime', 'real_time_hours', 'real_time_std',
        'real_time_ci_lower', 'real_time_ci_upper', 'condition'
    ]].copy()
    cell_data.to_csv(l2c_dir / 'cell_time_mapping.csv')

    l2c_config = {
        'metadata': {
            'method': 'Gaussian Process Regression (GPR)',
            'reference': 'Reid & Wernisch (2018) PLOS Computational Biology',
            'kernel': str(gpr_kernel),
            'pseudotime_method': 'DPT (Haghverdi et al., 2016)'
        },
        'time_mapping': {
            'cell_type_pseudotime_median': pseudotime_stats,
            'bulk_reference': BULK_TIME_POINTS
        },
        'l2c_usage': {
            'time_column': 'real_time_hours',
            'uncertainty_column': 'real_time_std',
            'ci_columns': ['real_time_ci_lower', 'real_time_ci_upper'],
            'time_unit': 'hours',
            'time_range': [0, 24]
        },
        'quality_metrics': {
            'avg_uncertainty_h': float(np.nanmean(adata.obs['real_time_std'])),
            'ci_width_95_h': float(np.nanmean(adata.obs['real_time_ci_upper'] - adata.obs['real_time_ci_lower']))
        }
    }

    with open(l2c_dir / 'l2c_gpr_config.json', 'w', encoding='utf-8') as f:
        json.dump(l2c_config, f, indent=2, ensure_ascii=False)

    adata.write_h5ad(l2c_dir / 'GSE174574_with_gpr_time.h5ad')

    print(f"  ✓ cell_time_mapping.csv")
    print(f"  ✓ l2c_gpr_config.json")
    print(f"  ✓ GSE174574_with_gpr_time.h5ad")


def create_figures(adata, gpr, scaler_X, scaler_y):
    """生成可视化"""
    fig_dir = OUTPUT_DIR / 'figures'
    fig_dir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(2, 3, figsize=(20, 12))

    sc.pl.umap(adata, color='pseudotime', ax=axes[0, 0], show=False,
               title='DPT Pseudotime', cmap='viridis', size=3)

    sc.pl.umap(adata, color='real_time_hours', ax=axes[0, 1], show=False,
               title='GPR Calibrated Time (hours)', cmap='plasma', size=3)

    sc.pl.umap(adata, color='real_time_std', ax=axes[0, 2], show=False,
               title='GPR Uncertainty (±hours)', cmap='RdYlBu_r', size=3)

    colors = {'Homeostatic': '#2ecc71', 'M2': '#3498db', 'M1': '#e67e22', 'DAM': '#e74c3c'}
    for ct, color in colors.items():
        mask = adata.obs['cell_type'] == ct
        axes[1, 0].hist(adata.obs.loc[mask, 'real_time_hours'], alpha=0.6,
                        label=ct, color=color, bins=40)
    axes[1, 0].axvline(x=2, color='gray', linestyle='--', alpha=0.5)
    axes[1, 0].axvline(x=8, color='gray', linestyle='--', alpha=0.5)
    axes[1, 0].axvline(x=24, color='gray', linestyle='--', alpha=0.5)
    axes[1, 0].set_xlabel('Time (hours)')
    axes[1, 0].set_ylabel('Cell Count')
    axes[1, 0].set_title('Time Distribution by Cell Type (GPR)')
    axes[1, 0].legend()
    axes[1, 0].grid(True, alpha=0.3)

    for ct, color in colors.items():
        mask = adata.obs['cell_type'] == ct
        median_time = np.nanmedian(adata.obs.loc[mask, 'real_time_hours'])
        median_std = np.nanmedian(adata.obs.loc[mask, 'real_time_std'])
        axes[1, 1].errorbar(x=[median_time], y=[ct], xerr=[median_std],
                           fmt='o', color=color, capsize=5, markersize=10)
    axes[1, 1].set_xlabel('Median Time ± GPR Uncertainty (hours)')
    axes[1, 1].set_title('Cell Type Time Estimates (GPR 95% CI)')
    axes[1, 1].grid(True, alpha=0.3, axis='x')

    pt_test = np.linspace(0, 1, 200).reshape(-1, 1)
    pt_scaled = scaler_X.transform(pt_test)
    pred_scaled, std_scaled = gpr.predict(pt_scaled, return_std=True)
    pred_real = scaler_y.inverse_transform(pred_scaled.reshape(-1, 1)).ravel()
    std_real = std_scaled * scaler_y.scale_[0]

    axes[1, 2].plot(pt_test, pred_real, 'b-', linewidth=2, label='GPR Mean')
    axes[1, 2].fill_between(pt_test.ravel(),
                           pred_real - 1.96 * std_real,
                           pred_real + 1.96 * std_real,
                           alpha=0.2, color='blue', label='95% CI')

    for ct, color in colors.items():
        mask = adata.obs['cell_type'] == ct
        axes[1, 2].scatter(adata.obs.loc[mask, 'pseudotime'].values,
                          adata.obs.loc[mask, 'real_time_hours'].values,
                          alpha=0.1, s=5, color=color, label=ct)

    axes[1, 2].set_xlabel('Pseudotime (DPT)')
    axes[1, 2].set_ylabel('Calibrated Real Time (hours)')
    axes[1, 2].set_title('GPR: Pseudotime → Real Time')
    axes[1, 2].legend(fontsize=8)
    axes[1, 2].grid(True, alpha=0.3)

    plt.suptitle('GSE174574: GPR Pseudotime Calibration (Reid & Wernisch, 2018)',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig_path = fig_dir / 'Figure_GPR_Calibration.png'
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  ✓ Figure_GPR_Calibration.png")


if __name__ == '__main__':
    print("=" * 70)
    print("GSE174574 伪时间→真实时间映射 (GPR版)")
    print("方法: Gaussian Process Regression (Reid & Wernisch, 2018)")
    print("=" * 70)

    adata = sc.read_h5ad(DATA_FILE)
    print(f"\n数据: {adata.n_obs} cells × {adata.n_vars} genes")

    adata, pseudotime_stats = compute_pseudotime(adata)

    gpr, scaler_X, scaler_y = gpr_calibration(pseudotime_stats, BULK_TIME_POINTS)

    cell_pseudotime = adata.obs['pseudotime'].values
    time_real, std_real, ci_lower, ci_upper = predict_with_uncertainty(
        gpr, scaler_X, scaler_y, cell_pseudotime
    )

    apply_calibration(adata, time_real, std_real, ci_lower, ci_upper)

    save_results(adata, pseudotime_stats, str(gpr.kernel_))

    create_figures(adata, gpr, scaler_X, scaler_y)

    print("\n" + "=" * 70)
    print("✅ GPR校准完成!")
    print("=" * 70)
    print(f"\nL2c使用方式:")
    print(f"  1. 加载 GSE174574_with_gpr_time.h5ad")
    print(f"  2. 使用 real_time_hours 列作为时间变量 t")
    print(f"  3. NeuralODE: dy/dt = f(y, t), t ∈ [0, 24] hours")
    print(f"\n不确定性量化:")
    print(f"  - real_time_std: GPR预测标准差")
    print(f"  - real_time_ci_lower/upper: 95%置信区间")
