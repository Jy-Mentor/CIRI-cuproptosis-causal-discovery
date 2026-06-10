"""
主流程编排模块
===============
负责组织 ML 管道执行流程，协调各子模块完成：
  1. GPU 环境检测
  2. 数据加载与预处理
  3. 构建特征工程×分类器笛卡尔积组合
  4. 药物可靶向性预测 (DT Task)
  5. 疾病相关性预测 (DG Task)
  6. 桥梁得分计算 (P_drug × P_disease)
  7. RRF 倒数排名融合
  8. 结果输出与 GAT 对比
"""

import os
import time
import numpy as np
import pandas as pd

from sklearn.preprocessing import StandardScaler

from .config import PipelineConfig
from .utils import log, print_header, compute_sample_weight
from .gpu_detector import detect_gpu_environment
from .data_loader import DataLoader
from .feature_engineering import FeatureEngineeringRegistry, RawFeatureStrategy
from .classifiers import ClassifierRegistry
from .trainer import ModelTrainer
from .ensemble import EnsembleIntegrator


def main():
    """主入口函数 — 执行完整的多算法桥梁靶点预测流程"""
    start_time = time.time()

    # ══════════════════════════════════════════════════════════════════════
    # 初始化
    # ══════════════════════════════════════════════════════════════════════
    config = PipelineConfig.create_default()
    config.ensure_output_dir()

    log("═" * 70)
    log("  多算法融合桥梁靶点预测系统 v3.1 (模块化重构)")
    log("  特征工程 × 基础分类器 = 论文级笛卡尔积组合")
    log(f"  随机种子: {config.training.seed}, {config.training.n_folds}折CV, "
        f"RRF_k={config.training.rrf_k}")
    log(f"  并行核心: {config.training.n_jobs}, GPU加速: {config.gpu.enabled}")
    log("═" * 70)

    # ── GPU 检测 ──
    config.gpu = detect_gpu_environment(
        gpu_enabled=config.gpu.enabled,
        n_jobs=config.training.n_jobs,
    )

    # ── 数据加载 ──
    data_loader = DataLoader(config.paths)
    data = data_loader.load_data()

    n_total = len(data.gene_index)
    n_unknown = data.unknown_mask.sum()
    n_dt = data.y_dt.sum()
    n_dg = data.y_dg.sum()
    log(f"  总基因: {n_total}, DT+={n_dt}, DG+={n_dg}, Unknown={n_unknown}")

    # ── 初始化各模块 ──
    fe_registry = FeatureEngineeringRegistry.create_default()
    clf_registry = ClassifierRegistry.create_default()
    trainer = ModelTrainer(config)
    integrator = EnsembleIntegrator(config.training)

    # ══════════════════════════════════════════════════════════════════════
    # STEP 2: 构建模型组合
    # ══════════════════════════════════════════════════════════════════════
    print_header("STEP 2: 构建 特征工程 × 分类器 笛卡尔积组合")
    model_combos, fe_map, clf_map = trainer.build_model_combinations(fe_registry, clf_registry)
    model_names = [m[0] for m in model_combos]
    log(f"  共 {len(model_combos)} 种模型组合")

    # ══════════════════════════════════════════════════════════════════════
    # STEP 3: DT Task (并行)
    # ══════════════════════════════════════════════════════════════════════
    print_header("STEP 3: 药物可靶向性预测 (DT Task) — 并行")
    drug_probas, dt_metrics = trainer.run_cv_parallel(
        model_combos, data.X, data.y_dt, data.unknown_mask, "DT",
    )

    # ══════════════════════════════════════════════════════════════════════
    # STEP 4: DG Task (并行)
    # ══════════════════════════════════════════════════════════════════════
    print_header("STEP 4: 疾病相关性预测 (DG Task) — 并行")
    disease_probas, dg_metrics = trainer.run_cv_parallel(
        model_combos, data.X, data.y_dg, data.unknown_mask, "DG",
    )

    # ══════════════════════════════════════════════════════════════════════
    # STEP 5: 桥梁得分
    # ══════════════════════════════════════════════════════════════════════
    print_header("STEP 5: 计算桥梁得分 (P_drug × P_disease)")
    bridge_scores = integrator.compute_bridge_scores(drug_probas, disease_probas)
    log(f"  有效组合: {len(bridge_scores)}/{len(model_combos)}")

    # ══════════════════════════════════════════════════════════════════════
    # STEP 6: RRF 集成
    # ══════════════════════════════════════════════════════════════════════
    print_header("STEP 6: RRF 倒数排名融合")
    rrf_scores, rrf_ranks = integrator.rrf_integrate(
        bridge_scores, n_total, data.unknown_mask,
    )

    # ══════════════════════════════════════════════════════════════════════
    # STEP 7: 生成输出文件
    # ══════════════════════════════════════════════════════════════════════
    print_header("STEP 7: 生成输出文件")

    _save_results(
        config=config,
        data=data,
        model_names=model_names,
        drug_probas=drug_probas,
        disease_probas=disease_probas,
        bridge_scores=bridge_scores,
        rrf_scores=rrf_scores,
        rrf_ranks=rrf_ranks,
        dt_metrics=dt_metrics,
        dg_metrics=dg_metrics,
        integrator=integrator,
        fe_map=fe_map,
        clf_map=clf_map,
    )

    # ══════════════════════════════════════════════════════════════════════
    # STEP 8: GAT 对比
    # ══════════════════════════════════════════════════════════════════════
    print_header("STEP 8: GAT 桥梁基因对比")

    df_unknown = _build_unknown_df(
        data.gene_index, data.unknown_mask,
        model_names, drug_probas, disease_probas, bridge_scores,
        rrf_scores, rrf_ranks,
    )
    integrator.compare_with_gat(df_unknown, config.paths.gat_bridge_path)

    # ══════════════════════════════════════════════════════════════════════
    # 完成
    # ══════════════════════════════════════════════════════════════════════
    elapsed = time.time() - start_time
    log("\n" + "=" * 70)
    log(f"[OK] 完成! 总耗时: {elapsed / 60:.1f} 分钟")
    log(f"  输出目录: {config.paths.output_dir}")
    log(f"  模型组合数: {len(model_combos)}")
    log(f"  每组合 {config.training.n_folds}折CV × 2任务 = "
        f"{len(model_combos) * config.training.n_folds * 2} 次训练")
    log("=" * 70)


# ═══════════════════════════════════════════════════════════════════════════════
# 内部辅助函数
# ═══════════════════════════════════════════════════════════════════════════════

def _build_unknown_df(gene_index, unknown_mask, model_names,
                      drug_probas, disease_probas, bridge_scores,
                      rrf_scores, rrf_ranks):
    """构建包含所有未知基因完整信息的 DataFrame"""
    df_all = pd.DataFrame({'gene_symbol': gene_index})
    for model in model_names:
        if model in drug_probas:
            df_all[f'{model}_drug_prob'] = drug_probas[model]
        if model in disease_probas:
            df_all[f'{model}_disease_prob'] = disease_probas[model]
        if model in bridge_scores:
            df_all[f'{model}_bridge_score'] = bridge_scores[model]
    df_all['RRF_score'] = rrf_scores
    df_all['final_rank'] = rrf_ranks
    return df_all[unknown_mask].sort_values('RRF_score', ascending=False)


def _save_results(config, data, model_names,
                  drug_probas, disease_probas, bridge_scores,
                  rrf_scores, rrf_ranks,
                  dt_metrics, dg_metrics,
                  integrator, fe_map, clf_map):
    """保存所有输出文件"""
    output_dir = config.paths.output_dir
    df_unknown = _build_unknown_df(
        data.gene_index, data.unknown_mask,
        model_names, drug_probas, disease_probas, bridge_scores,
        rrf_scores, rrf_ranks,
    )

    # 7a. 所有未知基因详情
    unknown_path = os.path.join(output_dir, "ml_bridge_genes_all.csv")
    df_unknown.to_csv(unknown_path, index=False, encoding='utf-8-sig')
    log(f"  所有未知基因: {unknown_path} ({len(df_unknown)} genes)")

    # 7b. Top-20
    top20 = df_unknown.head(20)
    top20_path = os.path.join(output_dir, "top20_bridge_genes_ml.csv")
    top20.to_csv(top20_path, index=False, encoding='utf-8-sig')
    log(f"  Top-20: {top20_path}")

    log("\n" + "=" * 70)
    log("TOP-10 桥梁基因 (RRF 融合)")
    log("=" * 70)
    for i, (_, row) in enumerate(top20.head(10).iterrows()):
        log(f"  {int(row['final_rank']):>4d}. {row['gene_symbol']:<12s}  "
            f"RRF={row['RRF_score']:.6f}")

    # 7c. 算法指标
    df_metrics = integrator.build_metrics_df(model_names, dt_metrics, dg_metrics)
    metrics_path = os.path.join(output_dir, "algorithm_metrics.csv")
    df_metrics.to_csv(metrics_path, index=False, encoding='utf-8-sig')
    log(f"  算法指标: {metrics_path}")

    # 7d. 各特征工程策略最佳模型
    log("\n--- 各特征工程策略的 DT 最佳模型 ---")
    fe_names = sorted(fe_map.keys())
    metrics_rows = df_metrics.to_dict('records') if not df_metrics.empty else []
    for fe in fe_names:
        fe_models = [
            r for r in metrics_rows
            if r['model'].startswith(fe + "__") and r['task'] == 'DT'
        ]
        if fe_models:
            best = max(fe_models, key=lambda x: x['auroc'])
            log(f"  {fe}: {best['model']} AUROC={best['auroc']:.4f}")

    # 7e. 特征重要性提取
    _extract_and_save_importance(
        config, data, model_names, dt_metrics, dg_metrics,
        metrics_rows, fe_map, clf_map,
    )


def _extract_and_save_importance(config, data, model_names,
                                 dt_metrics, dg_metrics,
                                 metrics_rows, fe_map, clf_map):
    """提取最佳模型的特征重要性"""
    from ml_bridge_target_pipeline.feature_engineering import (
        RawFeatureStrategy, SUPERVISED_FE_TYPES,
    )
    from ml_bridge_target_pipeline.classifiers import (
        GBStrategy, NBStrategy,  # requires_sample_weight
    )

    log("\n--- 特征重要性 (raw 特征空间, 最佳模型) ---")
    importance_dfs = []

    if not metrics_rows:
        return

    dt_best_row = max([r for r in metrics_rows if r['task'] == 'DT'],
                      key=lambda x: x['auroc'])
    dg_best_row = max([r for r in metrics_rows if r['task'] == 'DG'],
                      key=lambda x: x['auroc'])
    best_models = {'DT': dt_best_row['model'], 'DG': dg_best_row['model']}
    log(f"  DT最佳: {best_models['DT']} | DG最佳: {best_models['DG']}")

    from ml_bridge_target_pipeline.ensemble import EnsembleIntegrator

    for task_name, best_model_name in best_models.items():
        try:
            fe_name, clf_name = best_model_name.split("__", 1)
            y_task = data.y_dt if task_name == 'DT' else data.y_dg

            if len(np.unique(y_task)) < 2:
                continue

            # 标准化
            scaler = StandardScaler()
            X_s = scaler.fit_transform(data.X)

            fe_strategy = fe_map.get(fe_name)
            clf_strategy = clf_map.get(clf_name)
            if fe_strategy is None or clf_strategy is None:
                log(f"    [WARN] {best_model_name}: 未找到 strategy")
                continue

            # 特征工程
            fe_obj = fe_strategy.create_transformer(config.training.seed)
            if fe_obj is not None:
                if fe_strategy.is_supervised:
                    fe_obj.fit(X_s, y_task)
                else:
                    fe_obj.fit(X_s)
                X_fe = fe_obj.transform(X_s)
                if isinstance(X_fe, tuple):
                    X_fe = X_fe[0]
            else:
                X_fe = X_s

            # 训练分类器
            clf = clf_strategy.create_classifier(config.training.seed, config.gpu)
            if clf_strategy.requires_sample_weight:
                clf.fit(X_fe, y_task, sample_weight=compute_sample_weight(y_task))
            else:
                clf.fit(X_fe, y_task)

            # 提取重要性
            if isinstance(fe_strategy, RawFeatureStrategy):
                imp_df = EnsembleIntegrator.extract_importance(
                    best_model_name, clf, data.feature_cols,
                )
                if imp_df is not None:
                    imp_df['task'] = task_name
                    importance_dfs.append(imp_df)
                    log(f"    {best_model_name}: Top-5 = {list(imp_df['feature'].head(5))}")
            else:
                log(f"    {best_model_name}: 特征变换空间({fe_name}), 重要性跳过")

        except Exception as e:
            log(f"    [WARN] {best_model_name} 重要性提取失败: {e}")

    if importance_dfs:
        imp_all = pd.concat(importance_dfs, ignore_index=True)
        imp_path = os.path.join(config.paths.output_dir, "feature_importance.csv")
        imp_all.to_csv(imp_path, index=False, encoding='utf-8-sig')
        log(f"  特征重要性: {imp_path}")


if __name__ == "__main__":
    main()