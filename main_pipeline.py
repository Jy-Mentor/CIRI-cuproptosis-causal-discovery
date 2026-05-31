# -*- coding: utf-8 -*-
"""
CEHG-RNP 2.0 主流水线 (v3.0)
==============================

BCP × Cuproptosis × CIRI 靶点筛选系统

集成改进:
- v3.0: 并行化DAG执行引擎 (pipeline_engine)
- v3.0: 严格错误处理 + 状态机 (StageStatus)
- v3.0: 智能体调度集成 (agent_dispatch)
- v3.0: 内存数据共享 (data_sharing)
- v3.0: MCP工具集成
- v3.0: 输出验证 (StageDataValidator)

使用方法:
    python main_pipeline.py                  # 默认并行+宽松模式
    python main_pipeline.py --sequential     # 顺序模式
    python main_pipeline.py --strict         # 严格模式(失败即停)
    python main_pipeline.py --stage stage4   # 仅运行指定阶段
    python main_pipeline.py --dry-run        # 仅分析依赖关系
"""

import os
import sys
import argparse
import logging
import time
import json

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from config import BASE_DIR, RESULTS_DIR
from scripts.utils import setup_logger, ensure_dir
from scripts.pipeline_engine import (
    PipelineDAG, PipelineExecutor, create_default_dag,
    StageStatus, StageConfig
)
from scripts.agent_dispatch import init_agent_dispatch, PipelineAgentIntegrator
from scripts.data_sharing import get_shared_manager, get_validator

MAIN_DIR = os.path.join(RESULTS_DIR, "main_pipeline")
ensure_dir(MAIN_DIR)

logger = setup_logger("main_v3", os.path.join(MAIN_DIR, "main_pipeline_v3.log"))


def parse_args():
    parser = argparse.ArgumentParser(
        description="CEHG-RNP 2.0 主流水线 (v3.0 - 并行DAG引擎)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  python main_pipeline.py                  # 默认: 并行 + 宽松
  python main_pipeline.py --sequential     # 顺序执行
  python main_pipeline.py --strict         # 严格模式 (失败即停)
  python main_pipeline.py --dry-run        # 分析依赖,不执行
  python main_pipeline.py --stage stage4   # 仅运行stage4
        """
    )
    parser.add_argument("--sequential", action="store_true",
                       help="顺序执行 (禁用并行)")
    parser.add_argument("--strict", action="store_true",
                       help="严格模式 (任何阶段失败立即停止)")
    parser.add_argument("--dry-run", action="store_true",
                       help="干运行 (仅分析依赖关系)")
    parser.add_argument("--stage", type=str,
                       help="仅运行指定阶段 (如 stage4_seed_wgcna)")
    parser.add_argument("--workers", type=int, default=4,
                       help="最大并行工作线程数 (默认: 4)")
    parser.add_argument("--no-agents", action="store_true",
                       help="禁用智能体调度")
    parser.add_argument("--no-sharing", action="store_true",
                       help="禁用内存数据共享")
    return parser.parse_args()


def print_dag_info(dag: PipelineDAG):
    """打印DAG依赖关系分析"""
    levels = dag.get_dependency_order()

    logger.info("")
    logger.info("=" * 70)
    logger.info("DAG 依赖关系分析")
    logger.info("=" * 70)

    for i, level in enumerate(levels):
        logger.info(f"\n层级 {i+1} ({'并行' if len(level)>1 else '单独'}):")
        for name in level:
            config = dag.get_stage_config(name)
            if config:
                deps_str = ", ".join(config.dependencies) if config.dependencies else "无"
                logger.info(f"  >> {name}: {config.description}")
                logger.info(f"    脚本: {config.script}")
                logger.info(f"    依赖: [{deps_str}]")
                logger.info(f"    超时: {config.timeout}s, 重试: {config.max_retries}")

    # 关键路径分析
    level_sizes = [len(l) for l in levels]
    logger.info(f"\n层级数: {len(levels)}")
    logger.info(f"最大并行度: {max(level_sizes) if level_sizes else 1}")
    logger.info(f"无环验证: {'[OK] 通过' if dag.validate_no_cycles() else '[FAIL] 有环!'}")


def run_single_stage(dag: PipelineDAG, stage_name: str, args):
    """运行单个指定阶段"""
    config = dag.get_stage_config(stage_name)
    if config is None:
        logger.error(f"阶段不存在: {stage_name}")
        logger.info(f"可用阶段: {list(dag.stages.keys())}")
        return

    logger.info(f"运行单个阶段: {stage_name} - {config.description}")

    # 检查依赖是否满足
    for dep in config.dependencies:
        dep_dir = os.path.join(RESULTS_DIR, dep)
        if not os.path.exists(dep_dir):
            logger.warning(f"上游阶段 {dep} 输出目录不存在，可能缺少必要数据")

    executor = PipelineExecutor(
        dag=dag,
        base_dir=BASE_DIR,
        results_dir=RESULTS_DIR,
        logger=logger,
        mode="sequential",
        strict=args.strict,
        max_workers=1
    )
    executor._execute_single(stage_name)


def main():
    args = parse_args()

    logger.info("=" * 70)
    logger.info("CEHG-RNP 2.0 主流水线 (v3.0)")
    logger.info("BCP × Cuproptosis × CIRI 靶点筛选系统")
    logger.info("=" * 70)
    logger.info(f"模式: {'顺序' if args.sequential else '并行 (DAG)'}")
    logger.info(f"策略: {'严格' if args.strict else '宽松'}")
    logger.info(f"工作线程: {args.workers}")
    logger.info(f"智能体: {'启用' if not args.no_agents else '禁用'}")
    logger.info(f"数据共享: {'启用' if not args.no_sharing else '禁用'}")

    # ---- 1. 构建DAG ----
    dag = create_default_dag()

    # ---- 2. 分析依赖 ----
    print_dag_info(dag)

    if args.dry_run:
        logger.info("\n干运行完成，未执行任何阶段。")
        return

    # ---- 3. 单独阶段模式 ----
    if args.stage:
        run_single_stage(dag, args.stage, args)
        return

    # ---- 4. 初始化共享数据管理器 ----
    if not args.no_sharing:
        shared_mgr = get_shared_manager(RESULTS_DIR)
        logger.info(f"共享数据管理器: 已启用")
    else:
        shared_mgr = None

    # ---- 5. 初始化智能体调度 ----
    if not args.no_agents:
        agent_integrator = init_agent_dispatch(RESULTS_DIR)
    else:
        agent_integrator = PipelineAgentIntegrator()

    # ---- 6. 初始化输出验证器 ----
    validator = get_validator(RESULTS_DIR)

    # ---- 7. 创建执行器并运行 ----
    mode = "sequential" if args.sequential else "parallel"

    executor = PipelineExecutor(
        dag=dag,
        base_dir=BASE_DIR,
        results_dir=RESULTS_DIR,
        logger=logger,
        mode=mode,
        strict=args.strict,
        max_workers=args.workers
    )

    # 注册智能体回调
    if not args.no_agents:
        original_worker = executor._execute_single_worker

        def worker_with_agent(name):
            original_worker(name)
            state = dag.stages[name]
            result = executor.results.get(name)
            success = result is not None and result.is_success
            agent_integrator.on_stage_complete(name, success)

        executor._execute_single_worker = worker_with_agent

    # ---- 8. 执行管道 ----
    results = executor.run()

    # ---- 9. 验证输出 ----
    logger.info("\n" + "=" * 70)
    logger.info("阶段输出验证")
    logger.info("=" * 70)

    required_files_map = {
        "stage1_rma_degs": ["limma_degs.csv", "sample_annotations.csv"],
        "stage2_single_cell": ["sc_adata.h5ad", "cell_annotations.csv"],
        "stage3_enrichment": ["go_enrichment.csv", "kegg_enrichment.csv"],
        "stage4_seed_wgcna": ["wgcna_modules.csv", "seed_pool_genes.txt"],
        "stage5_ppi_mcode": ["ppi_topology.json", "node_degree_ranking.csv"],
        "stage6_sctenifold_knockout": ["gene_perturbation_scores.csv"],
        "stage7_ml_shap": ["gene_shap_importance.csv", "ml_model_performance.csv"],
        "stage8_final_targets": ["core_targets.csv", "tier1_targets.csv", "final_report.txt"],
        "brain_coexpression": ["brain_coexpression_features.csv", "feature_dimensions.json"],
    }

    for stage_name, required_files in required_files_map.items():
        is_valid, errors = validator.validate_stage_output(stage_name, required_files)
        status = "[OK]" if is_valid else "[FAIL]"
        logger.info(f"  {status} {stage_name}: {'通过' if is_valid else f'失败 ({errors})'}")

    # ---- 10. 智能体调度完成 ----
    success_count = sum(1 for r in results.values() if r.is_success)
    agent_integrator.on_pipeline_complete(success_count, len(dag.stages))

    agent_summary = agent_integrator.get_integration_summary()
    if agent_summary:
        logger.info(f"\n智能体调用摘要: {len(agent_summary)} 次调用")
        for call in agent_summary:
            logger.info(f"  - {call.get('agent', 'N/A')}: {str(call.get('query', ''))[:80]}")

    # ---- 11. 共享数据管理器统计 ----
    if shared_mgr:
        stats = shared_mgr.get_stats()
        logger.info(f"\n共享数据管理器统计:")
        logger.info(f"  缓存条目: {stats['cache_size']}")
        logger.info(f"  缓存命中率: {stats['cache_hit_rate']:.1%}")
        logger.info(f"  内存映射: {stats['memmap_count']} ({stats['memmap_total_mb']:.1f}MB)")

    # ---- 12. 最终检查 ----
    final_targets = os.path.join(RESULTS_DIR, "stage8_final_targets", "core_targets.csv")
    if os.path.exists(final_targets):
        import pandas as pd
        try:
            targets = pd.read_csv(final_targets)
            logger.info(f"\n[TARGET] 最终核心靶点 ({len(targets)} 个):")
            for _, row in targets.head(10).iterrows():
                score = row.get('Comprehensive', row.get('Score', 'N/A'))
                gene = row.get('Gene', 'N/A')
                tier = row.get('Tier', '')
                logger.info(f"  {gene}: {score} [{tier}]")
        except Exception as e:
            logger.warning(f"读取core_targets失败: {e}")

    logger.info("")
    logger.info("=" * 70)
    logger.info("流水线执行完成!")
    logger.info("=" * 70)


if __name__ == "__main__":
    main()
