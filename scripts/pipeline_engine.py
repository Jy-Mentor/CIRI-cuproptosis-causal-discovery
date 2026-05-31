# -*- coding: utf-8 -*-
"""
CEHG-RNP 2.0 并行化管道执行引擎
=================================

功能:
- DAG依赖图调度 (拓扑排序 + 并行执行)
- 阶段状态机 (PENDING → RUNNING → SUCCESS/FAILED → SKIPPED)
- 严格模式 vs 宽松模式
- 阶段超时控制
- 自动重试机制
- 执行统计与性能分析

参考:
- Airflow DAG调度模型
- concurrent.futures ThreadPoolExecutor/ProcessPoolExecutor

版本: v1.0 | 日期: 2026-05-28
"""

import os
import sys
import time
import json
import logging
import threading
import traceback
from pathlib import Path
from typing import Dict, List, Optional, Callable, Any, Tuple, Set
from enum import Enum
from dataclasses import dataclass, field
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from collections import defaultdict, deque

import numpy as np
import pandas as pd


class StageStatus(Enum):
    PENDING = "PENDING"
    RUNNING = "RUNNING"
    SUCCESS = "SUCCESS"
    FAILED = "FAILED"
    SKIPPED = "SKIPPED"
    TIMEOUT = "TIMEOUT"


@dataclass
class StageConfig:
    name: str
    script: str
    description: str
    dependencies: List[str] = field(default_factory=list)
    timeout: int = 1800
    max_retries: int = 1
    is_r_script: bool = False
    required_outputs: List[str] = field(default_factory=list)
    priority: int = 0


@dataclass
class StageResult:
    stage_name: str
    status: StageStatus
    start_time: float = 0.0
    end_time: float = 0.0
    elapsed: float = 0.0
    return_code: int = -1
    stdout: str = ""
    stderr: str = ""
    error_message: str = ""
    retry_count: int = 0

    @property
    def is_success(self) -> bool:
        return self.status == StageStatus.SUCCESS

    @property
    def is_failed(self) -> bool:
        return self.status in (StageStatus.FAILED, StageStatus.TIMEOUT)


class StageStateMachine:
    """阶段状态机，管理单个阶段的执行状态"""

    def __init__(self, config: StageConfig):
        self.config = config
        self.status = StageStatus.PENDING
        self.result: Optional[StageResult] = None
        self._lock = threading.Lock()

    def can_start(self, completed_stages: Set[str]) -> bool:
        """检查所有依赖是否已完成"""
        return all(dep in completed_stages for dep in self.config.dependencies)

    def transition(self, new_status: StageStatus) -> None:
        """状态转换"""
        with self._lock:
            self.status = new_status

    def record_result(self, result: StageResult) -> None:
        """记录执行结果"""
        with self._lock:
            self.result = result
            self.status = result.status


class PipelineDAG:
    """管道DAG调度器"""

    def __init__(self):
        self.stages: Dict[str, StageStateMachine] = {}
        self._adjacency: Dict[str, List[str]] = defaultdict(list)
        self._reverse_adjacency: Dict[str, List[str]] = defaultdict(list)

    def add_stage(self, config: StageConfig) -> None:
        """添加阶段到DAG"""
        name = config.name
        self.stages[name] = StageStateMachine(config)

        for dep in config.dependencies:
            self._adjacency[dep].append(name)
            self._reverse_adjacency[name].append(dep)

        if name not in self._adjacency:
            self._adjacency[name] = []

    def get_ready_stages(self, completed: Set[str]) -> List[str]:
        """获取所有依赖已满足的准备就绪阶段"""
        ready = []
        for name, state in self.stages.items():
            if state.status == StageStatus.PENDING:
                if state.can_start(completed):
                    ready.append(name)
        return ready

    def get_dependency_order(self) -> List[List[str]]:
        """
        拓扑排序获取执行层级
        返回: [[level_0_stages], [level_1_stages], ...]
        """
        in_degree = {name: len(deps) for name, deps in self._reverse_adjacency.items()}
        for name in self.stages:
            if name not in in_degree:
                in_degree[name] = 0

        queue = deque([name for name, deg in in_degree.items() if deg == 0])
        levels = []
        visited = set()

        while queue:
            level = []
            for _ in range(len(queue)):
                node = queue.popleft()
                if node in visited:
                    continue
                visited.add(node)
                level.append(node)

                for neighbor in self._adjacency[node]:
                    in_degree[neighbor] -= 1
                    if in_degree[neighbor] == 0:
                        queue.append(neighbor)

            if level:
                levels.append(level)

        for name in self.stages:
            if name not in visited:
                levels.append([name])

        return levels

    def validate_no_cycles(self) -> bool:
        """检测DAG是否有环"""
        visited = set()
        rec_stack = set()

        def has_cycle(node: str) -> bool:
            visited.add(node)
            rec_stack.add(node)
            for neighbor in self._adjacency.get(node, []):
                if neighbor not in visited:
                    if has_cycle(neighbor):
                        return True
                elif neighbor in rec_stack:
                    return True
            rec_stack.discard(node)
            return False

        for node in self.stages:
            if node not in visited:
                if has_cycle(node):
                    return False
        return True

    def get_downstream(self, stage_name: str) -> Set[str]:
        """获取某个阶段的所有下游阶段"""
        downstream = set()
        queue = deque([stage_name])
        while queue:
            node = queue.popleft()
            for neighbor in self._adjacency.get(node, []):
                if neighbor not in downstream:
                    downstream.add(neighbor)
                    queue.append(neighbor)
        return downstream

    def get_stage_config(self, name: str) -> Optional[StageConfig]:
        state = self.stages.get(name)
        return state.config if state else None


class PipelineExecutor:
    """管道执行器"""

    def __init__(
        self,
        dag: PipelineDAG,
        base_dir: str,
        results_dir: str,
        logger: logging.Logger,
        mode: str = "parallel",
        strict: bool = True,
        max_workers: int = 4
    ):
        self.dag = dag
        self.base_dir = base_dir
        self.results_dir = results_dir
        self.logger = logger
        self.mode = mode
        self.strict = strict
        self.max_workers = max_workers

        self.results: Dict[str, StageResult] = {}
        self.completed_stages: Set[str] = set()
        self.failed_stages: Set[str] = set()
        self._lock = threading.Lock()

        # 执行统计
        self.total_start_time = 0.0
        self.total_end_time = 0.0

    def _run_subprocess(self, config: StageConfig) -> StageResult:
        """通过子进程运行阶段"""
        import subprocess

        script_path = os.path.join(self.base_dir, "scripts", config.script)

        if not os.path.exists(script_path):
            return StageResult(
                stage_name=config.name,
                status=StageStatus.FAILED,
                error_message=f"脚本不存在: {script_path}"
            )

        start_time = time.time()
        result = StageResult(stage_name=config.name, status=StageStatus.RUNNING, start_time=start_time)

        try:
            if config.is_r_script:
                rscript_path = r"C:\R\R-4.5.2\bin\Rscript.exe"
                if not os.path.exists(rscript_path):
                    rscript_path = "Rscript"
                cmd = [rscript_path, script_path]
            else:
                cmd = [sys.executable, script_path]

            proc = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                encoding='utf-8',
                errors='replace',
                timeout=config.timeout,
                cwd=self.base_dir
            )

            result.return_code = proc.returncode
            result.stdout = proc.stdout[-5000:] if proc.stdout else ""
            result.stderr = proc.stderr[-5000:] if proc.stderr else ""
            result.end_time = time.time()
            result.elapsed = result.end_time - start_time

            if proc.returncode == 0:
                result.status = StageStatus.SUCCESS
            else:
                result.status = StageStatus.FAILED
                result.error_message = f"返回码: {proc.returncode}"

        except subprocess.TimeoutExpired:
            result.end_time = time.time()
            result.elapsed = result.end_time - start_time
            result.status = StageStatus.TIMEOUT
            result.error_message = f"超时 ({config.timeout}s)"

        except Exception as e:
            result.end_time = time.time()
            result.elapsed = result.end_time - start_time
            result.status = StageStatus.FAILED
            result.error_message = f"异常: {str(e)}"

        return result

    def _run_stage_with_retry(self, config: StageConfig) -> StageResult:
        """带重试机制的阶段执行"""
        last_result = None

        for attempt in range(config.max_retries + 1):
            if attempt > 0:
                self.logger.info(f"  [重试 {attempt}/{config.max_retries}] {config.name}")
                time.sleep(2)

            result = self._run_subprocess(config)
            result.retry_count = attempt

            if result.is_success:
                return result

            last_result = result

        return last_result

    def _execute_level(self, stage_names: List[str]) -> None:
        """并行执行同一层级的所有阶段"""
        if not stage_names:
            return

        if self.mode == "sequential" or len(stage_names) == 1:
            for name in stage_names:
                self._execute_single(name)
        else:
            with ThreadPoolExecutor(max_workers=min(self.max_workers, len(stage_names))) as executor:
                futures = {
                    executor.submit(self._execute_single_worker, name): name
                    for name in stage_names
                }
                for future in as_completed(futures):
                    name = futures[future]
                    try:
                        future.result()
                    except Exception as e:
                        self.logger.error(f"  {name} 执行异常: {e}")

    def _execute_single_worker(self, name: str) -> None:
        """工作线程中执行单个阶段"""
        config = self.dag.get_stage_config(name)
        if config is None:
            return

        state = self.dag.stages[name]

        self.logger.info("")
        self.logger.info("=" * 70)
        self.logger.info(f">> {name}: {config.description}")
        self.logger.info("=" * 70)

        state.transition(StageStatus.RUNNING)
        result = self._run_stage_with_retry(config)
        state.record_result(result)

        with self._lock:
            self.results[name] = result

        if result.is_success:
            with self._lock:
                self.completed_stages.add(name)
            self.logger.info(f"[OK] {name} 完成 ({result.elapsed:.1f}s)")
        else:
            with self._lock:
                self.failed_stages.add(name)

            if self.strict:
                self.logger.error(f"[FAIL] {name} {result.status.value} ({result.elapsed:.1f}s)")
                if result.error_message:
                    self.logger.error(f"  错误: {result.error_message}")
                # 标记所有下游阶段为SKIPPED
                downstream = self.dag.get_downstream(name)
                for ds in downstream:
                    ds_state = self.dag.stages[ds]
                    if ds_state.status == StageStatus.PENDING:
                        ds_state.transition(StageStatus.SKIPPED)
                        ds_state.record_result(StageResult(
                            stage_name=ds,
                            status=StageStatus.SKIPPED,
                            error_message=f"上游 {name} 失败"
                        ))
                        with self._lock:
                            self.results[ds] = ds_state.result
            else:
                self.logger.warning(f"[FAIL] {name} {result.status.value} ({result.elapsed:.1f}s) - 继续执行")
                with self._lock:
                    self.completed_stages.add(name)

        if result.stdout:
            lines = result.stdout.strip().split('\n')
            for line in lines[-5:]:
                self.logger.info(f"  {line}")

    def _execute_single(self, name: str) -> None:
        """顺序执行单个阶段 (非并行模式)"""
        self._execute_single_worker(name)

    def run(self) -> Dict[str, StageResult]:
        """执行完整管道"""
        if not self.dag.validate_no_cycles():
            self.logger.error("DAG中存在循环依赖!")
            return {}

        levels = self.dag.get_dependency_order()

        self.logger.info("=" * 70)
        self.logger.info("CEHG-RNP 2.0 管道执行引擎 (v1.0)")
        self.logger.info(f"模式: {'并行 (max_workers=' + str(self.max_workers) + ')' if self.mode == 'parallel' else '顺序'}")
        self.logger.info(f"策略: {'严格 (失败即停)' if self.strict else '宽松 (继续执行)'}")
        self.logger.info(f"执行层级: {len(levels)}")
        self.logger.info("=" * 70)

        for i, level in enumerate(levels):
            self.logger.info(f"\n--- 层级 {i+1}/{len(levels)}: {level} ---")
            for name in level:
                config = self.dag.get_stage_config(name)
                if config:
                    deps_str = ", ".join(config.dependencies) if config.dependencies else "无"
                    self.logger.info(f"  {name}: 依赖=[{deps_str}]")

        self.total_start_time = time.time()

        for i, level in enumerate(levels):
            if self.strict and self.failed_stages:
                self.logger.warning(f"\n跳过层级 {i+1} - 上游阶段失败")
                continue

            active_stages = [
                name for name in level
                if self.dag.stages[name].status == StageStatus.PENDING
                and all(dep not in self.failed_stages for dep in self.dag._reverse_adjacency.get(name, []))
            ]

            if not active_stages:
                continue

            self.logger.info(f"\n{'=' * 70}")
            self.logger.info(f"执行层级 {i+1}/{len(levels)} ({len(active_stages)} 个阶段)")
            self.logger.info(f"{'=' * 70}")

            self._execute_level(active_stages)

        self.total_end_time = time.time()
        total_elapsed = self.total_end_time - self.total_start_time

        self._print_summary(total_elapsed)
        self._save_report(total_elapsed)

        return self.results

    def _print_summary(self, total_elapsed: float) -> None:
        """打印执行摘要"""
        self.logger.info("")
        self.logger.info("=" * 70)
        self.logger.info("管道执行完成!")
        self.logger.info(f"总耗时: {total_elapsed:.1f}s ({total_elapsed/60:.1f}min)")
        self.logger.info("=" * 70)

        success_count = 0
        for name, state in self.dag.stages.items():
            status = state.status
            if status == StageStatus.SUCCESS:
                success_count += 1

            result = self.results.get(name)
            elapsed_str = f"{result.elapsed:.1f}s" if result and result.elapsed > 0 else "N/A"

            icon = {
                StageStatus.SUCCESS: "[OK]",
                StageStatus.FAILED: "[FAIL]",
                StageStatus.TIMEOUT: "[TIME]",
                StageStatus.SKIPPED: "[SKIP]",
                StageStatus.PENDING: "[PEND]",
                StageStatus.RUNNING: "[RUN]",
            }.get(status, "[???]")

            self.logger.info(f"  {icon} {name}: {status.value} ({elapsed_str})")

        self.logger.info(f"\n成功: {success_count}/{len(self.dag.stages)}")

        # 估算并行化收益
        serial_time = sum(
            (r.elapsed if r else 0) for r in self.results.values()
        )
        if serial_time > 0:
            speedup = serial_time / max(total_elapsed, 1)
            self.logger.info(f"并行加速比: {speedup:.2f}x (串行估算: {serial_time:.1f}s)")

    def _save_report(self, total_elapsed: float) -> None:
        """保存执行报告"""
        report_dir = os.path.join(self.results_dir, "pipeline_report")
        os.makedirs(report_dir, exist_ok=True)

        report = {
            "pipeline_version": "CEHG-RNP 2.0",
            "execution_mode": self.mode,
            "strict_mode": self.strict,
            "total_time_s": round(total_elapsed, 1),
            "stages": {}
        }

        for name, state in self.dag.stages.items():
            result = self.results.get(name)
            report["stages"][name] = {
                "status": state.status.value,
                "elapsed_s": round(result.elapsed, 1) if result else 0,
                "return_code": result.return_code if result else -1,
                "retry_count": result.retry_count if result else 0,
                "error": result.error_message[:200] if result and result.error_message else ""
            }

        report_path = os.path.join(report_dir, "execution_report.json")
        with open(report_path, 'w', encoding='utf-8') as f:
            json.dump(report, f, indent=2, ensure_ascii=False)

        self.logger.info(f"\n执行报告已保存: {report_path}")


def create_default_dag() -> PipelineDAG:
    """
    创建默认的CEHG-RNP 2.0管道DAG

    依赖分析:
    - brain_coexpression: 无依赖，独立下载脑共表达数据
    - Stage1 (RMA+limma): 无依赖，产生DEGs
    - Stage2 (SingleCell): 无依赖（独立加载GSE210986），但与Stage1并行无意义（共享I/O）
    - Stage3 (Enrichment): 依赖Stage1 DEGs
    - Stage4 (WGCNA): 依赖Stage1 DEGs + Stage2 markers → 依赖[1,2]
    - Stage5 (PPI): 依赖Stage1 DEGs + Stage3 enrichment → 依赖[1,3]
    - Stage6 (scTenifoldKnk): 依赖Stage2 single cell data → 依赖[2]
    - Stage7 (ML+SHAP): 依赖Stage4 seed pool + Stage2 data → 依赖[2,4]
    - Stage8 (Final): 依赖Stage5,6,7 → 依赖[5,6,7]

    可并行组:
    - Level 0: [brain_coexpression, Stage1] (可并行)
    - Level 1: [Stage2, Stage3] (可并行)
    - Level 2: [Stage4, Stage5, Stage6] (可并行)
    - Level 3: [Stage7]
    - Level 4: [Stage8]
    """
    dag = PipelineDAG()

    stages_config = [
        StageConfig(
            name="brain_coexpression",
            script="integrate_brain_coexpression.py",
            description="脑共表达特征整合 (PsychENCODE + Harmonizome)",
            dependencies=[],
            timeout=1800,
            required_outputs=["brain_coexpression_features.csv", "feature_dimensions.json"],
            priority=10
        ),
        StageConfig(
            name="stage1_rma_degs",
            script="stage1_rma_degs.R",
            description="RMA标准化 + limma差异分析",
            dependencies=[],
            timeout=900,
            is_r_script=True,
            required_outputs=["limma_degs.csv", "sample_annotations.csv"],
            priority=10
        ),
        StageConfig(
            name="stage2_single_cell",
            script="stage2_single_cell.py",
            description="单细胞RNA-seq分析 (Scanpy)",
            dependencies=[],
            timeout=1800,
            required_outputs=["sc_adata.h5ad", "cell_annotations.csv"],
            priority=9
        ),
        StageConfig(
            name="stage3_enrichment",
            script="stage3_enrichment.py",
            description="DEG功能富集分析 (GO/KEGG/GSEA)",
            dependencies=["stage1_rma_degs"],
            timeout=600,
            required_outputs=["go_enrichment.csv", "kegg_enrichment.csv"],
            priority=8
        ),
        StageConfig(
            name="stage4_seed_wgcna",
            script="stage4_seed_wgcna.py",
            description="三层种子池 + WGCNA",
            dependencies=["stage1_rma_degs", "stage2_single_cell"],
            timeout=1200,
            required_outputs=["wgcna_modules.csv", "seed_pool_genes.txt"],
            priority=7
        ),
        StageConfig(
            name="stage5_ppi_mcode",
            script="stage5_ppi_mcode.py",
            description="STRING PPI + MCODE模块识别",
            dependencies=["stage1_rma_degs", "stage3_enrichment"],
            timeout=900,
            required_outputs=["ppi_topology.json", "node_degree_ranking.csv"],
            priority=6
        ),
        StageConfig(
            name="stage6_sctenifold_knockout",
            script="stage6_sctenifold_knockout.py",
            description="scTenifoldKnk 虚拟敲除 (单细胞GRN)",
            dependencies=["stage2_single_cell"],
            timeout=2400,
            required_outputs=["gene_perturbation_scores.csv"],
            priority=5
        ),
        StageConfig(
            name="stage7_ml_shap",
            script="stage7_ml_shap.py",
            description="机器学习集成 + SHAP",
            dependencies=["stage2_single_cell", "stage4_seed_wgcna"],
            timeout=1200,
            required_outputs=["gene_shap_importance.csv"],
            priority=4
        ),
        StageConfig(
            name="stage8_final_targets",
            script="stage8_final_targets.py",
            description="分层筛选 + 最终核心靶点",
            dependencies=["stage5_ppi_mcode", "stage6_sctenifold_knockout", "stage7_ml_shap"],
            timeout=600,
            required_outputs=["core_targets.csv", "tier1_targets.csv"],
            priority=3
        ),
    ]

    for config in stages_config:
        dag.add_stage(config)

    return dag
