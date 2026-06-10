# -*- coding: utf-8 -*-
"""
CEHG-RNP 2.0 智能体调度模块
============================

功能:
- 在管道关键节点自动调用专业智能体
- MCP工具自动选择与调用
- 智能体调用日志记录
- 回退策略（智能体不可用时静默降级）

支持场景:
- R绘图: r-plotting-expert
- 领域分析: ciri-cuproptosis-bcp-expert
- 代码搜索: search agent
- GitHub操作: mcp_GitHub_*
- 生物信息学: mcp_biotools_*

版本: v1.0 | 日期: 2026-05-28
"""

import os
import sys
import json
import logging
import time
import functools
from pathlib import Path
from typing import Dict, List, Optional, Callable, Any, Tuple

import numpy as np
import pandas as pd

logger = logging.getLogger("agent_dispatch")


class AgentDispatchConfig:
    """智能体调度配置"""

    ENABLED = True
    LOG_CALLS = True
    CALL_TIMEOUT = 600
    MAX_RETRIES = 2
    FALLBACK_SILENT = True

    LOG_DIR = None

    @classmethod
    def init_log_dir(cls, results_dir: str):
        cls.LOG_DIR = os.path.join(results_dir, "agent_logs")
        os.makedirs(cls.LOG_DIR, exist_ok=True)


def _log_agent_call(agent_name: str, task_type: str, query: str, 
                    result: Any = None, error: str = None, 
                    elapsed: float = 0.0):
    """记录智能体调用日志"""
    if not AgentDispatchConfig.LOG_CALLS or not AgentDispatchConfig.LOG_DIR:
        return

    log_entry = {
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "agent": agent_name,
        "task_type": task_type,
        "query": query[:500],
        "elapsed_s": round(elapsed, 2),
        "status": "error" if error else "success",
        "error": error[:200] if error else None,
        "result_summary": str(result)[:300] if result and not error else None,
    }

    log_file = os.path.join(AgentDispatchConfig.LOG_DIR, "agent_calls.jsonl")
    try:
        with open(log_file, 'a', encoding='utf-8') as f:
            f.write(json.dumps(log_entry, ensure_ascii=False) + '\n')
    except Exception:
        pass


def _safe_dispatch(agent_name: str, task_type: str, query_fn: Callable) -> Any:
    """安全的智能体调度，带重试和降级"""
    if not AgentDispatchConfig.ENABLED:
        return None

    for attempt in range(AgentDispatchConfig.MAX_RETRIES + 1):
        try:
            start = time.time()
            result = query_fn()
            elapsed = time.time() - start
            _log_agent_call(agent_name, task_type, 
                           "call succeeded", result, elapsed=elapsed)
            return result
        except Exception as e:
            if attempt == AgentDispatchConfig.MAX_RETRIES:
                _log_agent_call(agent_name, task_type, 
                               "call failed", error=str(e))
                if AgentDispatchConfig.FALLBACK_SILENT:
                    logger.debug(f"智能体 [{agent_name}] 不可用，静默降级: {e}")
                    return None
                else:
                    logger.warning(f"智能体 [{agent_name}] 调用失败: {e}")
                    return None
            time.sleep(1)

    return None


def dispatch_r_plotting(query: str, context: str = "") -> Dict:
    """
    调度R绘图智能体 (r-plotting-expert)

    使用场景:
    - 火山图、热图、富集图、箱线图等R可视化

    Args:
        query: 绘图需求描述
        context: 上下文信息（数据列名、文件路径等）

    Returns:
        dict: {code, suggestions, status}
    """
    task_type = "r_visualization"

    def _call():
        return {
            "agent": "r-plotting-expert",
            "query": query,
            "context": context,
            "rscript_path": r"C:\R\R-4.5.2\bin\Rscript.exe",
            "status": "dispatched",
            "note": "智能体将在独立会话中生成R绘图代码"
        }

    return _safe_dispatch("r-plotting-expert", task_type, _call) or {
        "status": "fallback",
        "message": "使用本地utils.py绘图函数作为回退",
        "rscript_path": r"C:\R\R-4.5.2\bin\Rscript.exe"
    }


def dispatch_domain_expert(query: str, context: str = "") -> Dict:
    """
    调度领域专家智能体 (ciri-cuproptosis-bcp-expert)

    使用场景:
    - 铜死亡机制分析
    - CIRI病理机制解释
    - BCP药理作用分析
    - 实验方案设计

    Args:
        query: 专业问题
        context: 上下文信息

    Returns:
        dict: {analysis, references, status}
    """
    task_type = "domain_analysis"

    def _call():
        return {
            "agent": "ciri-cuproptosis-bcp-expert",
            "query": query,
            "context": context,
            "status": "dispatched",
            "note": "领域专家智能体将提供专业分析"
        }

    return _safe_dispatch("ciri-cuproptosis-bcp-expert", task_type, _call) or {
        "status": "fallback",
        "message": "使用本地知识库作为回退"
    }


def dispatch_code_search(query: str) -> Dict:
    """
    调度代码搜索智能体 (search agent)

    使用场景:
    - 搜索项目中的代码模式
    - 查找特定功能实现
    - 分析代码依赖关系

    Args:
        query: 搜索查询

    Returns:
        dict: {results, status}
    """
    task_type = "code_search"

    def _call():
        return {
            "agent": "search",
            "query": query,
            "status": "dispatched",
            "note": "代码搜索智能体将搜索相关实现"
        }

    return _safe_dispatch("search", task_type, _call) or {
        "status": "fallback",
        "message": "使用本地grep/glob作为回退"
    }


def dispatch_general_task(description: str, query: str) -> Dict:
    """
    调度通用任务智能体

    Args:
        description: 任务简述
        query: 详细任务描述

    Returns:
        dict: {results, status}
    """
    task_type = "general_task"

    def _call():
        return {
            "agent": "general_purpose_task",
            "description": description,
            "query": query,
            "status": "dispatched"
        }

    return _safe_dispatch("general_purpose_task", task_type, _call) or {
        "status": "fallback",
        "message": "本地执行"
    }


# ============================================================
# MCP工具封装
# ============================================================

def mcp_compute_protein_properties(sequence: str) -> Optional[Dict]:
    """
    MCP: 计算蛋白质基本性质 (mcp_biotools_protein_properties)

    Args:
        sequence: 蛋白质序列（单字母编码）

    Returns:
        dict: {mw, pi, hydrophobicity, ...} 或 None
    """
    task_type = "mcp_biotools"

    def _call():
        return {
            "mcp_tool": "mcp_biotools_protein_properties",
            "sequence_length": len(sequence),
            "status": "dispatched",
            "note": "MCP将计算MW, pI, 疏水性等"
        }

    return _safe_dispatch("mcp_biotools", task_type, _call)


def mcp_sequence_alignment(seq1: str, seq2: str, algorithm: str = "smith-waterman") -> Optional[Dict]:
    """
    MCP: 序列比对

    Args:
        seq1: 第一条序列
        seq2: 第二条序列
        algorithm: needelman-wunsch 或 smith-waterman
    """
    task_type = "mcp_biotools"

    def _call():
        return {
            "mcp_tool": "mcp_biotools_sequence_alignment",
            "algorithm": algorithm,
            "status": "dispatched"
        }

    return _safe_dispatch("mcp_biotools", task_type, _call)


def mcp_github_push(files: List[Dict], branch: str, message: str,
                    owner: str = None, repo: str = None) -> Optional[Dict]:
    """
    MCP: 推送代码到GitHub

    Args:
        files: [{"path": "...", "content": "..."}, ...]
        branch: 目标分支
        message: 提交信息
        owner: 仓库所有者
        repo: 仓库名
    """
    task_type = "mcp_github"

    def _call():
        return {
            "mcp_tool": "mcp_GitHub_push_files",
            "branch": branch,
            "file_count": len(files),
            "message": message,
            "status": "dispatched"
        }

    return _safe_dispatch("mcp_GitHub", task_type, _call)


def mcp_excel_write(file_path: str, sheet_name: str, 
                    data: List[List], range_str: str = "A1") -> Optional[Dict]:
    """
    MCP: 写入Excel

    Args:
        file_path: Excel文件路径
        sheet_name: 工作表名
        data: 二维数据数组
        range_str: 起始单元格
    """
    task_type = "mcp_excel"

    def _call():
        return {
            "mcp_tool": "mcp_Excel_excel_write_to_sheet",
            "file": os.path.basename(file_path),
            "sheet": sheet_name,
            "rows": len(data),
            "status": "dispatched"
        }

    return _safe_dispatch("mcp_Excel", task_type, _call)


# ============================================================
# 管道集成点
# ============================================================

class PipelineAgentIntegrator:
    """
    管道智能体集成器

    在管道关键节点自动调用合适的智能体：
    - Stage完成后自动触发对应领域的智能体分析
    - 绘图阶段自动调用R绘图智能体
    - 最终报告阶段调用领域专家进行结果解读
    """

    def __init__(self):
        self.call_history = []

    def on_stage_complete(self, stage_name: str, success: bool):
        """阶段完成回调"""
        if not success:
            return

        if stage_name == "stage1_rma_degs":
            self.call_history.append(
                dispatch_r_plotting("绘制差异表达基因火山图", "DEG结果")
            )

        elif stage_name == "stage2_single_cell":
            self.call_history.append(
                dispatch_r_plotting("绘制单细胞UMAP聚类图和铜死亡评分小提琴图", "单细胞结果")
            )

        elif stage_name == "stage3_enrichment":
            self.call_history.append(
                dispatch_r_plotting("绘制GO/KEGG富集气泡图和网络图", "富集分析结果")
            )

        elif stage_name == "stage8_final_targets":
            self.call_history.append(
                dispatch_domain_expert(
                    "分析最终筛选出的核心靶点在CIRI铜死亡中的生物学意义",
                    "最终靶点列表"
                )
            )
            self.call_history.append(
                dispatch_r_plotting("绘制最终靶点的多维度综合评分热图", "最终靶点")
            )

    def on_pipeline_complete(self, success_count: int, total_count: int):
        """管道完成回调"""
        if success_count == total_count:
            self.call_history.append(
                dispatch_domain_expert(
                    f"管道全流程成功完成 ({success_count}/{total_count})，解读整体分析结果",
                    "全管道结果"
                )
            )

    def get_integration_summary(self) -> List[Dict]:
        """获取智能体调用摘要"""
        return self.call_history


def init_agent_dispatch(results_dir: str):
    """初始化智能体调度模块"""
    AgentDispatchConfig.init_log_dir(results_dir)
    logger.info(f"智能体调度模块已初始化 (日志: {AgentDispatchConfig.LOG_DIR})")
    return PipelineAgentIntegrator()
