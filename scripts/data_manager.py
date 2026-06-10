# -*- coding: utf-8 -*-
"""
StageDataManager - 跨阶段数据访问管理器
=======================================
统一管理铜死亡生物信息学分析管线中各stage之间的数据加载和访问

功能:
- 统一加载接口（JSON、CSV等）
- 自动缓存机制避免重复加载
- 标准化错误处理和日志记录
- 列名变体自动识别
- 数据格式验证

版本: v1.0 | 日期: 2026-05-21
"""

import os
import json
import logging
from typing import Dict, List, Optional, Any, Tuple, Union
from pathlib import Path

import numpy as np
import pandas as pd


class StageDataManager:
    """
    跨阶段数据管理器
    
    封装所有跨stage数据加载逻辑，提供统一的接口来加载、缓存和管理
    来自不同分析阶段的数据文件。
    
    使用示例:
    ```python
    dm = StageDataManager()
    
    # 加载JSON文件
    topology = dm.load_json('stage5_string_ppi', 'ppi_topology.json')
    
    # 加载CSV并转换为字典
    shap_dict = dm.load_csv_as_dict(
        'stage7_ml_shap',
        'gene_shap_importance.csv',
        key_col='Gene',
        value_col='SHAP_importance'
    )
    
    # 批量加载一个stage的多个文件
    results = dm.load_stage('stage5_string_ppi', {
        'topology': ('json', {}),
        'degree_ranking': ('csv_dict', {'key_col': 'Gene', 'value_col': 'Degree'})
    })
    ```
    """
    
    # 列名变体映射表
    COLUMN_VARIANTS = {
        'Gene': ['Gene', 'gene', 'GeneSymbol', 'gene_symbol', 'Gene.symbol', 'SYMBOL', 'gene_name'],
        'gene': ['Gene', 'gene', 'GeneSymbol', 'gene_symbol', 'Gene.symbol', 'SYMBOL', 'gene_name'],
        'SHAP_importance': ['SHAP_importance', 'shap_importance', 'importance', 'SHAP', 'shap'],
        'perturbation_score': ['perturbation_score', 'score', 'perturbation', 'pert_score'],
        'Degree': ['Degree', 'degree', 'node_degree', 'Degree_ranking'],
        'Module': ['Module', 'module', 'ModuleColor', 'moduleColor', 'color', 'module_color'],
        'P_Value': ['P_Value', 'p_value', 'PValue', 'pvalue', 'padj', 'adj.P.Val', 'FDR'],
    }
    
    def __init__(self, results_dir: str = None, cache_enabled: bool = True):
        """
        初始化数据管理器
        
        Args:
            results_dir: 结果目录根路径，默认为config.RESULTS_DIR
            cache_enabled: 是否启用缓存，默认启用
        """
        if results_dir is None:
            from config import RESULTS_DIR
            self.results_dir = RESULTS_DIR
        else:
            self.results_dir = results_dir
        
        self.cache_enabled = cache_enabled
        self._cache = {}
        self.logger = logging.getLogger(__name__)
        
        if not self.logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter(
                '%(asctime)s | %(levelname)-8s | [DataManager] %(message)s',
                datefmt='%Y-%m-%d %H:%M:%S'
            )
            handler.setFormatter(formatter)
            self.logger.addHandler(handler)
            self.logger.setLevel(logging.INFO)
        
        self.logger.debug(f"StageDataManager initialized with results_dir: {self.results_dir}")
    
    def _get_filepath(self, stage: str, filename: str) -> str:
        """构建完整的文件路径"""
        return os.path.join(self.results_dir, stage, filename)
    
    def _find_column(self, df: pd.DataFrame, target: str) -> Optional[str]:
        """
        查找DataFrame中的列名，支持常见变体
        
        Args:
            df: DataFrame对象
            target: 目标列名
            
        Returns:
            实际列名，如果找不到则返回None
        """
        if target in df.columns:
            return target
        
        variants = self.COLUMN_VARIANTS.get(target, [target])
        for variant in variants:
            if variant in df.columns:
                if variant != target:
                    self.logger.debug(f"列名映射: {target} -> {variant}")
                return variant
        
        return None
    
    def load_json(self, stage: str, filename: str, default=None) -> Any:
        """
        加载JSON文件
        
        Args:
            stage: stage目录名称，如'stage5_string_ppi'
            filename: 文件名，如'ppi_topology.json'
            default: 文件不存在或加载失败时的默认值
            
        Returns:
            JSON解析后的数据对象
        """
        cache_key = f"json:{stage}:{filename}"
        if self.cache_enabled and cache_key in self._cache:
            return self._cache[cache_key]
        
        filepath = self._get_filepath(stage, filename)
        
        if not os.path.exists(filepath):
            self.logger.warning(f"文件不存在: {filepath}")
            if self.cache_enabled:
                self._cache[cache_key] = default
            return default
        
        try:
            with open(filepath, 'r', encoding='utf-8') as f:
                data = json.load(f)
            
            self.logger.info(f"JSON加载成功: {filename} ({stage})")
            
            if self.cache_enabled:
                self._cache[cache_key] = data
            
            return data
        except Exception as e:
            self.logger.error(f"JSON加载失败 {filepath}: {e}")
            if self.cache_enabled:
                self._cache[cache_key] = default
            return default
    
    def load_csv(self, stage: str, filename: str, **kwargs) -> Optional[pd.DataFrame]:
        """
        加载CSV文件为DataFrame
        
        Args:
            stage: stage目录名称
            filename: 文件名
            **kwargs: 传递给pd.read_csv的其他参数
            
        Returns:
            DataFrame对象，失败时返回None
        """
        cache_key = f"csv:{stage}:{filename}"
        if self.cache_enabled and cache_key in self._cache:
            return self._cache[cache_key]
        
        filepath = self._get_filepath(stage, filename)
        
        if not os.path.exists(filepath):
            self.logger.warning(f"文件不存在: {filepath}")
            return None
        
        try:
            df = pd.read_csv(filepath, **kwargs)
            self.logger.info(f"CSV加载成功: {filename} ({stage}), shape={df.shape}")
            
            if self.cache_enabled:
                self._cache[cache_key] = df
            
            return df
        except Exception as e:
            self.logger.error(f"CSV加载失败 {filepath}: {e}")
            return None
    
    def load_csv_as_dict(self, stage: str, filename: str,
                        key_col: str, value_col: str,
                        key_upper: bool = True, default=None) -> Dict:
        """
        加载CSV文件并转换为字典
        
        适用于将基因名映射到分数、模块等场景
        
        Args:
            stage: stage目录名称
            filename: 文件名
            key_col: 用作字典键的列名
            value_col: 用作字典值的列名
            key_upper: 是否将键转换为大写
            default: 失败时的默认值
            
        Returns:
            字典 {key_col: value_col}
        """
        cache_key = f"csv_dict:{stage}:{filename}:{key_col}:{value_col}"
        if self.cache_enabled and cache_key in self._cache:
            return self._cache[cache_key]
        
        default_dict = default if default is not None else {}
        
        df = self.load_csv(stage, filename)
        if df is None or df.empty:
            self.logger.warning(f"CSV为空或加载失败: {filename}")
            if self.cache_enabled:
                self._cache[cache_key] = default_dict
            return default_dict
        
        actual_key_col = self._find_column(df, key_col)
        actual_value_col = self._find_column(df, value_col)
        
        if actual_key_col is None:
            self.logger.error(f"找不到键列 {key_col}，可用列: {list(df.columns)}")
            if self.cache_enabled:
                self._cache[cache_key] = default_dict
            return default_dict
        
        if actual_value_col is None:
            self.logger.error(f"找不到值列 {value_col}，可用列: {list(df.columns)}")
            if self.cache_enabled:
                self._cache[cache_key] = default_dict
            return default_dict
        
        result = {}
        for _, row in df.iterrows():
            key = str(row[actual_key_col]).upper().strip() if key_upper else str(row[actual_key_col]).strip()
            if key and key != 'nan':
                result[key] = row[actual_value_col]
        
        self.logger.info(f"CSV字典转换成功: {filename} ({len(result)} 键值对)")
        
        if self.cache_enabled:
            self._cache[cache_key] = result
        
        return result
    
    def load_stage(self, stage: str, file_configs: Dict[str, Tuple]) -> Dict:
        """
        批量加载一个stage的多个文件
        
        Args:
            stage: stage目录名称
            file_configs: 文件配置字典
                {
                    'key_in_result': ('load_method', {method_kwargs})
                }
                
                支持的load_method:
                - 'json': 加载JSON文件，需要filename参数
                - 'csv': 加载CSV为DataFrame，需要filename参数
                - 'csv_dict': 加载CSV并转换为字典，需要filename, key_col, value_col参数
                
        Returns:
            包含所有加载数据的字典
            
        示例:
        ```python
        results = dm.load_stage('stage5_string_ppi', {
            'topology': ('json', {'filename': 'ppi_topology.json'}),
            'degree_ranking': ('csv_dict', {
                'filename': 'node_degree_ranking.csv',
                'key_col': 'Gene',
                'value_col': 'Degree'
            })
        })
        ```
        """
        results = {}
        
        for result_key, (method, kwargs) in file_configs.items():
            filename = kwargs.get('filename', result_key)
            
            if method == 'json':
                default = kwargs.get('default', {})
                results[result_key] = self.load_json(stage, filename, default)
                
            elif method == 'csv':
                read_kwargs = {k: v for k, v in kwargs.items() if k != 'filename'}
                df = self.load_csv(stage, filename, **read_kwargs)
                results[result_key] = df
                
            elif method == 'csv_dict':
                key_col = kwargs.get('key_col', 'Gene')
                value_col = kwargs.get('value_col', 'value')
                key_upper = kwargs.get('key_upper', True)
                default = kwargs.get('default', {})
                results[result_key] = self.load_csv_as_dict(
                    stage, filename, key_col, value_col, key_upper, default
                )
            else:
                self.logger.error(f"不支持的加载方法: {method}")
                results[result_key] = None
        
        return results
    
    def validate_stage_output(self, stage: str, required_files: List[str]) -> bool:
        """
        验证某个stage的输出文件是否存在
        
        Args:
            stage: stage目录名称
            required_files: 必需的文件名列表
            
        Returns:
            如果所有文件都存在则返回True
        """
        missing = []
        for filename in required_files:
            filepath = self._get_filepath(stage, filename)
            if not os.path.exists(filepath):
                missing.append(filename)
        
        if missing:
            self.logger.warning(f"Stage {stage} 缺少文件: {missing}")
            return False
        
        self.logger.info(f"Stage {stage} 验证通过: {len(required_files)} 文件")
        return True
    
    def clear_cache(self, stage: str = None):
        """
        清空缓存
        
        Args:
            stage: 如果指定，只清空该stage的缓存；否则清空所有缓存
        """
        if stage is None:
            self._cache.clear()
            self.logger.info("缓存已清空")
        else:
            keys_to_remove = [k for k in self._cache if f":{stage}:" in k]
            for key in keys_to_remove:
                del self._cache[key]
            self.logger.info(f"已清空 {stage} 的缓存 ({len(keys_to_remove)} 条目)")
    
    def get_cache_info(self) -> Dict:
        """获取缓存统计信息"""
        cache_size = len(self._cache)
        cache_by_stage = {}
        
        for key in self._cache:
            parts = key.split(':')
            if len(parts) >= 3:
                stage = parts[1]
                cache_by_stage[stage] = cache_by_stage.get(stage, 0) + 1
        
        return {
            'total_entries': cache_size,
            'by_stage': cache_by_stage,
            'enabled': self.cache_enabled
        }