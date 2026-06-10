"""
多算法融合桥梁靶点预测系统 (v3.x — 模块化重构)
==================================================
网络药理学第一阶段：药物-靶点-疾病桥梁基因发现

模块结构:
  - config.py:              配置管理
  - utils.py:               工具函数 (日志, 样本权重)
  - data_loader.py:         数据加载与预处理
  - feature_engineering.py: 特征工程策略
  - classifiers.py:         分类器构建
  - trainer.py:             模型训练与交叉验证
  - ensemble.py:            模型集成与排名
  - main.py:                主流程编排

版本: v3.1 (模块化重构)
"""

__version__ = "3.1.0"