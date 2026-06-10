# -*- coding: utf-8 -*-
"""HGT 基因-通路关联预测模型 - 模块化包。

hgt_pathway_discovery/
├── __init__.py          # 包初始化
├── config.yaml          # YAML 配置文件
├── config.py            # 配置加载器 (dataclass)
├── utils.py             # 工具函数 (日志/种子/标准化/降维/WandB)
├── data_loader.py       # 数据加载模块
├── build_graph.py       # 异构图构建 (增强版)
├── model.py             # HGT 模型 (增强版架构)
├── pretrain.py          # GAE 预训练模块
├── train.py             # 训练模块 (增强策略 + NeighborLoader)
├── inference.py         # 推理模块 (MC Dropout/集成)
└── main.py              # 主入口
"""

from .config import load_config, Config
from .utils import set_seed, setup_logging, setup_tensorboard, try_compile, setup_wandb
from .data_loader import load_all_data
from .build_graph import build_hetero_graph, normalize_node_features, apply_cached_normalization
from .model import HGTModel, focal_bce_loss, NegEdgeSampler
from .pretrain import pretrain_gae, HeteroGAE
from .train import (
    train_fold, train_final, cross_validate, train_fold_minibatch,
    remove_edges_from_data, PlattScaler,
)
from .inference import (
    predict_bridge_pathways,
    ensemble_predict_bridge_pathways,
    evaluate_inductive_subset,
    mc_dropout_predict_per_gene,
)