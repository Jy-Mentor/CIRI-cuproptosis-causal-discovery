"""
特征工程模块
=============
封装所有特征工程策略，提供统一接口。

支持的策略:
  - RAW: 原始特征
  - PCA: n_components=10/20/30/50/80/100
  - Lasso: max_features=20/30/50/80
  - PLS: n_components=5/10/15
  - KBest: k=50/100
"""

from abc import ABC, abstractmethod
from typing import Any, Dict, Optional, Type

from sklearn.decomposition import PCA
from sklearn.cross_decomposition import PLSRegression
from sklearn.feature_selection import SelectFromModel, SelectKBest, f_classif
from sklearn.linear_model import LassoCV


class FeatureEngineeringStrategy(ABC):
    """特征工程策略基类"""

    @abstractmethod
    def get_name(self) -> str:
        """获取策略唯一名称"""
        pass

    @abstractmethod
    def create_transformer(self, random_state: int):
        """创建转换器实例"""
        pass

    @property
    def is_supervised(self) -> bool:
        """是否需要 y 标签 (有监督特征工程)"""
        return False


class RawFeatureStrategy(FeatureEngineeringStrategy):
    """原始特征 — 不做任何变换"""

    def get_name(self) -> str:
        return "raw"

    def create_transformer(self, random_state: int):
        return None

    @property
    def is_supervised(self) -> bool:
        return False


class PCAStrategy(FeatureEngineeringStrategy):
    """PCA 降维"""

    def __init__(self, n_components: int):
        self.n_components = n_components

    def get_name(self) -> str:
        return f"pca_{self.n_components}"

    def create_transformer(self, random_state: int):
        return PCA(n_components=self.n_components, random_state=random_state)


class LassoStrategy(FeatureEngineeringStrategy):
    """LassoCV 特征选择"""

    def __init__(self, max_features: int):
        self.max_features = max_features

    def get_name(self) -> str:
        return f"lasso_{self.max_features}"

    def create_transformer(self, random_state: int):
        return SelectFromModel(
            LassoCV(cv=2, n_alphas=20, random_state=random_state,
                    max_iter=2000, n_jobs=-1),
            max_features=self.max_features,
        )

    @property
    def is_supervised(self) -> bool:
        return True


class PLSStrategy(FeatureEngineeringStrategy):
    """PLS 偏最小二乘"""

    def __init__(self, n_components: int):
        self.n_components = n_components

    def get_name(self) -> str:
        return f"pls_{self.n_components}"

    def create_transformer(self, random_state: int):
        return PLSRegression(n_components=self.n_components, scale=False)

    @property
    def is_supervised(self) -> bool:
        return True


class KBestStrategy(FeatureEngineeringStrategy):
    """ANOVA F-test 特征选择"""

    def __init__(self, k: int):
        self.k = k

    def get_name(self) -> str:
        return f"kbest_{self.k}"

    def create_transformer(self, random_state: int):
        return SelectKBest(f_classif, k=self.k)

    @property
    def is_supervised(self) -> bool:
        return True


# 定义需要 y 的有监督特征工程类型
SUPERVISED_FE_TYPES: set = {PLSRegression, SelectFromModel, SelectKBest}


class FeatureEngineeringRegistry:
    """特征工程策略注册表"""

    def __init__(self):
        self._strategies: Dict[str, FeatureEngineeringStrategy] = {}

    def register(self, strategy: FeatureEngineeringStrategy):
        """注册特征工程策略"""
        self._strategies[strategy.get_name()] = strategy

    def get_strategy(self, name: str) -> Optional[FeatureEngineeringStrategy]:
        """按名称获取策略"""
        return self._strategies.get(name)

    def get_all_strategies(self) -> Dict[str, FeatureEngineeringStrategy]:
        """获取所有已注册策略"""
        return dict(self._strategies)

    def get_names(self):
        """获取所有策略名称"""
        return list(self._strategies.keys())

    @classmethod
    def create_default(cls) -> 'FeatureEngineeringRegistry':
        """创建包含所有默认策略的注册表"""
        registry = cls()
        registry.register(RawFeatureStrategy())

        for n in [10, 20, 30, 50, 80, 100]:
            registry.register(PCAStrategy(n))
        for n in [20, 30, 50, 80]:
            registry.register(LassoStrategy(n))
        for n in [5, 10, 15]:
            registry.register(PLSStrategy(n))
        for k in [50, 100]:
            registry.register(KBestStrategy(k))

        return registry