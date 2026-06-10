"""
分类器模块
===========
封装所有基础分类器的构建逻辑，提供统一接口。

支持的分类器:
  - L1_LR: L1正则化逻辑回归
  - L2_LR: L2正则化逻辑回归
  - ElasticNet_LR: ElasticNet逻辑回归 (SGD)
  - RF: 随机森林
  - GB: 梯度提升
  - NB: 高斯朴素贝叶斯
  - SVC: 支持向量机 + 概率标定
  - ExtraTrees: 极端随机树
  - PAC: 被动攻击分类器 + 概率标定
  - XGBoost: (可选)
  - LightGBM: (可选)
"""

from abc import ABC, abstractmethod
from typing import Dict, Optional, Tuple

from sklearn.linear_model import LogisticRegression, SGDClassifier, PassiveAggressiveClassifier
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier, ExtraTreesClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC
from sklearn.calibration import CalibratedClassifierCV

from .config import GPUConfig


class ClassifierStrategy(ABC):
    """分类器策略基类"""

    @abstractmethod
    def get_name(self) -> str:
        """获取分类器唯一名称"""
        pass

    @abstractmethod
    def create_classifier(self, random_state: int, gpu_config: Optional[GPUConfig] = None):
        """创建分类器实例"""
        pass

    @property
    def supports_class_weight(self) -> bool:
        """是否支持 class_weight 参数"""
        return True

    @property
    def requires_sample_weight(self) -> bool:
        """是否需要手动传递 sample_weight"""
        return False


class L1LRStrategy(ClassifierStrategy):
    """L1正则化逻辑回归"""

    def get_name(self) -> str:
        return "L1_LR"

    def create_classifier(self, random_state: int, gpu_config=None):
        return LogisticRegression(
            penalty='l1', solver='liblinear', C=0.1,
            class_weight='balanced', max_iter=5000, random_state=random_state,
        )


class L2LRStrategy(ClassifierStrategy):
    """L2正则化逻辑回归"""

    def get_name(self) -> str:
        return "L2_LR"

    def create_classifier(self, random_state: int, gpu_config=None):
        return LogisticRegression(
            penalty='l2', C=1.0,
            class_weight='balanced', max_iter=5000, random_state=random_state,
        )


class ElasticNetLRStrategy(ClassifierStrategy):
    """ElasticNet 逻辑回归"""

    def get_name(self) -> str:
        return "ElasticNet_LR"

    def create_classifier(self, random_state: int, gpu_config=None):
        return SGDClassifier(
            loss='log_loss', penalty='elasticnet', alpha=0.001, l1_ratio=0.5,
            class_weight='balanced', max_iter=2000, random_state=random_state,
        )


class RFStrategy(ClassifierStrategy):
    """随机森林"""

    def get_name(self) -> str:
        return "RF"

    def create_classifier(self, random_state: int, gpu_config=None):
        return RandomForestClassifier(
            n_estimators=200, class_weight='balanced',
            random_state=random_state, n_jobs=1,  # n_jobs=1 避免嵌套并行
        )


class GBStrategy(ClassifierStrategy):
    """梯度提升"""

    def get_name(self) -> str:
        return "GB"

    def create_classifier(self, random_state: int, gpu_config=None):
        return GradientBoostingClassifier(
            n_estimators=100, learning_rate=0.1, random_state=random_state,
        )

    @property
    def supports_class_weight(self) -> bool:
        return False

    @property
    def requires_sample_weight(self) -> bool:
        return True


class NBStrategy(ClassifierStrategy):
    """高斯朴素贝叶斯"""

    def get_name(self) -> str:
        return "NB"

    def create_classifier(self, random_state: int, gpu_config=None):
        return GaussianNB()

    @property
    def supports_class_weight(self) -> bool:
        return False

    @property
    def requires_sample_weight(self) -> bool:
        return True


class SVCStrategy(ClassifierStrategy):
    """SVM + CalibratedClassifierCV (概率标定)"""

    def get_name(self) -> str:
        return "SVC"

    def create_classifier(self, random_state: int, gpu_config=None):
        base_svc = SVC(C=1.0, kernel='rbf', class_weight='balanced',
                       random_state=random_state, max_iter=5000, probability=False)
        return CalibratedClassifierCV(base_svc, cv=3, method='sigmoid')


class ExtraTreesStrategy(ClassifierStrategy):
    """极端随机树"""

    def get_name(self) -> str:
        return "ExtraTrees"

    def create_classifier(self, random_state: int, gpu_config=None):
        return ExtraTreesClassifier(
            n_estimators=200, class_weight='balanced',
            random_state=random_state, n_jobs=1,
        )


class PACStrategy(ClassifierStrategy):
    """被动攻击 + CalibratedClassifierCV"""

    def get_name(self) -> str:
        return "PAC"

    def create_classifier(self, random_state: int, gpu_config=None):
        base_pac = PassiveAggressiveClassifier(
            C=0.1, class_weight='balanced', max_iter=2000, random_state=random_state,
        )
        return CalibratedClassifierCV(base_pac, cv=3, method='sigmoid')


class XGBoostStrategy(ClassifierStrategy):
    """XGBoost — 支持 GPU 加速"""

    def get_name(self) -> str:
        return "XGBoost"

    def create_classifier(self, random_state: int, gpu_config=None):
        import xgboost as xgb
        params = {
            'n_estimators': 200, 'learning_rate': 0.1,
            'eval_metric': 'logloss', 'random_state': random_state,
            'verbosity': 0, 'n_jobs': 1,
        }
        if gpu_config is not None:
            params.update(gpu_config.xgb_params)
        return xgb.XGBClassifier(**params)


class LightGBMStrategy(ClassifierStrategy):
    """LightGBM — 支持 GPU 加速"""

    def get_name(self) -> str:
        return "LightGBM"

    def create_classifier(self, random_state: int, gpu_config=None):
        import lightgbm as lgb
        params = {
            'n_estimators': 200, 'learning_rate': 0.1,
            'class_weight': 'balanced', 'random_state': random_state,
            'verbose': -1,
        }
        if gpu_config is not None:
            params.update(gpu_config.lgb_params)
        return lgb.LGBMClassifier(**params)


class ClassifierRegistry:
    """分类器策略注册表"""

    def __init__(self):
        self._strategies: Dict[str, ClassifierStrategy] = {}

    def register(self, strategy: ClassifierStrategy):
        """注册分类器策略"""
        self._strategies[strategy.get_name()] = strategy

    def get_strategy(self, name: str) -> Optional[ClassifierStrategy]:
        """按名称获取策略"""
        return self._strategies.get(name)

    def get_all_strategies(self) -> Dict[str, ClassifierStrategy]:
        """获取所有已注册策略"""
        return dict(self._strategies)

    def get_names(self):
        """获取所有策略名称"""
        return list(self._strategies.keys())

    @classmethod
    def create_default(cls) -> 'ClassifierRegistry':
        """创建包含所有可用分类器的注册表"""
        registry = cls()

        # 基础分类器 (始终可用)
        registry.register(L1LRStrategy())
        registry.register(L2LRStrategy())
        registry.register(ElasticNetLRStrategy())
        registry.register(RFStrategy())
        registry.register(GBStrategy())
        registry.register(NBStrategy())
        registry.register(SVCStrategy())
        registry.register(ExtraTreesStrategy())
        registry.register(PACStrategy())

        # XGBoost (可选)
        try:
            import xgboost  # noqa
            registry.register(XGBoostStrategy())
        except ImportError:
            pass

        # LightGBM (可选)
        try:
            import lightgbm  # noqa
            registry.register(LightGBMStrategy())
        except ImportError:
            pass

        return registry