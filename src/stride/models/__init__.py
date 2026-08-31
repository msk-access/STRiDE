"""Bundled model resources and model factory for STRiDE.

Provides helper functions to locate the default trained model, MSI site-list files,
and factory functions for creating trainers and predictors.
"""

import importlib.resources
import logging

from .svm import SVMPredictor, SVMTrainer

logger = logging.getLogger(__name__)


def get_default_model_path() -> str:
    """Return the absolute path to the bundled default MSI model.

    The default model is an SGD classifier pipeline saved as a joblib file.
    Users can override this by passing ``--model-joblib`` on the CLI.
    """
    path = str(importlib.resources.files("stride.models") / "msi_sgd_v1.joblib")
    logger.debug("Default model path: %s", path)
    return path


def get_default_sites_path() -> str:
    """Return the absolute path to the bundled 170-site MSI locus list.

    Users can override this by passing ``--site-list`` on the CLI.
    """
    path = str(importlib.resources.files("stride.data") / "msi_sites_170.txt")
    logger.debug("Default sites path: %s", path)
    return path


def get_trainer(method: str = "svm", **kwargs):
    """Factory function to retrieve a Trainer instance by method name.

    Supported methods: 'svm', 'tabpfn', 'tabpfn_access_only', 'tabpfn_access_impact'
    """
    method_key = method.lower().replace("-", "_")
    if method_key == "svm":
        return SVMTrainer(**kwargs)
    elif method_key.startswith("tabpfn"):
        try:
            from .tabpfn import TabPFNTrainer

            variant = kwargs.pop("variant", None)
            if method_key == "tabpfn_access_only":
                variant = "access_only"
            elif method_key == "tabpfn_access_impact":
                variant = "access_impact"

            return TabPFNTrainer(variant=variant, **kwargs) if variant else TabPFNTrainer(**kwargs)
        except ImportError as err:
            raise ImportError(
                "TabPFN trainer requires optional dependencies. Install via: pip install '.[tabpfn]'"
            ) from err
    else:
        raise ValueError(
            f"Unknown model method '{method}'. Supported methods: 'svm', 'tabpfn', 'tabpfn_access_only', 'tabpfn_access_impact'"
        )


def get_predictor(method: str = "svm", **kwargs):
    """Factory function to retrieve a Predictor instance by method name.

    Supported methods: 'svm', 'tabpfn', 'tabpfn_access_only', 'tabpfn_access_impact'
    """
    method_key = method.lower().replace("-", "_")
    if method_key == "svm":
        return SVMPredictor(**kwargs)
    elif method_key.startswith("tabpfn"):
        try:
            from .tabpfn import TabPFNPredictor

            variant = kwargs.pop("variant", None)
            if method_key == "tabpfn_access_only":
                variant = "access_only"
            elif method_key == "tabpfn_access_impact":
                variant = "access_impact"

            return TabPFNPredictor(variant=variant, **kwargs) if variant else TabPFNPredictor(**kwargs)
        except ImportError as err:
            raise ImportError(
                "TabPFN predictor requires optional dependencies. Install via: pip install '.[tabpfn]'"
            ) from err
    else:
        raise ValueError(
            f"Unknown model method '{method}'. Supported methods: 'svm', 'tabpfn', 'tabpfn_access_only', 'tabpfn_access_impact'"
        )


__all__ = [
    "get_default_model_path",
    "get_default_sites_path",
    "get_trainer",
    "get_predictor",
    "SVMTrainer",
    "SVMPredictor",
]
