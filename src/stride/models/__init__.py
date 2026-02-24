"""Bundled model resources for STRiDE.

Provides helper functions to locate the default trained model and
MSI site-list files shipped with the package.
"""

import importlib.resources
import logging

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
