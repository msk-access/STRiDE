"""STRiDE — MSI prediction pipeline for MSK-ACCESS.

Extracts microsatellite instability features from paired tumor/normal BAMs
and predicts MSI/MSS status using a trained machine learning model.
"""

__version__ = "0.1.0"

from .models import SVMPredictor, SVMTrainer, get_predictor, get_trainer

__all__ = [
    "__version__",
    "get_trainer",
    "get_predictor",
    "SVMTrainer",
    "SVMPredictor",
]
