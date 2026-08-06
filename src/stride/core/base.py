from abc import ABC, abstractmethod
import pandas as pd
from typing import List, Dict, Any, Union, Optional
from pathlib import Path

class BaseTrainer(ABC):
    @abstractmethod
    def train_cv_refit(
        self,
        impact_records: pd.DataFrame,
        access_records: pd.DataFrame,
        test_records: pd.DataFrame,
        metric_pool: List[str],
        outdir: Path,
        **kwargs
    ) -> Dict[str, Any]:
        """
        Train the model using the standard workflow:
        1. Construct matrices for Impact, Access, and Test data.
        2. Fit imputer / preprocessing.
        3. Cross-Validation on Access cohort (with Impact as pretraining/fixed training).
        4. Select classification threshold to meet specificity constraints.
        5. Final refit on combination of Impact & Access training data.
        6. Evaluation on Test data.
        7. Save all trained model components and metadata.
        """
        pass


class BasePredictor(ABC):
    @abstractmethod
    def predict_sample(self, tsv_path: Union[str, Path]) -> Dict[str, Any]:
        """
        Evaluate a single genomic sample TSV file.
        """
        pass

    @abstractmethod
    def predict_batch(
        self,
        tsv_paths: List[Union[str, Path]],
        output_tsv: Optional[Union[str, Path]] = None
    ) -> pd.DataFrame:
        """
        Evaluate multiple sample TSV files in batch.
        """
        pass
