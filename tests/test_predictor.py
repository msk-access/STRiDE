"""Tests for stride.predictor — feature extraction, model loading, prediction."""

import os

import pandas as pd

from stride.predictor import (
    extract_all_features_from_tsv,
    map_msi_label,
    raise_on_failed_predictions,
    write_one_output_per_sample,
)

# ---------------------------------------------------------------------------
# extract_all_features_from_tsv
# ---------------------------------------------------------------------------


class TestExtractFeatures:
    """Verify per-site feature flattening from TSV."""

    def test_valid_tsv(self, sample_feature_tsv):
        feats = extract_all_features_from_tsv(sample_feature_tsv)
        assert feats is not None
        # Should have site-prefixed feature keys
        assert "chr1_100_p_value" in feats
        assert "chr2_200_entropy_diff" in feats

    def test_empty_file_returns_none(self, empty_feature_tsv):
        result = extract_all_features_from_tsv(empty_feature_tsv)
        assert result is None

    def test_malformed_file_returns_none(self, tmp_dir):
        path = os.path.join(tmp_dir, "malformed.tsv")
        with open(path, "w") as f:
            f.write("col_a\tcol_b\n")
            f.write("1\t2\n")
        result = extract_all_features_from_tsv(path)
        assert result is None


# ---------------------------------------------------------------------------
# map_msi_label
# ---------------------------------------------------------------------------


class TestMapMsiLabel:
    """Verify integer-to-label mapping for MSI status."""

    def test_positive_prediction(self):
        assert map_msi_label(1) == "MSI"

    def test_negative_prediction(self):
        assert map_msi_label(0) == "NA"

    def test_unexpected_value_returns_na(self):
        """Any non-1 value should map to NA."""
        assert map_msi_label(-1) == "NA"
        assert map_msi_label(2) == "NA"


# ---------------------------------------------------------------------------
# write_one_output_per_sample
# ---------------------------------------------------------------------------


class TestWriteOutput:
    """Verify per-sample prediction output writing with MAF-aligned columns."""

    def test_writes_files(self, tmp_dir):
        df = pd.DataFrame(
            {
                "Tumor_Sample_Barcode": ["S1", "S2"],
                "Matched_Norm_Sample_Barcode": ["N1", "N2"],
                "MSI_class_predicted": ["NA", "MSI"],
                "msi_score": [0.123, 0.987],
            }
        )
        paths = write_one_output_per_sample(df, tmp_dir)
        assert len(paths) == 2
        for p in paths:
            assert os.path.isfile(p)
            assert p.endswith("_msi.txt")

    def test_output_content(self, tmp_dir):
        df = pd.DataFrame(
            {
                "Tumor_Sample_Barcode": ["TestSample"],
                "Matched_Norm_Sample_Barcode": ["NormalSample"],
                "MSI_class_predicted": ["MSI"],
                "msi_score": [0.5],
            }
        )
        paths = write_one_output_per_sample(df, tmp_dir)
        result = pd.read_csv(paths[0], sep="\t")
        assert result["MSI_class_predicted"].iloc[0] == "MSI"
        assert result["Tumor_Sample_Barcode"].iloc[0] == "TestSample"
        assert result["Matched_Norm_Sample_Barcode"].iloc[0] == "NormalSample"
        assert result["msi_score"].iloc[0] == 0.5

    def test_output_na_status(self, tmp_dir):
        """Verify NA status is written correctly."""
        df = pd.DataFrame(
            {
                "Tumor_Sample_Barcode": ["Sample01"],
                "Matched_Norm_Sample_Barcode": [""],
                "MSI_class_predicted": ["NA"],
                "msi_score": [-0.3],
            }
        )
        paths = write_one_output_per_sample(df, tmp_dir)
        # keep_default_na=False prevents pandas from parsing "NA" as NaN
        result = pd.read_csv(paths[0], sep="\t", keep_default_na=False)
        assert result["MSI_class_predicted"].iloc[0] == "NA"

    def test_creates_output_dir(self, tmp_dir):
        nested = os.path.join(tmp_dir, "nested", "dir")
        df = pd.DataFrame(
            {
                "Tumor_Sample_Barcode": ["S1"],
                "Matched_Norm_Sample_Barcode": [""],
                "MSI_class_predicted": ["NA"],
                "msi_score": [0.1],
            }
        )
        paths = write_one_output_per_sample(df, nested)
        assert os.path.isdir(nested)
        assert len(paths) == 1


# ---------------------------------------------------------------------------
# get_predictor factory
# ---------------------------------------------------------------------------


class TestGetPredictorFactory:
    """Verify get_predictor creates appropriate predictor instances."""

    def test_get_svm_predictor(self):
        from stride.models import SVMPredictor, get_predictor

        pred = get_predictor("svm")
        assert isinstance(pred, SVMPredictor)

    def test_get_tabpfn_predictor_variants(self, monkeypatch):
        import pytest

        from stride.models import get_predictor

        model_file = "/models/SPECIAL_k2__ed__te_0.6975.joblib"
        try:
            from stride.models.tabpfn.predictor import TabPFNPredictor

            # Mock _load_components so we don't require heavy torch/tabpfn weights in fast unit test
            monkeypatch.setattr(TabPFNPredictor, "_load_components", lambda self: None)

            pred_default = get_predictor("tabpfn", model_path=model_file)
            assert isinstance(pred_default, TabPFNPredictor)
            assert pred_default.threshold == 0.6975

            pred_acc = get_predictor("tabpfn_access_only", model_path=model_file)
            assert isinstance(pred_acc, TabPFNPredictor)
            assert pred_acc.variant == "access_only"

            pred_imp = get_predictor("tabpfn_access_impact", model_path=model_file)
            assert isinstance(pred_imp, TabPFNPredictor)
            assert pred_imp.variant == "access_impact"
        except ImportError:
            with pytest.raises(ImportError, match="TabPFN"):
                get_predictor("tabpfn", model_path=model_file)

    def test_tabpfn_predictor_requires_threshold(self, monkeypatch):
        import pytest

        predictor_mod = pytest.importorskip("stride.models.tabpfn.predictor")
        monkeypatch.setattr(predictor_mod.TabPFNPredictor, "_load_components", lambda self: None)

        with pytest.raises(ValueError, match="--threshold"):
            predictor_mod.TabPFNPredictor(model_path="/models/my_model.joblib")

        pred = predictor_mod.TabPFNPredictor(model_path="/models/my_model.joblib", threshold=0.6)
        assert pred.threshold == 0.6


# ---------------------------------------------------------------------------
# TabPFN registry
# ---------------------------------------------------------------------------


class TestTabPFNRegistry:
    """Verify model-file resolution and calibrated threshold lookup."""

    def test_missing_model_names_expected_file(self):
        import pytest

        from stride.models.tabpfn.registry import resolve_tabpfn_model

        with pytest.raises(FileNotFoundError, match="SPECIAL_k2__ed__te_0.6975.joblib"):
            resolve_tabpfn_model("ao_top1", None)
        with pytest.raises(FileNotFoundError, match="SPECIAL_k3__wd__te__ne_0.6853.joblib"):
            resolve_tabpfn_model(None, None, cohort="access_impact")

    def test_threshold_from_manifest_filename(self):
        from stride.models.tabpfn.registry import resolve_tabpfn_model

        path, thr, info = resolve_tabpfn_model(
            None, "/anywhere/SPECIAL_k3__wd__te__ne_0.6853.joblib"
        )
        assert str(path) == "/anywhere/SPECIAL_k3__wd__te__ne_0.6853.joblib"
        assert thr == 0.6853
        assert info.model_id == "ai_top1"

    def test_threshold_fallbacks(self):
        from stride.models.tabpfn.registry import BUNDLED_MODEL_THRESHOLD, resolve_tabpfn_model

        assert resolve_tabpfn_model(None, "/m/tabpfn_finetuned.joblib")[1] == BUNDLED_MODEL_THRESHOLD
        assert resolve_tabpfn_model(None, "/m/custom_0.42.joblib")[1] == 0.42
        assert resolve_tabpfn_model(None, "/m/custom.joblib")[1] is None


# ---------------------------------------------------------------------------
# raise_on_failed_predictions
# ---------------------------------------------------------------------------


class TestRaiseOnFailedPredictions:
    """Verify per-sample prediction failures stop the run."""

    def test_all_succeeded(self):
        df = pd.DataFrame({"sample_id": ["S1", "S2"], "prediction": ["MSI", "MSS"]})
        raise_on_failed_predictions(df)

    def test_failure_raises_with_sample_and_message(self):
        import pytest

        df = pd.DataFrame(
            {"sample_id": ["S1", "S2"], "prediction": ["MSS", "ERROR: model exploded"]}
        )
        with pytest.raises(RuntimeError, match="S2: ERROR: model exploded"):
            raise_on_failed_predictions(df)
