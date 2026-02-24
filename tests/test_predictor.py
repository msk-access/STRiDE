"""Tests for stride.predictor — feature extraction, model loading, prediction."""

import os

import numpy as np
import pandas as pd
import pytest

from stride.predictor import extract_all_features_from_tsv, write_one_output_per_sample


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
# write_one_output_per_sample
# ---------------------------------------------------------------------------

class TestWriteOutput:
    """Verify per-sample prediction output writing."""

    def test_writes_files(self, tmp_dir):
        df = pd.DataFrame({
            "Sample_ID": ["S1", "S2"],
            "MSI_class_predicted": [0, 1],
            "msi_score": [0.123, 0.987],
        })
        paths = write_one_output_per_sample(df, tmp_dir)
        assert len(paths) == 2
        for p in paths:
            assert os.path.isfile(p)
            assert p.endswith("_msi.txt")

    def test_output_content(self, tmp_dir):
        df = pd.DataFrame({
            "Sample_ID": ["TestSample"],
            "MSI_class_predicted": [1],
            "msi_score": [0.5],
        })
        paths = write_one_output_per_sample(df, tmp_dir)
        result = pd.read_csv(paths[0], sep="\t")
        assert result["MSI_class_predicted"].iloc[0] == 1
        assert result["msi_score"].iloc[0] == 0.5

    def test_creates_output_dir(self, tmp_dir):
        nested = os.path.join(tmp_dir, "nested", "dir")
        df = pd.DataFrame({
            "Sample_ID": ["S1"],
            "MSI_class_predicted": [0],
            "msi_score": [0.1],
        })
        paths = write_one_output_per_sample(df, nested)
        assert os.path.isdir(nested)
        assert len(paths) == 1
