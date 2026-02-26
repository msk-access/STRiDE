"""Tests for stride.predictor — feature extraction, model loading, prediction."""

import os

import pandas as pd

from stride.predictor import (
    extract_all_features_from_tsv,
    map_msi_label,
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
