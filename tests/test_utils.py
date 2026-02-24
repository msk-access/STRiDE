"""Tests for stride.utils — logging, filename helpers, sample-list parsing."""

import logging
import os

import pytest

from stride.utils import read_samples_list, safe_name, setup_logging, strip_ext

# ---------------------------------------------------------------------------
# setup_logging
# ---------------------------------------------------------------------------


class TestSetupLogging:
    """Verify Rich-based logging configuration."""

    def test_default_level_is_info(self):
        setup_logging(verbose=False)
        assert logging.getLogger().level == logging.INFO

    def test_verbose_sets_debug(self):
        setup_logging(verbose=True)
        assert logging.getLogger().level == logging.DEBUG

    def test_idempotent(self):
        """Calling setup_logging twice should not crash (force=True)."""
        setup_logging(verbose=False)
        setup_logging(verbose=True)
        assert logging.getLogger().level == logging.DEBUG


# ---------------------------------------------------------------------------
# safe_name
# ---------------------------------------------------------------------------


class TestSafeName:
    """Verify filename sanitisation."""

    def test_passthrough_alphanum(self):
        assert safe_name("sample_01") == "sample_01"

    def test_replaces_spaces(self):
        assert safe_name("my sample") == "my_sample"

    def test_replaces_slashes(self):
        assert safe_name("path/to/file") == "path_to_file"

    def test_preserves_dots_hyphens(self):
        assert safe_name("v1.2-beta") == "v1.2-beta"

    def test_empty_string(self):
        assert safe_name("") == ""


# ---------------------------------------------------------------------------
# strip_ext
# ---------------------------------------------------------------------------


class TestStripExt:
    """Verify alignment extension stripping."""

    def test_bam_extension(self):
        assert strip_ext("/data/tumor.bam") == "tumor"

    def test_cram_extension(self):
        assert strip_ext("normal.cram") == "normal"

    def test_sam_extension(self):
        assert strip_ext("reads.sam") == "reads"

    def test_non_alignment_fallback(self):
        assert strip_ext("features.tsv") == "features"

    def test_no_extension(self):
        assert strip_ext("filename") == "filename"

    def test_full_path(self):
        assert strip_ext("/long/path/to/sample_tumor.bam") == "sample_tumor"


# ---------------------------------------------------------------------------
# read_samples_list
# ---------------------------------------------------------------------------


class TestReadSamplesList:
    """Verify CSV/TSV sample manifest parsing."""

    def test_tsv_with_standard_headers(self, tmp_dir):
        path = os.path.join(tmp_dir, "samples.tsv")
        with open(path, "w") as f:
            f.write("sample_id\ttumor_bam\tnormal_bam\n")
            f.write("S1\t/data/S1_T.bam\t/data/S1_N.bam\n")
            f.write("S2\t/data/S2_T.bam\t/data/S2_N.bam\n")

        result = read_samples_list(path)
        assert len(result) == 2
        assert result[0]["sample_id"] == "S1"
        assert result[1]["tumor_bam"] == "/data/S2_T.bam"

    def test_csv_with_alias_headers(self, tmp_dir):
        path = os.path.join(tmp_dir, "samples.csv")
        with open(path, "w") as f:
            f.write("sample,tumor,normal\n")
            f.write("P1,t1.bam,n1.bam\n")

        result = read_samples_list(path)
        assert len(result) == 1
        assert result[0]["sample_id"] == "P1"

    def test_skips_nan_rows(self, tmp_dir):
        path = os.path.join(tmp_dir, "samples.tsv")
        with open(path, "w") as f:
            f.write("sample_id\ttumor_bam\tnormal_bam\n")
            f.write("S1\tt.bam\tn.bam\n")
            f.write("\t\t\n")  # blank row

        result = read_samples_list(path)
        assert len(result) == 1

    def test_raises_on_missing_columns(self, tmp_dir):
        path = os.path.join(tmp_dir, "bad.tsv")
        with open(path, "w") as f:
            f.write("name\tpath\n")
            f.write("X\ty\n")

        with pytest.raises(ValueError, match="sample_id"):
            read_samples_list(path)
