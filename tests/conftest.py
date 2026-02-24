"""Shared fixtures for STRiDE test suite."""

import os
import tempfile
from unittest.mock import MagicMock

import pytest


@pytest.fixture
def tmp_dir():
    """Provide a temporary directory that is cleaned up after the test."""
    with tempfile.TemporaryDirectory() as d:
        yield d


@pytest.fixture
def sample_sites_file(tmp_dir):
    """Create a minimal MSI sites file for testing."""
    path = os.path.join(tmp_dir, "test_sites.txt")
    with open(path, "w") as f:
        # Header matching the real site list format
        f.write("chrom\tstart\tunit_len\tref_len\ttimes\tend\tname\trepeat_unit\tleft_flank\tright_flank\n")
        f.write("chr1\t100\t1\t1\t15\t115\tsite1\tA\tGCTAGCTA\tTCGATCGA\n")
        f.write("chr2\t200\t2\t2\t10\t220\tsite2\tAT\tACGTACGT\tTGCATGCA\n")
    return path


@pytest.fixture
def sample_feature_tsv(tmp_dir):
    """Create a minimal feature TSV for testing predictor functions."""
    path = os.path.join(tmp_dir, "msi_features_test_sample.tsv")
    with open(path, "w") as f:
        f.write("chrom\tstart\tp_value\tentropy_diff\tl1_distance\n")
        f.write("chr1\t100\t0.05\t0.3\t0.15\n")
        f.write("chr2\t200\t0.01\t0.7\t0.45\n")
    return path


@pytest.fixture
def empty_feature_tsv(tmp_dir):
    """Create an empty TSV for edge-case testing."""
    path = os.path.join(tmp_dir, "empty.tsv")
    open(path, "w").close()
    return path
