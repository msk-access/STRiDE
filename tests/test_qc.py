import os

import pytest

from stride.qc import generate_report, is_qc_available, parse_feature_tsv


@pytest.fixture
def dummy_feature_tsv(tmp_path):
    tsv_content = (
        "chrom\tstart\trepeat_unit\trepeat_count\tl1_distance\tl2_distance\t"
        "wasserstein_distance\tp_value\t"
        "tumor_mapq_mean\tnormal_mapq_mean\ttumor_total_coverage\tnormal_total_coverage\t"
        "entropy_diff\ttumor_entropy\tnormal_entropy\t"
        "tumor_bq_mean\t"
        "tumor_norm_freqs\tnormal_norm_freqs\ttumor_freqs\tnormal_freqs\n"
        "chr1\t1000\tA\t10\t0.5\t0.35\t0.003\t0.001\t"
        "60.0\t60.0\t100\t100\t"
        "0.1\t2.5\t2.4\t"
        "33.0\t"
        "0.1 0.2\t0.1 0.2\t10 20\t10 20\n"
        "chr2\t2000\tT\t15\t1.5\t0.90\t0.007\t0.0001\t"
        "30.0\t55.0\t200\t150\t"
        "0.5\t3.1\t2.6\t"
        "32.0\t"
        "0.5 0.5\t0.2 0.8\t100 100\t30 120\n"
    )
    p = tmp_path / "dummy_features.tsv"
    p.write_text(tsv_content)
    return str(p)


@pytest.mark.skipif(not is_qc_available(), reason="QC dependencies not installed")
def test_parse_feature_tsv(dummy_feature_tsv):
    df = parse_feature_tsv(dummy_feature_tsv)
    assert len(df) == 2
    assert "locus" in df.columns
    assert df.iloc[0]["locus"] == "chr1:1000 (10xA)"

    # Check if numpy arrays were parsed correctly
    assert len(df.iloc[0]["tumor_norm_freqs"]) == 2
    assert df.iloc[0]["tumor_norm_freqs"][0] == 0.1
    assert df.iloc[1]["normal_freqs"][1] == 120.0


@pytest.mark.skipif(not is_qc_available(), reason="QC dependencies not installed")
def test_generate_html_report(dummy_feature_tsv, tmp_path):
    out_html = str(tmp_path / "test_report.html")
    prediction = {"msi_status": "MSI", "msi_score": 0.95}

    generate_report(dummy_feature_tsv, out_html, format="html", prediction_result=prediction)

    assert os.path.exists(out_html)
    content = open(out_html).read()
    assert "STRiDE MSI Quality Control Report" in content
    assert "Site Explorer" in content
    assert "stride-table" in content  # Tabulator data table
