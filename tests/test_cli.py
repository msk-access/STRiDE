"""Tests for stride CLI — version, help, sub-command structure."""

import re

from typer.testing import CliRunner

from stride.cli import app

runner = CliRunner()


def _strip(text: str) -> str:
    """Remove ANSI escape codes for clean string matching."""
    return re.sub(r"\x1b\[[0-9;]*m", "", text)


class TestCLIBasics:
    """Verify top-level CLI behaviour."""

    def test_version(self):
        result = runner.invoke(app, ["--version"])
        assert result.exit_code == 0
        assert "stride 0.1.0" in _strip(result.output)

    def test_help(self):
        result = runner.invoke(app, ["--help"])
        assert result.exit_code == 0
        assert "features" in _strip(result.output)
        assert "predict" in _strip(result.output)
        assert "run" in _strip(result.output)

    def test_no_args_shows_help(self):
        """Typer with no_args_is_help returns exit code 0 or 2 and shows Usage."""
        result = runner.invoke(app, [])
        assert "Usage" in _strip(result.output)


class TestFeaturesCommand:
    """Verify the 'features' sub-command interface."""

    def test_help(self):
        result = runner.invoke(app, ["features", "--help"])
        assert result.exit_code == 0
        out = _strip(result.output)
        assert "--tumor-bam" in out
        assert "--normal-bam" in out
        assert "--site-list" in out
        assert "--verbose" in out

    def test_requires_tumor_and_normal(self):
        result = runner.invoke(app, ["features"])
        assert result.exit_code != 0  # Missing required --tumor-bam


class TestPredictCommand:
    """Verify the 'predict' sub-command interface."""

    def test_help(self):
        result = runner.invoke(app, ["predict", "--help"])
        assert result.exit_code == 0
        out = _strip(result.output)
        assert "--model-joblib" in out
        assert "--features-dir" in out
        assert "--out-dir" in out
        assert "--matched-norm-sample-barco" in out

    def test_requires_out_dir(self):
        result = runner.invoke(app, ["predict"])
        assert result.exit_code != 0  # Missing required --out-dir


class TestRunCommand:
    """Verify the 'run' sub-command interface."""

    def test_help(self):
        result = runner.invoke(app, ["run", "--help"])
        assert result.exit_code == 0
        out = _strip(result.output)
        assert "--tumor-bam" in out
        assert "--samples-list" in out
        assert "--model-joblib" in out
        assert "--delete-features" in out
        assert "Matched normal sample" in out

    def test_requires_out_dir(self):
        result = runner.invoke(app, ["run"])
        assert result.exit_code != 0  # Missing required --out-dir
