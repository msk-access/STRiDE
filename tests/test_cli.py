"""Tests for stride CLI — version, help, sub-command structure."""

from typer.testing import CliRunner

from stride.cli import app

runner = CliRunner()


class TestCLIBasics:
    """Verify top-level CLI behaviour."""

    def test_version(self):
        result = runner.invoke(app, ["--version"])
        assert result.exit_code == 0
        assert "stride 0.1.0" in result.output

    def test_help(self):
        result = runner.invoke(app, ["--help"])
        assert result.exit_code == 0
        assert "features" in result.output
        assert "predict" in result.output
        assert "run" in result.output

    def test_no_args_shows_help(self):
        """Typer with no_args_is_help returns exit code 0 or 2 and shows Usage."""
        result = runner.invoke(app, [])
        assert "Usage" in result.output


class TestFeaturesCommand:
    """Verify the 'features' sub-command interface."""

    def test_help(self):
        result = runner.invoke(app, ["features", "--help"])
        assert result.exit_code == 0
        assert "--tumor-bam" in result.output
        assert "--normal-bam" in result.output
        assert "--site-list" in result.output
        assert "--verbose" in result.output

    def test_requires_tumor_and_normal(self):
        result = runner.invoke(app, ["features"])
        assert result.exit_code != 0  # Missing required --tumor-bam


class TestPredictCommand:
    """Verify the 'predict' sub-command interface."""

    def test_help(self):
        result = runner.invoke(app, ["predict", "--help"])
        assert result.exit_code == 0
        assert "--model-joblib" in result.output
        assert "--features-dir" in result.output
        assert "--out-dir" in result.output

    def test_requires_out_dir(self):
        result = runner.invoke(app, ["predict"])
        assert result.exit_code != 0  # Missing required --out-dir


class TestRunCommand:
    """Verify the 'run' sub-command interface."""

    def test_help(self):
        result = runner.invoke(app, ["run", "--help"])
        assert result.exit_code == 0
        assert "--tumor-bam" in result.output
        assert "--samples-list" in result.output
        assert "--model-joblib" in result.output
        assert "--delete-features" in result.output

    def test_requires_out_dir(self):
        result = runner.invoke(app, ["run"])
        assert result.exit_code != 0  # Missing required --out-dir
