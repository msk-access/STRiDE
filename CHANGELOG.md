# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Typer CLI with sub-commands (`features`, `predict`, `run`)
- Rich logging with `RichHandler` and progress bars
- Bundled default MSI model and 170-site locus list
- `pyproject.toml` entry points, metadata, and dev/docs extras
- `.gitignore` and `CHANGELOG.md`
- Context manager support for `MSIProfileGenerator` (fixes BAM file handle leak)
- **QC HTML report** — interactive, self-contained report (`stride qc`)
  - Three-tab layout: Dashboard, Site Explorer, Data Table
  - Dashboard with 8 interactive Plotly cards (waterfalls, volcano plots, distance correlation, distributions, entropy scatter, quality violins, insert-size histogram)
  - Site Explorer: searchable combobox, ‹/› navigation, raw/normalized toggle, dynamic x-axis, quality badges, per-locus stats
  - Data Table powered by Tabulator.js with frozen columns, progress-bar formatters, filters, sort, CSV export
  - Light / dark theme toggle (persists across tabs)
  - `--generate-qc` flag for `stride run` to auto-generate reports


### Changed
- Migrated from `argparse` scripts to Typer CLI (all options, no positional args)
- Replaced `print()` statements with structured `logging` in `predictor.py`
- Moved `logging.basicConfig()` from module level to CLI layer

## [0.1.0] - 2026-02-23

### Added
- Initial STRiDE pipeline with MSI feature extraction and prediction
- `MSIProfileGenerator` for repeat counting and statistical features
- Batch and single-sample prediction modes
- 170 curated MSI loci site list
