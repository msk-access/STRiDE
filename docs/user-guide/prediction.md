# Prediction

## Overview

The `stride predict` command loads a trained model and predicts MSI-H vs MSS status from pre-computed feature TSVs.

## Usage

```bash
stride predict \
    --features-dir output/features/ \
    --out-dir predictions/
```

## Input Sources

You can provide feature TSVs via multiple methods:

```bash
# From a directory (recursive search for *.tsv)
stride predict --features-dir features/ --out-dir predictions/

# Explicit file list
stride predict --feature-files sample1.tsv --feature-files sample2.tsv --out-dir predictions/

# From a list file (one path per line)
stride predict --list-file feature_paths.txt --out-dir predictions/
```

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `--model-joblib` | bundled SGD model | Custom trained model path |
| `--features-dir` | — | Directory with feature TSVs |
| `--feature-files` | — | One or more explicit TSV paths |
| `--list-file` | — | File listing TSV paths |
| `--out-dir` | *required* | Output directory |
| `--verbose / -V` | off | Enable debug logging |

## Output

Each sample produces a `<sample_id>_msi.txt` file containing:

| Column | Description |
|--------|-------------|
| `Sample_ID` | Sample identifier |
| `MSI_class_predicted` | `0` = MSS, `1` = MSI-H |
| `msi_score` | Continuous decision score |
