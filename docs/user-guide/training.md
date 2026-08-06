# Training

## Overview

The `stride train` command fits a new **SVM** or **TabPFN** classification model on feature TSV cohorts (e.g. ACCESS or IMPACT cohorts with MSI and MSS samples).

## Usage

```bash
# Train an SVM model
stride train \
    --method svm \
    --access-msi-dir /path/to/access_msi/ \
    --access-mss-dir /path/to/access_mss/ \
    --out-dir trained_svm/

# Train a TabPFN model
stride train \
    --method tabpfn \
    --access-msi-dir /path/to/access_msi/ \
    --access-mss-dir /path/to/access_mss/ \
    --out-dir trained_tabpfn/
```

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `--method` / `--model` | `svm` | Model training algorithm (`svm` or `tabpfn`) |
| `--access-msi-dir` | *required* | Directory containing ACCESS MSI feature TSVs |
| `--access-mss-dir` | *required* | Directory containing ACCESS MSS feature TSVs |
| `--impact-msi-dir` | optional | Directory containing IMPACT MSI feature TSVs |
| `--impact-mss-dir` | optional | Directory containing IMPACT MSS feature TSVs |
| `--out-dir` | `trained_model` | Directory to save trained model artifacts |
| `--min-spec` | `0.95` | Target minimum specificity constraint for threshold selection |
| `--cv-folds` | `5` | Number of cross-validation folds |
| `--verbose / -V` | off | Enable debug logging |
