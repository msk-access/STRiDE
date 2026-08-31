# Prediction

## Overview

The `stride predict` command loads a trained model and predicts MSI status from pre-computed feature TSVs.

## Usage

```bash
# Predict using SVM model
stride predict \
    --model svm \
    --features-dir output/features/ \
    --out-dir predictions/

# Predict using TabPFN model
stride predict \
    --model tabpfn \
    --features-dir output/features/ \
    --out-dir predictions/
```

## Input Sources

You can provide feature TSVs via multiple methods:

```bash
# From a directory (recursive search for *.tsv)
stride predict --model svm --features-dir features/ --out-dir predictions/

# Explicit file list
stride predict --model svm --feature-files sample1.tsv --feature-files sample2.tsv --out-dir predictions/

# From a list file (one path per line)
stride predict --model svm --list-file feature_paths.txt --out-dir predictions/

# With matched normal barcode (MAF-standard)
stride predict --model svm --features-dir features/ --out-dir predictions/ \
    --matched-norm-sample-barcode P-0001234-N01
```

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `--model` | `svm` | Model architecture choice (`svm` or `tabpfn`) |
| `--model-path` / `--model-joblib` | bundled model | Custom trained model path |
| `--features-dir` | — | Directory with feature TSVs |
| `--feature-files` | — | One or more explicit TSV paths |
| `--list-file` | — | File listing TSV paths |
| `--out-dir` | *required* | Output directory |
| `--matched-norm-sample-barcode` | `""` | Matched normal sample barcode (MAF-standard) |
| `--verbose / -V` | off | Enable debug logging |


## Output

Each sample produces a `<tumor_barcode>_msi.txt` file containing:

| Column | Description |
|--------|-------------|
| `Tumor_Sample_Barcode` | Tumor sample identifier (MAF-standard) |
| `Matched_Norm_Sample_Barcode` | Matched normal identifier (empty if not provided) |
| `MSI_class_predicted` | `MSI` = Instable, `NA` = Not MSI |
| `msi_score` | Continuous decision score |
