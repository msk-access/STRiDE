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
| `--report / --no-report` | `True` | Generate interactive HTML performance report |
| `--verbose / -V` | off | Enable debug logging |

## Output

Training produces artifacts in the output directory:

```
trained_model/
├── model.joblib              # Trained classifier (ready to use with stride predict)
├── scaler.joblib             # Feature scaler (automatically included in model.joblib)
├── metadata.json             # Training configuration, data splits, classes
└── training_report.html      # Interactive performance report (default)
```

## Training Performance Report

By default, `stride train` generates an **interactive HTML report** (`training_report.html`) with comprehensive model evaluation metrics and visualizations.

### Report Sections

**Model Overview**
- Algorithm & hyperparameters
- Training data composition (samples per class)
- Cross-validation strategy
- Threshold optimization details

**Performance Metrics**
- Overall accuracy, precision, recall, F1-score
- Per-class metrics (MSI-H vs MSS)
- Confusion matrix (absolute counts + normalized)

**Visualizations**
- **ROC Curve**: Model discrimination ability across all decision thresholds
- **Precision-Recall Curve**: Trade-offs between precision and recall
- **Cross-Validation Scores**: Per-fold performance distribution
- **Threshold Optimization**: Specificity vs Sensitivity across thresholds
- **Class Distribution**: Training data composition (MSI vs MSS)

**Feature Analysis**
- Feature importance (for SVM: coefficients; for TabPFN: permutation importance)
- Feature statistics per class (mean, std, distribution)
- Top discriminative features

**Model Comparison** (when multiple folds exist)
- Fold-wise performance heatmap
- Standard deviation of metrics across folds
- Confidence intervals for predictions

### Example Report

```bash
stride train \
    --method svm \
    --access-msi-dir /data/msi_features/ \
    --access-mss-dir /data/mss_features/ \
    --out-dir trained_model/ \
    --report
```

Opens in browser:
```
trained_model/training_report.html
```

The report is **interactive** with:
- Downloadable charts
- Hover tooltips with exact values
- Filterable metric tables
- Export to PNG/SVG

### Disable Report Generation

To skip report generation (useful for batch training or resource constraints):

```bash
stride train \
    --method svm \
    --access-msi-dir /data/msi/ \
    --access-mss-dir /data/mss/ \
    --no-report
```

This saves computation time but still generates the trained model.
