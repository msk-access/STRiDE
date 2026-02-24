# Quick Start

## Single Sample

Run the full pipeline on one tumor/normal BAM pair:

```bash
stride run \
    --tumor-bam /path/to/tumor.bam \
    --normal-bam /path/to/normal.bam \
    --out-dir results/ \
    --verbose
```

This will:

1. Extract features from 170 MSI loci
2. Predict MSI-H vs MSS
3. Write prediction to `results/predictions/<sample>_msi.txt`

## Step-by-Step

### 1. Generate Features Only

```bash
stride features \
    --tumor-bam tumor.bam \
    --normal-bam normal.bam \
    --out-dir features_output/
```

### 2. Predict from Features

```bash
stride predict \
    --features-dir features_output/features/ \
    --out-dir predictions/
```

## Custom Model or Site List

Override the bundled defaults:

```bash
stride run \
    --tumor-bam tumor.bam \
    --normal-bam normal.bam \
    --model-joblib /path/to/custom_model.joblib \
    --site-list /path/to/custom_sites.txt \
    --out-dir results/
```
