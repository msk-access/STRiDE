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

### With QC Report

Generate an interactive QC dashboard for expert review:

```bash
stride run \
    --tumor-bam /path/to/tumor.bam \
    --normal-bam /path/to/normal.bam \
    --out-dir results/ \
    --qc
```

Output includes `results/qc/<sample>_qc.html` — interactive report with feature distributions, site explorer, and MSI verdict.

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

## Advanced Features

### Train a Custom Model

Train SVM or TabPFN models on labeled feature cohorts:

```bash
# Train an SVM model
stride train \
    --method svm \
    --access-msi-dir /path/to/msi_features/ \
    --access-mss-dir /path/to/mss_features/ \
    --out-dir my_svm_model/

# Train TabPFN model (state-of-the-art)
stride train \
    --method tabpfn \
    --access-msi-dir /path/to/msi_features/ \
    --access-mss-dir /path/to/mss_features/ \
    --out-dir my_tabpfn_model/

# Use your trained model
stride run \
    --tumor-bam tumor.bam \
    --normal-bam normal.bam \
    --model-joblib my_svm_model/model.joblib \
    --out-dir results/
```

See [Model Training Guide](../user-guide/training.md) for details.

### Generate QC Reports

Create interactive dashboards for expert review:

```bash
# Standalone QC report
stride qc \
    --feature-tsv features/msi_features_sample.tsv \
    --prediction predictions/sample_msi.txt \
    --output report.html

# Or generate during pipeline run
stride run --tumor-bam tumor.bam --normal-bam normal.bam --qc --out-dir results/
```

See [QC Report Guide](../user-guide/qc.md) for details.

### Nextflow on HPC

Process large batches on high-performance clusters:

```bash
nextflow run stride/main.nf -profile slurm \
    --samples_list samples.csv \
    --out_dir results/ \
    --max_cpus 16
```

See [Nextflow Guide](../user-guide/nextflow.md) for details.
