# STRiDE

**Microsatellite Instability prediction for MSK-ACCESS cfDNA sequencing.**

STRiDE extracts repeat-frequency features from paired tumor/normal BAMs at 170 curated microsatellite loci and classifies samples as **MSI** or **MSS** using trained machine learning models (**SVM** and **TabPFN**).

---

## Features

- **Multi-Model Support**: Predict and train using **SVM** or **TabPFN** classifiers (`--model svm` or `--model tabpfn`).
- **Feature Extraction**: Extract loci repeat frequency features and distance metrics (Wasserstein distance, L1, L2, JS divergence, entropy difference) directly from paired BAM files.
- **Interactive HTML QC Dashboards**: Generate self-contained, interactive HTML quality-control reports with Plotly and Tabulator.js data tables.
- **Typer CLI**: Simple, explicit command-line interface with sub-commands for every step of the workflow.

---

## Quick Installation

Within your Python / Conda environment:

```bash
git clone https://github.com/msk-access/STRiDE.git
cd STRiDE
pip install -e '.[all]'
```

*Optional extras:*
- `pip install -e '.[qc]'` — Install interactive HTML report generator dependencies.
- `pip install -e '.[tabpfn]'` — Install PyTorch & TabPFN dependencies.
- `pip install -e '.[all]'` — Install all core, QC, and TabPFN dependencies.

---

## CLI Usage Overview

Verify the installation:
```bash
stride --help
```

### 1. Extract Features from BAMs (`stride features`)
```bash
stride features \
    --tumor-bam sample_tumor.bam \
    --normal-bam sample_normal.bam \
    --out-dir output/
```

### 2. Predict MSI Status (`stride predict`)
Predict using the default or fine-tuned model (**SVM** or **TabPFN**):
```bash
# Predict using SVM model
stride predict --model svm --features-dir output/features/ --out-dir output/predictions/

# Predict using TabPFN model
stride predict --model tabpfn --features-dir output/features/ --out-dir output/predictions/
```

### 3. Train a New Model (`stride train`)
Train an SVM or TabPFN model on cohort TSV features:
```bash
stride train \
    --method svm \
    --access-msi-dir /path/to/access_msi \
    --access-mss-dir /path/to/access_mss \
    --out-dir trained_svm/
```

### 4. Generate Interactive HTML QC Report (`stride qc`)
Create an interactive HTML QC report for clinical review:
```bash
stride qc \
    --feature-tsv output/features/msi_features_sample.tsv \
    --prediction output/predictions/sample_prediction.txt \
    --output sample_qc_report.html
```

### 5. End-to-End Run (`stride run`)
Run feature extraction, prediction, and optional interactive QC in a single command:
```bash
stride run \
    --model svm \
    --tumor-bam sample_tumor.bam \
    --normal-bam sample_normal.bam \
    --out-dir output/ \
    --generate-qc
```

---

## Disclaimer

This pipeline is intended for research use within MSK-ACCESS and has not been validated for clinical deployment. Use and interpretation of results should follow institutional guidelines.

