# STRiDE

**Microsatellite Instability prediction for MSK-ACCESS cfDNA sequencing.**

STRiDE extracts repeat-frequency features from paired tumor/normal BAMs at 170 curated microsatellite loci and classifies samples as **MSI-H** or **MSS** using a trained machine learning model.

## Quick Start

```bash
pip install stride-msk

# Single-sample end-to-end
stride run \
    --tumor-bam /data/tumor.bam \
    --normal-bam /data/normal.bam \
    --out-dir results/

# Batch mode
stride run \
    --samples-list samples.csv \
    --out-dir results/
```

## Features

- 🧬 **170 curated MSI loci** — validated for cfDNA panel sequencing
- 📊 **36+ features per site** — statistical, distance, allele, and quality metrics
- 🤖 **Multiple classifiers** — SGD, SVM, TabPFN models; bundled default model
- 🎨 **Interactive QC dashboard** — visualize feature distributions, site explorer, expert review
- 🚀 **Nextflow DSL2 pipeline** — high-performance HPC execution with resource scheduling
- 🏋️ **Model training** — train custom SVM/TabPFN models on your cohorts
- 📋 **Rich CLI** — progress bars, colored tables, structured logging
- 🐳 **Docker support** — `ghcr.io/msk-access/stride`

## Disclaimer

!!! warning "Research Use Only"
    This pipeline is intended for research use within MSK-ACCESS and has not been validated for clinical deployment. Use and interpretation of results should follow institutional guidelines.
