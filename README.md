# STRiDE MSI Pipeline

STRiDE is a modular, end-to-end pipeline for Microsatellite Instability (MSI) calling from paired tumor/normal BAM files using MSK-ACCESS data.

This repository supports MSI feature generation from microsatellite loci, MSI/MSS prediction using a trained machine learning model, single-sample and batch execution, and generation of one prediction output file per sample.

## Overview

The pipeline is intentionally split into reusable core logic and user-facing scripts.

The src/ directory contains reusable MSI analysis and prediction code.  
The scripts/ directory contains command-line entry points that users actually run.

This separation improves maintainability, reproducibility, and scalability, and allows future extensions such as fragmentomics, cfRNA, or alternative specimen types without rewriting user interfaces.

## Repository Structure

.  
├── src/  
│   └── stride/  
│       ├── __init__.py  
│       ├── feature_generator.py   # MSI feature extraction from BAMs  
│       ├── predictor.py           # MSI/MSS prediction from feature TSVs  
│       ├── pipeline.py            # End-to-end orchestration  
│       └── utils.py               # Shared helpers  
├── scripts/  
│   ├── stride_features.py         # Feature generation only  
│   ├── stride_predict.py          # Prediction from precomputed features  
│   └── stride_run.py              # End-to-end (features + prediction)  
├── pyproject.toml  
├── environment.yml  
└── README.md  

## Installation

First, create the conda environment:

conda env create -f environment.yml  
conda activate STRiDE  

Then install the package in editable mode from the repository root:

pip install -e .  

This ensures that imports resolve cleanly, scripts run without needing to set PYTHONPATH, and code changes are picked up automatically.

## Inputs

Required inputs include a microsatellite site list (TSV), a tumor BAM, a normal BAM, and a trained MSI model saved as a .joblib file.

For batch mode, provide a TSV or CSV file containing at least the following columns:

sample_id    tumor_bam    normal_bam  

Example:

P001    /path/P001_tumor.bam    /path/P001_normal.bam  
P002    /path/P002_tumor.bam    /path/P002_normal.bam  

Column names are flexible; variations such as sample, tumor, and normal are also accepted.

## Usage

### Feature generation only

Generate a single feature TSV from tumor/normal BAMs:

python scripts/stride_features.py \
--site-list microsatellites.tsv \
--tumor-bam /path/sample_tumor.bam \
--normal-bam /path/sample_normal.bam \
--out-dir out  

This writes:

out/features/msi_features_<sample_id>.tsv  

### Prediction only (evaluate features)

Predict MSI/MSS from one or more feature TSVs:

python scripts/stride_predict.py \
--model-joblib model.joblib \
--features-dir out/features \
--out-dir out/predictions  

This writes one file per sample to:

out/predictions/<sample_id>_msi.txt  

### End-to-end single sample

This mode runs both feature generation and MSI prediction for a single sample.

python scripts/stride_run.py \
--site-list microsatellites.tsv \
--model-joblib model.joblib \
--tumor-bam /path/sample_tumor.bam \
--normal-bam /path/sample_normal.bam \
--out-dir out \
--keep-features  

Outputs include one prediction file per sample located at:

out/predictions/<sample_id>_msi.txt  

and, if --keep-features is specified, a feature file at:

out/features/msi_features_<sample_id>.tsv  

### End-to-end batch mode

This mode processes multiple samples defined in a samples list file.

python scripts/stride_run.py \
--site-list microsatellites.tsv \
--model-joblib model.joblib \
--samples-list samples.tsv \
--out-dir out \
--keep-features  

This produces one prediction file per sample in out/predictions/, and optionally retains per-sample feature TSV files in out/features/.

## Notes

By default, feature TSV files are deleted after prediction unless --keep-features is specified. Prediction outputs are always written as one file per sample.

The modular design allows feature generation, prediction, and orchestration to be used independently for validation, benchmarking, or future development.

## License and Disclaimer

This pipeline is intended for research use within MSK-ACCESS and has not been validated for clinical deployment. Use and interpretation of results should follow institutional guidelines.
