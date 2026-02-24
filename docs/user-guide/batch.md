# Batch Processing

## Overview

Process multiple samples in a single `stride run` invocation using a sample manifest.

## Sample Manifest Format

Create a CSV or TSV file with columns:

| Column | Aliases | Description |
|--------|---------|-------------|
| `sample_id` | `sample`, `id` | Unique sample identifier |
| `tumor_bam` | `tumor`, `tumor_path` | Path to tumor BAM |
| `normal_bam` | `normal`, `normal_path` | Path to normal BAM |

Example `samples.csv`:

```csv
sample_id,tumor_bam,normal_bam
PATIENT_001,/data/P001_T.bam,/data/P001_N.bam
PATIENT_002,/data/P002_T.bam,/data/P002_N.bam
PATIENT_003,/data/P003_T.bam,/data/P003_N.bam
```

## Usage

```bash
stride run \
    --samples-list samples.csv \
    --out-dir batch_results/ \
    --verbose
```

## Output Structure

```
batch_results/
├── features/
│   ├── msi_features_PATIENT_001.tsv
│   ├── msi_features_PATIENT_002.tsv
│   └── msi_features_PATIENT_003.tsv
└── predictions/
    ├── PATIENT_001_msi.txt
    ├── PATIENT_002_msi.txt
    └── PATIENT_003_msi.txt
```

## Cleanup

Remove intermediate feature TSVs after prediction:

```bash
stride run --samples-list samples.csv --out-dir results/ --delete-features
```
