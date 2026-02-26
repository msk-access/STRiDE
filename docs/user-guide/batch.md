# Batch Processing

## Overview

Process multiple samples in a single `stride run` invocation using a sample manifest.

## Sample Manifest Format

Create a CSV or TSV file with columns:

| Column | Aliases | Required | Description |
|--------|---------|:--------:|-------------|
| `sample_id` | `sample`, `id` | ✅ | Unique sample identifier (maps to `Tumor_Sample_Barcode` in output) |
| `tumor_bam` | `tumor`, `tumor_path` | ✅ | Path to tumor BAM |
| `normal_bam` | `normal`, `normal_path` | ✅ | Path to normal BAM |
| `matched_norm_sample_barcode` | `normal_barcode`, `normal_sample_barcode`, `matched_normal` | ❌ | Explicit matched-normal barcode. Defaults to normal BAM basename |

Example `samples.csv`:

```csv
sample_id,tumor_bam,normal_bam,matched_norm_sample_barcode
PATIENT_001,/data/P001_T.bam,/data/P001_N.bam,P001-N01-IM6
PATIENT_002,/data/P002_T.bam,/data/P002_N.bam,P002-N01-IM6
PATIENT_003,/data/P003_T.bam,/data/P003_N.bam,P003-N01-IM6
```

If `matched_norm_sample_barcode` is omitted, the normal barcode is auto-derived from the `normal_bam` filename:

```csv
sample_id,tumor_bam,normal_bam
PATIENT_001,/data/P001_T.bam,/data/P001_N.bam
PATIENT_002,/data/P002_T.bam,/data/P002_N.bam
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

Each `_msi.txt` prediction file contains MAF-aligned columns:
`Tumor_Sample_Barcode`, `Matched_Norm_Sample_Barcode`, `MSI_class_predicted`, `msi_score`.

## Cleanup

Remove intermediate feature TSVs after prediction:

```bash
stride run --samples-list samples.csv --out-dir results/ --delete-features
```
