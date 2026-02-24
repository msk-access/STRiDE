# Feature Generation

## Overview

The `stride features` command extracts 36+ statistical features from each of 170 curated microsatellite loci by comparing repeat-frequency distributions between tumor and normal BAMs.

## Usage

```bash
stride features \
    --tumor-bam /path/to/tumor.bam \
    --normal-bam /path/to/normal.bam \
    --out-dir output/ \
    --sample-id PATIENT_001
```

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `--tumor-bam` | *required* | Path to tumor BAM file |
| `--normal-bam` | *required* | Path to normal BAM file |
| `--out-dir` | `out` | Base output directory |
| `--site-list` | bundled 170-site list | Custom MSI loci file |
| `--output-tsv` | auto-generated | Explicit output TSV path |
| `--sample-id` | from BAM filename | Sample ID for output naming |
| `--min-coverage` | `20` | Minimum read coverage per site |
| `--max-repeat-bins` | `100` | Maximum repeat-count bins |
| `--verbose / -V` | off | Enable debug logging |

## Feature Categories

- **Statistical**: chi-squared p-value, Shannon entropy
- **Distance**: L1, L2, Wasserstein distance
- **Allele**: count differences at multiple thresholds
- **Quality**: mapping quality, base quality means
- **Insert size**: mean insert sizes (all/ref/alt reads)
