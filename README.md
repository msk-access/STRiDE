# STRiDE

Microsatellite Instability (MSI) prediction pipeline for MSK-ACCESS cfDNA sequencing.

STRiDE extracts repeat-frequency features from paired tumor/normal BAMs across 170 curated microsatellite loci, computes locus interaction values using ShapIQ, and classifies samples as **MSI** or **MSS** using fine-tuned machine learning models (**TabPFN** and **SVM**).

---

## Installation

### Prerequisites
- Python 3.10
- Conda / Micromamba

### Setup
```bash
# Clone the repository
git clone https://github.com/msk-access/STRiDE.git
cd STRiDE

# Create and activate environment
micromamba create -n stride python=3.10 -y
micromamba activate stride

# Install STRiDE with all dependencies (TabPFN, PyTorch, ShapIQ, QC reporting)
pip install -e '.[all]'

# Verify CLI
stride --help
```

---

## Nextflow Pipeline

For multi-sample runs or cluster execution, use the Nextflow workflow (`nextflow/main.nf`).

### 1. Prepare Sample Sheet (`samples.csv`)
Create a CSV file with sample IDs and paths to paired tumor/normal BAM files:

```csv
sample,tumor_bam,normal_bam
P-0080677-T01-XS1,/path/to/P-0080677-T01-XS1-T.bam,/path/to/P-0080677-T01-XS1-N.bam
P-0061710-T03-XS1,/path/to/P-0061710-T03-XS1-T.bam,/path/to/P-0061710-T03-XS1-N.bam
```

### 2. Run Locally
```bash
nextflow run nextflow/main.nf \
    --input samples.csv \
    --outdir results/ \
    --model tabpfn \
    --tabpfn_model ao_top1
```

### 3. Run on Slurm Cluster
```bash
nextflow run nextflow/main.nf \
    -profile slurm \
    --input samples.csv \
    --outdir results/ \
    --model tabpfn \
    --tabpfn_model ao_top1 \
    -resume
```

### Pipeline Parameters
| Parameter | Default | Description |
|:---|:---|:---|
| `--input` | *Required* | Path to sample sheet CSV (`sample,tumor_bam,normal_bam`) |
| `--outdir` | `results` | Output directory |
| `--model` | `tabpfn` | Model architecture: `tabpfn` or `svm` |
| `--tabpfn_model` | `ao_top1` | Pre-trained model preset (e.g. `ao_top1`, `ai_top1`) |
| `--threshold` | `auto` | Decision threshold (`auto` resolves calibrated threshold from manifest) |
| `--qc` | `true` | Generate interactive HTML interpretation dashboards |
| `--explain` | `true` | Compute ShapIQ locus attribution and driver loci |

---

## CLI Usage

### End-to-End Analysis (`stride run`)
Extracts features, runs prediction, and generates the interpretation dashboard in a single step:

```bash
stride run \
    --tumor-bam  sample_tumor.bam \
    --normal-bam sample_normal.bam \
    --sample-name SAMPLE_001 \
    --model tabpfn \
    --tabpfn-model ao_top1 \
    --threshold auto \
    --qc \
    --explain \
    --out-dir results/
```

### Modular Commands
- **Extract features only**:
  ```bash
  stride features --tumor-bam tumor.bam --normal-bam normal.bam --out-dir output/
  ```
- **Predict from features**:
  ```bash
  stride predict --model tabpfn --tabpfn-model ao_top1 --features-dir output/features/ --out-dir output/predictions/
  ```
- **List registered TabPFN models**:
  ```bash
  stride models
  ```
- **Generate standalone QC report**:
  ```bash
  stride qc --feature-tsv output/features/msi_features.tsv --prediction output/predictions/sample_prediction.tsv --output report.html
  ```

---

## Bundled TabPFN Models

| Model ID | Cohort | Rank | Features | Calibrated Threshold | Description |
|:---|:---|:---:|:---|:---:|:---|
| `ao_top1` *(Default)* | Access-Only | 1 | `entropy_diff, tumor_entropy` (ED + TE) | `0.6975` | Access-Only Top Model |
| `ao_top2` | Access-Only | 2 | `tumor_entropy, normal_entropy, n_alleles_diff_norm_6` | `0.5178` | Access-Only Rank 2 |
| `ao_top3` | Access-Only | 3 | `wasserstein_distance, tumor_entropy, normal_entropy` | `0.7213` | Access-Only Rank 3 |
| `ao_top4` | Access-Only | 4 | `wasserstein_distance, tumor_entropy` | `0.7213` | Access-Only Rank 4 |
| `ai_top1` | Access+Impact | 1 | `wasserstein_distance, tumor_entropy, normal_entropy` | `0.6853` | Access+Impact Top Model |
| `ai_top2` | Access+Impact | 2 | `wasserstein_distance, tumor_entropy, n_alleles_diff_norm_4, n_alleles_diff_norm_6` | `0.7319` | Access+Impact Rank 2 |
| `ai_top3` | Access+Impact | 3 | `wasserstein_distance, entropy_diff, tumor_entropy` | `0.7458` | Access+Impact Rank 3 |

---

## Output Files

Executing either Nextflow or `stride run` generates:

- `predictions/{sample}_prediction.tsv`: Final MSI call (`MSI` / `MSS`), probability score, model used, and decision threshold.
- `qc/{sample}_interpretation_reports.html`: Standalone interactive report containing:
  - Header with MSI status, probability score, and calibrated cutoff.
  - Model attribution section with ShapIQ waterfall plot, summary cards, and driver table.
  - Site Explorer ranking all 170 sites by interaction impact with embedded locus repeat distribution grid.
- `qc/{sample}_drivers.tsv`: Ranked tabular driver loci.

---

## Disclaimer

This pipeline is developed for research use within MSK-ACCESS.
