# Installation

## From PyPI

```bash
pip install stride
```

## Development Install

```bash
git clone https://github.com/msk-access/STRiDE.git
cd STRiDE
pip install -e ".[dev]"
```

## Docker

```bash
# Pull from GitHub Container Registry
docker pull ghcr.io/msk-access/stride:latest

# Run
docker run --rm -v /data:/data ghcr.io/msk-access/stride \
    run --tumor-bam /data/tumor.bam --normal-bam /data/normal.bam --out-dir /data/results

# Singularity (HPC)
singularity pull docker://ghcr.io/msk-access/stride:latest
singularity exec stride_latest.sif stride --version
```

## Requirements

- Python ≥ 3.9
- pysam ≥ 0.21.0
- scikit-learn, numpy, pandas, scipy
