#!/bin/bash
# ==============================================================================
# STRiDE Nextflow Pipeline — Cluster Job Submission Script (SLURM)
# ==============================================================================
#SBATCH --job-name=STRIDE_NF
#SBATCH --partition=cmobic_cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=stride_nf_%j.out
#SBATCH --error=stride_nf_%j.err

set -euo pipefail

echo "================================================================================"
echo "STRiDE Nextflow Pipeline — Head Job"
echo "Job ID       : ${SLURM_JOB_ID:-N/A}"
echo "Node         : $(hostname)"
echo "Submit Dir   : ${SLURM_SUBMIT_DIR:-$(pwd)}"
echo "Date / Time  : $(date)"
echo "================================================================================"

# ── 1. Environment Setup ──────────────────────────────────────────────────────
# Ensure Nextflow, Java, and Singularity are accessible
module purge || true
module load singularity || true
module load java || true
set +u
source ~/.bashrc 2>/dev/null || true
# If using a conda/mamba environment where nextflow & stride are installed:
micromamba activate "${CONDA_ENV_NAME:-tabpfn_msi2}" 2>/dev/null || \
conda activate "${CONDA_ENV_NAME:-tabpfn_msi2}" 2>/dev/null || true
set -u

# ── 2. Pipeline Parameters ────────────────────────────────────────────────────
# Override these via environment variables or CLI arguments
INPUT_CSV="${INPUT_CSV:-tests/data/qc/test_samplesheet.csv}"
OUTDIR="${OUTDIR:-results_nf_${SLURM_JOB_ID:-test}}"
MODEL="${MODEL:-svm}"                  # svm, tabpfn, tabpfn_access_only, tabpfn_access_impact
PROFILE="${PROFILE:-slurm}"            # slurm, local, docker, singularity
GENERATE_QC="${GENERATE_QC:-true}"
EXPLAIN="${EXPLAIN:-true}"
SHAPIQ_BUDGET="${SHAPIQ_BUDGET:-128}"

echo "Parameters:"
echo "  Input Sheet    : $INPUT_CSV"
echo "  Output Dir     : $OUTDIR"
echo "  Model          : $MODEL"
echo "  Profile        : $PROFILE"
echo "  Generate QC    : $GENERATE_QC"
echo "  Explain (Shap) : $EXPLAIN (budget: $SHAPIQ_BUDGET)"
echo "================================================================================"

# ── 3. Run Nextflow ───────────────────────────────────────────────────────────
nextflow run nextflow/main.nf \
    --input "$INPUT_CSV" \
    --outdir "$OUTDIR" \
    --model "$MODEL" \
    --generate_qc "$GENERATE_QC" \
    --explain "$EXPLAIN" \
    --shapiq_budget "$SHAPIQ_BUDGET" \
    -profile "$PROFILE" \
    -resume

echo "================================================================================"
echo "Nextflow Pipeline Finished at: $(date)"
echo "Output Directory: $OUTDIR"
echo "================================================================================"
