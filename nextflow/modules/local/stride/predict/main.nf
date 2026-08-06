// ──────────────────────────────────────────────────────────────────────────
// STRIDE_PREDICT — Standalone module for MSI prediction
//
// Wraps `stride predict` for users who have pre-computed feature TSVs
// and want to (re-)predict with a different model or parameters.
//
// Input:  Feature TSV, trained model .joblib
// Output: Prediction TXT with MAF-aligned columns
// ──────────────────────────────────────────────────────────────────────────

process STRIDE_PREDICT {
    tag "$meta.id"
    label 'process_low'

    publishDir "${params.outdir}/predictions", mode: params.publish_dir_mode

    container "ghcr.io/msk-access/stride:${params.stride_version}"

    input:
    tuple val(meta), path(features_tsv)
    path model_joblib

    output:
    tuple val(meta), path("*_msi.txt"),  emit: predictions
    tuple val(meta), path(features_tsv), emit: features     // pass-through for downstream
    path "versions.yml",                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Optional: matched-normal barcode for MAF-aligned output
    def norm_bc = meta.matched_norm_sample_barcode
        ? "--matched-norm-sample-barcode '${meta.matched_norm_sample_barcode}'"
        : ''

    // Model: pass only if real file
    def model_arg = model_joblib.name != 'NO_FILE' ? "--model-joblib ${model_joblib}" : ''

    """
    echo "── STRIDE_PREDICT ──────────────────────────────"
    echo "Sample:   ${prefix}"
    echo "Features: ${features_tsv}"
    echo "Model:    ${model_joblib}"
    echo "────────────────────────────────────────────────"

    stride predict \\
        --model ${params.model} \\
        --feature-files ${features_tsv} \\
        --out-dir . \\
        ${model_arg} \\
        ${norm_bc} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stride: \$(stride --version | sed 's/stride //')
    END_VERSIONS
    """
}
