// ──────────────────────────────────────────────────────────────────────────
// STRIDE_QC — Standalone module for QC report generation
//
// Wraps `stride qc` for users who want to (re-)generate QC reports
// from existing feature TSVs and prediction files.
//
// Input:  Feature TSV, optional prediction TXT
// Output: Interactive HTML QC report
// ──────────────────────────────────────────────────────────────────────────

process STRIDE_QC {
    tag "$meta.id"
    label 'process_low'

    publishDir "${params.outdir}/qc", mode: params.publish_dir_mode

    container "ghcr.io/msk-access/stride:${params.stride_version}"

    input:
    tuple val(meta), path(features_tsv), path(prediction_txt)

    output:
    tuple val(meta), path("*_qc.html"),    emit: qc_reports
    tuple val(meta), path("*_drivers.tsv"), emit: drivers, optional: true
    path "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Prediction file is optional — pass only if real file
    def pred_arg = prediction_txt.name != 'NO_FILE'
        ? "--prediction ${prediction_txt}"
        : ''

    // Explainability flag
    def explain_flag = params.explain ? '--explain' : '--no-explain'

    """
    echo "── STRIDE_QC ────────────────────────────────────"
    echo "Sample:     ${prefix}"
    echo "Features:   ${features_tsv}"
    echo "Prediction: ${prediction_txt}"
    echo "Explain:    ${params.explain}"
    echo "────────────────────────────────────────────────"

    stride qc \\
        --feature-tsv ${features_tsv} \\
        ${pred_arg} \\
        --output '${prefix}_qc.html' \\
        ${explain_flag} \\
        ${args}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stride: \$(stride --version | sed 's/stride //')
    END_VERSIONS
    """
}
