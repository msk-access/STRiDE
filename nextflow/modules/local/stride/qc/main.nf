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
    tuple val(meta), path("*_interpretation_reports.html"), emit: qc_reports
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

    def tabpfn_arg = (params.model.toString().toLowerCase().contains('tabpfn') && params.tabpfn_model) ? "--tabpfn-model '${params.tabpfn_model}'" : ''
    def thresh_arg = params.threshold ? "--threshold ${params.threshold}" : ''
    def model_arg  = (params.model_joblib && params.model_joblib.toString() != 'NO_FILE') ? "--model-joblib ${params.model_joblib}" : ''

    """
    echo "── STRIDE_QC ────────────────────────────────────"
    echo "Sample:     ${prefix}"
    echo "Features:   ${features_tsv}"
    echo "Prediction: ${prediction_txt}"
    echo "Explain:    ${params.explain}"
    echo "────────────────────────────────────────────────"

    stride qc \\
        --model ${params.model} \\
        ${tabpfn_arg} \\
        ${thresh_arg} \\
        ${model_arg} \\
        --feature-tsv ${features_tsv} \\
        ${pred_arg} \\
        --output '${prefix}_interpretation_reports.html' \\
        ${explain_flag} \\
        --shapiq-budget ${params.shapiq_budget} \\
        ${args}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stride: \$(stride --version | sed 's/stride //')
    END_VERSIONS
    """
}
