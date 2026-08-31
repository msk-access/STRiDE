// ──────────────────────────────────────────────────────────────────────────
// STRIDE_RUN — Primary module for end-to-end MSI prediction
//
// Wraps `stride run` which performs:
//   1. Feature extraction from tumor/normal BAM pairs (170 MSI loci)
//   2. MSI prediction using a trained SGD model
//   3. Optionally generates an interactive HTML QC report
//
// Inputs:  BAM pairs + index files, MSI site list, trained model
// Outputs: Feature TSV, prediction TXT, optional QC HTML report
//
// Follows nf-core module conventions:
//   - task.ext.when    for conditional execution
//   - task.ext.args    for extra CLI arguments
//   - task.ext.prefix  for output naming
//   - versions.yml     for software version tracking
// ──────────────────────────────────────────────────────────────────────────

process STRIDE_RUN {
    tag "$meta.id"
    label 'process_medium'

    // Output publishing — configurable via params
    publishDir "${params.outdir}/stride", mode: params.publish_dir_mode

    // Container — version-pinnable via params.stride_version
    container "ghcr.io/msk-access/stride:${params.stride_version}"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai),
                     path(normal_bam), path(normal_bai)
    path site_list
    path model_joblib

    output:
    tuple val(meta), path("features/*.tsv"),        emit: features
    tuple val(meta), path("predictions/*_msi.txt"), emit: predictions
    tuple val(meta), path("qc/*_qc.html"),          emit: qc_reports, optional: true
    tuple val(meta), path("qc/*_drivers.tsv"),      emit: drivers,    optional: true
    path "versions.yml",                            emit: versions

    // nf-core conditional execution guard
    when:
    task.ext.when == null || task.ext.when

    script:
    // Resolve process-level overrides (nf-core pattern)
    def args     = task.ext.args ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"

    // Optional: matched-normal sample barcode for MAF-aligned output
    def norm_bc  = meta.matched_norm_sample_barcode
        ? "--matched-norm-sample-barcode '${meta.matched_norm_sample_barcode}'"
        : ''

    // QC report generation is controlled by params.generate_qc
    def qc_flag  = params.generate_qc ? '--generate-qc' : ''

    // ShapIQ model explainability
    def explain_flag = params.explain ? '--explain' : '--no-explain'

    // Site list and model: pass only if real files (not NO_FILE sentinel)
    def site_arg  = site_list.name  != 'NO_FILE' ? "--site-list ${site_list}"    : ''
    def model_arg = model_joblib.name != 'NO_FILE' ? "--model-joblib ${model_joblib}" : ''

    """
    # Log inputs for debugging / Nextflow trace
    echo "── STRIDE_RUN ──────────────────────────────────"
    echo "Sample:     ${prefix}"
    echo "Tumor BAM:  ${tumor_bam}"
    echo "Normal BAM: ${normal_bam}"
    echo "Site list:  ${site_list}"
    echo "Model:      ${model_joblib}"
    echo "QC:         ${params.generate_qc}"
    echo "Explain:    ${params.explain}"
    echo "CPUs:       ${task.cpus}"
    echo "Memory:     ${task.memory}"
    echo "────────────────────────────────────────────────"

    stride run \\
        --model      ${params.model} \\
        --tumor-bam  ${tumor_bam} \\
        --normal-bam ${normal_bam} \\
        --out-dir    . \\
        --sample-id  '${prefix}' \\
        --min-coverage ${params.min_coverage} \\
        --max-repeat-bins ${params.max_repeat_bins} \\
        ${site_arg} \\
        ${model_arg} \\
        ${norm_bc} \\
        ${qc_flag} \\
        ${explain_flag} \\
        ${args}


    # nf-core version tracking
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stride: \$(stride --version | sed 's/stride //')
    END_VERSIONS
    """
}
