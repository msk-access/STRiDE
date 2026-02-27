// ──────────────────────────────────────────────────────────────────────────
// STRIDE_FEATURES — Standalone module for MSI feature extraction
//
// Wraps `stride features` for users who want to run feature extraction
// independently of prediction. Useful for:
//   - Pre-computing features for batch prediction
//   - Feature inspection before running the full pipeline
//
// Input:  Tumor/normal BAM pair with indexes, MSI site list
// Output: Feature TSV with per-locus statistics
// ──────────────────────────────────────────────────────────────────────────

process STRIDE_FEATURES {
    tag "$meta.id"
    label 'process_medium'

    publishDir "${params.outdir}/features", mode: params.publish_dir_mode

    container "ghcr.io/msk-access/stride:${params.stride_version}"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai),
                     path(normal_bam), path(normal_bai)
    path site_list

    output:
    tuple val(meta), path("features/*.tsv"), emit: features
    path "versions.yml",                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Site list: pass only if real file
    def site_arg = site_list.name != 'NO_FILE' ? "--site-list ${site_list}" : ''

    """
    echo "── STRIDE_FEATURES ─────────────────────────────"
    echo "Sample:     ${prefix}"
    echo "Tumor BAM:  ${tumor_bam}"
    echo "Normal BAM: ${normal_bam}"
    echo "────────────────────────────────────────────────"

    stride features \\
        --tumor-bam  ${tumor_bam} \\
        --normal-bam ${normal_bam} \\
        --out-dir    . \\
        --sample-id  '${prefix}' \\
        --min-coverage ${params.min_coverage} \\
        --max-repeat-bins ${params.max_repeat_bins} \\
        ${site_arg} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stride: \$(stride --version | sed 's/stride //')
    END_VERSIONS
    """
}
