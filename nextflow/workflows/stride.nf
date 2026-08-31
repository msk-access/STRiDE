// ──────────────────────────────────────────────────────────────────────────
// STRIDE Workflow — Named workflow for the STRiDE MSI pipeline
//
// Chains the STRIDE_RUN module and collects version information.
// This follows the py-gbcms pattern of separating workflow logic
// from the entry-point main.nf.
//
// Usage:
//   include { STRIDE } from './workflows/stride'
//   STRIDE ( ch_samples, ch_site_list, ch_model )
// ──────────────────────────────────────────────────────────────────────────

include { STRIDE_RUN } from '../modules/local/stride/run/main'

workflow STRIDE {
    take:
    ch_samples     // channel: [ val(meta), path(tumor_bam), path(tumor_bai),
                   //                       path(normal_bam), path(normal_bai) ]
    ch_site_list   // path: MSI site list TSV
    ch_model       // path: Trained model .joblib

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Run stride end-to-end (features → predict → optional QC)
    //
    STRIDE_RUN ( ch_samples, ch_site_list, ch_model )
    ch_versions = ch_versions.mix(STRIDE_RUN.out.versions)

    emit:
    features    = STRIDE_RUN.out.features     // [ val(meta), path(*.tsv) ]
    predictions = STRIDE_RUN.out.predictions  // [ val(meta), path(*_msi.txt) ]
    qc_reports  = STRIDE_RUN.out.qc_reports   // [ val(meta), path(*_qc.html) ]
    drivers     = STRIDE_RUN.out.drivers      // [ val(meta), path(*_drivers.tsv) ]
    versions    = ch_versions                 // path(versions.yml)
}

