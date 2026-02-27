#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    STRiDE — MSI Prediction Pipeline for MSK-ACCESS
========================================================================================
    Nextflow DSL2 pipeline wrapping the stride CLI for scalable,
    parallelized MSI prediction on HPC and cloud.

    Usage:
      nextflow run nextflow/main.nf \
          --input samples.csv \
          --outdir results/ \
          -profile docker

    Documentation: https://github.com/msk-access/STRiDE
========================================================================================
*/

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

if (!params.input) { exit 1, 'Sample sheet not specified! Use --input <path/to/samples.csv>' }

/*
========================================================================================
    IMPORT WORKFLOWS
========================================================================================
*/

include { STRIDE } from './workflows/stride'

/*
========================================================================================
    HELPER: Discover BAI index with both naming conventions
    Supports: sample.bam.bai and sample.bai
========================================================================================
*/

def discoverBai(String bamPath) {
    def bai1 = file("${bamPath}.bai")           // sample.bam.bai
    def bai2 = file(bamPath.replaceAll(/\.bam$/, '.bai'))  // sample.bai

    if (bai1.exists()) return bai1
    if (bai2.exists()) return bai2

    error "BAI index not found for ${bamPath}. Searched:\n  - ${bamPath}.bai\n  - ${bai2}"
}

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow {

    //
    // STEP 1: Parse sample sheet
    //
    // Expected CSV columns: sample_id, tumor_bam, normal_bam
    // Optional column:      matched_norm_sample_barcode
    //
    log.info """\
        ┌─────────────────────────────────────────────────┐
        │  S T R i D E   N e x t f l o w   P i p e l i n e │
        ├─────────────────────────────────────────────────┤
        │  Sample sheet : ${params.input}
        │  Output dir   : ${params.outdir}
        │  QC reports   : ${params.generate_qc}
        │  Container    : ghcr.io/msk-access/stride:${params.stride_version}
        └─────────────────────────────────────────────────┘
        """.stripIndent()

    Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, sep: ',', quote: '"')
        .map { row ->
            // Build meta map with sample ID and optional normal barcode
            def meta = [
                id: row.sample_id,
                matched_norm_sample_barcode:
                    row.containsKey('matched_norm_sample_barcode')
                    && row.matched_norm_sample_barcode
                    ? row.matched_norm_sample_barcode : ''
            ]

            def tumor_bam  = file(row.tumor_bam, checkIfExists: true)
            def normal_bam = file(row.normal_bam, checkIfExists: true)

            // Auto-discover BAI indexes (supports both .bam.bai and .bai)
            def tumor_bai  = discoverBai(row.tumor_bam)
            def normal_bai = discoverBai(row.normal_bam)

            log.debug "Parsed sample: ${meta.id} | T=${tumor_bam} | N=${normal_bam}"

            return [meta, tumor_bam, tumor_bai, normal_bam, normal_bai]
        }
        .set { ch_input }

    //
    // STEP 2: Resolve shared resources (site list + model)
    //
    // Uses bundled defaults from the STRiDE package if not specified
    //
    def ch_sites = params.site_list
        ? file(params.site_list, checkIfExists: true)
        : file("${projectDir}/../src/stride/data/msi_sites_170.tsv")

    def ch_model = params.model_joblib
        ? file(params.model_joblib, checkIfExists: true)
        : file("${projectDir}/../src/stride/data/msi_model.joblib")

    //
    // STEP 3: Run the pipeline
    //
    STRIDE ( ch_input, ch_sites, ch_model )
}
