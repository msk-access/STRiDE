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

def resolveInputFile(configuredPath, fallbackPaths) {
    if (configuredPath) return file(configuredPath, checkIfExists: true)

    return fallbackPaths.find { candidate -> candidate?.exists() } ?: file('NO_FILE')
}

def writePipelineOutputManifest(outdir) {
    if (!workflow.success) return

    def outdirPath = file(outdir)
    def manifestPath = outdirPath.resolve('manifest.json')
    def publishedFiles = []

    outdirPath.eachFileRecurse(groovy.io.FileType.FILES) { outputFile ->
        if (outputFile != manifestPath) {
            publishedFiles << [target: outputFile.toAbsolutePath().toString()]
        }
    }

    manifestPath.text = groovy.json.JsonOutput.toJson([published: publishedFiles])
}

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow {

    // Validate inputs
    if (!params.input) {
        error "Sample sheet not specified! Use --input <path/to/samples.csv>"
    }

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
            def norm_barcode = (row.containsKey('matched_norm_sample_barcode') && row.matched_norm_sample_barcode) ? row.matched_norm_sample_barcode : ''
            def meta = [
                id: row.sample_id,
                matched_norm_sample_barcode: norm_barcode
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
    def ch_sites = resolveInputFile(params.site_list, [
        file("${projectDir}/../src/stride/data/msi_sites_170.txt")
    ])

    def modelType = params.model.toString().toLowerCase()
    def defaultModels = []

    if (modelType.contains('tabpfn')) {
        // TabPFN preset model files are not bundled — the model must be passed explicitly
        if (!params.model_joblib) {
            error "--model_joblib is required for --model ${params.model}: pass the path to the TabPFN model .joblib (e.g. SPECIAL_k2__ed__te_0.6975.joblib for preset ao_top1; see `stride models`)"
        }
    } else {
        defaultModels << file("${projectDir}/../src/stride/models/svm/msi_sgd_v1.joblib")
        defaultModels << file("${projectDir}/../src/stride/models/msi_sgd_v1.joblib")
    }

    def ch_model = resolveInputFile(params.model_joblib, defaultModels)

    //
    // STEP 3: Run the pipeline
    //
    STRIDE ( ch_input, ch_sites, ch_model )
}

workflow.onComplete {
    writePipelineOutputManifest(params.outdir)
}
