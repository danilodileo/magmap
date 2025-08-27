#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/magmap
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/magmap
    Website: https://nf-co.re/magmap
    Slack  : https://nfcore.slack.com/channels/magmap
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MAGMAP                  } from './workflows/magmap'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_magmap_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_magmap_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_MAGMAP {

    take:
    samplesheet             // channel: samplesheet read in from --input
    genomeinfo              // channel: genome information sheet read in from --genomeinfo
    remote_genome_sources   // channel: NCBI-style genome summary files read in via --remote_genome_sources
    gtdb_metadata           // channel: GTDB metadata files
    gtdbtk_metadata         // channel: GTDB-Tk metadata files
    checkm_metadata         // channel: CheckM/CheckM2 metadata files

    main:

    //
    // WORKFLOW: Run pipeline
    //
    MAGMAP (
        samplesheet,
        genomeinfo,
        remote_genome_sources,
        gtdb_metadata,
        gtdbtk_metadata,
        checkm_metadata
    )
    emit:
    multiqc_report = MAGMAP.out.multiqc_report // channel: /path/to/multiqc_report.html
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.genomeinfo,
        params.remote_genome_sources,
        params.gtdb_metadata,
        params.gtdbtk_metadata,
        params.checkm_metadata
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_MAGMAP (
        PIPELINE_INITIALISATION.out.samplesheet,
        PIPELINE_INITIALISATION.out.genomeinfo,
        PIPELINE_INITIALISATION.out.remote_genome_sources,
        PIPELINE_INITIALISATION.out.gtdb_metadata,
        PIPELINE_INITIALISATION.out.gtdbtk_metadata,
        PIPELINE_INITIALISATION.out.checkm_metadata
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_MAGMAP.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
