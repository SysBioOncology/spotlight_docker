#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/spotlight
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/spotlight
    Website: https://nf-co.re/spotlight
    Slack  : https://nfcore.slack.com/channels/spotlight
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SPOTLIGHT  } from './workflows/spotlight'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_spotlight_pipeline'
include { PIPELINE_COMPLETION     } from'./subworkflows/local/utils_nfcore_spotlight_pipeline'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_SPOTLIGHT {

    // take:
    //     clinical_files_input
    //     path_codebook
    //     class_name
    //     clinical_file_out_file
    //     tumor_purity_threshold
    //     is_tcga,
    //     image_dir


    main:

    clinical_files_input = file(params.clinical_files_input, checkIfExists: true)
    path_codebook = file(params.path_codebook, checkIfExists:true)
    image_dir = file(params.image_dir , checkIfExists:true)
    checkpoint_path = file(params.checkpoint_path, checkIfExists:false)
    path_tissue_classes = file(params.path_tissue_classes)
    celltype_models = file(params.celltype_models)
    var_names_path = file(params.var_names_path)
    cell_types_path = file(params.cell_types_path)
    metadata_path = file(params.metadata_path)

    SPOTLIGHT (
        clinical_files_input,
        path_codebook,
        params.class_name,
        params.clinical_file_out_file,
        params.tumor_purity_threshold,
        params.is_tcga,
        image_dir,
        params.gradient_mag_filter,
        params.n_shards,
        params.bot_out,
        params.pred_out,
        params.model_name,
        checkpoint_path,
        params.slide_type,
        path_tissue_classes,
        celltype_models,
        var_names_path,
        params.prediction_mode,
        cell_types_path,
        params.n_outerfolds,
        params.out_prefix,
        params.abundance_threshold,
        params.shapiro_alpha,
        params.cutoff_path_length,
        params.n_clusters,
        params.max_dist,
        params.max_n_tiles_threshold,
        params.tile_size,
        params.overlap,
        metadata_path,
        params.merge_var,
        params.sheet_name

    )

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
    // PIPELINE_INITIALISATION (
    //     params.version,
    //     params.help,
    //     params.validate_params,
    //     params.monochrome_logs,
    //     args,
    //     params.outdir,
    //     params.input
    // )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_SPOTLIGHT (
        // PIPELINE_INITIALISATION.out.samplesheet
        // PIPELINE_INITIALISATION
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    // PIPELINE_COMPLETION (
    //     params.email,
    //     params.email_on_fail,
    //     params.plaintext_email,
    //     params.outdir,
    //     params.monochrome_logs,
    //     params.hook_url,
    // )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
