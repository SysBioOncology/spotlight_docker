#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SysBioOncology/spotlight_docker
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/SysBioOncology/spotlight_docker
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

include { SPOTLIGHT  } from './workflows/spotlight/main.nf'
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


    // main:

    clinical_files_input = file(params.clinical_files_input, checkIfExists: true)
    cell_types_path = file(params.cell_types_path)

    // Extracting histopatho features
    path_codebook = file(params.path_codebook, checkIfExists:true)
    image_dir = file(params.image_dir , checkIfExists:true)
    checkpoint_path = file(params.checkpoint_path, checkIfExists:false)
    path_tissue_classes = file(params.path_tissue_classes)

    // Immune deconv bulkRNAseq
    gene_exp_path = file(params.gene_exp_path)
    quantiseq_path = file(params.quantiseq_path)
    mcp_counter_path = file(params.mcp_counter_path)
    xcell_path = file(params.xcell_path)
    epic_path = file(params.epic_path)

    mcp_probesets = file(params.mcp_probesets)
    mcp_genes = file(params.mcp_genes)

    // Build multi-task cell type models
    clinical_file_path = file(params.clinical_file_path)
    thorsson_scores_path = file(params.thorsson_scores_path)
    estimate_scores_path = file(params.estimate_scores_path)
    absolute_tumor_purity_path = file(params.absolute_tumor_purity_path)
    gibbons_scores_path = file(params.gibbons_scores_path)
    mfp_gene_signatures_path = file(params.mfp_gene_signatures_path)
    bottleneck_features_path = file(params.bottleneck_features_path)
    var_names_path = file(params.var_names_path)
    target_features_path = file(params.target_features_path)

    // Predict tile-level cell type abundances
    celltype_models_path = file(params.celltype_models_path)

    // Spatial features
    metadata_path = file(params.metadata_path)


    tile_level_celltype_predictions = file(params.tile_level_celltype_predictions)

    SPOTLIGHT (
        // General parameters
        params.slide_type,
        params.out_prefix,

        // Immunedeconvolution
        params.is_tpm,
        params.deconv_tools,
        gene_exp_path,
        quantiseq_path,
        mcp_counter_path,
        xcell_path,
        epic_path,
        mcp_probesets, 
        mcp_genes,

        // Extracting histopatho features
        clinical_files_input,
        path_codebook,
        params.cancer_type,
        params.is_tumor,
        params.clinical_file_out_file,
        params.tumor_purity_threshold,
        params.is_tcga,
        image_dir,
        params.gradient_mag_filter,
        params.n_shards,
        params.bot_out_filename,
        params.pred_out_filename,
        params.model_name,
        checkpoint_path,
        path_tissue_classes,

        // Build multi-task cell type model(s)
        clinical_file_path,
        thorsson_scores_path,
        estimate_scores_path,
        absolute_tumor_purity_path,
        gibbons_scores_path,
        mfp_gene_signatures_path,
        bottleneck_features_path,
        var_names_path,
        target_features_path,
        params.model_cell_types,

        params.alpha_min,
        params.alpha_max,
        params.n_steps,
        params.n_outerfolds,
        params.n_innerfolds,
        params.n_tiles,
        params.split_level,

        // Predicing tile-level cell type quantification
        params.celltype_models_path,
        params.prediction_mode,
        cell_types_path,
        params.is_model_dir,

        // Spatial features
        tile_level_celltype_predictions,
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
