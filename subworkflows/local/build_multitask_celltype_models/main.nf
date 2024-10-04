// //
// // Subworkflow with functionality specific to the SysBioOncology/spotlight_docker pipeline
// //

// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */

include { BUILD_MULTITASK_CELLTYPE_MODEL} from '../../../modules/local/buildmultitaskcelltypemodel.nf'
include { PREPROCESSING_MULTITASK_MODEL_TARGET_FEATURES} from '../../../modules/local/preprocessingmultitaskmodeltargetfeatures.nf'

workflow BUILD_MULTITASK_CELLTYPE_MODELS {
    take: 
    clinical_file_path
    tpm_path
    thorsson_scores_path
    estimate_scores_path
    absolute_tumor_purity_path
    gibbons_scores_path
    mfp_gene_signatures_path
    immune_deconv_files
    //   Train multitask model
    bottleneck_features_path
    target_features
    model_cell_types
    alpha_min
    alpha_max
    n_steps
    n_outerfolds
    n_innerfolds
    n_tiles
    split_level
    slide_type
    
    main: 
    ch_model_cell_types     = Channel.of(model_cell_types.split(","))


    // Split into separate branches (expected 1 file per tool)
    ch_immune_deconv_files  = immune_deconv_files.branch {
        tool, filepath -> 
            quantiseq: tool == "quantiseq"
                return filepath
            mcp_counter: tool == "mcp_counter"
                return filepath
            xcell: tool == "xcell"
                return filepath
            epic: tool == "epic"
                return filepath
    }

    PREPROCESSING_MULTITASK_MODEL_TARGET_FEATURES(
        // Clinical file needs to have at least 'sample_submitter_id' and 'slide_submitter_id'
        clinical_file_path              = clinical_file_path,
        tpm_path                        = tpm_path, 
        thorsson_scores_path            = thorsson_scores_path, 
        estimate_scores_path            = estimate_scores_path,
        absolute_tumor_purity_path      = absolute_tumor_purity_path, 
        gibbons_scores_path             = gibbons_scores_path, 
        mfp_gene_signatures_path             = mfp_gene_signatures_path,
        mcp_counter_path                = ch_immune_deconv_files.mcp_counter,
        quantiseq_path                  = ch_immune_deconv_files.quantiseq,
        xcell_path                      = ch_immune_deconv_files.xcell,
        epic_path                       = ch_immune_deconv_files.epic,
    )
    // Close channels
    ch_var_names        = PREPROCESSING_MULTITASK_MODEL_TARGET_FEATURES.out.pkl.collect()
    ch_target_features  = PREPROCESSING_MULTITASK_MODEL_TARGET_FEATURES.out.csv.collect()

    ch_target_features  = ch_target_features.ifEmpty(target_features)

    BUILD_MULTITASK_CELLTYPE_MODEL(
        bottleneck_features_path    = bottleneck_features_path,
        ch_model_cell_types         = ch_model_cell_types, // channel
        alpha_min                   = alpha_min,
        alpha_max                   = alpha_max,
        n_steps                     = n_steps,
        n_outerfolds                = n_outerfolds,
        n_innerfolds                = n_innerfolds,
        n_tiles                     = n_tiles,
        split_level                 = split_level,
        slide_type                  = slide_type,
        var_names_path              = ch_var_names,
        target_features_path        = ch_target_features
    )

    ch_models = BUILD_MULTITASK_CELLTYPE_MODEL.out.pkl
        .map {
            it -> [ it[0].getParent().name, it  ] 
        }   

    emit: 
    var_names   = var_names_path
    models      = ch_models

}