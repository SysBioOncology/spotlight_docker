// //
// // Subworkflow with functionality specific to the SysBioOncology/spotlight_docker pipeline
// //

// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */

include { TRAIN_MULTITASK_CELLTYPE_MODEL} from '../../modules/local/trainmultitaskcelltypemodel.nf'
include { PREPROCESSING_MULTITASK_MODEL_TARGET_FEATURES} from '../../modules/local/preprocessingmultitaskmodeltargetfeatures.nf'

workflow TRAIN_MULTITASK_CELLTYPE_MODELS {
    take: 
    // Target features
    cancer_type
    quantiseq_path
    mpc_counter_path
    xcell_path
    epic_path

    // Train multitask model
    bottleneck_features_path
    category
    alpha_min
    alpha_max
    n_steps
    n_outerfolds
    n_innerfolds
    n_tilescancer_type
    split_level
    slide_type
    var_names_path
    target_features_path

    clinical_file_path
    tpm_path 
    thorsson_signatures_path
    estimate_signatures_path
    absolute_tumor_purity_path
    gibbons_signatures_path
    mcp_counter_path
    quantiseq_path
    xcell_path
    epic_path

        
    main: 

    PREPROCESSING_MULTITASK_MODEL_TARGET_FEATURES(
        clinical_file_path = clinical_file_path,
        tpm_path = tpm_path, 
        thorsson_signatures_path = thorsson_signatures_path, 
        estimate_signatures_path = estimate_signatures_path,
        absolute_tumor_purity_path = absolute_tumor_purity_path, 
        gibbons_signatures_path = gibbons_signatures_path, 
        mcp_counter_path = 
        quantiseq_path
        xcell_path
        epic_path

    )

    // TRAIN_MULTITASK_CELLTYPE_MODEL(
    //     bottleneck_features_path
    //     category
    //     alpha_min
    //     alpha_max
    //     n_steps
    //     n_outerfolds
    //     n_innerfolds
    //     n_tilescancer_type
    //     split_level
    //     slide_type
    //     var_names_path
    //     target_features_path
    // )


}