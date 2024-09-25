/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// include { paramsSummaryMap       } from 'plugin/nf-validation'
// include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_spotlight_pipeline'
include { EXTRACT_HISTOPATHO_FEATURES } from '../subworkflows/local/extract_histopatho_features.nf'
include { TF_LEARNING_CELLTYPE_QUANT} from '../subworkflows/local/tf_learning_celltyp_quant.nf'
include { DERIVE_SPATIAL_FEATURES } from '../subworkflows/local/derive_spatial_features.nf'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SPOTLIGHT {
    take:
        clinical_files_input
        path_codebook
        class_name
        out_file
        tumor_purity_threshold
        is_tcga
        image_dir
        gradient_mag_filter
        n_shards
        bot_out
        pred_out
        model_name
        checkpoint_path
        slide_type
        path_tissue_classes
        celltype_models
        var_names_path
        prediction_mode
        cell_types_path
        n_outerfolds
        out_prefix
        abundance_threshold
        shapiro_alpha
        cutoff_path_length
        n_clusters
        max_dist
        max_n_tiles_threshold
        tile_size
        overlap
        metadata_path
        merge_var
        sheet_name

    main:

    EXTRACT_HISTOPATHO_FEATURES(
        clinical_files_input,
        path_codebook,
        class_name,
        out_file,
        tumor_purity_threshold,
        is_tcga,
        image_dir,
        gradient_mag_filter,
        n_shards,
        bot_out,
        pred_out,
        model_name,
        checkpoint_path,
        slide_type,
        path_tissue_classes
    )

    TF_LEARNING_CELLTYPE_QUANT(
        features_input  = EXTRACT_HISTOPATHO_FEATURES.out.features,
        celltype_models = celltype_models,
        var_names_path = var_names_path,
        prediction_mode = prediction_mode,
        cell_types_path = cell_types_path,
        n_outerfolds =  n_outerfolds,
        slide_type = slide_type

    )

    DERIVE_SPATIAL_FEATURES(
        tile_level_cell_type_quantification = TF_LEARNING_CELLTYPE_QUANT.out.proba,
        out_prefix = out_prefix,
        slide_type = slide_type,
        abundance_threshold = abundance_threshold,
        cell_types = cell_types_path,
        shapiro_alpha = shapiro_alpha,
        cutoff_path_length = cutoff_path_length,
        n_clusters = n_clusters,
        max_dist = max_dist,
        max_n_tiles_threshold = max_n_tiles_threshold,
        tile_size = tile_size,
        overlap = overlap,
        metadata_path = metadata_path,
        is_tcga = is_tcga,
        merge_var = merge_var,
        sheet_name = sheet_name

    )
    // emit:
    //     ehf_txt = EXTRACT_HISTOPATHO_FEATURES.out.txt // channel: /path/to/createdclinicalfile.txt
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
