/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// include { paramsSummaryMap       } from 'plugin/nf-validation'

// Modules
include { CREATE_TPM_MATRIX                       } from '../../modules/local/createtpmmatrix.nf'
include { IMMUNEDECONV                            } from '../../modules/local/immunedeconv.nf'
// Subworkflows
include { EXTRACT_HISTOPATHO_FEATURES       } from '../../subworkflows/local/extract_histopatho_features/main.nf'
include { COMPUTE_SPATIAL_FEATURES          } from '../../subworkflows/local/compute_spatial_features/main.nf'
include { BUILD_MULTITASK_CELLTYPE_MODELS   } from '../../subworkflows/local/build_multitask_celltype_models/main.nf'
include { PREDICT_TILE_LEVEL_CELL_TYPE_ABUNDANCES } from '../../subworkflows/local/predict_tile_level_celltype_abundances/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SPOTLIGHT {
    take:
    // General parameters
    slide_type
    out_prefix

    // Immunedeconvolution
    is_tpm
    deconv_tools
    gene_exp_path
    quantiseq_path
    mcp_counter_path
    xcell_path
    epic_path
    mcp_probesets
    mcp_genes

    // Extracting histopatho features
    clinical_files_input
    path_codebook
    cancer_type
    is_tumor
    clinical_file_out_file
    tumor_purity_threshold
    is_tcga
    image_dir
    gradient_mag_filter
    n_shards
    bot_out_filename
    pred_out_filename
    model_name
    checkpoint_path
    path_tissue_classes

    // Build multi-task cell type model(s)
    clinical_file_path
    thorsson_scores_path
    estimate_scores_path
    absolute_tumor_purity_path
    gibbons_scores_path
    mfp_gene_signatures_path
    bottleneck_features_path
    var_names_path
    target_features_path
    model_cell_types

    alpha_min
    alpha_max
    n_steps
    n_outerfolds
    n_innerfolds
    n_tiles
    split_level

    // Predicing tile-level cell type quantification
    celltype_models_path
    prediction_mode
    cell_types_path
    is_model_dir

    // Spatial features
    tile_level_celltype_predictions
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
    // Get spotlight modules to run
    def spotlight_modules = params.spotlight_modules ? params.spotlight_modules.split(',').collect{ it.trim().toLowerCase() } : []

    // INITIALIZE 
    ch_tpm = !params.is_tpm 
    ? Channel.empty()
    : Channel.value(gene_exp_path) 

    // If histopathofeatures run use output of the module, else user has to give path
    ch_bottleneck_features = (spotlight_modules.contains("extracthistopatho") || bottleneck_features_path.name == "NO_NAME")
    ? Channel.empty() 
    : Channel.value(bottleneck_features_path)

    // If histopathofeatures is run, then use output, else user has to give
    ch_target_features = (spotlight_modules.contains("extracthistopatho") || target_features_path.name == "NO_NAME")
    ? Channel.empty() 
    : Channel.value(target_features_path)

    ch_clinical_file = (spotlight_modules.contains("extracthistopatho") || clinical_file_path.name == "NO_NAME")
    ? Channel.empty() 
    : Channel.value(clinical_file_path)

    ch_celltype_models = (spotlight_modules.contains("buildmodel") || file(celltype_models_path).name == "NO_NAME") 
    ? Channel.empty() 
    : channel.fromFilePairs("${file(celltype_models_path)}/**/{outer_models,x_train_scaler}.pkl", relative: false)
        .map { dummy, it -> [it[0].getParent().name, it ]}

    ch_var_names = (spotlight_modules.contains("buildmodel") || var_names_path.name == "NO_NAME") 
    ? Channel.empty()
    : Channel.value(var_names_path)

    ch_tile_level_celltype_predictions =  (spotlight_modules.contains('predicttiles')) || tile_level_celltype_predictions.name == "NO_NAME"
    ? Channel.empty()
    : Channel.value(tile_level_celltype_predictions)

    ch_immune_deconv_files = Channel.of(
        ["quantiseq", quantiseq_path],
        ["epic", epic_path],
        ["mcp_counter", mcp_counter_path],
        ["xcell", xcell_path])

    if (spotlight_modules.contains("extracthistopatho")) {
        EXTRACT_HISTOPATHO_FEATURES(
            clinical_files_input,
            path_codebook,
            cancer_type,
            is_tumor,
            out_prefix,
            tumor_purity_threshold,
            is_tcga,
            image_dir,
            gradient_mag_filter,
            n_shards,
            bot_out_filename,
            pred_out_filename,
            model_name,
            checkpoint_path,
            slide_type,
            path_tissue_classes
        )
        ch_bottleneck_features  = EXTRACT_HISTOPATHO_FEATURES.out.features
        ch_clinical_file        = EXTRACT_HISTOPATHO_FEATURES.out.clinical_file
    }

    if (spotlight_modules.contains("deconvbulk")) {
        ch_tpm.ifEmpty(gene_exp_path) | CREATE_TPM_MATRIX

        ch_immune_deconv_tmp = CREATE_TPM_MATRIX.out.txt
            .collect()
            .ifEmpty(gene_exp_path)
            .combine(Channel.of("quantiseq", "epic","mcp_counter", "xcell")) 
            // Reverse order
            .map { tpm_path, tool -> return [tool, tpm_path] }
            .combine(ch_immune_deconv_files, by: 0)
            .branch{
                tool, tpm_path, deconv_path ->  
                // Tools to use for deconvolution 
                invalid: deconv_path.name == "NO_FILE"
                    return [tool, tpm_path]
                // Tools already used
                valid: true
                    return [tool, deconv_path]
            }
        
        ch_tpm = CREATE_TPM_MATRIX.out.txt

        IMMUNEDECONV (
            ch_immune_deconv_tmp.invalid, 
            mcp_probesets                   = mcp_probesets, 
            mcp_genes                       = mcp_genes)
        
        ch_immune_deconv_files = IMMUNEDECONV.out.csv
            .mix(   ch_immune_deconv_tmp.valid  )

    }

    if (spotlight_modules.contains("buildmodel")) {
        
        BUILD_MULTITASK_CELLTYPE_MODELS(
            clinical_file_path                  = ch_clinical_file,
            tpm_path                            = ch_tpm,
            thorsson_scores_path                = thorsson_scores_path,
            estimate_scores_path                = estimate_scores_path,
            absolute_tumor_purity_path          = absolute_tumor_purity_path,
            gibbons_scores_path                 = gibbons_scores_path,    
            mfp_gene_signatures_path                 = mfp_gene_signatures_path,
            immune_deconv_files                 = ch_immune_deconv_files,
            bottleneck_features_path            = ch_bottleneck_features,
            target_features                     = ch_target_features,
            model_cell_types                    = model_cell_types,
            alpha_min                           = alpha_min,
            alpha_max                           = alpha_max,
            n_steps                             = n_steps,
            n_outerfolds                        = n_outerfolds,
            n_innerfolds                        = n_innerfolds,
            n_tiles                             = n_tiles,
            split_level                         = split_level,
            slide_type                          = slide_type
        )

        ch_var_names        = BUILD_MULTITASK_CELLTYPE_MODELS.out.var_names
        ch_celltype_models  = BUILD_MULTITASK_CELLTYPE_MODELS.out.models
    }

    if (spotlight_modules.contains('predicttiles')) {
        
        PREDICT_TILE_LEVEL_CELL_TYPE_ABUNDANCES(
            features_input          = ch_bottleneck_features,
            celltype_models_path    = ch_celltype_models,
            var_names_path          = ch_var_names,
            prediction_mode         = prediction_mode,
            cell_types_path         = cell_types_path,
            n_outerfolds            = n_outerfolds,
            slide_type              = slide_type,
            is_model_dir            = is_model_dir

        )
        ch_tile_level_celltype_predictions = PREDICT_TILE_LEVEL_CELL_TYPE_ABUNDANCES.out.proba
    }


    if (spotlight_modules.contains('computespatial')) {
        COMPUTE_SPATIAL_FEATURES(
            tile_level_cell_type_quantification = ch_tile_level_celltype_predictions,
            out_prefix                          = out_prefix,
            slide_type                          = slide_type,
            abundance_threshold                 = abundance_threshold,
            cell_types                          = cell_types_path,
            shapiro_alpha                       = shapiro_alpha,
            cutoff_path_length                  = cutoff_path_length,
            n_clusters                          = n_clusters,
            max_dist                            = max_dist,
            max_n_tiles_threshold               = max_n_tiles_threshold,
            tile_size                           = tile_size,
            overlap                             = overlap,
            metadata_path                       = metadata_path,
            is_tcga                             = is_tcga,
            merge_var                           = merge_var,
            sheet_name                          = sheet_name
        ) 
    }

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
