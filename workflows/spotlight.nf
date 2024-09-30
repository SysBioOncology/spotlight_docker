/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// include { paramsSummaryMap       } from 'plugin/nf-validation'

// Modules
include {   CREATE_TPM_MATRIX   } from '../modules/local/createtpmmatrix.nf'
include {   IMMUNEDECONV        } from '../modules/local/immunedeconv.nf'

// Subworkflows
include { EXTRACT_HISTOPATHO_FEATURES } from '../subworkflows/local/extract_histopatho_features.nf'
include { PREDICT_CELLTYPE_QUANTIFICATION_TILES} from '../subworkflows/local/predictcelltypequantificationtiles.nf'
include { DERIVE_SPATIAL_FEATURES } from '../subworkflows/local/derive_spatial_features.nf'
include { BUILD_MULTITASK_CELLTYPE_MODELS } from '../subworkflows/local/build_multi_task_celltype_models.nf'

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
    mpc_counter_path
    xcell_path
    epic_path

    // Extracting histopatho features
    clinical_files_input
    path_codebook
    class_name
    clinical_file_out_file
    tumor_purity_threshold
    is_tcga
    image_dir
    gradient_mag_filter
    n_shards
    bot_out
    pred_out
    model_name
    checkpoint_path
    path_tissue_classes

    // Predicing tile-level cell type quantification
    celltype_models
    var_names_path
    prediction_mode
    cell_types_path

    // Build multi-task cell type model(s)


    // Spatial features
    n_outerfolds
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


    // Initialize 
    // If histopathofeatures run use output of the module, else user has to give path
    histopatho_features_path = spotlight_modules.contains("extracthistopathofeatures") 
    ? Channel.empty() 
    : file(params.histopatho_features_path)

    // If histopathofeatures is run, then use output, else user has to give
    tile_level_cell_type_quantification_path = spotlight_modules.contains("extracthistopathofeatures")
    ? Channel.empty() 
    : file(params.tile_level_cell_type_quantification_path)

    // ch_immune_deconv_files =  Channel.of(  // channel: [tool name, csv]
    //         ["quantiseq", quantiseq_path],
    //         ["epic", epic_path],
    //         ["mpc_counter", mpc_counter_path],
    //         ["xcell", xcell_path],
    //     )

    if (spotlight_modules.contains("immunedeconv")) {
        CREATE_TPM_MATRIX(gene_exp_path)
        tpm_path = is_tpm ? Channel.fromPath(gene_exp_path) : Channel.empty()

        CREATE_TPM_MATRIX.out.txt.set { tpm_path }
        
        // Immune deconvolution
        Channel.of(                                 // channel: [tool name, csv]
            ["quantiseq", quantiseq_path],
            ["epic", epic_path],
            ["mpc_counter", mpc_counter_path],
            ["xcell", xcell_path],
        )
        .branch{
            tool, filepath ->  
            // Tools to use for deconvolution 
            invalid: deconv_tools.contains(tool) && (filepath.name == "NO_FILE")
                return [tool, tpm_path]
            // Tools already used
            valid: true
                return [tool, filepath]
        }.set {
            ch_immune_deconv_files
        }

        IMMUNEDECONV (
            ch_immune_deconv_files.invalid
        )
        
        
        IMMUNEDECONV.out.csv
        .collect()
        .mix(ch_immune_deconv_files.valid)
        .toSortedList()
        .set { immune_deconv_files }

        immune_deconv_files.view()

    }


    // if (spotlight_modules.contains("")) {
    //     BUILD_MULTITASK_CELLTYPE_MODELS(
    //         cancer_type = class_name,
    //         immune_deconv_files,
    //         clinical_file_path,
    //         thorsson_signatures_path,
    //         estimate_signatures_path,
    //         absolute_tumor_purity_path,
    //         gibbons_signatures_path,
    //         tpm_path




    //     )
    // }

    if (spotlight_modules.contains("extracthistopathofeatures")) {
        EXTRACT_HISTOPATHO_FEATURES(
            clinical_files_input,
            path_codebook,
            class_name,
            out_prefix,
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
        EXTRACT_HISTOPATHO_FEATURES.out.features.set { histopatho_features_path }
    }

    if (spotlight_modules.contains('predictcelltypequantificationtiles')) {

        PREDICT_CELLTYPE_QUANTIFICATION_TILES(
            features_input          = histopatho_features_path,
            celltype_models,
            var_names_path,
            prediction_mode,
            cell_types_path,
            n_outerfolds,
            slide_type,
        )

        PREDICT_CELLTYPE_QUANTIFICATION_TILES.out.proba.set {tile_level_cell_type_quantification_path}
    }


    if (spotlight_modules.contains('derivespatialfeatures')) {
        DERIVE_SPATIAL_FEATURES(
            tile_level_cell_type_quantification = tile_level_cell_type_quantification_path,
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
            merge_var,
            sheet_name
        ) 
    }

    // emit:
    //     ehf_txt = EXTRACT_HISTOPATHO_FEATURES.out.txt // channel: /path/to/createdclinicalfile.txt
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
