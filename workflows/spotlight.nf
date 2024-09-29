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
// include { CELLTYPE_QUANTIFICATION_BULKRNASEQ } from '../subworkflows/local/celltype_quantification_bulkrnaseq.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SPOTLIGHT {
    take:
        // Subworkflows 
        // skip_celltype_quantification_bulkrnaeq
        // skip_build_multitask_cell_type_model
        // skip_extract_histopatho_features 
        // skip_PREDICT_CELLTYPE_QUANTIFICATION_TILES 
        // skip_derive_spatial_features
        // Parameters CELLTYPE_QUANTIFICATION_BULKRNASEQ
        is_tpm
        deconv_tools
        gene_exp_path
        quantiseq_path
        mpc_counter_path
        xcell_path
        epic_path
        

        // clinical_files_input
        // path_codebook
        // class_name
        // out_file
        // tumor_purity_threshold
        // is_tcga
        // image_dir
        // gradient_mag_filter
        // n_shards
        // bot_out
        // pred_out
        // model_name
        // checkpoint_path
        // slide_type
        // path_tissue_classes
        // celltype_models
        // var_names_path
        // prediction_mode
        // cell_types_path
        // n_outerfolds
        // out_prefix
        // abundance_threshold
        // shapiro_alpha
        // cutoff_path_length
        // n_clusters
        // max_dist
        // max_n_tiles_threshold
        // tile_size
        // overlap
        // metadata_path
        // merge_var
        // sheet_name

    main:
    // Get spotlight modules to run
    def spotlight_modules = params.spotlight_modules ? params.spotlight_modules.split(',').collect{ it.trim().toLowerCase() } : []


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
        

        // // Immune deconvolution
        // Channel.of(                                 // channel: [tool name, csv]
        //     ["quantiseq", quantiseq_path],
        //     ["epic", epic_path],
        //     ["mpc_counter", mpc_counter_path],
        //     ["xcell", xcell_path],
        // )
        // .branch{
        //     tool, filepath ->  
        //     // Tools to use for deconvolution 
        //     invalid: deconv_tools.contains(tool) && (filepath.name == "NO_FILE")
        //         return [tool, tpm_path]
        //     // Tools already used
        //     valid: true
        //         return [tool, filepath]
        // }.set {
        //     ch_immune_deconv_files
        // }

        // IMMUNEDECONV (
        //     ch_immune_deconv_files.invalid
        // )
        
        
        // IMMUNEDECONV.out.csv
        // .collect()
        // .mix(ch_immune_deconv_files.valid)
        // .toSortedList()
        // .set { immune_deconv_files }


        // immune_deconv_files.view()
    }

    if (spotlight_modules.contains("extracthistopathofeatures")) {
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
    }

    if (spotlight_modules.contains('predictcelltypequantificationtiles')) {

        PREDICT_CELLTYPE_QUANTIFICATION_TILES(
            features_input  = EXTRACT_HISTOPATHO_FEATURES.out.features,
            celltype_models = celltype_models,
            var_names_path = var_names_path,
            prediction_mode = prediction_mode,
            cell_types_path = cell_types_path,
            n_outerfolds =  n_outerfolds,
            slide_type = slide_type
        )
    }

    if (spotlight_modules.contains('derivespatialfeatures')) {
        DERIVE_SPATIAL_FEATURES(
            tile_level_cell_type_quantification = PREDICT_CELLTYPE_QUANTIFICATION_TILES.out.proba,
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
    }

    // emit:
    //     ehf_txt = EXTRACT_HISTOPATHO_FEATURES.out.txt // channel: /path/to/createdclinicalfile.txt
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
