//
// Subworkflow with functionality specific to the SysBioOncology/spotlight_docker pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { PREDICT_TILE_LEVEL_CELL_TYPE_ABUNDANCE  } from '../../../modules/local/predicttilelevelcelltypeabundance.nf'
include { COMBINE_TILE_LEVEL_CELLTYPE_ABUNDANCE   } from '../../../modules/local/combinetilelevelcelltypeabundance.nf'

workflow PREDICT_TILE_LEVEL_CELL_TYPE_ABUNDANCES {
    take: 
        features_input         
        celltype_models_path   
        var_names_path       
        prediction_mode        
        cell_types_path       
        n_outerfolds         
        slide_type            
        is_model_dir          
    main: 


    PREDICT_TILE_LEVEL_CELL_TYPE_ABUNDANCE(
        features_input          = features_input,
        celltype_models_path    = celltype_models_path,
        var_names_path          = var_names_path,
        prediction_mode         = prediction_mode,
        cell_types_path         = cell_types_path,
        n_outerfolds            = n_outerfolds,
        slide_type              = slide_type,
        is_model_dir            = is_model_dir
    )
    ch_tile_level_celltype_predictions = PREDICT_TILE_LEVEL_CELL_TYPE_ABUNDANCE.out.csv.collect()
    

    COMBINE_TILE_LEVEL_CELLTYPE_ABUNDANCE( 
        features_input                  = features_input,
        tile_level_celltype_predictions = ch_tile_level_celltype_predictions,
        var_names_path                  = var_names_path,
        prediction_mode                 = prediction_mode,
        cell_types_path                 = cell_types_path,
        n_outerfolds                    = n_outerfolds,
        slide_type                      = slide_type

    )

    ch_target_features = COMBINE_TILE_LEVEL_CELLTYPE_ABUNDANCE.out.proba_csv.collect()


    emit: 
        proba = ch_target_features


}