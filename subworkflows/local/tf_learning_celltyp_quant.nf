//
// Subworkflow with functionality specific to the nf-core/spotlight pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { TILE_LEVEL_CELL_TYPE_QUANT} from '../../modules/local/tilelevelcelltypequant.nf'


workflow TF_LEARNING_CELLTYPE_QUANT {
    take:
    features_input
    celltype_models
    var_names_path
    prediction_mode
    cell_types_path
    n_outerfolds
    slide_type

    main:
    TILE_LEVEL_CELL_TYPE_QUANT(
        features_input,
        celltype_models,
        var_names_path,
        prediction_mode,
        cell_types_path,
        n_outerfolds,
        slide_type
    )

    emit:
        proba = TILE_LEVEL_CELL_TYPE_QUANT.out.proba_csv
        zscore = TILE_LEVEL_CELL_TYPE_QUANT.out.zscore_csv


}
