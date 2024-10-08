process COMBINE_TILE_LEVEL_CELLTYPE_ABUNDANCE {
    label 'processing_low'
    label 'tf_learning_celltyp_quant'

    input:
    path features_input
    path tile_level_celltype_predictions, stageAs: "tile_level_predictions/*"
    path var_names_path
    val prediction_mode
    path cell_types_path
    val n_outerfolds
    val slide_type

    output:
    path "${prediction_mode}_tile_predictions_proba.csv", emit: proba_csv
    path "${prediction_mode}_tile_predictions_zscores.csv", emit: zscores_csv

    script:
    cell_types_arg = cell_types_path.name == "NO_FILE" ? "": "--cell_types ${cell_types_path}"
    features_input_arg = slide_type == "FFPE" ? "" : "--features_input ${features_input}"
    
    """
    combine_tile_level_celltype_abundance.py \
        --tile_predictions_input_dir tile_level_predictions \
        --prediction_mode ${prediction_mode} \
        --var_names_path ${var_names_path} \
        --slide_type ${slide_type} \
        --n_outerfolds ${n_outerfolds} ${cell_types_arg} ${features_input_arg}
    """

    stub: 
    """
    touch "${prediction_mode}_tile_predictions_zscores.csv"
    touch "${prediction_mode}_tile_predictions_proba.csv"

    """
}
