process PREDICT_TILE_LEVEL_CELL_TYPE_ABUNDANCE {
    label 'processing_low'
    label 'tf_learning_celltyp_quant'

    input:
    path features_input
    tuple val(cell_type), path(celltype_models)
    path var_names_path
    val prediction_mode
    path cell_types_path
    val n_outerfolds
    val slide_type
    val is_model_dir

    output:
    path "${prediction_mode}_${cell_type}_tile_predictions_zscores.csv", emit: csv

    script:
    features_input_arg = slide_type == "FFPE" ? "" : "--features_input ${features_input}"
    
    """
    predict_tile_level_celltype_abundance.py \
        --models_dir \$PWD \
        --cell_type ${cell_type} \
        --prediction_mode ${prediction_mode} \
        --var_names_path ${var_names_path} \
        --slide_type ${slide_type} \
        --n_outerfolds ${n_outerfolds} ${features_input_arg} \
        --is_model_dir ${is_model_dir}
    """

    stub: 
    """
    touch "${prediction_mode}_${cell_type}_tile_predictions_zscores.csv"
    """
}
