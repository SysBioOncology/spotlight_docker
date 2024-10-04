process BUILD_MULTITASK_CELLTYPE_MODEL {
    input: 
        path bottleneck_features_path
        val cell_type
        val alpha_min
        val alpha_max
        val n_steps
        val n_outerfolds
        val n_innerfolds
        val n_tiles
        val split_level
        val slide_type
        path var_names_path
        path target_features_path

    output:
        tuple path("${cell_type}/cv_outer_splits.pkl"),
        path("${cell_type}/total_tile_selection.pkl"),
        path("${cell_type}/outer_models.pkl"),
        path("${cell_type}/x_train_scaler.pkl"),
        path("${cell_type}/y_train_scaler.pkl"),
        path("${cell_type}/outer_scores_slides_train.pkl"),
        path("${cell_type}/outer_scores_slides_test.pkl"),
        path("${cell_type}/outer_scores_tiles_train.pkl"),
        path("${cell_type}/outer_scores_tiles_test.pkl"), emit: pkl
    
    script: 
    """
    build_multitask_celltype_model.py \\
        --bottleneck_features_path ${bottleneck_features_path} \\
        --category ${cell_type} \\
        --alpha_min ${alpha_min} \\
        --alpha_max ${alpha_max} \\
        --n_steps ${n_steps} \\
        --n_outerfolds ${n_outerfolds} \\
        --n_innerfolds ${n_innerfolds} \\
        --n_tiles ${n_tiles} \\
        --split_level ${split_level} \\
        --slide_type ${slide_type} \\
        --n_cores ${task.cpus} \\
        --var_names_path ${var_names_path} \\
        --target_features_path ${target_features_path}
    """

    stub: 
    """
    mkdir -p ${cell_type}
    touch "${cell_type}/cv_outer_splits.pkl"
    touch "${cell_type}/total_tile_selection.pkl"
    touch "${cell_type}/outer_models.pkl"
    touch "${cell_type}/x_train_scaler.pkl"
    touch "${cell_type}/y_train_scaler.pkl"
    touch "${cell_type}/outer_scores_slides_train.pkl"
    touch "${cell_type}/outer_scores_slides_test.pkl"
    touch "${cell_type}/outer_scores_tiles_train.pkl"
    touch "${cell_type}/outer_scores_tiles_test.pkl"
    """
}