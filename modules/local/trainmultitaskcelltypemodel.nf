// process TRAIN_MULTITASK_CELLTYPE_MODEL {
//     input: 
//         path bottleneck_features_path
//         val category
//         val alpha_min
//         val alpha_max
//         val n_steps
//         val n_outerfolds
//         val n_innerfolds
//         val n_tiles
//         val split_level
//         val slide_type
//         path var_names_path
//         path target_features_path

//     output:
//         path "${category}/cv_outer_splits.pkl"
//         path "${category}/total_tile_selection.pkl"
//         path "${category}/outer_models.pkl"
//         path "${category}/x_train_scaler.pkl"
//         path "${category}/y_train_scaler.pkl"
//         path "${category}/outer_scores_slides_train.pkl"
//         path "${category}/outer_scores_slides_test.pkl"
//         path "${category}/outer_scores_tiles_train.pkl"
//         path "${category}/outer_scores_tiles_test.pkl"
    
//     script: 
//     """
//     train_multitask_celltype_model.py \\
//         --bottleneck_features_path ${bottleneck_features_path} \\
//         --category ${category} \\
//         --alpha_min ${alpha_min} \\
//         --alpha_max ${alpha_max} \\
//         --n_steps ${n_steps} \\
//         --n_outerfolds ${n_outerfolds} \\
//         --n_innerfolds ${n_innerfolds} \\
//         --n_tiles ${n_tiles} \\
//         --split_level ${split_level} \\
//         --slide_type ${slide_type} \\
//         --n_jobs ${task.cpus} \\
//         --var_names_path ${var_names_path} \\
//         --target_features_path ${target_features}
//     """

// }