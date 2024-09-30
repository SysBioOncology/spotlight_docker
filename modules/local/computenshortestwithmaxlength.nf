process COMPUTE_N_SHORTEST_PATHS_WITH_MAX_LENGTH {
    label 'compute_spatial_features'
    label 'spatial_network_features'

    input:
        path tile_quantification_path
        path cell_types
        path graphs_path
        val cutoff_path_length
        val slide_type
        val out_prefix

    output:
        path "${prefix}_features_shortest_paths_thresholded.csv"
        path "${prefix}_features_shortest_paths_thresholded_wide.csv", emit: csv
        path "${prefix}_graphs.pkl", optional: true

    script:
    prefix = out_prefix != "dummy" ? "${out_prefix}${slide_type}" : "${slide_type}"
    graphs_path_arg = graphs_path.name == "NO_FILE" ? "" : "--graphs_path ${graphs_path}"
    cell_types_arg = cell_types.name == "NO_FILE" ? "": "--cell_types ${cell_types}"
    """
    compute_nshortest_with_max_length.py \\
        --tile_quantification_path ${tile_quantification_path} \\
        --n_cores ${task.cpus} \\
        --cutoff_path_length ${cutoff_path_length} \\
        --prefix ${prefix} ${graphs_path_arg} ${cell_types_arg}
    """

    stub: 
    prefix = out_prefix != "dummy" ? "${out_prefix}${slide_type}" : "${slide_type}"

    """
    touch "${prefix}_features_shortest_paths_thresholded.csv"
    touch "${prefix}_features_shortest_paths_thresholded_wide.csv"
    touch "${prefix}_graphs.pkl"
    """
}

