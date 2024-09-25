process COMPUTE_COLOCALIZATION {
    label 'compute_spatial_features'
    label 'spatial_network_features'

    input:
        path tile_quantification_path
        path cell_types
        path graphs_path
        val abundance_threshold
        val slide_type
        val out_prefix

    output:
        path  "${prefix}_features_coloc_fraction_wide.csv", emit: csv
        path  "${prefix}_features_coloc_fraction.csv"
        path "${prefix}_graphs.pkl", optional: true

    script:
    prefix = out_prefix != "dummy" ? "${out_prefix}${slide_type}" : "${slide_type}"
    graphs_path_arg = graphs_path.name == "NO_FILE" ? "" : "--graphs_path ${graphs_path}"
    cell_types_arg = cell_types.name == "NO_FILE" ? "": "--cell_types ${cell_types}"
    """
    compute_colocalization.py \\
        --tile_quantification_path ${tile_quantification_path} \\
        --abundance_threshold ${abundance_threshold} \\
        --n_cores ${task.cpus} \\
        --prefix ${prefix} ${graphs_path_arg} ${cell_types_arg}
    """
}
