process COMPUTE_NODE_DEGREE_WITH_ES {
    label 'compute_spatial_features'
    label 'spatial_network_features'

    input:
        path tile_quantification_path
        path cell_types
        path graphs_path
        val shapiro_alpha
        val slide_type
        val out_prefix

    output:
        path "${prefix}_features_ND_ES.csv", emit: csv
        path "${prefix}_features_ND_sims.csv"
        path "${prefix}_features_ND.csv"
        path "${prefix}_features_ND_sim_assignments.pkl"
        path "${prefix}_shapiro_tests.csv"
        path "${prefix}_graphs.pkl", optional: true

    script:
    prefix = out_prefix != "dummy" ? "${out_prefix}${slide_type}" : "${slide_type}"
    graphs_path_arg = graphs_path.name == "NO_FILE" ? "" : "--graphs_path ${graphs_path}"
    cell_types_arg = cell_types.name == "NO_FILE" ? "": "--cell_types ${cell_types}"
    """
    compute_node_degree_with_es.py \\
        --tile_quantification_path ${tile_quantification_path} \\
        --n_cores ${task.cpus} \\
        --shapiro_alpha ${shapiro_alpha} \\
        --prefix ${prefix} ${graphs_path_arg} ${cell_types_arg}
    """

    stub: 
        prefix = out_prefix != "dummy" ? "${out_prefix}${slide_type}" : "${slide_type}"

    """
    touch "${prefix}_features_ND_ES.csv"
    touch "${prefix}_features_ND_sims.csv"
    touch "${prefix}_features_ND.csv"
    touch "${prefix}_features_ND_sim_assignments.pkl"
    touch "${prefix}_shapiro_tests.csv"
    touch "${prefix}_graphs.pkl"
    """
}
