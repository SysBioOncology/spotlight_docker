process GENERATE_GRAPHS {
    // Add label for time and memory
    label 'spatial_features'
    label 'compute_spatial_features'

    input:
    path tile_quantification_path
    val out_prefix
    val slide_type

    output:
    path "${prefix}_graphs.pkl", emit: pkl

    script:
    prefix = out_prefix != "dummy" ? "${out_prefix}${slide_type}" : "${slide_type}"
    """
    generate_graphs.py \\
        --tile_quantification_path ${tile_quantification_path} \\
        --n_cores ${task.cpus} \\
        --prefix ${prefix} \\
        --slide_type ${slide_type}
    """

    stub: 
    """
    touch "${prefix}_graphs.pkl"
    """
}
