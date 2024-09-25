process COMPUTE_NCLUSTERS {
    label 'spatial_clustering_features'
    label 'compute_spatial_features'

    input:
        tuple path(tiles), path(tiles_labeled)
        path cell_types
        val slide_type
        val out_prefix

    output:
        path "${prefix}_nclusters_wide.csv", emit: csv

    script:
    prefix = out_prefix != "dummy" ? "${out_prefix}${slide_type}" : "${slide_type}"
    cell_types_arg = cell_types.name == "NO_FILE" ? "": "--cell_types ${cell_types}"
    """
    compute_nclusters.py \\
        --all_slide_clusters_characterized ${tiles_labeled} \\
        --prefix ${prefix} ${cell_types_arg}
    """
}
