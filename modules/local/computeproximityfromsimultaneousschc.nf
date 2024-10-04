process COMPUTE_PROXIMITY_FROM_SIMULTANEOUS_SCHC {
    label 'spatial_clustering_features'
    label 'compute_spatial_features'

    input:
        tuple path(tiles), path(tiles_labeled)
        path cell_types
        val n_clusters
        val max_dist
        val max_n_tiles_threshold
        val tile_size
        val overlap
        val slide_type
        val out_prefix

    output:
        path "${prefix}_features_clust_all_schc_prox_wide.csv", emit: csv

    script:
    prefix = out_prefix != "dummy" ? "${out_prefix}${slide_type}" : "${slide_type}"
    cell_types_arg = cell_types.name == "NO_FILE" ? "": "--cell_types ${cell_types}"
    max_dist_arg = max_dist == "dummy" ? "" : "--max_dist ${max_dist}"
    """
    compute_proximity_from_simultaneous_schc.py \\
        --slide_clusters_characterized ${tiles_labeled} \\
        --tiles_schc ${tiles} \\
        --n_cores ${task.cpus} \\
        --max_n_tiles_threshold ${max_n_tiles_threshold} \\
        --n_clusters ${n_clusters} \\
        --tile_size ${tile_size} \\
        --overlap ${overlap} \\
        --prefix ${prefix} ${cell_types_arg} ${max_dist_arg}
    """

    stub: 
    prefix = out_prefix != "dummy" ? "${out_prefix}${slide_type}" : "${slide_type}"

    """
    touch "${prefix}_features_clust_all_schc_prox_wide.csv"
    """
}
