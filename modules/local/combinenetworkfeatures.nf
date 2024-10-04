process COMBINE_NETWORK_FEATURES {
    label 'spatial_features'
    label 'process_single'

    input:
        path all_largest_cc_sizes_wide
        path shortest_paths_wide
        path colocalization_wide
        val slide_type
        val out_prefix


    output:
        path  "${prefix}_all_graph_features.csv", emit: csv

    script:
    prefix = out_prefix != "dummy" ? "${out_prefix}${slide_type}" : "${slide_type}"

    """
    combine_network_features.py \\
        --all_largest_cc_sizes_wide ${all_largest_cc_sizes_wide} \\
        --shortest_paths_wide ${shortest_paths_wide} \\
        --colocalization_wide ${colocalization_wide} \\
        --prefix ${prefix}
    """


    stub: 
    prefix = out_prefix != "dummy" ? "${out_prefix}${slide_type}" : "${slide_type}"

    """
    touch "${prefix}_all_graph_features.csv"
    """
}
