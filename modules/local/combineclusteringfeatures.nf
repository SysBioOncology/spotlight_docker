process COMBINE_CLUSTERING_FEATURES {
    label 'spatial_features'
    label 'process_single'

    input:
        path frac_high_wide
        path num_clust_slide_wide
        path all_prox_df_wide
        path prox_indiv_schc_combined_wide
        val slide_type
        val out_prefix

    output:
        path "${prefix}_clustering_features.csv", emit: csv

    script:
    prefix = out_prefix != "dummy" ? "${out_prefix}${slide_type}" : "${slide_type}"
    """
    combine_clustering_features.py \\
        --frac_high_wide ${frac_high_wide} \\
        --num_clust_slide_wide ${num_clust_slide_wide} \\
        --all_prox_df_wide ${all_prox_df_wide} \\
        --prox_indiv_schc_combined_wide ${prox_indiv_schc_combined_wide} \\
        --prefix ${prefix}
    """

    stub: 
    prefix = out_prefix != "dummy" ? "${out_prefix}${slide_type}" : "${slide_type}"

    """
    touch "${prefix}_clustering_features.csv"
    """
}
