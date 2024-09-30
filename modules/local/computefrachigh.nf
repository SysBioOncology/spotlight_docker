process COMPUTE_FRAC_HIGH {
    label 'spatial_clustering_features'
    label 'compute_spatial_features'

    input:
        tuple path(tiles), path(tiles_labeled)
        val slide_type
        val out_prefix

    output:
        path "${prefix}_frac_high_wide.csv", emit: csv

    script:
    prefix = out_prefix != "dummy" ? "${out_prefix}${slide_type}" : "${slide_type}"
    """
    compute_frac_high.py \\
        --slide_indiv_clusters_labeled ${tiles_labeled} \\
        --prefix ${prefix}
    """

    stub: 
    prefix = out_prefix != "dummy" ? "${out_prefix}${slide_type}" : "${slide_type}"

    """
    touch "${prefix}_frac_high_wide.csv"
    """
}
