process COMPUTE_PROXIMITY_FROM_INDIV_SCHC_COMBINE {
    label 'spatial_clustering_features'
    label 'compute_spatial_features'

    input:
    path prox_between
    path prox_within
    val slide_type
    val out_prefix

    output:
    path "${prefix}_features_clust_indiv_schc_prox.csv", emit: csv

    script:
    prefix = out_prefix != "dummy" ? "${out_prefix}${slide_type}" : "${slide_type}"
    """
    compute_proximity_from_indiv_schc_combine.py \\
        --prox_between ${prox_between} \\
        --prox_within ${prox_within} \\
        --prefix ${prefix}
    """
}
