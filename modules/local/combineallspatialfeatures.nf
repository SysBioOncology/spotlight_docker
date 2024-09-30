process COMBINE_ALL_SPATIAL_FEATURES {
    label 'spatial_features'
    label 'process_single'
    input:
        path graph_features
        path clustering_features
        path metadata_path
        val is_tcga
        val merge_var
        val sheet_name
        val slide_type
        val out_prefix

    output:
        path "${prefix}_all_features_combined.csv", emit: csv

    script:
    is_tcga_to_int= is_tcga ? 1 : 0
    prefix = out_prefix != "dummy" ? "${out_prefix}${slide_type}" : "${slide_type}"

    metadata_arg = metadata_path.name == "NO_FILE" ? "" : "--metadata_path ${metadata_path}"
    sheet_name_arg = sheet_name == "dummy" ? "" : "--sheet_name ${sheet_name}"
    """

    combine_all_spatial_features.py \\
        --graph_features ${graph_features} \\
        --clustering_features ${clustering_features} \\
        --is_tcga ${is_tcga_to_int} \\
        --merge_var ${merge_var} \\
        --prefix ${prefix} \\
        ${metadata_arg} ${sheet_name_arg}


    """

    stub: 
    prefix = out_prefix != "dummy" ? "${out_prefix}${slide_type}" : "${slide_type}"

    """
    touch "${prefix}_all_features_combined.csv"
    """
}
