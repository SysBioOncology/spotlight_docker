def determine_out_prefix(slide_type) {
    def res = ""
    switch (slide_type) {
        case { it == "FF" }:
            res = "predictions.txt"
            break
        case { it == "FFPE" }:
            res = "predictions-*.parquet"
            break
        default:
            res = "*"
            break
    }
    return res
}

// # ----------------------------------------------------- #
// # ---- Post-processing of predictions and futures ----- #
// # ----------------------------------------------------- #
process POST_PROCESSING_PREDICTIONS {
    label 'process_medium'
    label "extract_histo_patho_features"

    input:
    path path_codebook
    path path_tissue_classes
    path pred_train_file
    val cancer_type
    val slide_type

    output:
    path out_prefix, emit: txt_parquet
    path "ok.txt"

    script:
    out_prefix = determine_out_prefix(slide_type)
    """

    post_process_predictions.py \
        --slide_type ${slide_type} \
        --path_codebook "${path_codebook}" \
        --cancer_type ${cancer_type} \
        --pred_train_file ${pred_train_file} \
        --path_tissue_classes "${path_tissue_classes}" && touch \$PWD/ok.txt
    """

    stub: 
    out_prefix = determine_out_prefix(slide_type)
    """
    touch ${out_prefix}
    touch "ok.txt"
    """
}
