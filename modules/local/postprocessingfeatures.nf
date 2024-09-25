def determine_out_prefix(slide_type) {
    def res = ""
    switch (slide_type) {
        case { it == "FF" }:
            res = "features.txt"
            break
        case { it == "FFPE" }:
            res = "features-*.parquet"
            break
        default:
            res = "*"
            break
    }
    return res
}

process POST_PROCESSING_FEATURES {
    label 'process_medium'
    label "extract_histo_patho_features"


    input:
    path bot_train_file
    val slide_type
    val is_tcga

    output:
    path out_prefix, emit: txt_parquet
    path "ok.txt"

    script:
    is_tcga_numeric = is_tcga ? 1 : 0
    out_prefix = determine_out_prefix(slide_type)
    """
    post_process_features.py \
        --bot_train_file ${bot_train_file} \
        --slide_type ${slide_type} \
        --is_tcga ${is_tcga_numeric} && touch ok.txt
    """
}

