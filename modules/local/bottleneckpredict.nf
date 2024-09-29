
// # ------------------------------------------------------ #
// # ---- Compute predictions and bottlenecks features ---- #
// # ------------------------------------------------------ #

process BOTTLENECK_PREDICT {
    label 'mem_64G'
    label 'time_4h'
    label 'extract_histo_patho_features'

    input:
    val bot_out
    val pred_out
    path tf_records
    val model_name
    val checkpoint_path

    output:
    path "${bot_out}.txt", emit: bot_txt
    path "${pred_out}.txt", emit: pred_txt
    path "ok.txt"

    script:
    """
    bottleneck_predict.py \
        --num_classes=42 \
        --bot_out ${bot_out}.txt \
        --pred_out ${pred_out}.txt \
        --model_name ${model_name} \
        --checkpoint_path ${checkpoint_path} \
        --file_dir \$PWD/ && touch \$PWD/ok.txt

    """

    stub: 
    """
    touch ${bot_out}.txt
    touch ${pred_out}.txt
    touch ok.txt
    """
}
