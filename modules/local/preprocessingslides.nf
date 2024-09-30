
process PREPROCESSING_SLIDES {
    label 'process_medium'
    label "extract_histo_patho_features"

    input:
    path file_info_train
    val n_shards

    output:
    path "images_train_*-of-*.tfrecord", emit: tfrecords

    script:
    """
    pre_processing.py \
        --file_info_train ${file_info_train} \
        --N_shards ${n_shards}\
    """


    stub: 
    """

    touch "images_train_00001-of-00${n_shards}.tfrecord"
    touch "images_train_00${n_shards}-of-00${n_shards}.tfrecord"

    """

}
