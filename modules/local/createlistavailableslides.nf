process CREATE_LIST_AVAIL_SLIDES {
    label 'process_single'
    label "extract_histo_patho_features"

    input:
    path clinical_file_path
    path image_dir

    output:
    path "avail_slides_for_img.csv", emit: csv

    script:
    """
    create_list_avail_img_for_tiling.py \
        --slides_folder ${image_dir} \
        --clinical_file_path ${clinical_file_path}

    """

    stub: 
    """
    echo slide_id,slide_filename > "avail_slides_for_img.csv"
    echo slide_id1,slide1.svs >> "avail_slides_for_img.csv"
    echo slide_id2,slide2.svs >> "avail_slides_for_img.csv"

    """

}
