
process TILING_SINGLE_SLIDE {
    label 'process_medium'
    label 'extract_histo_patho_features'

    input:
    tuple val(slide_id), val(filename_slide), path(slide_path)
    val gradient_mag_filter

    output:
    path "${slide_id}_*.jpg", emit: jpg

    script:
    """
    create_tiles_from_slides.py \\
        --slide_path ${slide_path} \\
        --gradient_mag_filter ${gradient_mag_filter} \\
        --filename_slide ${filename_slide}

    """

    stub: 
    """
    touch ${slide_id}_tile_1.jpg
    touch ${slide_id}_tile_1.jpg

    """
}
