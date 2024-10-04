process CREATE_CLINICAL_FILE {
    label 'process_single'
    label "extract_histo_patho_features"

    input:
    path clinical_files_input
    val class_name
    val out_prefix
    path path_codebook
    val tumor_purity_threshold
    val is_tcga
    path image_dir
    val slide_type

    output:
    path "${out_prefix}.txt", emit: txt

    script:
    def args   = task.ext.args   ?: ''
    out_prefix = task.ext.prefix ? task.ext.prefix : out_prefix

    if ((is_tcga && slide_type == "FF")) {
        """
        create_clinical_file.py \\
            --class_name ${class_name} \\
            --clinical_files_input ${clinical_files_input} \\
            --tumor_purity_threshold ${tumor_purity_threshold} \\
            --path_codebook ${path_codebook}
        """
    } else {
        template_txt="${projectDir}/assets/tmp_clinical_file.txt"
        list_txt =  "list_images.txt"
        """
        ls ${image_dir} | tee ${list_txt}
        awk -v a=81 -v b="${class_name}" -v c=41 'FNR==NR{print; next}{split(\$1, tmp, "."); OFS="\t"; print tmp[1], tmp[1], \$1, a, b, c}' ${template_txt} ${list_txt} > ${out_prefix}.txt

        """
    }

    stub: 
    """
    touch ${out_prefix}.txt
    """
}
