process CREATE_CLINICAL_FILE {
    label 'process_single'
    label "extract_histo_patho_features"

    input:
    path clinical_files_input
    val class_name
    val out_file
    path path_codebook
    val tumor_purity_threshold
    val is_tcga
    path image_dir

    output:
    path "${out_file}.txt", emit: txt

    script:
    def args   = task.ext.args   ?: ''
    def missing_clinical_file = clinical_files_input.name == 'NO_FILE'
    if (is_tcga | !missing_clinical_file ) {
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
        awk -v a=81 -v b="${class_name}" -v c=41 'FNR==NR{print; next}{split(\$1, tmp, "."); OFS="\t"; print tmp[1], tmp[1], \$1, a, b, c}' ${template_txt} ${list_txt} > ${out_file}.txt

        """
    }
}
