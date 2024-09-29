process PREPROCESSING_MULTITASK_MODEL_TARGET_FEATURES {
    input: 
    path clinical_file_path
    path tpm_path
    path thorsson_signatures_path
    path estimate_signatures_path
    path absolute_tumor_purity_path
    path gibbons_signatures_path
    path mcp_counter_path
    path quantiseq_path
    path xcell_path
    path epic_path

    output:
        path "ensembled_selected_tasks.csv", emit: csv
        path "task_selection_names.pkl", emit: pkl

    script: 
    """
    preprocessing_multitask_model_target_features \\
        --clinical_file_path ${clinical_file_path} \\
        --tpm_path ${tpm_path} \\
        --thorsson_signatures_path ${thorsson_signatures_path} \\
        --estimate_signatures_path ${estimate_signatures_path} \\
        --absolute_tumor_purity_path ${absolute_tumor_purity_path} \\
        --gibbons_signatures_path ${gibbons_signatures_path} \\
        --mcp_counter_path ${mcp_counter_path} \\
        --quantiseq_path ${quantiseq_path} \\
        --xcell_path ${xcell_path} \\
        --epic_path ${epic_path} \\


    """

    stub: 
    """

    touch "ensembled_selected_tasks.csv"
    touch "task_selection_names.pkl"
    """
}