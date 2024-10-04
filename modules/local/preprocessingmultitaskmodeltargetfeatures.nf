process PREPROCESSING_MULTITASK_MODEL_TARGET_FEATURES {
    input: 
        path clinical_file_path
        path tpm_path
        path thorsson_scores_path
        path estimate_scores_path
        path absolute_tumor_purity_path
        path gibbons_scores_path
        path mfp_gene_signatures_path
        path mcp_counter_path
        path quantiseq_path
        path xcell_path
        path epic_path

    output:
        path "ensembled_selected_tasks.csv", emit: csv
        path "task_selection_names.pkl", emit: pkl

    script: 
    """
    preprocessing_multitask_model_target_features.py \\
        --clinical_file_path ${clinical_file_path} \\
        --tpm_path ${tpm_path} \\
        --thorsson_scores_path ${thorsson_scores_path} \\
        --estimate_scores_path ${estimate_scores_path} \\
        --absolute_tumor_purity_path ${absolute_tumor_purity_path} \\
        --gibbons_scores_path ${gibbons_scores_path} \\
        --mcp_counter_path ${mcp_counter_path} \\
        --quantiseq_path ${quantiseq_path} \\
        --xcell_path ${xcell_path} \\
        --epic_path ${epic_path} \\
        --mfp_gene_signatures_path ${mfp_gene_signatures_path}
    """
    
    stub: 
    """

    touch "ensembled_selected_tasks.csv"
    touch "task_selection_names.pkl"
    """
}