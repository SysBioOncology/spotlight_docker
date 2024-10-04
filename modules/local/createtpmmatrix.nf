process CREATE_TPM_MATRIX {
    label 'rcontainer'

    input: 
        path gene_exp_path

    output:
        path "tpm.txt", emit: txt

    // Run only if 'params.is_tpm'
    when: (task.ext.when || task.ext.when == null) && (gene_exp_path.name != "NO_FILE")

    script: 
    """

    /usr/local/bin/_entrypoint.sh Rscript ${projectDir}/bin/create_tpm_matrix.R --gene_exp_path ${gene_exp_path}
    """

    stub: 
    """
    touch tpm.txt
    """
}