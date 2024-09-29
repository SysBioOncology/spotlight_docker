process IMMUNEDECONV {
    label 'rcontainer'

    containerOptions {"--mount type=bind,src=\$PWD,dst=\$PWD --mount type=tmpfs,destination=/tmpdir -w \$PWD" }

    input: 
        tuple val(tool), path(tpm_path)

    output: 
        tuple val(tool), path("${tool}.csv"), emit: csv

    script: 
    """
    immunedeconv.R \\
        --tpm_path ${tpm_path} \\
        --tool ${tool}
    """    
    

    stub:
    """
    touch "${tool}.csv"
    """

}