//
// Subworkflow with functionality specific to the SysBioOncology/spotlight_docker pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CREATE_TPM_MATRIX } from '../../modules/local/createtpmmatrix.nf'
include { IMMUNEDECONV as QUANTISEQ } from '../../modules/local/immunedeconv.nf'
include { IMMUNEDECONV as MCP_COUNTER } from '../../modules/local/immunedeconv.nf'
include { IMMUNEDECONV as XCELL } from '../../modules/local/immunedeconv.nf'
include { IMMUNEDECONV as EPIC } from '../../modules/local/immunedeconv.nf'

workflow CELLTYPE_QUANTIFICATION_BULKRNASEQ {
    take: 
    gene_exp_path

    main: 
    


    // CREATE_TPM_MATRIX(gene_exp_path)
    // tpm_path = is_tpm ? Channel.value(gene_exp_path) : CREATE_TPM_MATRIX.out.txt
    
    // // Immune deconvolution
    // Channel.of(                                 // channel: [tool name, csv]
    //     ["quantiseq", quantiseq_path],
    //     ["mpc_counter", mpc_counter_path],
    //     ["xcell", xcell_path],
    //     ["epic", epic_path]
    // )
    // .branch{
    //     tool, filepath ->  
    //     // Tools to run 
    //     invalid: deconv_tools.contains(tool) && (filepath.name == "NO_FILE")
    //         return [tool, tpm_path]
    //     // Tools already used
    //     valid: true
    //         return [tool, filepath]
    // }.set {
    //     ch_immune_deconv_files
    // }

    // IMMUNEDECONV (
    //     ch_immune_deconv_files.invalid
    // ).collect()


    // emit: 
    //     tpm_path        = tpm_path
    //     quantiseq       = QUANTISEQ.out == null ? quantiseq_path : QUANTISEQ.out.csv
    //     mcp_counter     = MCP_COUNTER.out == null ? mcp_counter_path : MCP_COUNTER.out.csv
    //     xcell           = XCELL.out == null ? xcell_path : XCELL.out.csv
    //     epic            = EPIC.out == null ? epic_path : EPIC.out.csv
}