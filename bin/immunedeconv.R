#!/usr/bin/env Rscript
# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))


pacman::p_load(GaitiLabUtils, immunedeconv, install = FALSE)
pacman::p_load(glue, data.table, tidyverse, stringr, install = FALSE)

parser <- setup_default_argparser(
    description = "Compute cell fractions with immunedeconv using tpm from bulkRNAseq data", default_output = "."
)
parser$add_argument("--tpm_path", type = "character", default = NULL, help = "Full path of input file incl. extension")
parser$add_argument("--tool", type = "character", default = "quantiseq", help = "Tool to use for quantification (default=quantiseq)")
parser$add_argument("--probesets", type = "character", default = "assets/mcp_counter/probesets.txt", help = "Path to probesets.txt")
parser$add_argument("--genes", type = "character", default = "assets/mcp_counter/genes.txt", help = "Path to genes.txt")
args <- parser$parse_args()

# args <- list()
# args$tpm_path <- "output_test_tcga/create/tpm.txt"
# args$probesets <- "assets/mcp_counter/probesets.txt"
# args$genes <- "assets/mcp_counter/genes.txt"
# args$output_dir <- "/Users/joankant/Desktop/gaitigroup/Users/Joan/spotlight_docker/output_test_tcga/immunedeconv"


# # Set up logging
logr <- init_logging(log_level = args$log_level)

log_info("Create output directory...")
GaitiLabUtils::create_dir(args$output_dir)

log_info(getwd())
log_info(list.files(getwd()))
# immunedeconv offers the following methods:
# - quantiseq
# - mcp_counter
# - xcell
# - epic

log_info("Load tpm data")
tpm <- data.table::fread(args$tpm_path,
    check.names = FALSE,
    sep = "\t",
    header = TRUE,
) %>%
    data.frame(row.names = 1) %>%
    data.matrix()

create_dir(args$output_dir)

## 1) quanTIseq
if (args$tool == "quantiseq") {
    log_info("Running quanTIseq...")
    tpm_quantiseq <- tpm
    rownames(tpm_quantiseq) <- toupper(rownames(tpm_quantiseq))
    cell_fractions <- immunedeconv::deconvolute(tpm_quantiseq, "quantiseq", tumor = TRUE, arrays = FALSE, scale_mrna = TRUE) %>% column_to_rownames("cell_type")

    write.csv(cell_fractions, file.path(args$output_dir, "quantiseq.csv"), row.names = TRUE)
} else if (args$tool == "mcp_counter") {
    ## 2) MCP Counter
    log_info("Running MCP Counter...")
    user_probesets <- read.table(args$probesets, sep = "\t", stringsAsFactors = FALSE, colClasses = "character")
    user_genesets <- read.table(args$genes, sep = "\t", stringsAsFactors = FALSE, header = TRUE, colClasses = "character", check.names = FALSE)

    cell_fractions <- immunedeconv::deconvolute(tpm, "mcp_counter",
        probesets = user_probesets,
        genes = user_genesets
    ) %>% column_to_rownames("cell_type")


    write.csv(cell_fractions, file.path(args$output_dir, "mcp_counter.csv"), row.names = TRUE)
} else if (args$tool == "xcell") {
    ## 3) XCell
    log_info("Running XCell...")
    cell_fractions <- immunedeconv::deconvolute(tpm, "xcell") %>% column_to_rownames("cell_type")

    write.csv(cell_fractions, file.path(args$output_dir, "xcell.csv"), row.names = TRUE)
} else if (args$tool == "epic") {
    ## 4) EPIC
    log_info("Running EPIC...")
    cell_fractions <- immunedeconv::deconvolute_epic(tpm, tumor = TRUE, scale_mrna = TRUE)

    write.csv(cell_fractions, file.path(args$output_dir, "epic.csv"), row.names = TRUE)
} else {
    log_info("No valid tool given, please choose: 'quantiseq', 'mcp', 'xcell' or 'epic'...")
}

log_info("COMPLETED!")
