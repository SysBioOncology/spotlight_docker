#!/usr/bin/env RScript
# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

pacman::p_load(GaitiLabUtils, install = FALSE)

# Load libraries
pacman::p_load(glue, data.table, tidyverse, stringr, install = FALSE)

# Define input arguments when running from bash
parser <- setup_default_argparser(
    description = "Create TPM matrix", default_output = "."
)
parser$add_argument("--gene_exp_path", type = "character", help = "Path to gene expression")

args <- parser$parse_args()

# Set up logging
logr <- init_logging(log_level = args$log_level)

log_info("Create output directory...")
GaitiLabUtils::create_dir(args$output_dir)

log_info("Load gene-expression...")
gene_exp <- data.table::fread(
    args$gene_exp_path,
    check.names = FALSE,
    sep = "\t",
    header = TRUE
) %>% data.frame(row.names = 1)

gene_exp_mat <- gene_exp[, gene_exp["gene_id", ] == "scaled_estimate"]
gene_exp_mat <- gene_exp_mat[-1, ]

# Extract gene symbols from Hybridization Ref (rownames)
genes <- str_split(rownames(gene_exp_mat), "\\|", simplify = TRUE)[, 1]

# Remove rows without gene symbols
gene_exp_mat <- gene_exp_mat[genes != "?", ]

# Convert all columns to numeric data type
gene_exp_mat <- gene_exp_mat %>%
    mutate_if(is.character, as.numeric) %>%
    mutate(gene = str_split(rownames(gene_exp_mat), "\\|", simplify = TRUE)[, 1]) %>%
    remove_rownames() %>%
    group_by(gene) %>%
    summarize(across(where(is.numeric), mean))

log_info("Save gene exp...")
write.table(gene_exp_mat,
    file = file.path(args$output_dir, "tpm.txt"),
    row.names = FALSE,
    col.names = TRUE, sep = "\t"
)
log_info("COMPLETED!")
