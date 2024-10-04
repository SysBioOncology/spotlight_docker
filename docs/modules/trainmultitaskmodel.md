# Building a multi-task model for a cancer type using TCGA bulkRNAseq data

## Set up  

1. Download TCGA bulkRNAseq data via the [Firehose Tool](https://gdac.broadinstitute.org) from the BROAD Institute, the files required are: “illuminahiseq_rnaseqv2-RSEM_genes”.
2. Unzip the downloaded file (.tar.gz)
3. Download the signatures/published scores, see table below.

| Parameter                  | Reference                                                 | Additional info                                                                                                                                                |
| -------------------------- | --------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `absolute_tumor_purity_path` | https://gdc.cancer.gov/about-data/publications/panimmune  | Download the 'Score for 160 Genes Signatures in Tumor Samples' or use [direct link]( https://api.gdc.cancer.gov/data/80a82092-161d-4615-9d96-e858f113618d)        |
| `estimate_scores_path` | https://bioinformatics.mdanderson.org/estimate/index.html | Download the relevant file for the cancer type of interest, use the RNA-seqV2 column on the page.                                                                                                    |
| `gibbons_scores_path` | https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5503821/     | Download the 'Supp Datafile S1.' or use the (direct link)[https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5503821/bin/NIHMS840944-supplement-Supp_Datafile_S1.xlsx] |
| `thorsson_scores_path` | https://gdc.cancer.gov/about-data/publications/panimmune  | Download the 'ABSOLUTE purity/ploidy file', or use [direct link](https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5)                           |
