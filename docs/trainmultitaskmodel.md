# Building a multi-task model for a cancer type using TCGA bulkRNAseq data

## Prerequisites 

1. Download TCGA bulkRNAseq data via the [Firehose Tool](https://gdac.broadinstitute.org) from the BROAD Institute, the files required are: “illuminahiseq_rnaseqv2-RSEM_genes”.
2. Unzip the downloaded file (.tar.gz)
3. Download the signatures/published scores. 

Scores and clinical data retrieved from publications
| file name | short description | publication | source | 
| --------- | ------------ | -------------- | -------------- |
| [ImmuneSubtypes_Thorsson.xlsx](data/ImmuneSubtypes_Thorsson.xlsx) | Six immune subtypes | https://doi.org/10.1016/j.immuni.2018.03.023 | https://www.cell.com/cms/10.1016/j.immuni.2018.03.023/attachment/1b63d4bc-af31-4a23-99bb-7ca23c7b4e0a/mmc2.xlsx * |
| [MFP.xlsx](data/MFP.xlsx) | Four TME subtypes (Molecular Functional Portrait MFP) | https://doi.org/10.1016/j.ccell.2021.04.014 | paper supplementary Bagaev_mmc6.xlsx (not open access) | 
| [Fges.tsv](data/Fges.tsv) | 29 functional gene expression signatures (Fges) used to derive the 4 subtypes | https://doi.org/10.1016/j.ccell.2021.04.014 | paper supplementary signatures-tcga.tsv (not open access) + | 
| [ESTIMATE.xlsx](data/ESTIMATE.xlsx) | Estimation of STromal and Immune cells in MAlignant Tumours using Expression data (ESTIMATE) | https://doi.org/10.1038/ncomms3612 | Supplementay Data 2 in https://www.nature.com/articles/ncomms3612#Sec21 |
| [Scores_160_Signatures.tsv](data/Scores_160_Signatures.tsv) | TILs transcriptomic score ("LIexpression_score" in the table) | https://doi.org/10.1016/j.immuni.2018.03.023 | Scores for 160 Genes Signatures in Tumor Samples (Scores_160_Signatures.tsv.gz) https://gdc.cancer.gov/about-data/publications/panimmune |

* All the material related to the Immune Landscape of Cancer publication are available at: https://gdc.cancer.gov/about-data/publications/panimmune
* Additional supplementary information of the Bagaev paper, i.e. information on the genes used to derive Fges and description of the Fges is available in the folder [Bagaev](data/Bagaev/)
