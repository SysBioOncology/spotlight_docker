# Integrating histopathology and transcriptomics for spatial profiling of the tumor microenvironment: a melanoma case study

> Pipeline has been successfully tested with fresh frozen (FF) slides only. Currently, only the TF models for SKCM (melanoma) are readily available in /data/TF_models.

Our pipeline (SPoTLIghT) to derive spatial graph-based interpretable features from H&E (fresh-frozen, FF) tissue slides is packaged as a docker or a singularity/apptainer container.

The pipeline comprises the following steps:
1. Extract 1,536 histopathological features from Inception V4 model.
2. Predict tile level abundances for the different cell types.
3. Compute spatial features.

See also the figures below.

![](src/spotlight_a.jpg)
![](src/spotlight_b.jpg)

## Run SPoTLIghT

1. Pull the Docker container:

```bash
docker pull joank23/spotlight:latest
```

Alternatively you can use Singularity/Apptainer (HPCs):

```bash
# 1. save docker as tar or tar.gz (compressed)
docker save joank23/spotlight -o spotlight.tar.gz
# 2. build apptainer (.sif) from docker (.tar)
apptainer build spotlight.sif docker-archive:spotlight.tar.gz
```

> Please rename your images file names, so they only include "-", to follow the same sample coding used by the TCGA.

1. Download retrained models to extract the histopathological features, available from Fu et al., Nat Cancer, 2020 ([Retrained_Inception_v4](https://www.ebi.ac.uk/biostudies/bioimages/studies/S-BSST292)). 
Once you unzip the folder, extract the files to the `data/checkpoint/Retrained_Inception_v4/` folder.
2. If a TCGA dataset is used, please download metadata (i.e. "biospecimen -> TSV", unzip and keep slide.tsv), then rename `slide.tsv` to `clinical_file_TCGA_{cancer_type_abbrev}` such as `clinical_file_TCGA_SKCM.tsv` and copy to `/data`. Example dataset TCGA-SKCM can be downloaded [here](https://portal.gdc.cancer.gov/projects/TCGA-SKCM). For non-TCGA datasets, please omit this step.
3. Setup your paths and variables in `run_pipeline.sh`
4. Set a config ensuring compatibility with available resources, you can use `custom.config` as a template. (see `nextflow.config` for all parameters, if a parameter is 'assets/NO_FILE' or 'dummy', they are optional parameters, if not used please leave as is)
5. Run the Nextflow Pipeline as follows:

```bash

nextflow run . -profile apptainer -c "${your_config_file}"

````

## Output documentation

SPoTLIghT generates the following output directory structure:

```bash
{outdir}
├── 1_extract_histopatho_features
│   ├── avail_slides_for_img.csv
│   ├── features-0.parquet
│   ├── file_info_train.txt
│   ├── generated_clinical_file.txt
│   ├── ok.txt
│   ├── predictions-0.parquet
│   ├── process_train
│   │   ├── images_train_00001-of-00320.tfrecord
│   │   ├── images_train_00002-of-00320.tfrecord 
│   │   ├── images_train_00004-of-00320.tfrecord
│   └── tiles
│       ├── xenium-skin-panel_10165_10165.jpg
│       ├── xenium-skin-panel_10165_10627.jpg
│       ├── xenium-skin-panel_10165_11089.jpg
├── 2_tile_level_quantification
│   ├── test_tile_predictions_proba.csv
│   └── test_tile_predictions_zscores.csv
├── 3_spatial_features
│   ├── clustering_features
│   │   ├── FFPE_all_schc_clusters_labeled.csv
│   │   ├── FFPE_all_schc_tiles.csv
│   │   ├── FFPE_all_schc_tiles_raw.csv
│   │   ├── FFPE_features_clust_all_schc_prox_wide.csv
│   │   ├── FFPE_features_clust_indiv_schc_prox_between.csv
│   │   ├── FFPE_features_clust_indiv_schc_prox.csv
│   │   ├── FFPE_features_clust_indiv_schc_prox_within.csv
│   │   ├── FFPE_frac_high_wide.csv
│   │   ├── FFPE_graphs.pkl
│   │   ├── FFPE_indiv_schc_clusters_labeled.csv
│   │   ├── FFPE_indiv_schc_tiles.csv
│   │   ├── FFPE_indiv_schc_tiles_raw.csv
│   │   └── FFPE_nclusters_wide.csv
│   ├── FFPE_all_features_combined.csv
│   ├── FFPE_all_graph_features.csv
│   ├── FFPE_clustering_features.csv
│   ├── FFPE_graphs.pkl
│   └── network_features
│       ├── FFPE_features_coloc_fraction.csv
│       ├── FFPE_features_coloc_fraction_wide.csv
│       ├── FFPE_features_lcc_fraction_wide.csv
│       ├── FFPE_features_ND.csv
│       ├── FFPE_features_ND_ES.csv
│       ├── FFPE_features_ND_sim_assignments.pkl
│       ├── FFPE_features_ND_sims.csv
│       ├── FFPE_features_shortest_paths_thresholded.csv
│       ├── FFPE_features_shortest_paths_thresholded_wide.csv
│       ├── FFPE_graphs.pkl
│       └── FFPE_shapiro_tests.csv
├── bottleneck
│   ├── bot_train.txt
│   ├── ok.txt
│   └── pred_train.txt
└── pipeline_info
    ├── execution_report_2024-09-23_21-07-41.html
    ├── execution_timeline_2024-09-23_21-07-41.html
    ├── execution_trace_2024-09-23_21-07-41.txt
    └── pipeline_dag_2024-09-23_21-07-41.html
```
