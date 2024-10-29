# Integrating histopathology and transcriptomics for spatial tumor microenvironment profiling in a melanoma case study

> Pipeline has been successfully tested with fresh frozen (FF) slides only. Currently, only the TF models for SKCM (melanoma) are readily available in /data/TF_models.

Our pipeline (SPoTLIghT) to derive spatial graph-based interpretable features from H&E (fresh-frozen, FF) tissue slides is available as a Nextflow pipeline.

The pipeline comprises the following steps:
1. Extract 1,536 histopathological features from Inception V4 model.
2. Predict tile level abundances for the different cell types.
3. Compute spatial features.

See also the figures below.

![](src/spotlight_a.jpg)
![](src/spotlight_b.jpg)

## Required software

* Docker: `27.2.0`
* Apptainer: `1.0.2`
* Nextflow: `24.04.4 build 5917`

> These were the versions used for testing the pipeline.

## Run SPoTLIghT

*Use the SKCM multi-task models to infer spatial features*

1. Create apptainer/singularity container from Docker image:

```bash
# 1. save docker as tar or tar.gz (compressed)
docker save joank23/spotlight -o spotlight.tar.gz
# 2. build apptainer (.sif) from docker (.tar)
apptainer build spotlight.sif docker-archive:spotlight.tar.gz

# 1. save docker as tar or tar.gz (compressed)
docker save joank23/immunedeconvr -o immunedeconvr.tar.gz
# 2. build apptainer (.sif) from docker (.tar)
apptainer build immunedeconvr.sif docker-archive:immunedeconvr.tar.gz

```

> Please rename your images file names, so they only include "-", to follow the same sample coding used by the TCGA.

2. Download retrained models to extract the histopathological features, available from Fu et al., Nat Cancer, 2020 ([Retrained_Inception_v4](https://www.ebi.ac.uk/biostudies/bioimages/studies/S-BSST292)). 
Once you unzip the folder, extract the files to the `data/checkpoint/Retrained_Inception_v4/` folder.
3. If a TCGA dataset is used, please download metadata (i.e. "biospecimen -> TSV", unzip and keep slide.tsv), then rename `slide.tsv` to `clinical_file_TCGA_{cancer_type_abbrev}` such as `clinical_file_TCGA_SKCM.tsv` and copy to `/data`. Example dataset TCGA-SKCM can be downloaded [here](https://portal.gdc.cancer.gov/projects/TCGA-SKCM). For non-TCGA datasets, please omit this step.
4. Setup your paths and variables in `run_pipeline.sh`
5. Set a config ensuring compatibility with available resources, you can use `nf-custom.config` as a template. (see also `nextflow.config` for default settings). Set also the paths to the containers.
6. Please set parameters in 'nf-params.yml', if a parameter is 'assets/NO_FILE' or 'dummy', they are optional parameters, if not used please leave as is (see also `nextflow.config` for default settings)
7. Run the Nextflow Pipeline as follows:

```bash

nextflow run . -profile apptainer -c "${your_config_file}" -outdir ${your_output_directory} -params-file "nf-params.yml"

````

> For more information, incl. building your own multi-task model for a TCGA dataset (FF slides), please read the [docs](docs/)

## Output documentation

SPoTLIghT generates the following output directory structure for a run with FFPE slides:

> Of note, when using FF slides, the directories `features_format_parquet` and `predictions_format_parquet` won't be created, instead, the files 'features.txt' and 'predictions.txt' are created.

```bash
{outdir}
├── 1_extract_histopatho_features
│   ├── avail_slides_for_img.csv
│   ├── bot_train.txt
│   ├── features_format_parquet
│   │   ├── features-0.parquet
│   │   ├── features-1.parquet
│   ├── file_info_train.txt
│   ├── generated_clinical_file.txt
│   ├── pred_train.txt
│   ├── predictions_format_parquet
│   │   ├── predictions-0.parquet
│   │   ├── predictions-0.parquet
│   ├── process_train
│   │   ├── images_train_00001-of-00320.tfrecord
│   │   ├── images_train_00002-of-00320.tfrecord 
│   │   ├── images_train_00004-of-00320.tfrecord
│   └── tiles
│       ├── xenium-skin-panel_10165_10165.jpg
│       ├── xenium-skin-panel_10165_10627.jpg
│       ├── xenium-skin-panel_10165_11089.jpg
├── 2_deconv_bulk_rnaseq
│   ├── epic.csv
│   ├── mcp_counter.csv
│   ├── quantiseq.csv
│   ├── tpm.txt
│   └── xcell.csv
├── 3_build_multi_task_celltype_model
│   ├── CAFs
│   │   ├── cv_outer_splits.pkl
│   │   ├── outer_models.pkl
│   │   ├── outer_scores_slides_test.pkl
│   │   ├── outer_scores_slides_train.pkl
│   │   ├── outer_scores_tiles_test.pkl
│   │   ├── outer_scores_tiles_train.pkl
│   │   ├── total_tile_selection.pkl
│   │   ├── x_train_scaler.pkl
│   │   └── y_train_scaler.pkl
│   ├── T_cells
│   │   ├── cv_outer_splits.pkl
│   │   ├── outer_models.pkl
│   │   ├── outer_scores_slides_test.pkl
│   │   ├── outer_scores_slides_train.pkl
│   │   ├── outer_scores_tiles_test.pkl
│   │   ├── outer_scores_tiles_train.pkl
│   │   ├── total_tile_selection.pkl
│   │   ├── x_train_scaler.pkl
│   │   └── y_train_scaler.pkl
│   ├── endothelial_cells
│   │   ├── cv_outer_splits.pkl
│   │   ├── outer_models.pkl
│   │   ├── outer_scores_slides_test.pkl
│   │   ├── outer_scores_slides_train.pkl
│   │   ├── outer_scores_tiles_test.pkl
│   │   ├── outer_scores_tiles_train.pkl
│   │   ├── total_tile_selection.pkl
│   │   ├── x_train_scaler.pkl
│   │   └── y_train_scaler.pkl
│   ├── ensembled_selected_tasks.csv
│   ├── task_selection_names.pkl
│   └── tumor_purity
│       ├── cv_outer_splits.pkl
│       ├── outer_models.pkl
│       ├── outer_scores_slides_test.pkl
│       ├── outer_scores_slides_train.pkl
│       ├── outer_scores_tiles_test.pkl
│       ├── outer_scores_tiles_train.pkl
│       ├── total_tile_selection.pkl
│       ├── x_train_scaler.pkl
│       └── y_train_scaler.pkl
├── 4_tile_level_quantification
│   ├── test_tile_predictions_proba.csv
│   └── test_tile_predictions_zscores.csv
├── 5_spatial_features
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
└── pipeline_info
    ├── execution_report_2024-09-23_21-07-41.html
    ├── execution_timeline_2024-09-23_21-07-41.html
    ├── execution_trace_2024-09-23_21-07-41.txt
    └── pipeline_dag_2024-09-23_21-07-41.html
```
