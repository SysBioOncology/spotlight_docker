# SysBioOncology/spotlight_docker: Output

## Introduction

This document describes the output produced by the pipeline. 

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [SysBioOncology/spotlight\_docker: Output](#sysbiooncologyspotlight_docker-output)
  - [Introduction](#introduction)
  - [Pipeline overview](#pipeline-overview)
    - [Extract histopathological features](#extract-histopathological-features)
    - [Deconvolution of bulkRNAseq data](#deconvolution-of-bulkrnaseq-data)
    - [Building a multi-task cell type model to predict cell type abundances on a tile-level](#building-a-multi-task-cell-type-model-to-predict-cell-type-abundances-on-a-tile-level)
    - [Predicting tile-level cell type abundances using the multi-task models](#predicting-tile-level-cell-type-abundances-using-the-multi-task-models)
    - [Compute spatial features using the tile-level cell type abundances](#compute-spatial-features-using-the-tile-level-cell-type-abundances)
    - [Pipeline information](#pipeline-information)

### Extract histopathological features

<details markdown="1">
<summary>Output files</summary>

* `1_extract_histopatho_features/`
  + `avail_slides_for_img.csv`
  + `bot_train.txt`
  + `features_format_parquet`
    - `features-0.parquet`
    - `features-1.parquet`
  + `file_info_train.txt`
  + `generated_clinical_file.txt`
  + `pred_train.txt`
  + `predictions_format_parquet`
    - `predictions-0.parquet`
    - `predictions-0.parquet`
  + `process_train/`
    - `images_train_00001-of-00320.tfrecord`
    - `images_train_00002-of-00320.tfrecord`
    - `images_train_00004-of-00320.tfrecord`
  + `tiles/`
    - `xenium-skin-panel_10165_10165.jpg`
    - `xenium-skin-panel_10165_10627.jpg`
    - `xenium-skin-panel_10165_11089.jpg`

</details>

### Deconvolution of bulkRNAseq data

<details markdown="1">
<summary>Output files</summary>

* `2_deconv_bulk_rnaseq/`
  + `epic.csv`
  + `mcp_counter.csv`
  + `quantiseq.csv`

  + `tpm.txt`

  + `xcell.csv`

</details>

### Building a multi-task cell type model to predict cell type abundances on a tile-level 

<details markdown="1">
<summary>Output files</summary>

* `3_build_multi_task_celltype_model/`
  + `CAFs/`
    - `cv_outer_splits.pkl`
    - `outer_models.pkl`
    - `outer_scores_slides_test.pkl`
    - `outer_scores_slides_train.pkl`
    - `outer_scores_tiles_test.pkl`
    - `outer_scores_tiles_train.pkl`
    - `total_tile_selection.pkl`
    - `x_train_scaler.pkl`
    - `y_train_scaler.pkl`
  + `T_cells/`
      - `cv_outer_splits.pkl`
    - `outer_models.pkl`
    - `outer_scores_slides_test.pkl`
    - `outer_scores_slides_train.pkl`
    - `outer_scores_tiles_test.pkl`
    - `outer_scores_tiles_train.pkl`
    - `total_tile_selection.pkl`
    - `x_train_scaler.pkl`
    - `y_train_scaler.pkl`
  + `endothelial_cells/`
    - `cv_outer_splits.pkl`
    - `outer_models.pkl`
    - `outer_scores_slides_test.pkl`
    - `outer_scores_slides_train.pkl`
    - `outer_scores_tiles_test.pkl`
    - `outer_scores_tiles_train.pkl`
    - `total_tile_selection.pkl`
    - `x_train_scaler.pkl`
    - `y_train_scaler.pkl`
  + `ensembled_selected_tasks.csv`
  + `task_selection_names.pkl`
  + `tumor_purity/`
    - `cv_outer_splits.pkl`
    - `outer_models.pkl`
    - `outer_scores_slides_test.pkl`
    - `outer_scores_slides_train.pkl`
    - `outer_scores_tiles_test.pkl`
    - `outer_scores_tiles_train.pkl`
    - `total_tile_selection.pkl`
    - `x_train_scaler.pkl`
    - `y_train_scaler.pkl`
</details>

### Predicting tile-level cell type abundances using the multi-task models

<details markdown="1">
<summary>Output files</summary>

* `4_tile_level_quantification/`
  + `test_tile_predictions_proba.csv`
  + `test_tile_predictions_zscores.csv`

</details>

### Compute spatial features using the tile-level cell type abundances

<details markdown="1">
<summary>Output files</summary>

* `5_spatial_features/`
  + `clustering_features/`
    - `FFPE_all_schc_clusters_labeled.csv`

    - `FFPE_all_schc_tiles.csv`

    - `FFPE_all_schc_tiles_raw.csv`

    - `FFPE_features_clust_all_schc_prox_wide.csv`

    - `FFPE_features_clust_indiv_schc_prox_between.csv`

    - `FFPE_features_clust_indiv_schc_prox.csv`

    - `FFPE_features_clust_indiv_schc_prox_within.csv`

    - `FFPE_frac_high_wide.csv`

    - `FFPE_graphs.pkl`

    - `FFPE_indiv_schc_clusters_labeled.csv`

    - `FFPE_indiv_schc_tiles.csv`

    - `FFPE_indiv_schc_tiles_raw.csv`

    - `FFPE_nclusters_wide.csv`

  + `FFPE_all_features_combined.csv`

  + `FFPE_all_graph_features.csv`

  + `FFPE_clustering_features.csv`

  + `FFPE_graphs.pkl`

  + `network_features/`
    - `FFPE_features_coloc_fraction.csv`

    - `FFPE_features_coloc_fraction_wide.csv`

    - `FFPE_features_lcc_fraction_wide.csv`

    - `FFPE_features_ND.csv`

    - `FFPE_features_ND_ES.csv`

    - `FFPE_features_ND_sim_assignments.pkl`

    - `FFPE_features_ND_sims.csv`

    - `FFPE_features_shortest_paths_thresholded.csv`

    - `FFPE_features_shortest_paths_thresholded_wide.csv`

    - `FFPE_graphs.pkl`

    - `FFPE_shapiro_tests.csv`

</details>

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

* `pipeline_info/`
  + Reports generated by Nextflow: `execution_report.html`,                                                                  `execution_timeline.html`,  `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  + Reports generated by the pipeline: `pipeline_report.html`,  `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
