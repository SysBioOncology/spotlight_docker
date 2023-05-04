#!/bin/bash
repo_path="/home/olapuent/Desktop/spatial_localization/repo_manuscript/spotlight"
cd $repo_path/run_scripts

# ---- Setup virtual environment ---- #

bash setup_virtualenv.sh

# ---- Run full pipeline (step by step) ---- #

bash 1_extract_histopatho_features.sh

bash 2_TF_training.sh

bash 3_tile_level_cell_type_quantification.sh

bash 4_compute_spatial_features.sh
