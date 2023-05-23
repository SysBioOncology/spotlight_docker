#!/bin/bash

###############################
## Compute spatial features  ##
###############################

# Define path to spotlight repository
repo_path="/Users/joankant/Library/CloudStorage/OneDrive-TUEindhoven/spotlight"
# repo_path="/home/olapuent/Desktop/spatial_localization/repo_manuscript/spotlight"

# Define type of slide
slide_type="FF"

# Define path to output directory
#output_dir=$repo_path/tests/output/SKCM_${slide_type}
output_dir=$repo_path/tests/output/SKCM_FF_FFP

# Define path to cell-type quantification maps
#tile_quantification_path=$repo_path/tests/output/SKCM_${slide_type}/tcga_validation_tile_predictions_proba.csv
tile_quantification_path=$repo_path/tests/output/SKCM_FF_FFPE/tcga_train_validation_tile_predictions_proba.csv

# ---------------------------------- #
# ---- Compute all features -------- #
# ---------------------------------- #
run_mode=1
python $repo_path/Python/3_spatial_characterization/computing_features.py \
    --workflow_mode=$run_mode \
    --tile_quantification_path=$tile_quantification_path \
    --output_dir=$output_dir \
    --metadata_path=$output_dir/metadata.csv \
    # --slide_type=$slide_type \ # OPTIONAL BY DEFAULT FF
    # --cell_types=$cell_types \   # OPTIONAL
    #--graphs_path=$graphs_path  # OPTIONAL

# # ---------------------------------- #
# # ---- Compute network features ---- #
# # ---------------------------------- #
# workflow=2
# python $repo_path/Python/computing_features.py \
#     --workflow=$workflow \
#     --tile_quantification_path=$tile_quantification_path \
#     --output_dir=$output_dir
#     # --slide_type=$slide_type \ # OPTIONAL BY DEFAULT FF
#     # --cell_types=$cell_types \   # OPTIONAL
#     # --graphs_path=$graphs_path \ # OPTIONAL

# # ------------------------------------- #
# # ---- Compute clustering features ---- #
# # ------------------------------------- #
# workflow=3
# python $repo_path/Python/computing_features.py \
#     --workflow=$workflow \
#     --tile_quantification_path=$tile_quantification_path \
#     --output_dir=$output_dir \
#     # --slide_type=$slide_type \ # OPTIONAL BY DEFAULT FF
#     # --cell_types=$cell_types \   # OPTIONAL
#     # --graphs_path=$graphs_path \ # OPTIONAL


# # -------------------------- #
# # ---- Combine features ---- #
# # -------------------------- #
# workflow=4
# python $repo_path/Python/computing_features.py \
#     --workflow=$workflow \
#     --tile_quantification_path=$tile_quantification_path \
#     --output_dir=$output_dir \
#     # --slide_type=$slide_type \ # OPTIONAL BY DEFAULT FF
#     # --cell_types=$cell_types \   # OPTIONAL
#     # --graphs_path=$graphs_path \ # OPTIONAL
