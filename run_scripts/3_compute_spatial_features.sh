#!/bin/bash

###############################
## Compute spatial features  ##
###############################

# ----------------------------------- #
# --------- Setup file paths -------- #
# ----------------------------------- #

# General setup
repo_dir=$(pwd)

# command line rguments
echo "Cell-type quantification: $1"
echo "Slide type: $2";
echo "Folder output: $3";

# Define input features
tile_quantification_path=$1

# Define type of slide
slide_type=$2

# Define output directory
output_dir=$3

# ---------------------------------- #
# ---- Compute all features -------- #
# ---------------------------------- #
run_mode=1
python $repo_path/Python/3_spatial_characterization/computing_features.py \
    --workflow_mode=$run_mode \
    --tile_quantification_path=$tile_quantification_path \
    --output_dir=$output_dir/3_spatial_features \
    --metadata_path=$output_dir/3_spatial_features/metadata.csv \
    --slide_type=$slide_type # OPTIONAL BY DEFAULT FF
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
