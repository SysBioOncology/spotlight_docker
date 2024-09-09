#!/bin/bash

#####################################################################
## Compute cell-type quantification from transfer learning models  ##
#####################################################################

# ----------------------------------- #
# --------- Setup file paths -------- #
# ----------------------------------- #

# General setup
repo_dir="/project"

# command line arguments
echo "Slide type: $1";

# Define type of slide
slide_type=$1

# Fixed dir
output_dir=${repo_dir}/output
histopatho_features_dir=${output_dir}/1_histopathological_features

# Transfer Learning trained models directory (default: use of FF here)
models_dir=${repo_dir}/data/TF_models/SKCM_FF
var_names_path=${repo_dir}/Python/2_train_multitask_models/task_selection_names.pkl

# Compute predictions using models learned from unseen folds
prediction_mode="test" # (tcga_validation, tcga_train_validation)

echo "Prediction mode: $prediction_mode"

# ---------------------------------------------------- #
# ---- Predict cell type abundances on tile level ---- #
# ---------------------------------------------------- #

# For now, we use models trained on FF slides

python ${repo_dir}/Python/2_train_multitask_models/tile_level_cell_type_quantification.py \
    --models_dir $models_dir \
    --output_dir "$output_dir/2_tile_level_quantification" \
    --histopatho_features_dir $histopatho_features_dir \
    --prediction_mode $prediction_mode \
    --var_names_path $var_names_path \
    --slide_type $slide_type
