#!/bin/bash

#####################################################################
## Compute cell-type quantification from transfer learning models  ##
#####################################################################


# ----------------------------------- #
# --------- Setup file paths -------- #
# ----------------------------------- #

# General setup
repo_dir=$(pwd)

# command line rguments
echo "Features folder: $1"
echo "Slide type: $2";
echo "Folder output: $3";

# Define type of slide
features_dir=$1

# Define type of slide
slide_type=$2

# Define output directory
output_dir=$3

# Transfer Learning trained models directory
models_dir=$repo_dir/data/TF_models/SKCM_${slide_type}

# Histopathological features
if [[ $slide_type == "FF" ]];then
    histopatho_features_path=$output_dir/features.txt
elif [[ $slide_type == "FFPE" ]];then
    histopatho_features_path=$output_dir/features_format_parquet/
fi

# Task signatures
var_names_path=$repo_dir/Python/2_train_multitask_models/task_selection_names.pkl

# Compute predictions using models learned from unseen folds
prediction_mode="test" # (tcga_validation, tcga_train_validation)

# ---------------------------------------------------- #
# ---- Predict cell type abundances on tile level ---- #
# ---------------------------------------------------- #

# For now, we use models trained on FF slides

python $repo_dir/Python/2_train_multitask_models/tile_level_cell_type_quantification.py \
    --models_dir=$models_dir \
    --output_dir=$output_dir \
    --histopatho_features_path=$histopatho_features_path \
    --prediction_mode=$prediction_mode \
    --var_names_path=$var_names_path \
    --slide_type=$slide_type
