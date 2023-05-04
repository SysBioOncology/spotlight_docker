#!/bin/bash

#####################################################################
## Compute cell-type quantification from transfer learning models  ##
#####################################################################

# Define path to spotlight repository
#repo_path="/Users/joankant/Library/CloudStorage/OneDrive-TUEindhoven/spotlight"
repo_path="/home/olapuent/Desktop/spatial_localization/repo_manuscript/spotlight"
run_scripts_dir=$repo_path/run_scripts

# ----------------------------------- #
# ---- Setup virtual environment ---- #
# ----------------------------------- #
# python3 -m pip install virtualenv (virtualenv needs to be installed)
python3 -m virtualenv spotlight_env

# for bash
$run_scripts_dir/spotlight_env/bin/activate
#source $run_scripts_dir/spotlight_env/bin/activate
$run_scripts_dir/spotlight_env/bin/pip install -r $repo_path/requirements.txt

# ----------------------------------- #
# --------- Setup file paths -------- #
# ----------------------------------- #
models_dir=$repo_path/analysis/SKCM_FF/output/model_runs
output_dir=$repo_path/tests/output/SKCM_FF

# Molecular Immune Phenotype (Bagaev et al.)
MFP_path=$repo_path/data/MFP.xlsx

# Define type of slide
slide_type="FFPE"

# Path to trained models
#models_dir=$repo_path/analysis/SKCM_${slide_type}/output/model_runs
models_dir=$repo_path/analysis/SKCM_FF/output/model_runs

# Path to save output
#output_dir=$repo_path/tests/output/SKCM_${slide_type}
output_dir=$repo_path/tests/output/SKCM_FF_FFPE


# Histopathological features
if [[ $slide_type == "FF" ]];then
    histopatho_features_path=$repo_path/analysis/SKCM_${slide_type}/output/features.txt
elif [[ $slide_type == "FFPE" ]];then
    histopatho_features_path=$repo_path/analysis/SKCM_${slide_type}/output/features_format_parquet/
fi

# Compute predictions using models learned from unseen folds
#prediction_mode="tcga_validation"
prediction_mode="tcga_train_validation"

# Task signatures
var_names_path=$repo_path/Python/2_train_multitask_models/task_selection_names.pkl

# ---------------------------------------------------- #
# ---- Predict cell type abundances on tile level ---- #
# ---------------------------------------------------- #
echo $slide_type

python $repo_path/Python/2_train_multitask_models/tile_level_cell_type_quantification.py \
    --models_dir=$models_dir \
    --output_dir=$output_dir \
    --MFP_path=$MFP_path \
    --histopatho_features_path=$histopatho_features_path \
    --prediction_mode=$prediction_mode \
    --var_names_path=$var_names_path \
    --slide_type=$slide_type
