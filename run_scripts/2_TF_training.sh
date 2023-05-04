#!/bin/bash

# Args:
# cur_dir: current directory
# target_features_path: path pointing to file containing quantification of the cell types
# histopatho_features_path: path pointing to the file containing the histopathological features
# output_folder: path pointing ot a folder where the trained models can be stored (within the script for each cell type a subfolder is automatically created)
# alpha_min: min. value for hyperparameter alpha
# alpha_max: max. value for hyperparameter alpha
# n_steps: number of steps for hyperparameter grid

# Load virtual environment
venv_dir="/home/olapuent/Desktop/spatial_localization/repo_manuscript/spotlight/spotlight_venv"
source $venv_dir/bin/activate

# General setup
repo_dir=/Users/joankant/Library/CloudStorage/OneDrive-TUEindhoven/spotlight
output_dir=$repo_dir/output
slide_type=FFPE # FF



# ---------------------------------------------------- #
# ---- Cell type quantification with immunedeconv ---- #
# ---------------------------------------------------- #
tpm_path=""
# Computing cell type quantification with immunedeconv: EPIC, quanTIseq, xCell, MPC Counter
Rscript $repo_dir/R/immunedeconv.R \
    --tpm_path=$tpm_path \
    --output_dir=$output_dir

# ------------------------------------------------------------ #
# ---- Process and create file with tasks for TF learning ---- #
# ------------------------------------------------------------ #
cancer_type="SKCM"
clinical_file_path="/home/olapuent/Desktop/spatial_localization/repo_manuscript/spotlight/data/clinical/SKCM/clinical_data_SKCM.csv"
python $repo_dir/Python/2_train_multitask_models/processing_transcriptomics.py \
    --cancer_type=$cancer_type \
    --slide_type=$slide_type \
    --tpm_path=$tpm_path \
    --clinical_file_path=$clinical_file_path \
    --output_dir=$output_dir

# --------------------------------------------- #
# ---- TF learning: Multi-task Lasso Model ---- #
# --------------------------------------------- #
alpha_min=-4
alpha_max=-1
n_steps=40
cell_types_path="/Users/joankant/Library/CloudStorage/OneDrive-TUEindhoven/spotlight/data/params/cell_types.csv"

cat $cell_types_path | while read celltype; do
	echo $celltype
	python $repo_dir/Python/2_train_multitask_models/run_TF_pipeline.py \
	--output_dir=$output_dir \
	--category=$celltype \
	--alpha_min=$alpha_min \
	--alpha_max=$alpha_max \
	--n_steps=$n_steps \
	--slide_type=$slide_type
done