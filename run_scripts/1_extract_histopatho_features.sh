#!/bin/bash

# Adjusted pipeline from PC-CHiP workflow:
# Fu, Y., Jung, A.W., Torne, R.V. et al. Pan-cancer computational histopathology reveals mutations, tumor composition and prognosis. Nat Cancer 1, 800–810 (2020).>

# To execute the pipeline, define: input_dir, output_dir, cancer_type, class_name, checkpoint_path and TCGA_clinical_files.
# cancertype = TCGA_abbreviation_tumor/normal e.g. TCGA_COAD_tumor (see https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations)
# class_name = e.g. COAD_T (see codebook.txt)

# Define type of slide (Fresh-Frozen [FF] vs Formalin-Fixed Paraffin-Embedded [FFPE])

# General setup
repo_dir=$(pwd)

# ommand line rguments
echo "Folder images: $1";
echo "Folder clinical file: $2";
echo "Model checkpoints: $3";
echo "Slide type: $4";
echo "Folder output: $5";

slides_dir=$1
clinical_files_dir=$2                                                                           
checkpoint_path=$3
slide_type=$4
output_dir=$5

# ---------------------------------- #
# ---- create new clinical file ---- #
# ---------------------------------- #

########################
## For TCGA datasets: ##
########################

# Multi-class: e.g. LUAD and LUSC
# - clinical_files_folder: Directory expected to be in this format: {clinical_files_dir}/clinical_file_TCGA_{cancer_type}.tsv
# - class_names: txt or csv file with on each line a class name (TCGA abbreviation)

# Single-class: e.g. SKCM
# - clinical_files_folder: Give the full path to the clinical file
# - class_names: Give single the class name (TCGA abbreviation)

tumor_purity_threshold=80
class_name=SKCM_T

# Create a filtered clinical file based on thresholded tumor purity
#python $repo_dir/Python/1_extract_histopathological_features/myslim/create_clinical_file.py \
#   --class_names $class_name \
#   --clinical_files_dir $clinical_files_dir \
#   --output_dir $output_dir \
#   --tumor_purity_threshold=$tumor_purity_threshold # (OPTIONAL: by default 80)

# Formatted clinical file available at:
## FF
#clin_file_path="/home/olapuent/Desktop/spatial_localization/repo_manuscript/spotlight/analysis/SKCM_FF/output/generated_clinical_file.txt"
## FFPE (made based on FF generated clinical file, so all slides passed the filtering)
# clin_file_path="/home/olapuent/Desktop/spatial_localization/repo_manuscript/spotlight/analysis/SKCM_FFPE/output/FFPE_generated_clinical_file.txt"

########################
## Any other datasets: ##
########################

# TODO: add code for other datasets

# --------------------------------------------------------- #
# ---- image tiling and image conversion to TF records ---- #
# --------------------------------------------------------- #

python $repo_dir/Python/1_extract_histopathological_features/pre_processing.py \
    --slides_folder=$slides_dir \
    --output_folder=$output_dir \
    --clinical_file_path=$clinical_files_dir/generated_clinical_file.txt

# ------------------------------------------------------ #
# ---- Compute predictions and bottlenecks features ---- #
# ------------------------------------------------------ #

# Compute predictions and bottlenecks features using the Retrained_Inception_v4 checkpoints
model_name="inception_v4"
python $repo_dir/Python/1_extract_histopathological_features/myslim/bottleneck_predict.py \
    --num_classes=42 \
    --bot_out=$output_dir/bot.train.txt \
    --pred_out=$output_dir/pred.train.txt \
    --model_name=$model_name \
    --checkpoint_path=$checkpoint_path \

# ----------------------------------------------------- #
# ---- Post-processing of predictions and futures ----- #
# ----------------------------------------------------- #

# Transform bottleneck features, add dummy variable for tissue type for each tile, save predictions in seperate files
# (= input for pipeline part 2)

## FFPE
python $repo_dir/Python/1_extract_histopathological_features/post_processing.py \
    --output_dir=$output_dir \
    --slide_type=$slide_type

# outputs two files: $output_dir/features $output_dir/predictions
