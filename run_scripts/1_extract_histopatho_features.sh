#!/bin/bash

# Adjusted pipeline from PC-CHiP workflow:
# Fu, Y., Jung, A.W., Torne, R.V. et al. Pan-cancer computational histopathology reveals mutations, tumor composition and prognosis. Nat Cancer 1, 800–810 (2020).>

# To execute the pipeline, define: input_dir, output_dir, cancer_type, class_name, checkpoint_path and TCGA_clinical_files.
# cancertype = TCGA_abbreviation_tumor/normal e.g. TCGA_COAD_tumor (see https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations)
# class_name = e.g. COAD_T (see codebook.txt)

# Define type of slide (Fresh-Frozen [FF] vs Formalin-Fixed Paraffin-Embedded [FFPE])

# General setup
repo_dir=$(pwd)

# command line rguments
echo "Folder images: $1";
echo "Model checkpoints: $2";
echo "Slide type: $3";
echo "Folder output: $4";

slides_dir=$1
checkpoint_path=$2
slide_type=$3
output_dir=$4
clinical_file=$output_dir/1_histopathological_features/tmp_clinical_file.txt

# ---------------------------------- #
# ---- create new clinical file ---- #
# ---------------------------------- #

ls $slides_dir | tee list_images.txt
awk -v a=81 -v b="SKCM_T" -v c=41 'FNR==NR{print; next}{split($1, tmp, "."); OFS="\t"; print tmp[1], tmp[1], $1, a, b, c}' $clinical_file list_images.txt > $output_dir/1_histopathological_features/final_clinical_file.txt

# --------------------------------------------------------- #
# ---- image tiling and image conversion to TF records ---- #
# --------------------------------------------------------- #

python $repo_dir/Python/1_extract_histopathological_features/pre_processing.py \
    --slides_folder=$slides_dir \
    --output_folder=$output_dir/1_histopathological_features \
    --clinical_file_path=$output_dir/1_histopathological_features/final_clinical_file.txt

# ------------------------------------------------------ #
# ---- Compute predictions and bottlenecks features ---- #
# ------------------------------------------------------ #

# Compute predictions and bottlenecks features using the Retrained_Inception_v4 checkpoints
model_name="inception_v4"
python $repo_dir/Python/1_extract_histopathological_features/myslim/bottleneck_predict.py \
    --num_classes=42 \
    --bot_out=$output_dir/1_histopathological_features/bot_train.txt \
    --pred_out=$output_dir/1_histopathological_features/pred_train.txt \
    --model_name=$model_name \
    --checkpoint_path=$checkpoint_path \
    --file_dir=$output_dir/1_histopathological_features/process_train

# ----------------------------------------------------- #
# ---- Post-processing of predictions and futures ----- #
# ----------------------------------------------------- #

# Transform bottleneck features, add dummy variable for tissue type for each tile, save predictions in seperate files
# (= input for pipeline part 2)

python $repo_dir/Python/1_extract_histopathological_features/post_processing.py \
    --output_dir=$output_dir/1_histopathological_features \
    --slide_type=$slide_type

# outputs two files: $output_dir/features $output_dir/predictions
