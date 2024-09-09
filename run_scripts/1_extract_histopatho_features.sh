#!/bin/bash

# Adjusted pipeline from PC-CHiP workflow:
# Fu, Y., Jung, A.W., Torne, R.V. et al. Pan-cancer computational histopathology reveals mutations, tumor composition and prognosis. Nat Cancer 1, 800–810 (2020).>

# To execute the pipeline, define: input_dir, output_dir, cancer_type, class_name, checkpoint_path and TCGA_clinical_files.
# cancertype = TCGA_abbreviation_tumor/normal e.g. TCGA_COAD_tumor (see https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations)
# class_name = e.g. COAD_T (see codebook.txt)

# Define type of slide (Fresh-Frozen [FF] vs Formalin-Fixed Paraffin-Embedded [FFPE])

# General setup
repo_dir="/project"

echo "REPO_DIR: ${repo_dir}"
echo "Model checkpoints: $1";
echo "Dir w/ clinical files: $2";

echo "Slide type: $3";
echo "Class names: $4";
echo "Tumor purity threshold: $5";
echo "Model name: $6";

# ---- Relative to /project ---- #
# User input
checkpoint_path=${repo_dir}/$1
clinical_files_dir=${repo_dir}/$2

# Fixed dir
slides_dir=${repo_dir}/images
output_dir=${repo_dir}/output

# Fixed files
path_codebook=${repo_dir}/Python/1_extract_histopathological_features/codebook.txt
path_tissue_classes=${repo_dir}/Python/1_extract_histopathological_features/tissue_classes.csv

# ---- Parameters ---- #
slide_type=$3
class_names=$4
tumor_purity_threshold=$5
model_name=$6
is_tcga=$7 # true or false

# ---------------------------------- #
# ---- create new clinical file ---- #
# ---------------------------------- #

if [ "$is_tcga" = true ]
then
    python $repo_dir/Python/1_extract_histopathological_features/myslim/create_clinical_file.py \
        --class_names $class_names \
        --clinical_files_dir $clinical_files_dir \
        --tumor_purity_threshold $tumor_purity_threshold \
        --output_dir $output_dir/1_histopathological_features \
        --path_codebook ${path_codebook}

    clinical_file=$output_dir/1_histopathological_features/generated_clinical_file.txt
else
    clinical_file=${repo_dir}/data/tmp_clinical_file.txt
    ls $slides_dir | tee ${output_dir}/list_images.txt
    awk -v a=81 -v b="${class_names}" -v c=41 'FNR==NR{print; next}{split($1, tmp, "."); OFS="\t"; print tmp[1], tmp[1], $1, a, b, c}' $clinical_file ${output_dir}/list_images.txt > $output_dir/1_histopathological_features/final_clinical_file.txt
    clinical_file=$output_dir/1_histopathological_features/final_clinical_file.txt
fi

# --------------------------------------------------------- #
# ---- image tiling and image conversion to TF records ---- #
# --------------------------------------------------------- #

python $repo_dir/Python/1_extract_histopathological_features/pre_processing.py \
    --slides_folder $slides_dir \
    --output_folder $output_dir/1_histopathological_features \
    --clinical_file_path $clinical_file

# ------------------------------------------------------ #
# ---- Compute predictions and bottlenecks features ---- #
# ------------------------------------------------------ #

# Compute predictions and bottlenecks features using the Retrained_Inception_v4 checkpoints
python $repo_dir/Python/1_extract_histopathological_features/myslim/bottleneck_predict.py \
    --num_classes=42 \
    --bot_out=$output_dir/1_histopathological_features/bot_train.txt \
    --pred_out=$output_dir/1_histopathological_features/pred_train.txt \
    --model_name $model_name \
    --checkpoint_path $checkpoint_path \
    --file_dir $output_dir/1_histopathological_features/process_train

# ----------------------------------------------------- #
# ---- Post-processing of predictions and futures ----- #
# ----------------------------------------------------- #

# Transform bottleneck features, add dummy variable for tissue type for each tile, save predictions in seperate files
# (= input for pipeline part 2)

python $repo_dir/Python/1_extract_histopathological_features/post_processing.py \
    --output_dir $output_dir/1_histopathological_features \
    --slide_type $slide_type \
    --path_codebook $path_codebook \
    --path_tissue_classes $path_tissue_classes \

# # outputs two files: $output_dir/features $output_dir/predictions
