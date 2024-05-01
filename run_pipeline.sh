#!/usr/bin/env bash
#SBATCH -J spotlight-docker
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --partition=veryhimem
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=slurm_out/%x_%j.out
#SBATCH --error=slurm_out/%x_%j.out

module load apptainer

# Directory 'spotlight_docker'

work_dir="/path/to/spotlight_docker"
spotlight_sif="path/to/spotlight_sif"

# Define directories/files in container (mounted)

folder_images="/path/to/images_dir"
output_dir="/path/to/output_dir"

# Relative to docker, i.e. start with /data

checkpoint="/data/checkpoint/Retrained_Inception_v4/model.ckpt-100000"
clinical_files_dir="/data/path/to/clinical/TCGA/file.tsv"

# Remaining parameters (this configuration has been tested)
slide_type="FF"
tumor_purity_threshold=80
class_names="SKCM_T"
model_name="inception_v4"

echo "Create output directory: ${output_dir}..."
mkdir -p ${output_dir}

echo "Binding directories..."
# Bind directories + give r/o/w access (do not touch)
# Automatically binds the following
# - Included in repository: /data, /Python 
# - Defined by used: {folder_images},  {output_dir}
export APPTAINER_BINDPATH=${work_dir}/data/:/project/data:ro,${folder_images}:/project/images:ro,${output_dir}:/project/output:rw,${work_dir}/run_scripts:/project/run_scripts:ro,${work_dir}/Python:/project/Python:ro

echo "Run pipeline..."
echo "Extract histopathological features (1 out of 3)"
apptainer exec \
    --cleanenv \
    -c \
    ${spotlight_sif} \
    bash "/project/run_scripts/1_extract_histopatho_features.sh" ${checkpoint} ${clinical_files_dir} ${slide_type} ${class_names} ${tumor_purity_threshold} ${model_name}

echo "Tile level cell type quanitification (2 out of 3)"
apptainer exec \
    --cleanenv \
    -c \
    ${spotlight_sif} \
    bash "/project/run_scripts/2_tile_level_cell_type_quantification.sh" $slide_type

echo "Compute spatial features (3 out of 3)"
apptainer exec \
    --cleanenv \
    -c \
    ${spotlight_sif} \
    bash "/project/run_scripts/3_compute_spatial_features.sh" ${slide_type}

echo "COMPLETED!"