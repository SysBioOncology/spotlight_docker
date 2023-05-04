import argparse
import os
import git
import pandas as pd
import sys
from myslim.create_file_info_train import format_tile_data_structure
from myslim.create_tiles_from_slides import create_tiles_from_slides
from myslim.datasets.convert import _convert_dataset
REPO_DIR = git.Repo(os.getcwd(), search_parent_directories=True).working_tree_dir
sys.path.append(f"{REPO_DIR}/Python/libs")

def execute_preprocessing(slides_dir, output_dir, clinical_file_path, N_shards=320):
    """
    Execute several pre-processing steps necessary for extracting the histopathological features
    1. Create tiles from slides
    2. Construct file necessary for the deep learning architecture
    3. Convert images of tiles to TF records

    Args:
        slides_dir (str): path pointing to folder with all whole slide images (.svs files)
        output_dir (str): path pointing to folder for storing all created files by script
        clinical_file_path (str): path pointing to formatted clinical file (either generated or manually formatted)
        N_shards (int): default: 320
        checkpoint_path (str): path pointing to checkpoint to be used

    Returns:
        {output_dir}/tiles/{tile files}
        {output_dir}/file_info_train.txt file specifying data structure of the tiles required for inception architecture (to read the TF records)
        {output_dir}/process_train/{TFrecord file} files that store the data as a series of binary sequencies

    """
    full_output_path=f"{output_dir}/1_histopathological_features"
    if not os.path.exists(full_output_path):
        os.makedirs(full_output_path)
    # Create an empty folder for TF records if folder doesn't exist
    process_train_dir = f"{output_dir}/process_train"
    if not os.path.exists(process_train_dir):
        os.makedirs(process_train_dir)

    # Perform image tiling, only kept images of interest
    create_tiles_from_slides(slides_dir=slides_dir, output_dir=full_output_path, clinical_file_path=clinical_file_path)

    # File required for training
    format_tile_data_structure(
        slides_dir=slides_dir,
        output_dir=full_output_path,
        clinical_file_path=clinical_file_path,
    )

    # Convert tiles from jpg to TF record1
    file_info = pd.read_csv(f"{output_dir}/file_info_train.txt", sep="\t")
    training_filenames = list(file_info["tile_path"].values)
    training_classids = [int(id) for id in list(file_info["class_id"].values)]
    tps = [int(id) for id in list(file_info["percent_tumor_cells"].values)]
    Qs = list(file_info["jpeg_quality"].values)

    _convert_dataset(
        split_name="train",
        filenames=training_filenames,
        tps=tps,
        Qs=Qs,
        classids=training_classids,
        output_dir=process_train_dir,
        NUM_SHARDS=N_shards,
    )

    print("Finished converting dataset")
    print(f"The converted data is stored in the directory: {process_train_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--slides_dir", help="Set slides folder")
    parser.add_argument("--output_dir", help="Set output folder")
    parser.add_argument("--clinical_file_path", help="Set clinical file path")
    parser.add_argument("--N_shards", help="Number of shards", default=320)
    args = parser.parse_args()
    execute_preprocessing(
        slides_dir=args.slides_dir,
        output_dir=args.output_dir,
        clinical_file_path=args.clinical_file_path,
        N_shards=args.N_shards,
    )
