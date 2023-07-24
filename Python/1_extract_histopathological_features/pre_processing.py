import argparse
import os
import pandas as pd
import sys

from myslim.create_file_info_train import format_tile_data_structure
from myslim.create_tiles_from_slides import create_tiles_from_slides
from myslim.datasets.convert import _convert_dataset

sys.path.append(f"{os.path.dirname(os.getcwd())}/Python/libs")
REPO_DIR = os.path.dirname(os.getcwd())

def execute_preprocessing(slides_folder, output_folder, clinical_file_path, N_shards=320):
    """
    Execute several pre-processing steps necessary for extracting the histopathological features
    1. Create tiles from slides
    2. Construct file necessary for the deep learning architecture
    3. Convert images of tiles to TF records

    Args:
        slides_folder (str): path pointing to folder with all whole slide images (.svs files)
        output_folder (str): path pointing to folder for storing all created files by script
        clinical_file_path (str): path pointing to formatted clinical file (either generated or manually formatted)
        N_shards (int): default: 320
        checkpoint_path (str): path pointing to checkpoint to be used

    Returns:
        {output_folder}/tiles/{tile files}
        {output_folder}/file_info_train.txt file specifying data structure of the tiles required for inception architecture (to read the TF records)
        {output_folder}/process_train/{TFrecord file} files that store the data as a series of binary sequencies

    """

    # Create an empty folder for TF records if folder doesn't exist
    process_train_dir = f"{output_folder}/process_train"
    if not os.path.exists(process_train_dir):
        os.makedirs(process_train_dir)

    # Perform image tiling, only kept images of interest
    create_tiles_from_slides(slides_folder=slides_folder, output_folder=output_folder, clinical_file_path=clinical_file_path)

    # File required for training
    format_tile_data_structure(
        slides_folder=slides_folder,
        output_folder=output_folder,
        clinical_file_path=clinical_file_path
    )

    # Convert tiles from jpg to TF record1
    file_info = pd.read_csv(f"{output_folder}/file_info_train.txt", sep="\t")
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
    parser.add_argument("--slides_folder", help="Set slides folder")
    parser.add_argument("--output_folder", help="Set output folder")
    parser.add_argument("--clinical_file_path", help="Set clinical file path")
    parser.add_argument("--N_shards", help="Number of shards", default=320)
    args = parser.parse_args()
    execute_preprocessing(
        slides_folder=args.slides_folder,
        output_folder=args.output_folder,
        clinical_file_path=args.clinical_file_path,
        N_shards=args.N_shards,
    )
