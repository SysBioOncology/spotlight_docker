import argparse
import os
import os.path
import sys

import pandas as pd
import git
REPO_DIR= git.Repo('.', search_parent_directories=True).working_tree_dir
sys.path.append(f"{REPO_DIR}/Python/libs")

# trunk-ignore(flake8/E402)
import DL.utils as utils

# trunk-ignore(flake8/E402)
from openslide import OpenSlide


def format_tile_data_structure(slides_folder, output_dir, clinical_file_path):
    """
    Specifying the tile data structure required to store tiles as TFRecord files (used in convert.py)

    Args:
        slides_folder (str): path pointing to folder with all whole slide images (.svs files)
        output_dir (str): path pointing to folder for storing all created files by script
        clinical_file_path (str): path pointing to formatted clinical file (either generated or manually formatted)

    Returns:
        {output_dir}/file_info_train.txt containing the path to the individual tiles, class name, class id, percent of tumor cells and JPEG quality

    """
    tiles_folder = output_dir + "/tiles"

    clinical_file = pd.read_csv(clinical_file_path, sep="\t")
    clinical_file.dropna(how="all", inplace=True)
    clinical_file.drop_duplicates(inplace=True)
    clinical_file.drop_duplicates(subset="slide_submitter_id", inplace=True)

    # 2) Determine the paths paths of jpg tiles
    all_tile_names = os.listdir(tiles_folder)
    jpg_tile_names = []
    jpg_tile_paths = []

    for tile_name in all_tile_names:
        if "jpg" in tile_name:
            jpg_tile_names.append(tile_name)
            jpg_tile_paths.append(tiles_folder + "/" + tile_name)

    # 3) Get corresponding data from the clinical file based on the tile names
    jpg_tile_names_stripped = [
        utils.get_slide_submitter_id(jpg_tile_name) for jpg_tile_name in jpg_tile_names
    ]
    jpg_tile_names_df = pd.DataFrame(
        jpg_tile_names_stripped, columns=["slide_submitter_id"]
    )
    jpg_tiles_df = pd.merge(
        jpg_tile_names_df, clinical_file, on=["slide_submitter_id"], how="left"
    )

    # 4) Determine jpeg_quality of slides
    slide_quality = []
    for slide_name in jpg_tiles_df.image_file_name.unique():
        img = OpenSlide("{}/{}".format(slides_folder, slide_name))
        image_description = img.properties.values.__self__.get("tiff.ImageDescription").split("|")[0]
        image_description_split = image_description.split(" ")
        jpeg_quality = image_description_split[-1]
        slide_quality.append([slide_name, "RGB" + jpeg_quality])

    slide_quality_df = pd.DataFrame(
        slide_quality, columns=["image_file_name", "jpeg_quality"]
    )
    jpg_tiles_df = pd.merge(
        jpg_tiles_df, slide_quality_df, on=["image_file_name"], how="left"
    )
    jpg_tiles_df["tile_path"] = jpg_tile_paths

    # Create output dataframe
    output = jpg_tiles_df[
        ["tile_path", "class_name", "class_id", "jpeg_quality", "percent_tumor_cells"]
    ]
    output.to_csv(output_dir + "/file_info_train.txt", index=False, sep="\t")

    print("Finished creating the necessary file for training in the next step")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--slides_folder", help="Set slides folder")
    parser.add_argument("--output_dir", help="Set output folder")
    parser.add_argument("--clin_path", help="Set clinical file path")
    args = parser.parse_args()

    format_tile_data_structure(
        slides_folder=args.slides_folder,
        output_dir=args.output_dir,
        clinical_file_path=args.clin_path,
    )
