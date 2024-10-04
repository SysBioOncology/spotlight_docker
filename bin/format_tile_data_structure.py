#!/usr/bin/env python3

import argparse
import glob
import os
import os.path
import time
from argparse import ArgumentParser as AP
from os.path import abspath
from pathlib import Path

import DL.utils as utils
import pandas as pd
import tiffslide as openslide


def get_args():
    # Script description
    description = """Format tile data structure"""

    # Add parser
    parser = AP(
        description=description, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # Sections
    parser.add_argument("--slides_folder", help="Set slides folder", default="")
    parser.add_argument("--tiles_folder", help="Directory with the tiles", default="")
    parser.add_argument("--output_dir", help="Set output folder", default="")
    parser.add_argument("--clin_path", help="Set clinical file path", default="")
    parser.add_argument("--is_tcga", help="Set output folder", type=int)

    parser.add_argument("--version", action="version", version="0.1.0")
    arg = parser.parse_args()
    arg.output_dir = abspath(arg.output_dir)

    if (arg.output_dir != "") & (not os.path.isdir(arg.output_dir)):
        os.mkdir(arg.output_dir)
    return arg


def assess_tcga_slide_quality(slide_name, slides_folder):
    print("{}/{}".format(slides_folder, slide_name))

    img = openslide.OpenSlide("{}/{}".format(slides_folder, slide_name))
    image_description = str(img.properties["tiff.ImageDescription"]).split("|")[0]
    image_description_split = image_description.split(" ")
    jpeg_quality = image_description_split[-1]
    return [slide_name, "RGB" + jpeg_quality]


def format_tile_data_structure(
    slides_folder, tiles_folder, output_dir, clinical_file_path, is_tcga=True
):
    """
    Specifying the tile data structure required to store tiles as TFRecord files (used in convert.py)

    Args:
        slides_folder (str): path pointing to folder with all whole slide images (.svs files)
        output_dir (str): path pointing to folder for storing all created files by script
        clinical_file_path (str): path pointing to formatted clinical file (either generated or manually formatted)
        is_tcga (bool): default = TRUE

    Returns:
        {output_dir}/file_info_train.txt containing the path to the individual tiles, class name, class id, percent of tumor cells and JPEG quality

    """
    if tiles_folder == "":
        tiles_folder = output_dir

    clinical_file = pd.read_csv(clinical_file_path, sep="\t")
    clinical_file.dropna(how="all", inplace=True)
    clinical_file.drop_duplicates(inplace=True)
    clinical_file.drop_duplicates(subset="slide_submitter_id", inplace=True)

    # 2) Determine the paths paths of jpg tiles
    jpg_tile_names = glob.glob1(Path(args.tiles_folder), "*.jpg")
    jpg_tile_paths = [Path(tiles_folder, tile_name) for tile_name in jpg_tile_names]

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
    if is_tcga:
        for slide_name in jpg_tiles_df.image_file_name.unique():
            slide_quality.append(
                assess_tcga_slide_quality(
                    slide_name=slide_name, slides_folder=slides_folder
                )
            )
    else:
        jpeg_quality = 100  # assuming no loss
        slide_quality = [
            [slide_name, f"RGB{jpeg_quality}"]
            for slide_name in jpg_tiles_df.image_file_name.unique()
        ]

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
    return output


def main(args):
    output = format_tile_data_structure(
        slides_folder=args.slides_folder,
        tiles_folder=args.tiles_folder,
        output_dir=args.output_dir,
        clinical_file_path=args.clin_path,
        is_tcga=args.is_tcga,
    )
    output.to_csv(Path(args.output_dir, "file_info_train.txt"), index=False, sep="\t")

    print(
        "Finished creating the necessary file for computing the features in the next step"
    )


if __name__ == "__main__":
    args = get_args()
    st = time.time()
    main(args)
    rt = time.time() - st
    print(f"Script finished in {rt // 60:.0f}m {rt % 60:.0f}s")
