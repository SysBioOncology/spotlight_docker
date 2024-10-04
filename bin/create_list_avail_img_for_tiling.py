#!/usr/bin/env python3
import argparse
import os
import time
from argparse import ArgumentParser as AP
from os.path import abspath
from pathlib import Path

import pandas as pd


def get_args():
    # Script description
    description = """Creating list of available slide images that have to be tiled"""

    # Add parser
    parser = AP(
        description=description, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser = argparse.ArgumentParser()
    parser.add_argument("--slides_folder", help="Set slides folder", default=None)
    parser.add_argument("--output_dir", help="Set output folder", default="")
    parser.add_argument("--clinical_file_path", help="Set clinical file path")
    parser.add_argument("--version", action="version", version="0.1.0")
    arg = parser.parse_args()
    arg.output_dir = abspath(arg.output_dir)
    arg.slides_folder = abspath(arg.slides_folder)

    if (arg.output_dir != "") & (not os.path.isdir(arg.output_dir)):
        os.mkdir(arg.output_dir)
    return arg


def create_list_avail_img_for_tiling(slides_folder, clinical_file_path):
    """
    Create tiles from slides
    Dividing the whole slide images into tiles with a size of 512 x 512 pixels, with an overlap of 50 pixels at a magnification of 20x. In addition, remove blurred and non-informative tiles by using the weighted gradient magnitude.

    Source:
    Fu, Y., Jung, A. W., Torne, R. V., Gonzalez, S., Vöhringer, H., Shmatko, A., Yates, L. R., Jimenez-Linan, M., Moore, L., & Gerstung, M. (2020). Pan-cancer computational histopathology reveals mutations, tumor composition and prognosis. Nature Cancer, 1(8), 800–810. https://doi.org/10.1038/s43018-020-0085-8

    Args:
        slides_folder (str): path pointing to folder with all whole slide images (.svs files)
        clinical_file_path (str): path pointing to file with clinical file

    Returns:
        txt containing list of slides available for tiling
    """

    # Subset images of interest (present in generated clinical file)
    clinical_file = pd.read_csv(clinical_file_path, sep="\t", index_col=False)
    print(clinical_file)
    clinical_file.dropna(how="all", inplace=True)
    clinical_file.drop_duplicates(inplace=True)
    clinical_file.drop_duplicates(subset="slide_submitter_id", inplace=True)
    subset_images = clinical_file.image_file_name.tolist()
    print(subset_images)

    # Check if slides are among our data
    available_images = os.listdir(slides_folder)
    print(available_images)
    images_for_tiling = list(set(subset_images) & set(available_images))

    return pd.DataFrame(
        [[name.split(".")[0], name] for name in images_for_tiling],
        columns=["slide_id", "slide_filename"],
    )


def main(args):
    list_avail_img = create_list_avail_img_for_tiling(
        slides_folder=args.slides_folder, clinical_file_path=args.clinical_file_path
    )

    list_avail_img.to_csv(
        Path(args.output_dir, "avail_slides_for_img.csv"), index=False
    )
    print("Generated list of available images for tiling...")


if __name__ == "__main__":
    args = get_args()
    st = time.time()
    main(args)
    rt = time.time() - st
    print(f"Script finished in {rt // 60:.0f}m {rt % 60:.0f}s")
