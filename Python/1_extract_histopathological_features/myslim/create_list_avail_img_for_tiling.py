#!/usr/bin/python
import argparse
import DL.image as im
import tiffslide as openslide
import os
import sys

import numpy as np
import pandas as pd
from PIL import Image


def create_list_avail_img_for_tiling(slides_folder, clinical_file_path, output_folder):
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

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Subset images of interest (present in generated clinical file)
    clinical_file = pd.read_csv(clinical_file_path, sep="\t")
    clinical_file.dropna(how="all", inplace=True)
    clinical_file.drop_duplicates(inplace=True)
    clinical_file.drop_duplicates(subset="slide_submitter_id", inplace=True)
    subset_images = clinical_file.image_file_name.tolist()
    print(subset_images)

    # Check if slides are among our data
    available_images = os.listdir(slides_folder)
    images_for_tiling = list(set(subset_images) & set(available_images))

    pd.DataFrame([[name.split(".")[0], name] for name in images_for_tiling], columns=["slide_id", "slide_filename"]).to_csv(
        (f"{output_folder}/avail_slides_for_img.csv"), index=False)

    print("Generated list of available images for tiling...")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--slides_folder", help="Set slides folder")
    parser.add_argument("--output_folder", help="Set output folder")
    parser.add_argument("--clinical_file_path", help="Set clinical file path")
    args = parser.parse_args()
    create_list_avail_img_for_tiling(slides_folder=args.slides_folder,
                                     output_folder=args.output_folder, clinical_file_path=args.clinical_file_path)
