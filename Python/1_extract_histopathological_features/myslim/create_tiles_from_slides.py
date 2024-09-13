#!/usr/bin/python
import tiffslide as openslide
import os
import sys

import numpy as np
import pandas as pd
from PIL import Image

sys.path.append(f"{os.path.dirname(os.getcwd())}/Python/libs")
REPO_DIR = os.path.dirname(os.getcwd())

import DL.image as im


def create_tiles_from_slides(slides_folder, output_folder, clinical_file_path):
    """
    Create tiles from slides
    Dividing the whole slide images into tiles with a size of 512 x 512 pixels, with an overlap of 50 pixels at a magnification of 20x. In addition, remove blurred and non-informative tiles by using the weighted gradient magnitude.

    Source:
    Fu, Y., Jung, A. W., Torne, R. V., Gonzalez, S., Vöhringer, H., Shmatko, A., Yates, L. R., Jimenez-Linan, M., Moore, L., & Gerstung, M. (2020). Pan-cancer computational histopathology reveals mutations, tumor composition and prognosis. Nature Cancer, 1(8), 800–810. https://doi.org/10.1038/s43018-020-0085-8

    Args:
        slides_folder (str): path pointing to folder with all whole slide images (.svs files)
        output_folder (str): path pointing to folder for storing all created files by script (i.e. .jpg files for the created tiles)

    Returns:
        jpg files for the created tiles in the specified folder {output_folder}/tiles

    """

    # Create folder for storing the tiles if non-existent
    tiles_folder = "{}/tiles".format(output_folder)
    if not os.path.exists(tiles_folder):
        os.makedirs(tiles_folder)
        print(tiles_folder)

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

    print(len(images_for_tiling), 'images available:')
    counter = 1
    for slide_filename in images_for_tiling:
        # Accept different file types
        if slide_filename.endswith(('.svs', '.ndpi', '.tif')):
            print(counter, ':', slide_filename)
            slide = openslide.OpenSlide(
                "{}/{}".format(slides_folder, slide_filename))
            slide_name = slide_filename.split(".")[0]
            if (
                str(slide.properties["tiff.ImageDescription"]).find(
                    "AppMag = 40"
                )
                != -1
            ):
                region_size = 1024
                tile_size = 924
            else:
                region_size = 512
                tile_size = 462
            [width, height] = slide.dimensions
            for x_coord in range(1, width, tile_size):
                for y_coord in range(1, height, tile_size):
                    slide_region = slide.read_region(
                        location=(x_coord, y_coord),
                        level=0,
                        size=(region_size, region_size),
                    )
                    slide_region_converted = slide_region.convert("RGB")
                    tile = slide_region_converted.resize(
                        (512, 512), Image.ANTIALIAS)
                    grad = im.getGradientMagnitude(np.array(tile))
                    unique, counts = np.unique(grad, return_counts=True)
                    if counts[np.argwhere(unique <= 20)].sum() < 512 * 512 * 0.6:
                        tile.save(
                            "{}/{}_{}_{}.jpg".format(
                                tiles_folder, slide_name, x_coord, y_coord
                            ),
                            "JPEG",
                            optimize=True,
                            quality=94,
                        )
            counter = counter + 1

    print("Finished creating tiles from the given slides")


if __name__ == "__main__":
    create_tiles_from_slides(sys.argv[1], sys.argv[2])
