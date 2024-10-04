#!/usr/bin/env python3
import argparse
import glob
import os
import time
from argparse import ArgumentParser as AP
from os.path import abspath
from pathlib import Path

import DL.image as im
import numpy as np
import tiffslide as openslide
from PIL import Image


def get_args():
    # Script description
    description = """Creating tiles from a slide"""

    # Add parser
    parser = AP(
        description=description, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # Sections
    parser = argparse.ArgumentParser()
    parser.add_argument("--filename_slide", help="Name of slide", default="")
    parser.add_argument("--slides_folder", help="Set slides folder", default=None)
    parser.add_argument("--slide_path", help="Path to individual slide", default=None)
    parser.add_argument("--output_dir", help="Set output folder", default="")
    parser.add_argument("--clin_path", help="Set clinical file path", default=None)
    parser.add_argument(
        "--gradient_mag_filter", help="Threshold for filtering", default=20
    )
    parser.add_argument("--version", action="version", version="0.1.0")
    arg = parser.parse_args()
    arg.output_dir = abspath(arg.output_dir)

    if (arg.output_dir != "") & (not os.path.isdir(arg.output_dir)):
        arg.output_dir = Path(arg.output_dir, "tiles")
        os.mkdir(arg.output_dir)
    return arg


def create_tiles_from_slide(
    slide_filename, slides_folder, gradient_mag_filter=20, slide_path=None
):
    """
    Create tiles from a single slide
    Dividing the whole slide images into tiles with a size of 512 x 512 pixels, with an overlap of 50 pixels at a magnification of 20x. In addition, remove blurred and non-informative tiles by using the weighted gradient magnitude.

    Source:
    Fu, Y., Jung, A. W., Torne, R. V., Gonzalez, S., Vöhringer, H., Shmatko, A., Yates, L. R., Jimenez-Linan, M., Moore, L., & Gerstung, M. (2020). Pan-cancer computational histopathology reveals mutations, tumor composition and prognosis. Nature Cancer, 1(8), 800–810. https://doi.org/10.1038/s43018-020-0085-8

    Args:
        slide_filename (str): name of slide to use for creating tiles
        slides_folder (str): path pointing to folder with all whole slide images (.svs files)
        tiles (str): path pointing to folder for storing all created files by script (i.e. .jpg files for the created tiles)
        grad_mag_filter (int): remove tiles that are blurred or non-informative based on weighted gradient magnitude (default=20)

    Returns:
        jpg files for the created tiles in the specified folder {output_dir}/tiles

    """
    # Accept different file types
    if slide_filename.endswith((".svs", ".ndpi", ".tif")):
        if slide_path is not None:
            slide = openslide.OpenSlide(slide_path)
        else:
            slide = openslide.OpenSlide("{}/{}".format(slides_folder, slide_filename))
        slide_name = slide_filename.split(".")[0]
        if str(slide.properties["tiff.ImageDescription"]).find("AppMag = 40") != -1:
            region_size = 1024
            tile_size = 924
        else:
            region_size = 512
            tile_size = 462
        [width, height] = slide.dimensions

        list_of_tiles = []
        for x_coord in range(1, width, tile_size):
            for y_coord in range(1, height, tile_size):
                slide_region = slide.read_region(
                    location=(x_coord, y_coord),
                    level=0,
                    size=(region_size, region_size),
                )
                slide_region_converted = slide_region.convert("RGB")
                tile = slide_region_converted.resize((512, 512), Image.ANTIALIAS)
                grad = im.getGradientMagnitude(np.array(tile))
                unique, counts = np.unique(grad, return_counts=True)
                if (
                    counts[np.argwhere(unique <= int(gradient_mag_filter))].sum()
                    < 512 * 512 * 0.6
                ):
                    list_of_tiles.append((tile, slide_name, x_coord, y_coord))
        return list_of_tiles


def main(args):
    list_of_tiles = create_tiles_from_slide(
        slides_folder=args.slides_folder,
        gradient_mag_filter=args.gradient_mag_filter,
        slide_filename=args.filename_slide,
        slide_path=args.slide_path,
    )
    n_tiles = len(list_of_tiles)
    for tile in list_of_tiles:
        tile[0].save(
            "{}/{}_{}_{}.jpg".format(args.output_dir, tile[1], tile[2], tile[3]),
            "JPEG",
            optimize=True,
            quality=94,
        )
        # Check if all tiles were saved
    assert len(glob.glob1(Path(args.output_dir), "*.jpg")) == n_tiles


if __name__ == "__main__":
    args = get_args()
    st = time.time()
    main(args)
    rt = time.time() - st
    print(f"Script finished in {rt // 60:.0f}m {rt % 60:.0f}s")
