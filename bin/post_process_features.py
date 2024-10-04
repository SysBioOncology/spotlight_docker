#!/usr/bin/env python3
import argparse
import os
import time
from argparse import ArgumentParser as AP
from os.path import abspath
from pathlib import Path

import dask.dataframe as dd


import DL.utils as utils
import pandas as pd


def get_args():
    # Script description
    description = """Post processing features"""

    # Add parser
    parser = AP(
        description=description, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # Sections
    parser.add_argument(
        "--output_dir", help="Set output folder (default='.')", default="."
    )
    parser.add_argument(
        "--create_parquet_subdir",
        help="Whether to create a subdirectory called 'features_format_parquet' if slide_type == 'FFPE', default=False",
        default=False,
    )
    parser.add_argument("--slide_type", help="Type of tissue slide (FF or FFPE)")
    parser.add_argument(
        "--is_tcga", help="Is TCGA dataset, default=False", type=int, default=0
    )
    parser.add_argument("--bot_train_file", type=str, default=None, help="Txt file")
    parser.add_argument("--version", action="version", version="0.1.0")
    arg = parser.parse_args()

    if arg.bot_train_file is None:
        arg.bot_train_file = Path(arg.output_dir, "bot_train.txt")

    if arg.create_parquet_subdir:
        arg.output_dir = abspath(Path(arg.output_dir, "features_format_parquet"))

    if not os.path.isdir(arg.output_dir):
        os.mkdir(arg.output_dir)

    return arg


def handle_ff_slides(bot_train_file, is_tcga):
    features_raw = pd.read_csv(bot_train_file, sep="\t", header=None)
    # Extract the DL features (discard: col1 = tile paths, col2 = true class id)
    features = features_raw.iloc[:, 2:]
    features.columns = list(range(1536))
    # Add new column variables that define each tile
    features["tile_ID"] = [
        utils.get_tile_name(tile_path) for tile_path in features_raw.iloc[:, 0]
    ]
    features["Coord_X"] = [i[-2] for i in features["tile_ID"].str.split("_")]
    features["Coord_Y"] = [i[-1] for i in features["tile_ID"].str.split("_")]
    # FIX add sample_submitter_id and slide_submitter_id depending on is_tcga
    if is_tcga:
        features["sample_submitter_id"] = features["tile_ID"].str[0:16]
        features["slide_submitter_id"] = features["tile_ID"].str[0:23]
        features["Section"] = features["tile_ID"].str[20:23]
    else:
        features["sample_submitter_id"] = features["tile_ID"].str.split("_").str[0]
    return features


def handle_ffpe_slides(bot_train_file, is_tcga):
    features_raw = dd.read_csv(bot_train_file, sep="\t", header=None)
    features_raw["tile_ID"] = features_raw.iloc[:, 0]
    features_raw.tile_ID = features_raw.tile_ID.map(lambda x: x.split("/")[-1])
    features_raw["tile_ID"] = features_raw["tile_ID"].str.replace(".jpg'", "")
    features = features_raw.map_partitions(lambda df: df.drop(columns=[0, 1]))
    new_names = list(map(lambda x: str(x), list(range(1536))))
    new_names.append("tile_ID")
    features.columns = new_names
    # FIX add sample_submitter_id and slide_submitter_id depending on is_tcga
    if is_tcga:
        features["sample_submitter_id"] = features["tile_ID"].str[0:16]
        features["slide_submitter_id"] = features["tile_ID"].str[0:23]
        features["Section"] = features["tile_ID"].str[20:23]
    else:
        features["sample_submitter_id"] = features["tile_ID"].str.split("_").str[0]
    features["Coord_X"] = features["tile_ID"].str.split("_").str[1]
    features["Coord_Y"] = features["tile_ID"].str.split("_").str[-1]
    return features


def post_process_features(bot_train_file, slide_type="FF", is_tcga="TCGA"):
    """
    Format extracted histopathological features from bot.train.txt file generated by myslim/bottleneck_predict.py and extract the 1,536 features, tile names. Extract several variables from tile ID.

    Args:
        bot_train_file (txt)
        slide_type (str)
        is_tcga (bool)

    Returns:
            features (dataframe) contains the 1,536 features, followed by the sample_submitter_id, tile_ID, slide_submitter_id, Section, Coord_X and Coord_Y and in the rows the tiles
    """
    # Read histopathological computed features
    if slide_type == "FF":
        return handle_ff_slides(bot_train_file=bot_train_file, is_tcga=is_tcga)
    elif slide_type == "FFPE":
        return handle_ffpe_slides(bot_train_file=bot_train_file, is_tcga=is_tcga)
    else:
        raise Exception("Invalid `slide_type`, please choose 'FF' or 'FFPE' ")


def main(args):
    features = post_process_features(
        bot_train_file=args.bot_train_file,
        slide_type=args.slide_type,
        is_tcga=args.is_tcga,
    )
    if args.slide_type == "FF":
        #  Save features to .csv file
        features.to_csv(Path(args.output_dir, "features.txt"), sep="\t", header=True)
    elif args.slide_type == "FFPE":
        features.to_parquet(
            path=args.output_dir, compression="gzip", name_function=utils.name_function
        )
    print("Finished post-processing of features...")


if __name__ == "__main__":
    args = get_args()
    st = time.time()
    main(args)
    rt = time.time() - st
    print(f"Script finished in {rt // 60:.0f}m {rt % 60:.0f}s")
