#!/usr/bin/env python3
import os
import pandas as pd
import argparse
from os import path
from argparse import ArgumentParser as AP
from os.path import abspath
import time
from pathlib import Path


def get_args():
    # Script description
    description = """Combining all computed spatial features"""

    # Add parser
    parser = AP(description=description,
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output_dir", type=str,
                        help="Path to output folder to store generated files", required=False, default = "")
    parser.add_argument("--prefix", type=str,
                        help="Prefix for output file", default="")
    parser.add_argument("--graph_features", type=str,
                        help="Path to tab-separated file with the graph features", default="")

    parser.add_argument("--clustering_features", type=str,
                        help="Path to tab-separated file with the graph features", default="")
    parser.add_argument("--metadata_path", type=str,
                        help="Path to tab-separated file with metadata", default="")
    parser.add_argument("--is_tcga", type=int,
                        help="dataset is from TCGA (default: True)", default=True, required=False)
    parser.add_argument("--merge_var", type=str,
                        help="Variable to merge metadata and computed features on", default="slide_submitter_id")
    parser.add_argument("--sheet_name", type=str,
                        help="Name of sheet for merging in case a path to xls(x) file is given for metadata_path", default=None)

    parser.add_argument("--version", action="version", version="0.1.0")
    arg = parser.parse_args()
    arg.output_dir = abspath(arg.output_dir)

    if ((arg.output_dir != "") & (not os.path.isdir(arg.output_dir))):
        # Create an empty folder for TF records if folder doesn't exist
        arg.output_dir = Path(arg.output_dir,"process_train")
        os.mkdir(arg.output_dir)
    return arg


def combine_all_spatial_features(graph_features, clustering_features, metadata_path="", is_TCGA=False, merge_var="slide_submitter_id", sheet_name=None):
    """
    Combine network and clustering features into a single file. If metadata_path is not None, add the metadata as well, based on variable slide_submitter_id

    Args:
        output_dir (str): directory containing the graph and clustering features
        slide_type (str): slide type to identify correct files for merging, either "FF" or "FFPE" (default="FF")
        metadata_path (str): path to file containing metadata
        is_TCGA (bool): whether data is from TCGA
        merge_var (str): variable on which to merge (default: slide_submitter_id)

    """
    all_features_graph = pd.read_csv(graph_features, sep=",", index_col=False)
    all_features_clustering = pd.read_csv(clustering_features, sep=",", index_col=False)

    all_features_combined = pd.merge(
        all_features_graph, all_features_clustering)

    # Add additional identifiers for TCGA
    if is_TCGA:
        all_features_combined["TCGA_patient_ID"] = all_features_combined.slide_submitter_id.str[0:12]
        all_features_combined["TCGA_sample_ID"] = all_features_combined.slide_submitter_id.str[0:15]
        all_features_combined["sample_submitter_id"] = all_features_combined.slide_submitter_id.str[0:16]

    # Add metadata if available
    if path.isfile(metadata_path):
        file_extension = metadata_path.split(".")[-1]
        if (file_extension.startswith("xls")):
            if sheet_name is None:
                metadata = pd.read_excel(metadata_path)
        elif (file_extension == "txt") or (file_extension == "csv"):
            metadata = pd.read_csv(metadata_path, sep="\t")
        all_features_combined = pd.merge(
            all_features_combined, metadata, on=merge_var, how="left")
    return all_features_combined

def main(args):
    print("Post-processing: combining all features")
    combine_all_spatial_features(args.graph_features, args.clustering_features, args.metadata_path, args.is_tcga, args.merge_var, args.sheet_name).to_csv(
        Path(args.output_dir, f"{args.prefix}_all_features_combined.csv"), sep="\t", index=False)


if __name__ == "__main__":
    args = get_args()
    st = time.time()
    main(args)
    rt = time.time() - st
    print(f"Script finished in {rt // 60:.0f}m {rt % 60:.0f}s")
