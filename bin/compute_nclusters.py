#!/usr/bin/env python3
import os
import sys
import joblib
import pandas as pd
from joblib import Parallel, delayed
import argparse
from os import path

# Own modules
import features.clustering as clustering
import features.features as features
import features.graphs as graphs
import features.utils as utils
from model.constants import DEFAULT_SLIDE_TYPE, DEFAULT_CELL_TYPES, METADATA_COLS


import multiprocessing

from argparse import ArgumentParser as AP
from os.path import abspath
import time
from pathlib import Path



def get_args():
      # Script description
    description = """Compute Spatial Network Features: Compute Connectedness"""

    # Add parser
    parser = AP(description=description,
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--all_slide_clusters_characterized", type=str,
                        help="Path to csv file", required=True)
    parser.add_argument("--output_dir", type=str,
                        help="Path to output folder to store generated files", required=False, default = "")
    parser.add_argument("--slide_type", type=str,
                        help="Type of slides 'FFPE' or 'FF' used for naming generated files (default: 'FF')", default="FF")
    parser.add_argument("--cell_types_path", type=str,
                        help="Path to file with list of cell types (default: CAFs, endothelial_cells, T_cells, tumor_purity)", default=None)
    parser.add_argument("--graphs_path", type=str,
                        help="Path to pkl with generated graphs in case this was done before (OPTIONAL) if not specified, graphs will be generated", default=None)
    parser.add_argument("--prefix", type=str,
                        help="Prefix for output file", default="")
    parser.add_argument("--cutoff_path_length", type=int,
                        help="Max path length for proximity based on graphs", default=2, required=False)
    parser.add_argument("--shapiro_alpha", type=float,
                        help="Choose significance level alpha (default: 0.05)", default=0.05, required=False)
    parser.add_argument("--abundance_threshold", type=float,
                        help="Threshold for assigning cell types based on predicted probability (default: 0.5)", default=0.5, required=False)

    parser.add_argument("--n_clusters", type=int,
                        help="Number of clusters for SCHC (default: 8)", required=False, default=8)
    parser.add_argument("--n_cores", type = int, help = "Number of cores to use (parallelization)")

    parser.add_argument("--version", action="version", version="0.1.0")
    arg = parser.parse_args()
    arg.output_dir = abspath(arg.output_dir)

    if ((arg.output_dir != "") & (not os.path.isdir(arg.output_dir))):
        # Create an empty folder for TF records if folder doesn't exist
        arg.output_dir = Path(arg.output_dir,"process_train")
        os.mkdir(arg.output_dir)

        if arg.n_cores is None:
            arg.n_cores = multiprocessing.cpu_count()
    return arg



def compute_nclusters(all_slide_clusters_characterized, cell_types = DEFAULT_CELL_TYPES):

    if cell_types is None:
        cell_types = DEFAULT_CELL_TYPES

    # # Assign a cell type label based on the mean of all cluster means across all slides
    # all_slide_clusters_characterized = clustering.characterize_clusters(
    #     tiles_all_schc)
    all_slide_clusters_characterized = pd.read_csv(all_slide_clusters_characterized,
    sep = ",", header = 0)


    # Count the number of clusters per cell type for each slide
    num_clust_by_slide = features.n_clusters_per_cell_type(
        all_slide_clusters_characterized, cell_types=cell_types)

    num_clust_by_slide_sub = num_clust_by_slide.copy()
    num_clust_by_slide_sub = num_clust_by_slide_sub.drop(
        columns=["is_assigned", "n_clusters"])

    num_clust_slide_wide = num_clust_by_slide_sub.pivot(
        index=["slide_submitter_id"], columns=["cell_type"])["fraction"]
    new_cols = [('fraction {0} clusters'.format(col))
                for col in num_clust_slide_wide.columns]
    num_clust_slide_wide.columns = new_cols
    num_clust_slide_wide = num_clust_slide_wide.sort_index(
        axis="columns").reset_index()
    return num_clust_slide_wide


def main(args):
    num_clust_slide_wide = compute_nclusters(
          all_slide_clusters_characterized = args.all_slide_clusters_characterized,
          cell_types=args.cell_types_path)

    num_clust_slide_wide.to_csv(
        Path(args.output_dir, f"{args.prefix}_nclusters_wide.csv", sep = "\t", index = False))

if __name__ == "__main__":
    args = get_args()
    st = time.time()
    main(args)
    rt = time.time() - st
    print(f"Script finished in {rt // 60:.0f}m {rt % 60:.0f}s")
