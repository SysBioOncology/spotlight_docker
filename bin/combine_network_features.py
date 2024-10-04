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
    description = """Combine all network Features"""

    # Add parser
    parser = AP(
        description=description, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--all_largest_cc_sizes_wide", type=str, help="Path to csv file ", required=True
    )
    parser.add_argument(
        "--shortest_paths_wide", type=str, help="Path to csv file ", required=True
    )
    parser.add_argument(
        "--colocalization_wide", type=str, help="Path to csv file ", required=True
    )

    parser.add_argument(
        "--output_dir",
        type=str,
        help="Path to output folder to store generated files",
        required=False,
        default="",
    )

    parser.add_argument("--prefix", type=str, help="Prefix for output file", default="")

    parser.add_argument("--version", action="version", version="0.1.0")
    arg = parser.parse_args()
    arg.output_dir = abspath(arg.output_dir)

    if (arg.output_dir != "") & (not os.path.isdir(arg.output_dir)):
        # Create an empty folder for TF records if folder doesn't exist
        arg.output_dir = Path(arg.output_dir, "process_train")
        os.mkdir(arg.output_dir)

    return arg


def combine_network_features(
    all_largest_cc_sizes_wide, shortest_paths_wide, colocalization_wide
):
    """
    Compute network features
    1. effect sizes based on difference in node degree between simulated slides and actual slide
    2. fraction largest connected component
    3. number of shortest paths with a max length.

    Args:
        tile_quantification_path (str)
        output_dir (str)
        slide_type (str): type of slide either 'FF' or 'FFPE'
        cell_types (list): list of cell types
        graphs_path (str): path to pkl file with generated graphs [optional]
        abundance_threshold (float): threshold for assigning cell types to tiles based on the predicted probability (default=0.5)
        shapiro_alpha (float): significance level for shapiro tests for normality (default=0.05)
        cutoff_path_length (int): max. length of shortest paths (default=2)

    Returns:
        all_effect_sizes (DataFrame): dataframe containing the slide_submitter_id, center, neighbor, effect_size (Cohen's d), Tstat, pval, and the pair (string of center and neighbor)
        all_sims_nd (DataFrame): dataframe containing slide_submitter_id, center, neighbor, simulation_nr and degree (node degree)
        all_mean_nd_df (DataFrame): dataframe containing slide_submitter_id, center, neighbor, mean_sim (mean node degree across the N simulations), mean_obs
        all_largest_cc_sizes (DataFrame): dataframe containing slide_submitter_id, cell type and type_spec_frac (fraction of LCC w.r.t. all tiles for cell type)
        shortest_paths_slide (DataFrame): dataframe containing slide_submitter_id, source, target, pair and n_paths (number of shortest paths for a pair)
        all_dual_nodes_frac (DataFrame): dataframe containing slide_submitter_id, pair, counts (absolute) and frac

    """
    all_largest_cc_sizes_wide = pd.read_csv(
        all_largest_cc_sizes_wide, index_col=False, sep="\t"
    )
    shortest_paths_wide = pd.read_csv(shortest_paths_wide, index_col=False, sep="\t")
    colocalization_wide = pd.read_csv(colocalization_wide, index_col=False, sep="\t")

    all_features = pd.merge(all_largest_cc_sizes_wide, shortest_paths_wide)
    all_features = pd.merge(all_features, colocalization_wide)
    return all_features


def main(args):
    combine_network_features(
        args.all_largest_cc_sizes_wide,
        args.shortest_paths_wide,
        args.colocalization_wide,
    ).to_csv(
        Path(args.output_dir, f"{args.prefix}_all_graph_features.csv"), index=False
    )


if __name__ == "__main__":
    args = get_args()
    st = time.time()
    main(args)
    rt = time.time() - st
    print(f"Script finished in {rt // 60:.0f}m {rt % 60:.0f}s")
