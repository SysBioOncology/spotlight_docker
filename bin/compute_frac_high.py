#!/usr/bin/env python3
import os
import pandas as pd
import argparse

# Own modules
import features.features as features


from argparse import ArgumentParser as AP
from os.path import abspath
import time
from pathlib import Path

def get_args():
      # Script description
    description = """Compute Clustering Features"""

    # Add parser
    parser = AP(description=description,
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--slide_indiv_clusters_labeled", type=str,
                        help="Path to csv file", required=True)
    parser.add_argument("--output_dir", type=str,
                        help="Path to output folder to store generated files", required=False, default = "")
    parser.add_argument("--prefix", type=str,
                        help="Prefix for output file", default="")
    parser.add_argument("--version", action="version", version="0.1.0")
    arg = parser.parse_args()
    arg.output_dir = abspath(arg.output_dir)

    if ((arg.output_dir != "") & (not os.path.isdir(arg.output_dir))):
        # Create an empty folder for TF records if folder doesn't exist
        arg.output_dir = Path(arg.output_dir,"process_train")
        os.mkdir(arg.output_dir)

    return arg

def compute_frac_high(slide_indiv_clusters_labeled):
    """
    Compute clustering features:

    Args:
        slide_indiv_clusters_labeled: dataframe with the labeled clusters based on individual cell type SCHC

    Returns:
        frac_high_wide (DataFrame)

    """

    slide_indiv_clusters_labeled = pd.read_csv(slide_indiv_clusters_labeled)

    # Count the fraction of 'high' clusters
    frac_high = features.n_high_clusters(slide_indiv_clusters_labeled)

    frac_high_sub = frac_high[frac_high["is_high"]].copy()
    frac_high_sub = frac_high_sub.drop(
        columns=["is_high", "n_clusters", "n_total_clusters"])

    frac_high_wide = frac_high_sub.pivot(index=["slide_submitter_id"], columns=[
                                         "cell_type_map"])["fraction"]
    new_cols = [('fraction {0} clusters labeled high'.format(col))
                for col in frac_high_wide.columns]
    frac_high_wide.columns = new_cols
    frac_high_wide = frac_high_wide.sort_index(axis="columns").reset_index()
    return frac_high_wide

def main(args):
    compute_frac_high(
          slide_indiv_clusters_labeled = args.slide_indiv_clusters_labeled).to_csv(
        Path(args.output_dir, f"{args.prefix}_frac_high_wide.csv", index = False))

if __name__ == "__main__":
    args = get_args()
    st = time.time()
    main(args)
    rt = time.time() - st
    print(f"Script finished in {rt // 60:.0f}m {rt % 60:.0f}s")
