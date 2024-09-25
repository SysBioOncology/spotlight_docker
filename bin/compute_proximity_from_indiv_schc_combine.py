#!/usr/bin/env python3
import os
import pandas as pd
import argparse
import multiprocessing
from argparse import ArgumentParser as AP
from os.path import abspath
import time
from pathlib import Path

def get_args():
      # Script description
    description = """Compute Clustering Features: Proximity"""

    # Add parser
    parser = AP(description=description,
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--prox_within", type=str,
                        help="Path to csv file", required=True)
    parser.add_argument("--prox_between", type=str,
                        help="Path to csv file", required=True)
    parser.add_argument("--prefix", type=str,
                        help="Prefix for output file", default="")
    parser.add_argument("--output_dir", type=str,
                        help="Path to output folder to store generated files",
                        required=False, default = "")
    parser.add_argument("--version", action="version", version="0.1.0")
    arg = parser.parse_args()
    arg.output_dir = abspath(arg.output_dir)

    if ((arg.output_dir != "") & (not os.path.isdir(arg.output_dir))):
        # Create an empty folder for TF records if folder doesn't exist
        arg.output_dir = Path(arg.output_dir,"process_train")
        os.mkdir(arg.output_dir)

    return arg


def compute_proximity_from_indiv_schc_combine(prox_within, prox_between):

    results_schc_indiv_within = pd.read_csv(prox_within, header = 0, index_col=False)
    results_schc_indiv_between = pd.read_csv(prox_between, header = 0, index_col = False)

    # Concatenate within and between computed proximity values
    prox_indiv_schc_combined = pd.concat(
        [results_schc_indiv_within, results_schc_indiv_between])

    # Remove rows with a proximity of NaN
    prox_indiv_schc_combined = prox_indiv_schc_combined.dropna(axis=0)

    prox_indiv_schc_combined.comparison = prox_indiv_schc_combined.comparison.replace(dict(zip(['cluster1=True-cluster2=True', 'cluster1=True-cluster2=False',
                                                                                                'cluster1=False-cluster2=True', 'cluster1=False-cluster2=False'], ["high-high", "high-low", "low-high", "low-low"])))
    prox_indiv_schc_combined["pair (comparison)"] = [
        f"{pair} ({comp})" for pair, comp in prox_indiv_schc_combined[["pair", "comparison"]].to_numpy()]
    prox_indiv_schc_combined = prox_indiv_schc_combined.drop(
        axis=1, labels=["pair", "comparison"])


    prox_indiv_schc_combined_wide = prox_indiv_schc_combined.pivot(
        index=["slide_submitter_id"], columns=["pair (comparison)"])["proximity"]
    new_cols = [
        f'prox CC {col.replace("_", " ")}' for col in prox_indiv_schc_combined_wide.columns]
    prox_indiv_schc_combined_wide.columns = new_cols
    prox_indiv_schc_combined_wide = prox_indiv_schc_combined_wide.reset_index()
    return prox_indiv_schc_combined_wide

def main(args):
    compute_proximity_from_indiv_schc_combine(prox_within=args.prox_within, prox_between= args.prox_between).to_csv(
        Path(args.output_dir,
              f"{args.prefix}_features_clust_indiv_schc_prox.csv"), sep="\t",
        index=False)

if __name__ == "__main__":
    args = get_args()
    st = time.time()
    main(args)
    rt = time.time() - st
    print(f"Script finished in {rt // 60:.0f}m {rt % 60:.0f}s")
