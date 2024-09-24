#!/usr/bin/env python3
import os
import pandas as pd
import argparse
from argparse import ArgumentParser as AP
from os.path import abspath
import time
from pathlib import Path

def get_args():
      # Script description
    description = """Combine all clustering features"""

    # Add parser
    parser = AP(description=description,
                formatter_class=argparse.RawDescriptionHelpFormatter)


    parser.add_argument("--output_dir", type=str,
                        help="Path to output folder to store generated files", required=False, default = "")
    parser.add_argument("--prefix", type=str,
                        help="Prefix for output file", default="")

    parser.add_argument("--frac_high_wide", type=str,
                        help="Path to csv", default="")

    parser.add_argument("--num_clust_slide_wide", type=str,
                        help="Path to csv", default="")

    parser.add_argument("--all_prox_df_wide", type=str,
                        help="Path to csv", default="")
    parser.add_argument("--prox_indiv_schc_combined_wide", type=str,
                        help="Path to csv", default="")
    parser.add_argument("--version", action="version", version="0.1.0")
    arg = parser.parse_args()
    arg.output_dir = abspath(arg.output_dir)

    if ((arg.output_dir != "") & (not os.path.isdir(arg.output_dir))):
        # Create an empty folder for TF records if folder doesn't exist
        arg.output_dir = Path(arg.output_dir,"process_train")
        os.mkdir(arg.output_dir)
    return arg

def combine_clustering_features(frac_high_wide,
                                 num_clust_slide_wide,
                                all_prox_df_wide,
                                prox_indiv_schc_combined_wide):

    frac_high_wide = pd.read_csv(frac_high_wide, index_col = False, header = 0)
    num_clust_slide_wide = pd.read_csv(num_clust_slide_wide, index_col = False, header = 0)

    all_prox_df_wide =  pd.read_csv(all_prox_df_wide, index_col = False, header =0)
    prox_indiv_schc_combined_wide = pd.read_csv(prox_indiv_schc_combined_wide, index_col = False, header = 0, sep = "\t")


    # Store features
    all_features = pd.merge(frac_high_wide, num_clust_slide_wide, on=[
                            "slide_submitter_id"])
    all_features = pd.merge(all_features, all_prox_df_wide)
    all_features = pd.merge(all_features, prox_indiv_schc_combined_wide)
    # all_features = pd.merge(all_features, shape_feature_means_wide)

    return all_features

def main(args):

    combine_clustering_features(frac_high_wide = args.frac_high_wide,
                                num_clust_slide_wide = args.num_clust_slide_wide, all_prox_df_wide = args.all_prox_df_wide, prox_indiv_schc_combined_wide = args.prox_indiv_schc_combined_wide).to_csv(
        Path(args.output_dir, f"{args.prefix}_clustering_features.csv"), index=False)

if __name__ == "__main__":
    args = get_args()
    st = time.time()
    main(args)
    rt = time.time() - st
    print(f"Script finished in {rt // 60:.0f}m {rt % 60:.0f}s")
