#!/usr/bin/env python3
import os
import pandas as pd
from joblib import Parallel, delayed
import argparse
import features.features as features
from model.constants import DEFAULT_CELL_TYPES


import multiprocessing
from argparse import ArgumentParser as AP
from os.path import abspath
import time
from pathlib import Path

def get_args():
      # Script description
    description = """Compute Clustering features Features: Compute proximity (between clusters)"""

    # Add parser
    parser = AP(description=description,
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--slide_clusters", type=str,
                        help="Path to csv file", required=True)
    parser.add_argument("--tiles_schc", type=str,
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
    parser.add_argument("--max_dist", type=int,
                        help="Max dist", required=False, default= None)

    parser.add_argument("--max_n_tiles_threshold", type=int,
                        help="Max dist", required=False, default= 2)
    parser.add_argument("--tile_size", type=int,
                        help="Max dist", required=False, default= 512)
    parser.add_argument("--overlap", type=int,
                        help="Max dist", required=False, default= 50)

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


def compute_proximity_from_indiv_schc_between(
        slide_clusters, tiles_schc,
        cell_types=None,
        n_clusters=8, max_dist=None,
        max_n_tiles_threshold=2,
        tile_size=512,
        overlap=50,
        n_cores = multiprocessing.cpu_count()):

    if cell_types is None:
        cell_types = DEFAULT_CELL_TYPES

    all_slide_indiv_clusters = pd.read_csv(slide_clusters, sep = ",", header = 0)
    slide_submitter_ids = list(set(all_slide_indiv_clusters.slide_submitter_id))

    slide_indiv_clusters_labeled = pd.read_csv(tiles_schc,  sep = ",", header = 0)

    ##########################################################################
    # ---- Compute proximity features (individual cell type clustering) ---- #
    ##########################################################################

    # Computing proximity for clusters derived for each cell type individually
    # Between clusters
    slide_submitter_ids = list(set(slide_indiv_clusters_labeled.slide_submitter_id))
    results_schc_indiv = Parallel(n_jobs=n_cores)(delayed(features.compute_proximity_clusters_pairs)(all_slide_indiv_clusters,
                                                                                                      slide_submitter_id=id, method="individual_between",
                                                                                                       n_clusters=n_clusters,
                                                                                                       cell_types=cell_types,
                                                                                                       max_dist=max_dist,
                                                                                                       max_n_tiles_threshold=max_n_tiles_threshold,
                                                                                                        tile_size=tile_size,
                                                                                                        overlap=overlap) for id in slide_submitter_ids)
    prox_indiv_schc = pd.concat(results_schc_indiv)

    # Formatting
    prox_indiv_schc = pd.merge(prox_indiv_schc, slide_indiv_clusters_labeled, left_on=[
                               "slide_submitter_id", "cluster1_label", "cluster1"], right_on=["slide_submitter_id", "cell_type_map", "cluster_label"])
    prox_indiv_schc = prox_indiv_schc.drop(
        columns=["cell_type_map", "cluster_label"])
    prox_indiv_schc = prox_indiv_schc.rename(
        columns={"is_high": "cluster1_is_high"})
    prox_indiv_schc = pd.merge(prox_indiv_schc, slide_indiv_clusters_labeled, left_on=[
                               "slide_submitter_id", "cluster2_label", "cluster2"], right_on=["slide_submitter_id", "cell_type_map", "cluster_label"])
    prox_indiv_schc = prox_indiv_schc.rename(
        columns={"is_high": "cluster2_is_high"})
    prox_indiv_schc = prox_indiv_schc.drop(
        columns=["cell_type_map", "cluster_label"])

    # Order matters
    prox_indiv_schc["ordered_pair"] = [
        f"{i}-{j}" for i, j in prox_indiv_schc[["cluster1_label", "cluster2_label"]].to_numpy()]
    prox_indiv_schc["comparison"] = [
        f"cluster1={i}-cluster2={j}" for i, j in prox_indiv_schc[["cluster1_is_high", "cluster2_is_high"]].to_numpy()]

    # Post-processing
    results_schc_indiv = pd.concat(Parallel(n_jobs=n_cores)(delayed(features.post_processing_proximity)(
        prox_df=prox_indiv_schc, slide_submitter_id=id, method="individual_between") for id in slide_submitter_ids))

    return results_schc_indiv

def main(args):
    compute_proximity_from_indiv_schc_between(
        slide_clusters = args.slide_clusters,
        tiles_schc = args.tiles_schc,
        cell_types=args.cell_types_path,
        n_cores = args.n_cores,
        n_clusters=args.n_clusters,
        max_dist=args.max_dist,
        max_n_tiles_threshold=args.max_n_tiles_threshold,
        tile_size=args.tile_size,
        overlap=args.overlap).to_csv(
        Path(args.output_dir, f"{args.prefix}_features_clust_indiv_schc_prox_between.csv"),index=False)


if __name__ == "__main__":
    args = get_args()
    st = time.time()
    main(args)
    rt = time.time() - st
    print(f"Script finished in {rt // 60:.0f}m {rt % 60:.0f}s")
