#!/usr/bin/env python3
import argparse
import multiprocessing
import os
import time
from argparse import ArgumentParser as AP
from os.path import abspath
from pathlib import Path

import features.features as features
import pandas as pd
from joblib import Parallel, delayed
from model.constants import DEFAULT_CELL_TYPES


def get_args():
    # Script description
    description = """Compute Spatial Network Features: Compute Connectedness"""

    # Add parser
    parser = AP(
        description=description, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--slide_clusters_characterized",
        type=str,
        help="Path to csv file",
        required=True,
    )
    parser.add_argument(
        "--tiles_schc", type=str, help="Path to csv file", required=True
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        help="Path to output folder to store generated files",
        required=False,
        default="",
    )
    parser.add_argument(
        "--slide_type",
        type=str,
        help="Type of slides 'FFPE' or 'FF' used for naming generated files (default: 'FF')",
        default="FF",
    )
    parser.add_argument(
        "--cell_types_path",
        type=str,
        help="Path to file with list of cell types (default: CAFs, endothelial_cells, T_cells, tumor_purity)",
        default=None,
    )
    parser.add_argument(
        "--graphs_path",
        type=str,
        help="Path to pkl with generated graphs in case this was done before (OPTIONAL) if not specified, graphs will be generated",
        default=None,
    )
    parser.add_argument("--prefix", type=str, help="Prefix for output file", default="")
    parser.add_argument(
        "--cutoff_path_length",
        type=int,
        help="Max path length for proximity based on graphs",
        default=2,
        required=False,
    )
    parser.add_argument(
        "--shapiro_alpha",
        type=float,
        help="Choose significance level alpha (default: 0.05)",
        default=0.05,
        required=False,
    )
    parser.add_argument(
        "--abundance_threshold",
        type=float,
        help="Threshold for assigning cell types based on predicted probability (default: 0.5)",
        default=0.5,
        required=False,
    )
    parser.add_argument(
        "--max_dist", type=int, help="Max dist", required=False, default=None
    )
    parser.add_argument(
        "--max_n_tiles_threshold", type=int, help="Max dist", required=False, default=2
    )
    parser.add_argument(
        "--tile_size", type=int, help="Max dist", required=False, default=512
    )
    parser.add_argument(
        "--overlap", type=int, help="Max dist", required=False, default=50
    )
    parser.add_argument(
        "--n_clusters",
        type=int,
        help="Number of clusters for SCHC (default: 8)",
        required=False,
        default=8,
    )
    parser.add_argument(
        "--n_cores", type=int, help="Number of cores to use (parallelization)"
    )

    parser.add_argument("--version", action="version", version="0.1.0")
    arg = parser.parse_args()
    arg.output_dir = abspath(arg.output_dir)

    if (arg.output_dir != "") & (not os.path.isdir(arg.output_dir)):
        # Create an empty folder for TF records if folder doesn't exist
        arg.output_dir = Path(arg.output_dir, "process_train")
        os.mkdir(arg.output_dir)

        if arg.n_cores is None:
            arg.n_cores = multiprocessing.cpu_count()
    return arg


def compute_proximity_from_simultaneous_schc(
    slide_clusters_characterized,
    tiles_schc,
    cell_types=DEFAULT_CELL_TYPES,
    n_clusters=8,
    max_dist=None,
    max_n_tiles_threshold=2,
    tile_size=512,
    overlap=50,
    n_cores=multiprocessing.cpu_count(),
):
    all_slide_clusters_characterized = pd.read_csv(
        slide_clusters_characterized, sep=",", header=0, index_col=0
    )

    slide_submitter_ids = list(set(all_slide_clusters_characterized.slide_submitter_id))

    tiles_all_schc = pd.read_csv(tiles_schc, sep=",", header=0, index_col=0)

    # Computing proximity for clusters derived with all cell types simultaneously
    clusters_all_schc_long = all_slide_clusters_characterized.melt(
        id_vars=["slide_submitter_id", "cluster_label"],
        value_name="is_assigned",
        var_name="cell_type",
    )
    # remove all cell types that are not assigned to the cluster
    clusters_all_schc_long = clusters_all_schc_long[
        clusters_all_schc_long["is_assigned"]
    ]
    clusters_all_schc_long = clusters_all_schc_long.drop(columns="is_assigned")

    results_schc_all = Parallel(n_jobs=n_cores)(
        delayed(features.compute_proximity_clusters_pairs)(
            tiles=tiles_all_schc,
            slide_submitter_id=id,
            n_clusters=n_clusters,
            cell_types=cell_types,
            max_dist=max_dist,
            max_n_tiles_threshold=max_n_tiles_threshold,
            tile_size=tile_size,
            overlap=overlap,
            method="all",
        )
        for id in slide_submitter_ids
    )
    prox_all_schc = pd.concat(results_schc_all)

    # Label clusters (a number) with the assigned cell types
    prox_all_schc = pd.merge(
        prox_all_schc,
        clusters_all_schc_long,
        left_on=["slide_submitter_id", "cluster1"],
        right_on=["slide_submitter_id", "cluster_label"],
    )
    prox_all_schc = prox_all_schc.rename(columns={"cell_type": "cluster1_label"})
    prox_all_schc = prox_all_schc.drop(columns=["cluster_label"])

    prox_all_schc = pd.merge(
        prox_all_schc,
        clusters_all_schc_long,
        left_on=["slide_submitter_id", "cluster2"],
        right_on=["slide_submitter_id", "cluster_label"],
    )
    prox_all_schc = prox_all_schc.rename(columns={"cell_type": "cluster2_label"})

    # Order doesn't matter: x <->
    prox_all_schc["pair"] = [
        f"{sorted([i, j])[0]}-{sorted([i, j])[1]}"
        for i, j in prox_all_schc[["cluster1_label", "cluster2_label"]].to_numpy()
    ]
    prox_all_schc = prox_all_schc[
        (
            (prox_all_schc.cluster1 == prox_all_schc.cluster2)
            & (prox_all_schc.cluster2_label != prox_all_schc.cluster1_label)
        )
        | (prox_all_schc.cluster1 != prox_all_schc.cluster2)
    ]

    # slides = prox_all_schc[["MFP", "slide_submitter_id"]].drop_duplicates().to_numpy()
    slide_submitter_ids = list(set(prox_all_schc.slide_submitter_id))

    # Post Processing
    results_schc_all = Parallel(n_jobs=n_cores)(
        delayed(features.post_processing_proximity)(
            prox_df=prox_all_schc, slide_submitter_id=id, method="all"
        )
        for id in slide_submitter_ids
    )
    all_prox_df = pd.concat(results_schc_all)
    # Remove rows with a proximity of NaN
    all_prox_df = all_prox_df.dropna(axis=0)

    all_prox_df_wide = all_prox_df.pivot(
        index=["slide_submitter_id"], columns=["pair"]
    )["proximity"]
    new_cols = [
        f'prox CC {col.replace("_", " ")} clusters' for col in all_prox_df_wide.columns
    ]
    all_prox_df_wide.columns = new_cols
    all_prox_df_wide = all_prox_df_wide.reset_index()

    return all_prox_df_wide


def main(args):
    compute_proximity_from_simultaneous_schc(
        slide_clusters_characterized=args.slide_clusters_characterized,
        tiles_schc=args.tiles_schc,
        cell_types=args.cell_types_path,
        n_cores=args.n_cores,
        n_clusters=args.n_clusters,
        max_dist=args.max_dist,
        max_n_tiles_threshold=args.max_n_tiles_threshold,
        tile_size=args.tile_size,
        overlap=args.overlap,
    ).to_csv(
        Path(args.output_dir, f"{args.prefix}_features_clust_all_schc_prox_wide.csv"),
        index=False,
    )


if __name__ == "__main__":
    args = get_args()
    st = time.time()
    main(args)
    rt = time.time() - st
    print(f"Script finished in {rt // 60:.0f}m {rt % 60:.0f}s")
