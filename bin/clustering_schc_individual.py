#!/usr/bin/env python3
import argparse
import multiprocessing
import os
import time
from argparse import ArgumentParser as AP
from os.path import abspath
from pathlib import Path

# Own modules
import features.clustering as clustering
import features.graphs as graphs
import joblib
import pandas as pd
from joblib import Parallel, delayed
from model.constants import DEFAULT_CELL_TYPES


def get_args():
    # Script description
    description = """Spatial Clustering with SCHC (individual)"""

    # Add parser
    parser = AP(
        description=description, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--tile_quantification_path",
        type=str,
        help="Path to csv file with tile-level quantification (predictions)",
        required=True,
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


def clustering_schc_individual(
    tile_quantification_path: str,
    cell_types=None,
    graphs_path=None,
    n_cores=multiprocessing.cpu_count,
):
    if cell_types is None:
        cell_types = DEFAULT_CELL_TYPES

    predictions = pd.read_csv(tile_quantification_path, sep="\t")
    slide_submitter_ids = list(set(predictions.slide_submitter_id))

    #####################################
    # ---- Constructing the graphs ---- #
    #####################################

    if graphs_path is None:
        results = Parallel(n_jobs=n_cores)(
            delayed(graphs.construct_graph)(
                predictions=predictions, slide_submitter_id=id
            )
            for id in slide_submitter_ids
        )
        # Extract/format graphs
        all_graphs = {
            list(slide_graph.keys())[0]: list(slide_graph.values())[0]
            for slide_graph in results
        }

    else:
        all_graphs = joblib.load(graphs_path)

    ######################################################################################
    # ---- Fraction of highly abundant cell types (individual cell type clustering) ---- #
    ######################################################################################

    # Spatially Hierarchical Constrained Clustering with all quantification of all cell types for each individual cell type
    slide_indiv_clusters = Parallel(n_jobs=n_cores)(
        delayed(clustering.schc_individual)(predictions, all_graphs[id], id)
        for id in slide_submitter_ids
    )
    all_slide_indiv_clusters = pd.concat(slide_indiv_clusters, axis=0)

    # Add metadata
    all_slide_indiv_clusters = pd.merge(
        predictions, all_slide_indiv_clusters, on="tile_ID"
    )

    # Add abundance label 'high' or 'low' based on cluster means
    slide_indiv_clusters_labeled = clustering.label_cell_type_map_clusters(
        all_slide_indiv_clusters
    )

    all_slide_indiv_clusters_final = all_slide_indiv_clusters.drop(
        axis=1, columns=cell_types
    )  # drop the predicted probabilities

    return (
        all_slide_indiv_clusters,
        slide_indiv_clusters_labeled,
        all_slide_indiv_clusters_final,
        all_graphs,
    )


def main(args):
    (
        all_slide_indiv_clusters,
        slide_indiv_clusters_labeled,
        all_slide_indiv_clusters_final,
        all_graphs,
    ) = clustering_schc_individual(
        tile_quantification_path=args.tile_quantification_path,
        cell_types=args.cell_types_path,
        graphs_path=args.graphs_path,
        n_cores=args.n_cores,
    )

    all_slide_indiv_clusters.to_csv(
        Path(args.output_dir, f"{args.prefix}_indiv_schc_tiles_raw.csv", index=False)
    )

    all_slide_indiv_clusters_final.to_csv(
        Path(args.output_dir, f"{args.prefix}_indiv_schc_tiles.csv"), index=False
    )

    slide_indiv_clusters_labeled.to_csv(
        Path(
            args.output_dir,
            f"{args.prefix}_indiv_schc_clusters_labeled.csv",
            index=False,
        )
    )

    if args.graphs_path is None:
        joblib.dump(all_graphs, Path(args.output_dir, f"{args.prefix}_graphs.pkl"))


if __name__ == "__main__":
    args = get_args()
    st = time.time()
    main(args)
    rt = time.time() - st
    print(f"Script finished in {rt // 60:.0f}m {rt % 60:.0f}s")
