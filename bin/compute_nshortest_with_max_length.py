#!/usr/bin/env python3
import argparse
from argparse import ArgumentParser as AP
import os
import sys
import joblib
import pandas as pd
from joblib import Parallel, delayed
import argparse
from os import path
import multiprocessing
# Own modules
import features.clustering as clustering
import features.features as features
import features.graphs as graphs
import features.utils as utils
from model.constants import DEFAULT_SLIDE_TYPE, DEFAULT_CELL_TYPES, NUM_CORES, METADATA_COLS

from os.path import abspath
import time
from pathlib import Path


def get_args():
    # Script description
    description = """Compute Spatial Network Features: Compute Connectedness"""

    # Add parser
    parser = AP(description=description,
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tile_quantification_path", type=str,
                        help="Path to csv file with tile-level quantification (predictions)", required=True)
    parser.add_argument("--output_dir", type=str,
                        help="Path to output folder to store generated files", required=False, default="")
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
    parser.add_argument("--n_cores", type=int,
                        help="Number of cores to use (parallelization)")

    parser.add_argument("--version", action="version", version="0.1.0")
    arg = parser.parse_args()
    arg.output_dir = abspath(arg.output_dir)

    if ((arg.output_dir != "") & (not os.path.isdir(arg.output_dir))):
        # Create an empty folder for TF records if folder doesn't exist
        arg.output_dir = Path(arg.output_dir, "process_train")
        os.mkdir(arg.output_dir)

        if arg.n_cores is None:
            arg.n_cores = multiprocessing.cpu_count()
    return arg


def compute_n_shortest_paths_with_max_length(tile_quantification_path,
                                             cell_types=None,
                                             graphs_path=None,
                                             cutoff_path_length=2,
                                             n_cores=multiprocessing.cpu_count()):
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
    if cell_types is None:
        cell_types = DEFAULT_CELL_TYPES

    predictions = pd.read_csv(tile_quantification_path, sep="\t")
    slide_submitter_ids = list(set(predictions.slide_submitter_id))

    #####################################
    # ---- Constructing the graphs ---- #
    #####################################

    # TODO use 'generate_graphs.py' for this to replace
    if graphs_path is None:
        results = Parallel(n_jobs=n_cores)(
            delayed(graphs.construct_graph)(
                predictions=predictions, slide_submitter_id=id)
            for id in slide_submitter_ids
        )
        # Extract/format graphs
        all_graphs = {
            list(slide_graph.keys())[0]: list(slide_graph.values())[0]
            for slide_graph in results
        }
    else:
        all_graphs = joblib.load(graphs_path)

    #######################################################
    # ---- Compute N shortest paths with max. length ---- #
    #######################################################

    results = Parallel(n_jobs=n_cores)(
        delayed(features.compute_n_shortest_paths_max_length)(
            predictions=predictions, slide_submitter_id=id, graph=all_graphs[
                id], cutoff=cutoff_path_length
        )
        for id in slide_submitter_ids
    )
    # Formatting and count the number of shortest paths with max length
    all_shortest_paths_thresholded = pd.concat(results, axis=0)
    all_shortest_paths_thresholded["n_paths"] = 1
    proximity_graphs = (
        all_shortest_paths_thresholded.groupby(
            ["slide_submitter_id", "source", "target"]
        )
        .sum(numeric_only=True)
        .reset_index()
    )
    # Post-processing
    proximity_graphs["pair"] = [f"{source}-{target}" for source,
                                target in proximity_graphs[["source", "target"]].to_numpy()]
    proximity_graphs = proximity_graphs.drop(columns=["path_length"])

    # Additional formatting
    shortest_paths_wide = proximity_graphs.pivot(
        index=["slide_submitter_id"], columns="pair")["n_paths"]
    new_cols = [
        f'Prox graph {col.replace("_", " ")} clusters' for col in shortest_paths_wide.columns]
    shortest_paths_wide.columns = new_cols
    shortest_paths_wide = shortest_paths_wide.reset_index()

    return (proximity_graphs, shortest_paths_wide, all_graphs)


def main(args):
    proximity_graphs, shortest_paths_wide, all_graphs = compute_n_shortest_paths_with_max_length(
        tile_quantification_path=args.tile_quantification_path,
        cell_types=args.cell_types_path,
        graphs_path=args.graphs_path, cutoff_path_length=args.cutoff_path_length,
        n_cores=args.n_cores)

    proximity_graphs.to_csv(
        Path(args.output_dir,
             f"{args.prefix}_features_shortest_paths_thresholded.csv"),
        sep="\t",
        index=False)

    shortest_paths_wide.to_csv(
        Path(args.output_dir,
             f"{args.prefix}_features_shortest_paths_thresholded_wide.csv"),
        sep="\t",
        index=False)

    if (args.graphs_path is None):
        joblib.dump(all_graphs,
                    Path(args.output_dir,
                         f"{args.prefix}_graphs.pkl"))


if __name__ == "__main__":
    args = get_args()
    st = time.time()
    main(args)
    rt = time.time() - st
    print(f"Script finished in {rt // 60:.0f}m {rt % 60:.0f}s")
