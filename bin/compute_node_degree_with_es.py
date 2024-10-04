#!/usr/bin/env python3
import argparse
import multiprocessing
import os
import time
from argparse import ArgumentParser as AP
from os.path import abspath
from pathlib import Path

import features.features as features
import features.graphs as graphs
import features.utils as utils
import joblib
import pandas as pd
from joblib import Parallel, delayed
from model.constants import DEFAULT_CELL_TYPES


def get_args():
    # Script description
    description = """Spatial Network (graph) features: Compute 'mean_ND' and 'ND_effsize' features"""

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
        "--max_dist",
        type=int,
        help="Maximum distance between clusters",
        required=False,
        default=None,
    )
    parser.add_argument(
        "--max_n_tiles_threshold",
        type=int,
        help="Number of tiles for computing max. distance between two points in two different clusters",
        default=2,
        required=False,
    )
    parser.add_argument(
        "--tile_size",
        type=int,
        help="Size of tile (default: 512)",
        default=512,
        required=False,
    )
    parser.add_argument(
        "--overlap",
        type=int,
        help="Overlap of tiles (default: 50)",
        default=50,
        required=False,
    )

    parser.add_argument(
        "--metadata_path",
        type=str,
        help="Path to tab-separated file with metadata",
        default="",
    )
    parser.add_argument(
        "--is_TCGA",
        type=bool,
        help="dataset is from TCGA (default: True)",
        default=True,
        required=False,
    )
    parser.add_argument(
        "--merge_var",
        type=str,
        help="Variable to merge metadata and computed features on",
        default=None,
    )
    parser.add_argument(
        "--sheet_name",
        type=str,
        help="Name of sheet for merging in case a path to xls(x) file is given for metadata_path",
        default=None,
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


def compute_node_degree_with_es(
    tile_quantification_path,
    cell_types=None,
    graphs_path=None,
    shapiro_alpha=0.05,
    n_cores=multiprocessing.cpu_count(),
):
    """
    Compute network features
    'mean_ND': Average number of neighbor nodes of cell type B surrounding nodes of cell type A. Nodes are assigned a cell type label if the predicted probability of that node (tile) for the given cell type is higher than 0.5.
    'ND_effsize': Cohen's d measure of effect size computed comparing the mean_ND(A,B) with the null distribution obtained by recomputing the mean_ND randomly assigning the A or B cell type label to each node preserving the total number of cell type A and B nodes in the network. For a negative effect size, the true average mean_ND(A,B) is larger than the simulated average mean_ND(A,B) meaning that the two cell types in the actual slide are closer together compared to a random distribution of these two cell types. Vice versa for a positive effect size. Nodes are assigned a cell type label as described above.


    Args:
        tile_quantification_path (str): path to tile quantification path (csv)
        cell_types (list): list of cell types (or path to csv file with cell types)
        graphs_path (str): path to pkl file with generated graphs [optional]
        shapiro_alpha (float): significance level for shapiro tests for normality (default=0.05)

    Returns:
        all_effect_sizes (DataFrame): dataframe containing the slide_submitter_id, center, neighbor, effect_size (Cohen's d), Tstat, pval, and the pair (string of center and neighbor)
        all_sims_nd (DataFrame): dataframe containing slide_submitter_id, center, neighbor, simulation_nr and degree (node degree)
        all_mean_nd_df (DataFrame): dataframe containing slide_submitter_id, center, neighbor, mean_sim (mean node degree across the N simulations), mean_obs

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

    ###############################################
    # ---- Compute ES based on ND difference ---- #
    ###############################################
    nd_results = Parallel(n_jobs=n_cores)(
        delayed(features.node_degree_wrapper)(all_graphs[id], predictions, id)
        for id in slide_submitter_ids
    )
    nd_results = list(filter(lambda id: id is not None, nd_results))

    # Format results
    all_sims_nd = []
    all_mean_nd_df = []
    example_simulations = {}

    for sim_assignments, sim, mean_nd_df in nd_results:
        all_mean_nd_df.append(mean_nd_df)
        all_sims_nd.append(sim)
        example_simulations.update(sim_assignments)

    all_sims_nd = pd.concat(all_sims_nd, axis=0).reset_index()
    all_mean_nd_df = pd.concat(all_mean_nd_df).reset_index(drop=True)

    # Testing normality
    shapiro_tests = Parallel(n_jobs=n_cores)(
        delayed(utils.test_normality)(
            sims_nd=all_sims_nd,
            slide_submitter_id=id,
            alpha=shapiro_alpha,
            cell_types=cell_types,
        )
        for id in all_sims_nd.slide_submitter_id.unique()
    )
    all_shapiro_tests = pd.concat(shapiro_tests, axis=0)

    # Computing Cohen's d effect size and perform t-test
    effect_sizes = Parallel(n_jobs=n_cores)(
        delayed(features.compute_effect_size)(
            all_mean_nd_df, all_sims_nd, slide_submitter_id
        )
        for slide_submitter_id in all_sims_nd.slide_submitter_id.unique()
    )
    all_effect_sizes = pd.concat(effect_sizes, axis=0)
    all_effect_sizes["pair"] = [
        f"{c}-{n}" for c, n in all_effect_sizes[["center", "neighbor"]].to_numpy()
    ]

    return (
        all_effect_sizes,
        all_sims_nd,
        all_mean_nd_df,
        example_simulations,
        all_shapiro_tests,
        all_graphs,
    )


def main(args):
    (
        all_effect_sizes,
        all_sims_nd,
        all_mean_nd_df,
        example_simulations,
        all_shapiro_tests,
        all_graphs,
    ) = compute_node_degree_with_es(
        tile_quantification_path=args.tile_quantification_path,
        cell_types=args.cell_types_path,
        graphs_path=args.graphs_path,
        shapiro_alpha=args.shapiro_alpha,
        n_cores=args.n_cores,
    )

    all_effect_sizes.to_csv(
        Path(args.output_dir, f"{args.prefix}_features_ND_ES.csv"),
        sep="\t",
        index=False,
    )
    all_sims_nd.to_csv(
        Path(args.output_dir, f"{args.prefix}_features_ND_sims.csv"),
        sep="\t",
        index=False,
    )
    all_mean_nd_df.to_csv(
        Path(args.output_dir, f"{args.prefix}_features_ND.csv"), sep="\t", index=False
    )
    joblib.dump(
        example_simulations,
        Path(args.output_dir, f"{args.prefix}_features_ND_sim_assignments.pkl"),
    )
    all_shapiro_tests.to_csv(
        Path(args.output_dir, f"{args.prefix}_shapiro_tests.csv"), index=False, sep="\t"
    )

    if args.graphs_path is None:
        joblib.dump(all_graphs, Path(args.output_dir, f"{args.prefix}_graphs.pkl"))


if __name__ == "__main__":
    args = get_args()
    st = time.time()
    main(args)
    rt = time.time() - st
    print(f"Script finished in {rt // 60:.0f}m {rt % 60:.0f}s")
