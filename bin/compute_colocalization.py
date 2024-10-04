#!/usr/bin/env python3
import argparse
import multiprocessing
import os
import time
from argparse import ArgumentParser as AP
from os.path import abspath
from pathlib import Path

# Own modules
import features.features as features
import features.graphs as graphs
import features.utils as utils
import joblib
import pandas as pd
from joblib import Parallel, delayed
from model.constants import DEFAULT_CELL_TYPES


def get_args():
    # Script description
    description = """Compute Spatial Network Features: Compute co-localization"""

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
        "--abundance_threshold",
        type=float,
        help="Threshold for assigning cell types based on predicted probability (default: 0.5)",
        default=0.5,
        required=False,
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


def compute_colocalization(
    tile_quantification_path,
    cell_types=None,
    graphs_path=None,
    abundance_threshold=0.5,
    n_cores=multiprocessing.cpu_count(),
):
    """
    Compute network features: co-localization

    Args:
        tile_quantification_path (str)
        output_dir (str)
        cell_types (list): list of cell types
        graphs_path (str): path to pkl file with generated graphs [optional]
        abundance_threshold (float): threshold for assigning cell types to tiles based on the predicted probability (default=0.5)

    Returns:
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

    #######################################################
    # ---- Compute connectedness and co-localization ---- #
    #######################################################

    all_dual_nodes_frac = []
    for id in slide_submitter_ids:
        slide_data = utils.get_slide_data(predictions, id)
        node_cell_types = utils.assign_cell_types(
            slide_data=slide_data, cell_types=cell_types, threshold=abundance_threshold
        )

        dual_nodes_frac = features.compute_dual_node_fractions(
            node_cell_types, cell_types
        )
        dual_nodes_frac["slide_submitter_id"] = id
        all_dual_nodes_frac.append(dual_nodes_frac)

    all_dual_nodes_frac = pd.concat(all_dual_nodes_frac, axis=0)

    colocalization_wide = all_dual_nodes_frac.pivot(
        index=["slide_submitter_id"], columns="pair"
    )["frac"]
    new_cols = [
        f'Co-loc {col.replace("_", " ")} clusters'
        for col in colocalization_wide.columns
    ]
    colocalization_wide.columns = new_cols
    colocalization_wide = colocalization_wide.reset_index()

    return (all_dual_nodes_frac, colocalization_wide, all_graphs)


def main(args):
    all_dual_nodes_frac, colocalization_wide, all_graphs = compute_colocalization(
        tile_quantification_path=args.tile_quantification_path,
        cell_types=args.cell_types_path,
        graphs_path=args.graphs_path,
        abundance_threshold=args.abundance_threshold,
        n_cores=args.n_cores,
    )

    all_dual_nodes_frac.to_csv(
        Path(args.output_dir, f"{args.prefix}_features_coloc_fraction.csv"),
        sep="\t",
        index=False,
    )

    colocalization_wide.to_csv(
        Path(args.output_dir, f"{args.prefix}_features_coloc_fraction_wide.csv"),
        sep="\t",
        index=False,
    )

    if args.graphs_path is None:
        joblib.dump(all_graphs, Path(args.output_dir, f"{args.prefix}_graphs.pkl"))


if __name__ == "__main__":
    args = get_args()
    st = time.time()
    main(args)
    rt = time.time() - st
    print(f"Script finished in {rt // 60:.0f}m {rt % 60:.0f}s")
