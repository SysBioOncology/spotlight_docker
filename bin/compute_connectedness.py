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
    description = """Compute Spatial Network Features: Compute LCC"""

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


def compute_connectedness(
    tile_quantification_path,
    cell_types=None,
    graphs_path=None,
    abundance_threshold=0.5,
    n_cores=multiprocessing.cpu_count(),
):
    """
    Compute network features: LCC
    Normalized size of the largest connected component (LCC) for cell type A. This is defined as the largest set of nodes of cell type A
    connected with at least one path between every pair of nodes,
    divided by the total number of nodes of cell type A. Nodes are assigned a cell type label as described above.

    Args:
        tile_quantification_path (str)
        cell_types (list): list of cell types
        graphs_path (str): path to pkl file with generated graphs [optional]
        abundance_threshold (float): threshold for assigning cell types to tiles based on the predicted probability (default=0.5)

    Returns:
        all_largest_cc_sizes (DataFrame): dataframe containing slide_submitter_id, cell type and type_spec_frac (fraction of LCC w.r.t. all tiles for cell type)
        all_graphs (dict): dictionary with graphs [only if not existing]

    """
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

    #######################################################
    # ---- Compute connectedness and co-localization ---- #
    #######################################################

    all_largest_cc_sizes = []
    for id in slide_submitter_ids:
        slide_data = utils.get_slide_data(predictions, id)
        node_cell_types = utils.assign_cell_types(
            slide_data=slide_data, cell_types=cell_types, threshold=abundance_threshold
        )
        lcc = features.determine_lcc(
            graph=all_graphs[id],
            cell_type_assignments=node_cell_types,
            cell_types=cell_types,
        )
        lcc["slide_submitter_id"] = id
        all_largest_cc_sizes.append(lcc)

    all_largest_cc_sizes = pd.concat(all_largest_cc_sizes, axis=0)
    all_largest_cc_sizes = all_largest_cc_sizes.reset_index(drop=True)
    all_largest_cc_sizes_wide = all_largest_cc_sizes.pivot(
        index=["slide_submitter_id"], columns="cell_type"
    )["type_spec_frac"]
    new_cols = [
        f'LCC {col.replace("_", " ")} clusters'
        for col in all_largest_cc_sizes_wide.columns
    ]
    all_largest_cc_sizes_wide.columns = new_cols
    all_largest_cc_sizes_wide = all_largest_cc_sizes_wide.reset_index()

    return (all_largest_cc_sizes_wide, all_graphs)


def main(args):
    all_largest_cc_sizes_wide, all_graphs = compute_connectedness(
        tile_quantification_path=args.tile_quantification_path,
        cell_types=args.cell_types_path,
        graphs_path=args.graphs_path,
        abundance_threshold=args.abundance_threshold,
        n_cores=args.n_cores,
    )
    all_largest_cc_sizes_wide.to_csv(
        Path(args.output_dir, f"{args.prefix}_features_lcc_fraction_wide.csv"),
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
