#!/usr/bin/env python3
import argparse
import multiprocessing
from argparse import ArgumentParser as AP
import os
import joblib
import pandas as pd
from joblib import Parallel, delayed
import argparse

# Own modules
import features.graphs as graphs
from model.constants import DEFAULT_CELL_TYPES

from os.path import abspath
import time
from pathlib import Path

def get_args():
      # Script description
    description = """Generating graphs for computing spatial network features"""

    # Add parser
    parser = AP(description=description,
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tile_quantification_path", type=str,
                        help="Path to csv file with tile-level quantification (predictions)", required=True)
    parser.add_argument("--output_dir", type=str,
                        help="Path to output folder to store generated files", default = "")
    parser.add_argument("--slide_type", type=str,
                        help="Type of slides 'FFPE' or 'FF' used for naming generated files (default: 'FF')", default="FF")
    parser.add_argument("--cell_types_path", type=str,
                        help="Path to file with list of cell types (default: CAFs, endothelial_cells, T_cells, tumor_purity)", default=None)
    parser.add_argument("--prefix", type=str,
                        help="Prefix for output file", default="")
    parser.add_argument("--n_cores", type = int, help = "Number of cores to use (parallelization)")
    parser.add_argument("--version", action="version", version="0.1.0")
    arg = parser.parse_args()
    arg.output_dir = abspath(arg.output_dir)

    if ((arg.output_dir != "") & (not os.path.isdir(arg.output_dir))):
        os.mkdir(arg.output_dir)

    if arg.n_cores is None:
        arg.n_cores = multiprocessing.cpu_count()
    return arg

def generate_graphs(tile_quantification_path, cell_types=None, n_cores = multiprocessing.cpu_count()):
    """
    Generating graphs

    Args:
        tile_quantification_path (str)
        cell_types (list): list of cell types
        n_cores (int): Number of cores to use (parallelization)

    Returns:
    Graphs for all slides (dict)

    """
    if cell_types is None:
        cell_types = DEFAULT_CELL_TYPES

    predictions = pd.read_csv(tile_quantification_path, sep="\t")
    slide_submitter_ids = list(set(predictions.slide_submitter_id))

    #####################################
    # ---- Constructing the graphs ---- #
    #####################################

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
    return all_graphs

def main(args):
    all_graphs = generate_graphs(
        tile_quantification_path = args.tile_quantification_path,
        n_cores=args.n_cores)
    out_filepath =  Path(args.output_dir,
                     f"{args.prefix}_graphs.pkl")

    joblib.dump(all_graphs, out_filepath)
    print(f"Generated all graphs and stored in: {out_filepath}")

if __name__ == "__main__":
    args = get_args()
    st = time.time()
    main(args)
    rt = time.time() - st
    print(f"Script finished in {rt // 60:.0f}m {rt % 60:.0f}s")
