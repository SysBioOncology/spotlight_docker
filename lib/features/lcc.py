import sys
import os
import argparse
from argparse import ArgumentParser as AP
import networkx as nx
import pandas as pd
import time
from os.path import abspath

# Own modules
import features.features as features
import features.utils as utils
from model.constants import DEFAULT_CELL_TYPES

# Point to folder with custom imports
sys.path.append(f"{os.path.dirname(os.getcwd())}/Python/libs")

# def get_args():
#     # Script description
#     description = """Computing LCC"""

#     # Add parser
#     parser = AP(description=description,
#                 formatter_class=argparse.RawDescriptionHelpFormatter)

#     # Sections
#     parser.add_argument(
#         "--clinical_files_input",
#         help="Path to either a folder for multiple cancer types or single txt file.", required=False,
#         default=None
#     )
#     # TODO add arguments
#     parser.add_argument("--version", action="version", version="0.1.0")
#     arg = parser.parse_args()
#     arg.output = abspath(arg.output)
#     return arg

def determine_lcc(graph, cell_type_assignments, cell_types=None):
    """ Determine the fraction of the largest connected component (LCC) of a
    cell type w.r.t. to all nodes (tiles) of that cell type.
    1. Determine the number of nodes N in the LCC for the probability map of a
    cell type.
    2. Determine the total number of nodes (tiles) T for that cell type
    3. Determine the fraction of nodes that are connected: N/T

    Args:
        graph (Networkx Graph): graph representing the slide constructed with Networkx
        cell_type_assignments (DataFrame): Dataframe containing the cell type labels of the individual tiles indicated with booleans based on P > threshold
        cell_types (list): list of cell types
    """
    if cell_types is None:
        cell_types = DEFAULT_CELL_TYPES

    lcc = []
    for cell_type in cell_types:
        graph_temp = graph.copy()
        graph_temp.remove_nodes_from(
            list(cell_type_assignments[~cell_type_assignments[cell_type]].index)
        )
        if len(graph_temp.nodes()) > 0:
            # Get largest component
            # include only cell type specific tiles
            lcc_frac = len(max(nx.connected_components(graph_temp), key=len)) / len(
                graph_temp.nodes()
            )
            lcc.append([cell_type, lcc_frac])
    lcc = pd.DataFrame(lcc, columns=["cell_type", "type_spec_frac"])
    return lcc


def lcc_wrapper(id, slide_data, predictions, graph, cell_types, abundance_threshold):
    slide_data = utils.get_slide_data(predictions, id)
    node_cell_types = utils.assign_cell_types(
        slide_data=slide_data, cell_types=cell_types, threshold=abundance_threshold)
    lcc = features.determine_lcc(
        graph=graph, cell_type_assignments=node_cell_types, cell_types=cell_types
    )
    lcc["slide_submitter_id"] = id



# def main(args):
#     if not os.path.isdir(args.output_dir):
#         os.mkdir(args.output_dir)
#     lcc_wrapper(args.id, args.slide_data, args.predictions, args.graph, args.cell_types, args.abundance_threshold)

# if __name__ == "__main__":
#     args = get_args()
#     st = time.time()
#     main(args)
#     rt = time.time() - st
#     print(f"Script finished in {rt // 60:.0f}m {rt % 60:.0f}s")

