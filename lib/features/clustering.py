import sys
import os
import networkx as nx
import pandas as pd
from sklearn.cluster import AgglomerativeClustering
from sklearn.preprocessing import StandardScaler

# Point to folder with custom imports
sys.path.append(f"{os.path.dirname(os.getcwd())}/Python/libs")

from model.constants import *
import features.utils as utils

def schc_all(predictions, graph, slide_submitter_id, n_clusters=8, cell_types=None):
    """
    Spatially Constrained Hierarchical Clustering
    Using the adjacency matrix (connectivity matrix) from the graph of a slide as a condition for merging clusters.

    Args:
        predictions (DataFrame): dataframe containing the predicted abundances for the cell types
        graph (Networkx Graph): graph representing the slide constructed with Networkx
        slide_submitter_id (str): e.g. TCGA-77-6842-01A-01-BS1
        n_clusters (int): number of clusters to create
        cell_types (list): list of cell types, e.g. CAFs, T_cells, endothelial_cells, tumor purity

    Returns:
        slide_data (DataFrame): dataframe containing the predictions for the given slide_submitter_id for each tile and its cluster label (max. value=n_clusters)
    """
    if cell_types is None:
        cell_types = DEFAULT_CELL_TYPES
    slide_data = utils.get_slide_data(predictions, slide_submitter_id)
    var_to_select = cell_types + ["Coord_X", "Coord_Y"]
    scaled = StandardScaler().fit_transform(slide_data[var_to_select])

    # Clustering setup
    model = AgglomerativeClustering(
        linkage="ward",
        connectivity=nx.adjacency_matrix(graph),
        n_clusters=n_clusters,
    )
    # Fit algorithm to the data
    model.fit(scaled)
    # Add label to data for given slide_submitter_id
    slide_data["cluster_label"] = model.labels_
    return slide_data


# Fixed number of clusters
def schc_individual(
    predictions, graph, slide_submitter_id, n_clusters=8, cell_types=None
):
    """
       Spatially Constrained Hierarchical Clustering using one cell type at a time
       Using the adjacency matrix (connectivity matrix) from the graph of a slide as a condition for merging clusters.

    Args:
        predictions (DataFrame): dataframe containing the predicted abundances for the cell types
        graph (Networkx Graph): graph representing the slide constructed with Networkx
        slide_submitter_id (str): e.g. TCGA-77-6842-01A-01-BS1
        n_clusters (int): number of clusters to create
        cell_types (list): list of cell types, e.g. CAFs, T_cells, endothelial_cells, tumor purity

       Returns:
           slide_data (DataFrame): dataframe containing the tiles for the given slide and their cluster labels (max. value=n_clusters)
    """
    if cell_types is None:
        cell_types = DEFAULT_CELL_TYPES

    slide_data = utils.get_slide_data(predictions, slide_submitter_id)
    for cell_type in cell_types:
        var_to_select = [cell_type, "Coord_X", "Coord_Y"]
        scaled = StandardScaler().fit_transform(slide_data[var_to_select])

        model = AgglomerativeClustering(
            linkage="ward",
            connectivity=nx.adjacency_matrix(graph),
            n_clusters=n_clusters,
        )
        # Fit algorithm to the data
        model.fit(scaled)
        slide_data["{}_label".format(cell_type)] = model.labels_
    vars_to_select = ["{}_label".format(i) for i in cell_types]
    return slide_data[["tile_ID"] + vars_to_select]


def characterize_clusters(clusters, cell_types=None):
    """
    Determine the cell type labels for each cluster in a slide based on the overall mean of all cluster means across all slides
    (For simultaneous clustering)
    Args:
        clusters (DataFrame): Dataframe containing the predictions and the cluster id for each tile for each slide.
        cell_types (list): list of cell types

    Returns:
        (DataFrame): Dataframe containing slide_submitter_id and cluster_label as id variables, and a column for each cell type containing booleans, True if value > overall mean of cluster means, and otherwise False
    """
    if cell_types is None:
        cell_types = DEFAULT_CELL_TYPES

    # Average of prediction across all clusters for each cell type
    cluster_means = (
        clusters.groupby(["slide_submitter_id", "cluster_label"])
        .mean()[cell_types]
        .mean()
    )

    # Determine for each cluster the cell types
    return (
        clusters.groupby(["slide_submitter_id", "cluster_label"]).mean()[
            cell_types
        ]
        > cluster_means
    ).reset_index()


def label_cell_type_map_clusters(clusters, cell_types=None):
    """
    Label the clusters from the individuall cell type maps with 'high' or 'low'
    Labelling only when the cluster mean probablility (prediction) is greater then the average mean of the cell type prediction across all clusters

    Args:
        clusters (DataFrame): clusters for each cell type map (in long format) from the function schc_individual() should also contain the metadata
        cell_types (list): list of cell types

    Returns:
         Returns:
        (DataFrame): Long format dataframe containing slide_submitter_id and cluster_label as id variables, and a column containing the abundance label ('high' or 'low').
    """
    if cell_types is None:
        cell_types = DEFAULT_CELL_TYPES

    clusters = clusters[
        [ "slide_submitter_id"]
        + ["{}_label".format(i) for i in cell_types]
        + cell_types
    ]
    clusters_long = clusters.melt(
        id_vars=[ "slide_submitter_id"] + cell_types,
        value_vars=["{}_label".format(i) for i in cell_types],
        value_name="cluster_label",
        var_name="cell_type_map",
    )
    clusters_long["cell_type_map"] = clusters_long["cell_type_map"].replace(
        dict(zip(["{}_label".format(i) for i in cell_types], cell_types))
    )

    # Means of each cluster for each cell type for each slide
    slide_means = (
        clusters_long.groupby(
            ["slide_submitter_id", "cell_type_map", "cluster_label"]
        )
        .mean()
        .reset_index()
    )

    # Mean per cell type across clusters
    overall_cluster_means = (
        slide_means.groupby(["cell_type_map"]).mean()[cell_types].to_numpy().diagonal()
    )
    overall_cluster_means = dict(zip(cell_types, overall_cluster_means))

    # Labelling the clusters for each cell type
    all_labeled_slides = []
    for cell_type in cell_types:
        subset = slide_means.loc[slide_means["cell_type_map"] == cell_type].copy()
        subset.loc[:, "is_high"] = (
            subset.loc[:, cell_type] > overall_cluster_means[cell_type]
        )
        subset = subset.drop(columns=cell_types)
        all_labeled_slides.append(subset)

    return pd.concat(all_labeled_slides, axis=0)
