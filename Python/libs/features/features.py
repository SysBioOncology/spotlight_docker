"""
Compute different features from:
1. Graphs
2. Spatially Constrained Hierarchical Clustering (i.e. agglomerative clustering
with connectivity constraints)
"""
import itertools
import math
import statistics as statpy
import sys
import os
import networkx as nx
import numpy as np
import pandas as pd
import scipy.stats as stats
from scipy.spatial import ConvexHull
from sklearn.metrics import pairwise_distances_argmin_min

# Point to folder with custom imports
sys.path.append(f"{os.path.dirname(os.getcwd())}/Python/libs")

# Own modules
from model.constants import *
import features.utils as utils

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
    return pd.DataFrame(lcc, columns=["cell_type", "type_spec_frac"])


def compute_dual_node_fractions(cell_type_assignments, cell_types=None):
    """
    Co-localization of cell types within tiles
    1. Determine the number of tiles that are labeled with both cell types in
    a cell type pair
    2. Determine the number of tiles that are assigned with at least one of
    the cell types in the pair.
    3. Determine the fraction

    Args:
        cell_type_assignments (DataFrame): Dataframe containing the cell type labels of the individual tiles indicated with booleans based on P > threshold
        cell_types (list): list of cell types

    Returns:
        (DataFrame): dataframe containing the fraction and the number of tiles per cell type pair
    """
    if cell_types is None:
        cell_types = DEFAULT_CELL_TYPES
    cell_type_pairs = list(itertools.combinations(cell_types, 2))
    out = []
    for cell_type1, cell_type2 in cell_type_pairs:
        # Number of tiles that have both of the cell types in the pair
        num_overlap = sum(
            cell_type_assignments[[cell_type1, cell_type2]].sum(axis=1) == 2
        )

        # Number of tiles that are at least one of the cell types in the pair
        num_tiles = sum(cell_type_assignments[[cell_type1, cell_type2]].sum(axis=1) > 0)

        if num_tiles > 0:
            out.append(
                [f"{cell_type1}-{cell_type2}", num_overlap, num_overlap / num_tiles]
            )
    return pd.DataFrame(out, columns=["pair", "counts", "frac"])


def compute_n_shortest_paths_max_length(
    graph, cell_type_assignments=None, cell_types=None, cutoff=2, predictions=None, slide_submitter_id=None
):
    """
    Determine the number of shortest paths that have a path length of max N

    Args:
        graph (Networkx Graph): graph representing the slide constructed with Networkx
        cell_type_assignments (DataFrame): Dataframe containing the cell type labels of the individual tiles indicated with booleans based on P > threshold
        cell_types (list): list of cell types
        cutoff (int): max. path length

    Returns
        (DataFrame): dataframe containing the source and target cell types and the path length
    """
    if cell_types is None:
        cell_types = DEFAULT_CELL_TYPES
    cell_type_interactions = list(itertools.combinations(cell_types, 2)) + list(
        zip(cell_types, cell_types)
    )

    slide_data = utils.get_slide_data(predictions, slide_submitter_id)
    cell_type_assignments = utils.assign_cell_types(slide_data, cell_types)

    # Compute all shortest paths with max length of cutoff
    all_paths = nx.all_pairs_shortest_path(graph, cutoff=cutoff)
    all_paths = pd.DataFrame(
        all_paths,
        columns=["source", "target_paths"]
    )
    all_paths = all_paths.set_index("source")

    shortest_paths_diff = []
    if len(all_paths) > 0:
        temp = []
        # Extract shortest paths per source node
        for source in all_paths.index:
            target_paths_for_single_source = pd.DataFrame(
                pd.Series(all_paths.loc[source, "target_paths"]), columns=["path"]
            )
            # Compute path length from defined path
            target_paths_for_single_source["path_length"] = (
                target_paths_for_single_source["path"].str.len() - 1
            )
            # Only keep paths with a path length > 0
            paths_for_node = target_paths_for_single_source[
                target_paths_for_single_source["path_length"] != 0
            ].reset_index()
            paths_for_node = paths_for_node.rename(columns={"index": "target"})
            paths_for_node["source"] = source
            temp.append(paths_for_node)
        all_shortest_paths = pd.concat(temp, axis=0)
        # Select the shortest paths for each cell types (based on starting/source cell type)
        for source, target in cell_type_interactions:
            source_ids = list(
                cell_type_assignments[cell_type_assignments[source]].index
            )
            target_ids = list(
                cell_type_assignments[cell_type_assignments[target]].index
            )
            # Only select the source nodes (match node id with cell type)
            shortest_paths_subset = all_shortest_paths[
                (
                    all_shortest_paths["source"].isin(source_ids)
                    & all_shortest_paths["target"].isin(target_ids)
                    & all_shortest_paths["source"]
                )
            ]
            # Drop source, target columns with node ids (int)
            shortest_paths_subset = shortest_paths_subset.drop(
                columns=["source", "target"]
            )
            # Add metadata: cell types instead of node ids
            shortest_paths_subset["source"] = source
            shortest_paths_subset["target"] = target
            shortest_paths_diff.append(shortest_paths_subset)

        shortest_paths = pd.concat(shortest_paths_diff, axis=0)
        if slide_submitter_id is not None:
            shortest_paths["slide_submitter_id"] = slide_submitter_id
        else:
            shortest_paths["slide_submitter_id"] = np.unique(slide_data["slide_submitter_id"])[0]
        return shortest_paths
    return pd.DataFrame()


def n_clusters_per_cell_type(clusters_characterized, cell_types=None):
    """
    Determine the number of cluster per cell type when the clustering was done by using all cell types
    (simultaneous clustering)
    Args:
        clusters_characterized (DataFrame):
        cell_types (list): list of cell types

    Returns:
        (DataFrame): Dataframe containing the fraction of clusters per cell type per slide

    """
    if cell_types is None:
        cell_types = DEFAULT_CELL_TYPES

    clusters_characterized_long = clusters_characterized.melt(
        id_vars=[ "slide_submitter_id", "cluster_label"],
        value_vars=cell_types,
        var_name="cell_type",
        value_name="is_assigned",
    )

    num_clust_per_cell_type_slide = (
        clusters_characterized_long.groupby([ "slide_submitter_id", "cell_type"])
        .sum()
        .reset_index()
    )
    # Count the total number of clusters
    num_clust_per_cell_type_slide["n_clusters"] = (
        clusters_characterized_long.groupby([ "slide_submitter_id", "cell_type"])
        .count()
        .reset_index()["cluster_label"]
    )
    num_clust_per_cell_type_slide = num_clust_per_cell_type_slide.drop(
        columns=["cluster_label"]
    )
    num_clust_per_cell_type_slide["fraction"] = (
        num_clust_per_cell_type_slide["is_assigned"]
        / num_clust_per_cell_type_slide["n_clusters"]
    )
    return num_clust_per_cell_type_slide


def n_high_clusters(labeled_clusters):
    """
    Determine the fraction of clusters that are labeled as high for each cell type map, for each slide.

    Args:
        labeled_cluster (DataFrame): long dataframe containing all tiles for all slides with their cluster id and the assigned abundance label

    Returns:
        labeled_slides_grouped (DataFrame): Dataframe with slide_submitter_id and cell_type_map as ID variables, and the fraction of clusters labeled 'high' for each cell type map.
    """
    # Determine the number of clusters labeled as 'high'
    n_per_label = (
        labeled_clusters.groupby(
            [ "slide_submitter_id", "cell_type_map", "is_high"]
        )
        .count()
        .reset_index()
    )
    n_per_label = n_per_label.rename(columns={"cluster_label": "n_clusters"})
    # Determine the number of clussters per cell type map
    n_clusters = (
        n_per_label.groupby(["slide_submitter_id", "cell_type_map"])
        .sum()["n_clusters"]
        .reset_index()
    )
    n_clusters = n_clusters.rename(columns={"n_clusters": "n_total_clusters"})

    # Combine
    labeled_slides_grouped = pd.merge(
        n_per_label, n_clusters, on=[ "slide_submitter_id", "cell_type_map"]
    )
    # Determine the fraction of clusters labeled as high w.r.t the total number of clusters per cell type map (for each slide)
    labeled_slides_grouped["fraction"] = (
        labeled_slides_grouped["n_clusters"]
        / labeled_slides_grouped["n_total_clusters"]
    )
    return labeled_slides_grouped


def compute_proximity_clusters_pairs(
    tiles,
    slide_submitter_id,
    method="all",
    n_clusters=8,
    cell_types=None,
    max_dist=None,
    max_n_tiles_threshold=2,
    tile_size=512,
    overlap=50,
):
    """Determine the proximity of cluster pairs
    general method:
    1. Set a threshold of the max distance between clusters using the size of tiles
    2. Compute the minimum distances between cluster pairs
    3. Count the number of minimum distances <= threshold
    4. Compute fraction: (num distances <= threshold) / (size of largest cluster)

    _compute_proximity: determine the pairwise minimum distances between clusters
    _all: determine proximity using the clusters from SCHC with all cell type predictions (clusters labeled by cell type)
    _individual_between: determine proximity using clusters from SCHC between different cell type maps
    _individual_within: determine proximity using clusters from SCHC within each cell type map

    Args:
        tiles (DataFrame):
        slide_submitter_id (str): e.g. TCGA-77-6842-01A-01-BS1
        method (str): one of the following: all, individual_between, individual_within
        n_clusters (int)
        cell_types (list): list of cell types
        max_dist (int): max. distances between two points in two different clusters
        max_n_tiles_threshold (int): number of tiles for computing max. distance between two points in two different clusters
        tile_size (int):
        overlap (int): overlap of tiles

    Returns:
        (DataFrame): depending on the chosen method, but a dataframe containing the slide_submitter_id, cluster ids of the two compared clusters and the computed proximity.

    """

    if max_dist is None:
        tile_size = tile_size
        overlap = overlap
        actual_tile_size = tile_size - overlap
        total_X_dist = max_n_tiles_threshold * actual_tile_size
        max_dist = (total_X_dist**2 + total_X_dist**2) ** 0.5

    if method.startswith("individual") or (cell_types is None):
        cell_types = ["CAFs", "T_cells", "endothelial_cells", "tumor_purity"]

    def compute_proximity(coords1, coords2, max_dist):
        """
        Compute the proximity

        Args:
            coords1 (array/DataFrame): coordinates of tiles from a cluster A
            coords2 (array/DataFrame): coordinates of tiles from a cluster B
            max_dist (float): max. dist between two tiles from different clusters

        Returns:
            frac (float): computed proximity
        """
        if (len(coords1) > 0) & (len(coords2) > 0):

            if len(coords1) > len(coords2):
                _, dist = pairwise_distances_argmin_min(coords1, coords2)
            else:
                _, dist = pairwise_distances_argmin_min(coords2, coords1)
            total_counts = sum(dist <= max_dist)
            frac = total_counts / max([len(coords1), len(coords2)])
            return frac
        return None

    def _all(tiles, slide_submitter_id, n_clusters, max_dist):
        """
        Args:
            tiles (DataFrame): tiles for the slide(s) containing the at least the coordinates and the cluster labels (int).
            slide_submitter_id (str): e.g. TCGA-77-6842-01A-01-BS1
            n_clusters (int): number of clusters used for SCHC
            max_dist (float): max. dist between two tiles from different clusters

        Returns:
            (DataFrame): dataframe containing the slide_submitter_id,
            the cluster ids for which the proximity was computed and the actual computed proximity
        """
        tiles_sub = tiles[tiles.slide_submitter_id == slide_submitter_id]
        tiles = tiles_sub.drop_duplicates()

        # Check clusters that have multiple labels assigned
        multi_assigned = []
        for c in range(n_clusters):
            if sum(tiles.cluster_label == c) > 1:
                multi_assigned.append(c)

        cluster_pairs = list(itertools.combinations(range(8), 2)) + list(
            zip(multi_assigned, multi_assigned)
        )
        out = []
        for i, j in cluster_pairs:
            cluster1_tiles = tiles.loc[
                (tiles.slide_submitter_id == slide_submitter_id)
                & (tiles.cluster_label == i),
                ["Coord_X", "Coord_Y"],
            ]
            cluster2_tiles = tiles.loc[
                (tiles.slide_submitter_id == slide_submitter_id)
                & (tiles.cluster_label == j),
                ["Coord_X", "Coord_Y"],
            ]
            out.append(
                [
                    slide_submitter_id,
                    i,
                    j,
                    compute_proximity(cluster1_tiles, cluster2_tiles, max_dist),
                ]
            )
        return pd.DataFrame(
            out, columns=["slide_submitter_id", "cluster1", "cluster2", "proximity"]
        )

    def _individual_between(
        tiles, slide_submitter_id, n_clusters, cell_types, max_dist
    ):
        """
        Args:
            tiles (DataFrame): tiles for the slide(s) containing the at least the coordinates and the cluster labels (int).
            slide_submitter_id (str): e.g. TCGA-77-6842-01A-01-BS1
            n_clusters (int): number of clusters used for SCHC
            cell_types (list): list of cell types
            max_dist (float): max. dist between two tiles from different clusters

        Returns:
            (DataFrame): Dataframe containing the slide_submitter_id,  the cluster ids and the cell type labels for which the proximity was computed and the actual computed proximity
        """
        out = []
        cell_type_pairs = list(itertools.combinations(cell_types, 2))
        cluster_pairs = list(
            itertools.product(list(range(n_clusters)), list(range(n_clusters)))
        )
        for cell_type1, cell_type2 in cell_type_pairs:
            for i, j in cluster_pairs:
                cluster1_tiles = tiles.loc[
                    (tiles.slide_submitter_id == slide_submitter_id)
                    & (tiles[f"{cell_type1}_label"] == i),
                    ["Coord_X", "Coord_Y"],
                ]
                cluster2_tiles = tiles.loc[
                    (tiles.slide_submitter_id == slide_submitter_id)
                    & (tiles[f"{cell_type2}_label"] == j),
                    ["Coord_X", "Coord_Y"],
                ]
                out.append(
                    [
                        slide_submitter_id,
                        cell_type1,
                        i,
                        cell_type2,
                        j,
                        compute_proximity(cluster1_tiles, cluster2_tiles, max_dist),
                    ]
                )
        return pd.DataFrame(
            out,
            columns=[
                "slide_submitter_id",
                "cluster1_label",
                "cluster1",
                "cluster2_label",
                "cluster2",
                "proximity",
            ],
        )

    def _individual_within(tiles, slide_submitter_id, n_clusters, cell_types, max_dist):
        """
        Args:
            tiles (DataFrame): tiles for the slide(s) containing the at least the coordinates and the cluster labels (int).
            slide_submitter_id (str): e.g. TCGA-77-6842-01A-01-BS1
            n_clusters (int): number of clusters used for SCHC
            cell_types (list): list of cell types
            max_dist (float): max. dist between two tiles from different clusters

        Returns:
            (DataFrame): Dataframe containing the slide_submitter_id, the cluster ids for which the proximity was computed and the actual computed proximity
        """
        cluster_pairs = list(itertools.combinations(range(n_clusters), 2))
        out = []
        for cell_type in cell_types:
            for i, j in cluster_pairs:
                cluster1_tiles = tiles.loc[
                    (tiles.slide_submitter_id == slide_submitter_id)
                    & (tiles[f"{cell_type}_label"] == i),
                    ["Coord_X", "Coord_Y"],
                ]
                cluster2_tiles = tiles.loc[
                    (tiles.slide_submitter_id == slide_submitter_id)
                    & (tiles[f"{cell_type}_label"] == j),
                    ["Coord_X", "Coord_Y"],
                ]
                out.append(
                    [
                        slide_submitter_id,
                        cell_type,
                        i,
                        j,
                        compute_proximity(cluster1_tiles, cluster2_tiles, max_dist),
                    ]
                )
        return pd.DataFrame(
            out,
            columns=[
                "slide_submitter_id",
                "cell_type",
                "cluster1",
                "cluster2",
                "proximity",
            ],
        )

    if method == "all":
        return _all(
            tiles=tiles,
            slide_submitter_id=slide_submitter_id,
            n_clusters=n_clusters,
            max_dist=max_dist,
        )

    elif method == "individual_between":
        return _individual_between(
            tiles=tiles,
            slide_submitter_id=slide_submitter_id,
            n_clusters=n_clusters,
            max_dist=max_dist,
            cell_types=cell_types,
        )
    elif method == "individual_within":
        return _individual_within(
            tiles=tiles,
            slide_submitter_id=slide_submitter_id,
            n_clusters=n_clusters,
            max_dist=max_dist,
            cell_types=cell_types,
        )
    raise Exception( "Choose a valid method: 'all', 'individual_between' or 'individual_within'")

def post_processing_proximity(prox_df, slide_submitter_id, method="all"):
    """
    Select top-3 computed proximity values and take the mean

    Args:
        prox_df (DataFrame): dataframe containing the slide_submitter_id, pair and proximity for all possible cluster pairs
        slide_submitter_id (str): e.g. TCGA-77-6842-01A-01-BS1
        method (str): 'all' or 'individual_between', 'individual_within' depending on whether the clustering was done for all cell types simultaneously or one cell type at a time

    Returns:
        (DataFrame): dataframe containing  slide_submitter_id, cell type pair and computed proximity

    """
    if method == "all":
        out = []
        for pair in prox_df.pair.unique():
            prox_top3_mean = (
                prox_df.loc[
                    (prox_df.slide_submitter_id == slide_submitter_id)
                    & (prox_df.pair == pair),
                    "proximity",
                ]
                .sort_values(ascending=False)[:3]
                .mean()
            )
            out.append([ slide_submitter_id, pair, prox_top3_mean])
        return pd.DataFrame(
            out, columns=[ "slide_submitter_id", "pair", "proximity"]
        )

    elif method == "individual_between":
        out = []
        for comparison in prox_df.comparison.unique():
            for pair in prox_df.ordered_pair.unique():
                prox_top3_mean = (
                    prox_df.loc[
                        (prox_df.slide_submitter_id == slide_submitter_id)
                        & (prox_df.comparison == comparison)
                        & (prox_df.ordered_pair == pair),
                        "proximity",
                    ]
                    .sort_values(ascending=False)[:3]
                    .mean()
                )
                out.append(
                    [ slide_submitter_id, comparison, pair, prox_top3_mean])

        return pd.DataFrame(
            out,
            columns=[ "slide_submitter_id", "comparison", "pair", "proximity"],
        )
    elif method == "individual_within":
        out = []

        for comparison in prox_df.comparison.unique():
            for pair in prox_df.pair.unique():
                prox_top3_mean = (
                    prox_df.loc[
                        (prox_df.slide_submitter_id == slide_submitter_id)
                        & (prox_df.comparison == comparison)
                        & (prox_df.pair == pair),
                        "proximity",
                    ]
                    .sort_values(ascending=False)[:3]
                    .mean()
                )
                out.append(
                    [ slide_submitter_id, comparison, pair, prox_top3_mean]
                )

        return pd.DataFrame(
            out,
            columns=[ "slide_submitter_id", "comparison", "pair", "proximity"],
        )


def compute_shape_features(
    tiles,
    slide_submitter_id,
    tile_size=512,
    overlap=50,
    cell_types=None,
    method="all",
):
    """
    Compute solidity and roundness for each cluster using a convexhull

    Args:
        tiles (DataFrame): at least containing the slide submitter id, cluster label (int) and coordinates
        slide_submitter_id (str)
        method (str): 'individual' or 'all' dependent on whether the clustering was done using one cell type at a tie or using all cell types simultaneously
        tile_size (int)
        overlap (int): overlap of tiles
        cell_types (list): list of cell types only required for 'method=individual'

    Returns:
        (DataFrame): dataframe of slide_submitter_id, cluster_label and the computed solidity and roundness
    """
    if cell_types is None:
        cell_types = DEFAULT_CELL_TYPES
    slide_data = tiles[tiles.slide_submitter_id == slide_submitter_id]
    if (method == "all"):
        cluster_columns = [f"is_{cell_type}_cluster" for cell_type in cell_types]
        out = []
        # Find exclusive slides
        subset = slide_data[slide_data[cluster_columns].sum(axis=1) == 1]
        if len(subset) > 0:
            for cluster_label in slide_data.cluster_label.unique():
                slide_cluster = slide_data.loc[
                    (slide_data.cluster_label == cluster_label), ["Coord_X", "Coord_Y"]
                ].to_numpy()
                try:
                    real_tile_size = tile_size - overlap
                    ch = ConvexHull(slide_cluster)
                    area = len(slide_cluster) * (real_tile_size**2)
                    convex_area = ch.volume
                    convex_perimeter = ch.area
                    solidity = area / convex_area
                    roundness = (4 * area * math.pi) / (convex_perimeter**2)

                    # If ConvexHull is smaller than actual shape, fix at 1.
                    solidity = (solidity * (solidity <= 1)) or 1
                    roundness = (roundness * (roundness <= 1)) or 1

                    for cell_type in pd.Series(cell_types)[
                        subset.loc[subset.cluster_label == cluster_label, cluster_columns]
                        .iloc[0]
                        .tolist()
                    ]:
                        out.append(
                            [
                                slide_submitter_id,
                                cluster_label,
                                cell_type,
                                solidity,
                                roundness,
                            ]
                        )
                # trunk-ignore(flake8/E722)
                except:
                    print(f"Failed to create a convex hull for cluster={cluster_label} of slide {slide_submitter_id}")
    elif(method == "individual"):
        out = []
        for c in cell_types:
                cluster_col = f"{c}_label"
                for cluster_label in slide_data[cluster_col].unique():
                        slide_cluster = slide_data.loc[
                                (slide_data[cluster_col] == cluster_label), ["Coord_X", "Coord_Y"]
                        ].to_numpy()
                        try:
                                real_tile_size = tile_size - overlap
                                ch = ConvexHull(slide_cluster)
                                area = len(slide_cluster) * (real_tile_size**2)
                                convex_area = ch.volume
                                convex_perimeter = ch.area
                                solidity = area / convex_area
                                roundness = (4 * area * math.pi) / (convex_perimeter**2)

                                # If ConvexHull is smaller than actual shape, fix at 1.
                                solidity = (solidity * (solidity <= 1)) or 1
                                roundness = (roundness * (roundness <= 1)) or 1

                                out.append(
                                        [
                                        slide_submitter_id,
                                        cluster_label,
                                        c,
                                        solidity,
                                        roundness,
                                        ]
                                )
                        # trunk-ignore(flake8/E722)
                        except:
                                print(f"Failed to create a convex hull for cluster={cluster_label} of slide {slide_submitter_id}")

    if len(out) > 0:
        return pd.DataFrame(
            out,
            columns=[
                "slide_submitter_id",
                "cluster_label",
                "cell_type",
                "solidity",
                "roundness",
            ],
        )
    return pd.DataFrame(
        columns=[
            "slide_submitter_id",
            "cluster_label",
            "cell_type",
            "solidity",
            "roundness",
        ]
    )


def compute_node_degree(G, cell_type_assignments=None, cell_types=None):
    """
    Compute the node degree using the following approach:
    1. Determining the non-border nodes (nodes that have 8 neighbors)
    2. Number of neighbors per cell type for each source (center) node (regardless of the cell type of that node)
    3. Determine average degree w.r.t. source node cell type (nodes can be multiple cell types)

    Args:
        graph (Networkx Graph): graph representing the slide constructed with Networkx
        cell_type_assignments (DataFrame): Dataframe containing the cell type labels of the individual tiles indicated with booleans based on P > threshold
        cell_types (list): list of cell types

    Returns:
        (DataFrame): long dataframe containing the computed node degree for the center-neighbor pair
    """
    if cell_types is None:
        cell_types = DEFAULT_CELL_TYPES

    # 1. Determining the non-border nodes (nodes that have 8 neighbors)
    neighbors_per_node = [list(G.neighbors(node)) for node in G.nodes()]
    num_neighbors = pd.Series(
        [len(neighbor_set) for neighbor_set in neighbors_per_node], index=G.nodes()
    )
    source_nodes = num_neighbors.index[num_neighbors == 8]

    # Store
    neighbors = dict.fromkeys(source_nodes)
    node_degree = []

    # 2. Number of neighbors per cell type for each source (center) node (regardless of the cell type of that node)
    for node in source_nodes:
        num_neighbors_per_type = dict.fromkeys(cell_types)
        # a. Determine the neighbors (return indices) of the node (regardless of the center node's cell type as it could be multiple)
        # b. Check to which cell types those neighbors correspond
        # c. Add the number of neighbors per cell type = degree
        for cell_type in cell_types:
            num_neighbors_per_type[cell_type] = len(
                set(G.neighbors(node))
                & set(cell_type_assignments.index[cell_type_assignments[cell_type]])
            )
        neighbors[node] = num_neighbors_per_type
    neighbors_df = pd.DataFrame(neighbors).T

    # 3. Determine average degree w.r.t. source node cell type (nodes can be multiple cell types)
    for cell_type in cell_types:
        # Check if there are nodes/tiles that have the cell type
        node_ids = list(
            set(cell_type_assignments.index[cell_type_assignments[cell_type]])
            & set(source_nodes)
        )
        if len(node_ids) > 0:
            neighbors_subset = neighbors_df.loc[node_ids].copy()
            neighbors_subset.loc[:, "center"] = cell_type
            node_degree.append(neighbors_subset)

    # Transform in long format
    if len(node_degree) > 0:
        return pd.melt(
            pd.concat(node_degree, axis=0),
            id_vars=["center"],
            var_name="neighbor",
            value_name="degree",
        )
    else:
        return []


def simulate_assigning_cell_type_positions(
    G, cell_type_assignments, cell_types=None, num_repeats=100, seed=123
):
    """
    Creating simulated slides from true slides by randomly distributing the assigned tiles"

    Args:
        graph (Networkx Graph): graph representing the slide constructed with Networkx
        cell_type_assignments (DataFrame): Dataframe containing the cell type labels of the individual tiles indicated with booleans based on P > threshold
        cell_types (list): list of cell types
        num_repeats (int): number of simulations to create for each slide
        seed (int): for reproducibility

    Returns:
        (DataFrame): dataframe containing the computed node degrees for the different center-neighbor pairs for the number of repeats.
    """
    np.random.seed(seed=seed)
    if cell_types is None:
        cell_types =DEFAULT_CELL_TYPES

    # Setup
    total_per_type = cell_type_assignments[cell_types].sum(axis=0)
    num_tiles = len(G.nodes())
    node_degree_repeats = []
    for repeat in range(num_repeats):
        assigned_tiles = dict.fromkeys(cell_types)
        # Randomize the location of the cell type, total number of cell type tiles is fixed.
        for cell_type in cell_types:
            N_true = int(total_per_type[cell_type])
            N_false = num_tiles - N_true
            pool = N_true * [True] + N_false * [False]
            assigned_tiles[cell_type] = np.random.choice(
                pool,
                replace=False,
                size=num_tiles,
            )
        assigned_tiles = pd.DataFrame(assigned_tiles)
        # Compute the node degree for the different pair for the simulated slide
        node_degree = compute_node_degree(G, assigned_tiles, cell_types)
        if len(node_degree) > 0:
            node_degree["simulation_nr"] = repeat
            node_degree_repeats.append(node_degree)

    if len(node_degree_repeats) > 0:
        return pd.concat(node_degree_repeats)
    else:
        return []


def compute_effect_size(true_mean_nd, sims_nd_df, slide_submitter_id, cell_types=None):
    """
    Compute Cohen's d and do one-sample t-test based on computed effect size

    Args:
        true_mean_nd (DataFrame): computed mean node degree for each cell type pair
        sims_nd_df (DataFrame): computed node degree for each cell type pair for all simulations
        slide_submitter_id (str): e.g. TCGA-77-6842-01A-01-BS1
        cell_types (list): list of cell types

    Returns:
        (DataFrame): dataframe containing the computed effect size, T-statistic and p-value for the one-sample t-test for each center-neighbor pair for a slide
    """
    if cell_types is None:
        cell_types = DEFAULT_CELL_TYPES

    computed_effect_sizes = []
    center_neighbor_pairs = sorted(
        list(itertools.permutations(cell_types, 2)) + list(zip(cell_types, cell_types))
    )

    for center, neighbor in center_neighbor_pairs:
        # Select center-neighbor pair for a slide
        mask1 = (
            (sims_nd_df.slide_submitter_id == slide_submitter_id)
            & (sims_nd_df.center == center)
            & (sims_nd_df.neighbor == neighbor)
        )
        mask2 = (
            (true_mean_nd.slide_submitter_id == slide_submitter_id)
            & (true_mean_nd.center == center)
            & (true_mean_nd.neighbor == neighbor)
        )
        # Only if the slide has this pair, continue
        if (sum(mask1) > 0) & (sum(mask2) > 0):
            # Compute Cohen's d
            sample_sd = statpy.stdev(sims_nd_df.loc[mask1, "degree"])
            sample_mean = sims_nd_df.loc[mask1, "degree"].mean()
            theoretical_mean = true_mean_nd.loc[mask2, "mean_obs"].values[0]
            effect_size = utils.cohens_d(theoretical_mean, sample_mean, sample_sd)
            # Perform one-sample t-test with H1 based on effect size
            if effect_size > 0:
                Tstat, pval = stats.ttest_1samp(
                    sims_nd_df.loc[mask1, "degree"],
                    theoretical_mean,
                    alternative="greater",
                )
            else:
                Tstat, pval = stats.ttest_1samp(
                    sims_nd_df.loc[mask1, "degree"],
                    theoretical_mean,
                    alternative="less",
                )

            computed_effect_sizes.append(
                [slide_submitter_id, center, neighbor, effect_size, Tstat, pval]
            )
    return pd.DataFrame(
        computed_effect_sizes,
        columns=[
            "slide_submitter_id",
            "center",
            "neighbor",
            "effect_size",
            "Tstat",
            "pval",
        ],
    )


def node_degree_wrapper(graph, predictions, slide_submitter_id, cell_types=None, keep_example_simulations=True):
    """
    Compute Effect size based on difference between true node degree and node degree computed from simulated slides.

    NOTE: For documentation of used functions, please see the separate documentation.

    Args:
        graph (Networkx Graph): graph representing the slide constructed with Networkx
        predictions (DataFrame): dataframe containing the predicted abundances for the cell types
        slide_submitter_id (str): e.g. TCGA-77-6842-01A-01-BS1
        cell_types (list): list of cell types

    Returns:
        assigned_tiles (dict): dict with single key 'slide_submitter_id' and with a DataFrame as value containing a single simulation of a slide
        (DataFrame): dataframe containing the computed node degrees for each simulation for all cell type pairs
        mean_nd_df (DataFrame): dataframe containing the mean node degree across the N simulations and the true mean node degree for all ell type pairs

    """
    if cell_types is None:
        cell_types = DEFAULT_CELL_TYPES

    slide_data =utils.get_slide_data(predictions, slide_submitter_id)
    cell_type_assignments = utils.assign_cell_types(slide_data, cell_types=cell_types)
    center_neighbor_pairs = sorted(list(itertools.permutations(cell_types, 2)) + list(zip(cell_types, cell_types)))
    mean_nd = []

    # Compute true node degree
    node_degree_original = compute_node_degree(graph, cell_type_assignments, cell_types=cell_types)
    if (len(node_degree_original) > 0):
        # Simulate slides and compute their node degrees
        node_degree_sim = simulate_assigning_cell_type_positions(graph, cell_type_assignments, cell_types=cell_types)

        # Add slide submitter id
        node_degree_sim["slide_submitter_id"] = slide_submitter_id
        node_degree_original["slide_submitter_id"] = slide_submitter_id

        # Keep one simulation for visualization
        if keep_example_simulations:
            total_per_type = cell_type_assignments[cell_types].sum(axis=0)
            num_tiles = len(graph.nodes())
            assigned_tiles = dict.fromkeys(cell_types)
            for cell_type in cell_types:
                N_true = int(total_per_type[cell_type])
                N_false = num_tiles - N_true
                pool = N_true * [True] +  N_false * [False]
                assigned_tiles[cell_type] = np.random.choice(pool, replace=False, size=num_tiles)

            assigned_tiles = pd.DataFrame(assigned_tiles)
            assigned_tiles["slide_submitter_id"] = slide_submitter_id

        # Comparing node degree means
        node_degree_original_mean = node_degree_original.groupby(["center", "neighbor"]).mean(numeric_only=True).reset_index()
        node_degree_sim_mean = node_degree_sim.groupby(["center", "neighbor", "simulation_nr"]).mean(numeric_only=True)["degree"].reset_index()

        # Subset the mean node degree per cell type pair
        for center, neighbor in list(center_neighbor_pairs):
            pair_original = node_degree_original_mean.loc[(node_degree_original_mean["center"] == center) & (node_degree_original_mean["neighbor"] == neighbor), "degree"].to_list()
            pair_sim = node_degree_sim_mean.loc[(node_degree_sim_mean["center"] == center) & (node_degree_sim_mean["neighbor"] == neighbor), "degree"]
            if (len(pair_sim) > 0) & (len(pair_original)> 0):
                # Comparing mean node degrees from simulations (100 values) with the observed mean node degree (single value)
                # Compute fraction of times the mean node degree from simulations is greater
                mean_nd.append([slide_submitter_id, center, neighbor,np.mean(pair_sim), np.mean(pair_original)])


        mean_nd_df = pd.DataFrame(mean_nd, columns=["slide_submitter_id", "center", "neighbor",  "mean_sim", "mean_obs"])
        node_degree_sim_df = node_degree_sim.groupby(["slide_submitter_id", "center", "neighbor", "simulation_nr"])["degree"].mean().reset_index(name="degree")

        if keep_example_simulations:
            return [assigned_tiles, node_degree_sim_df, mean_nd_df]
        else:
            return [node_degree_sim_df, mean_nd_df]
