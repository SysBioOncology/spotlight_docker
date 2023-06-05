import os
import sys
import joblib
import pandas as pd
from joblib import Parallel, delayed
import argparse
from os import path

sys.path.append(f"{os.path.dirname(os.getcwd())}/Python/libs")

# Own modules
import features.clustering as clustering
import features.features as features
import features.graphs as graphs
import features.utils as utils
from model.constants import DEFAULT_SLIDE_TYPE, DEFAULT_CELL_TYPES, NUM_CORES, METADATA_COLS

def compute_network_features(tile_quantification_path, output_dir, slide_type=DEFAULT_SLIDE_TYPE, cell_types=None, graphs_path=None,
abundance_threshold=0.5, shapiro_alpha=0.05, cutoff_path_length=2):
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

    if graphs_path is None:
        results = Parallel(n_jobs=NUM_CORES)(
            delayed(graphs.construct_graph)(predictions=predictions, slide_submitter_id=id)
            for id in slide_submitter_ids
        )
        # Extract/format graphs
        all_graphs = {
            list(slide_graph.keys())[0]: list(slide_graph.values())[0]
            for slide_graph in results
        }
        joblib.dump(all_graphs, f"{output_dir}/{slide_type}_graphs.pkl")
    else:
        all_graphs = joblib.load(graphs_path)

    #######################################################
    # ---- Compute connectedness and co-localization ---- #
    #######################################################

    all_largest_cc_sizes = []
    all_dual_nodes_frac = []
    for id in slide_submitter_ids:
        slide_data = utils.get_slide_data(predictions, id)
        node_cell_types = utils.assign_cell_types(slide_data=slide_data, cell_types=cell_types, threshold=abundance_threshold)
        lcc = features.determine_lcc(
            graph=all_graphs[id], cell_type_assignments=node_cell_types, cell_types=cell_types
        )
        lcc["slide_submitter_id"] = id
        all_largest_cc_sizes.append(lcc)

        dual_nodes_frac = features.compute_dual_node_fractions(node_cell_types, cell_types)
        dual_nodes_frac["slide_submitter_id"] = id
        all_dual_nodes_frac.append(dual_nodes_frac)

    all_largest_cc_sizes = pd.concat(all_largest_cc_sizes, axis=0)
    all_dual_nodes_frac = pd.concat(all_dual_nodes_frac, axis=0)

    #######################################################
    # ---- Compute N shortest paths with max. length ---- #
    #######################################################

    results = Parallel(n_jobs=NUM_CORES)(
        delayed(features.compute_n_shortest_paths_max_length)(
            predictions=predictions, slide_submitter_id=id, graph=all_graphs[id], cutoff=cutoff_path_length
        )
        for id in slide_submitter_ids
    )
    # Formatting and count the number of shortest paths with max length
Ã    all_shortest_paths_thresholded = pd.concat(results, axis=0)
    all_shortest_paths_thresholded["n_paths"] = 1
    proximity_graphs = (
        all_shortest_paths_thresholded.groupby(
            ["slide_submitter_id", "source", "target"]
        )
        .sum(numeric_only=True)
        .reset_index()
    )
    # Post-processing
    proximity_graphs["pair"] = [f"{source}-{target}" for source, target in proximity_graphs[["source", "target"]].to_numpy()]
    proximity_graphs = proximity_graphs.drop(columns=["path_length"])

    # ---- Formatting ---- #
    all_largest_cc_sizes = all_largest_cc_sizes.reset_index(drop=True)
    all_largest_cc_sizes_wide = all_largest_cc_sizes.pivot(index=["slide_submitter_id"], columns="cell_type")["type_spec_frac"]
    new_cols = [f'LCC {col.replace("_", " ")} clusters' for col in all_largest_cc_sizes_wide.columns]
    all_largest_cc_sizes_wide.columns = new_cols
    all_largest_cc_sizes_wide = all_largest_cc_sizes_wide.reset_index()

    shortest_paths_wide = proximity_graphs.pivot(index=["slide_submitter_id"], columns="pair")["n_paths"]
    new_cols = [f'Prox graph {col.replace("_", " ")} clusters' for col in shortest_paths_wide.columns]
    shortest_paths_wide.columns = new_cols
    shortest_paths_wide = shortest_paths_wide.reset_index()

    colocalization_wide = all_dual_nodes_frac.pivot(index=["slide_submitter_id"], columns="pair")["frac"]
    new_cols = [f'Co-loc {col.replace("_", " ")} clusters' for col in colocalization_wide.columns]
    colocalization_wide.columns = new_cols
    colocalization_wide = colocalization_wide.reset_index()

    all_features = pd.merge(all_largest_cc_sizes_wide, shortest_paths_wide)
    all_features = pd.merge(all_features, colocalization_wide)

    # ---- Save to file ---- #
    all_largest_cc_sizes_wide.to_csv(f"{output_dir}/{slide_type}_features_lcc_fraction.csv", sep="\t", index=False)
    proximity_graphs.to_csv(f"{output_dir}/{slide_type}_features_shortest_paths_thresholded.csv", sep="\t", index=False)
    all_dual_nodes_frac.to_csv(f"{output_dir}/{slide_type}_features_coloc_fraction.csv", sep="\t", index=False)
    all_features.to_csv(f"{output_dir}/{slide_type}_all_graph_features.csv", sep="\t", index=False)

    ###############################################
    # ---- Compute ES based on ND difference ---- #
    ###############################################
    # Remove one slide for which node degree could not be resolved (no node with 8 neighbours)
    #problematic_slide = 'TCGA-D3-A2JE-06A-01-TS1' # just 63 tiles
    #filtered_slides = list(filter(lambda id: id != problematic_slide, slide_submitter_ids))
    nd_results = Parallel(n_jobs=NUM_CORES)(delayed(features.node_degree_wrapper)(all_graphs[id], predictions, id) for id in slide_submitter_ids)
    nd_results = list(filter(lambda id: id != None, nd_results))

    # Format results
    all_sims_nd = []
    all_mean_nd_df = []
    example_simulations = {}

    for sim_assignments, sim, mean_nd_df in nd_results:
        all_mean_nd_df.append(mean_nd_df)
        all_sims_nd.append(sim)
        example_simulations.update(sim_assignments)

    all_sims_nd = pd.concat(all_sims_nd, axis=0).reset_index()
    all_mean_nd_df =pd.concat(all_mean_nd_df).reset_index(drop=True)

    # Testing normality
    shapiro_tests = Parallel(n_jobs=NUM_CORES)(delayed(utils.test_normality)(sims_nd=all_sims_nd, slide_submitter_id=id, alpha=shapiro_alpha, cell_types=cell_types) for id in all_sims_nd.slide_submitter_id.unique())
    all_shapiro_tests = pd.concat(shapiro_tests, axis=0)

    # Computing Cohen's d effect size and perform t-test
    effect_sizes = Parallel(n_jobs=NUM_CORES)(delayed(features.compute_effect_size)(all_mean_nd_df, all_sims_nd, slide_submitter_id) for slide_submitter_id in all_sims_nd.slide_submitter_id.unique())
    all_effect_sizes = pd.concat(effect_sizes, axis=0)
    all_effect_sizes["pair"] = [f"{c}-{n}" for c, n in all_effect_sizes[["center", "neighbor"]].to_numpy()]

    # ---- Save to file ---- #
    all_effect_sizes.to_csv(
        f"{output_dir}/{slide_type}_features_ND_ES.csv", sep="\t", index=False)
    all_sims_nd.to_csv(
        f"{output_dir}/{slide_type}_features_ND_sims.csv", sep="\t", index=False)
    all_mean_nd_df.to_csv(
        f"{output_dir}/{slide_type}_features_ND.csv", sep="\t", index=False)
    joblib.dump(example_simulations,
        f"{output_dir}/{slide_type}_features_ND_sim_assignments.pkl")
    all_shapiro_tests.to_csv(f"{output_dir}/{slide_type}_shapiro_tests.csv", index=False, sep="\t")

def compute_clustering_features(
    tile_quantification_path, output_dir, slide_type=DEFAULT_SLIDE_TYPE, cell_types=None, graphs_path=None, n_clusters=8, max_dist=None, max_n_tiles_threshold=2, tile_size=512, overlap=50):
    output_dir = f"{output_dir}/4_spatial_features"
    if not path.exists(output_dir):
        os.mkdir(output_dir)
    if cell_types is None:
        cell_types = DEFAULT_CELL_TYPES

    predictions = pd.read_csv(tile_quantification_path, sep="\t")
    predictions = predictions[cell_types + METADATA_COLS]

    slide_submitter_ids = list(set(predictions.slide_submitter_id))

    #####################################
    # ---- Constructing the graphs ---- #
    #####################################

    if graphs_path is None:
        results = Parallel(n_jobs=NUM_CORES)(
            delayed(graphs.construct_graph)(predictions=predictions, slide_submitter_id=id)
            for id in slide_submitter_ids
        )
        # Extract/format graphs
        all_graphs = {
            list(slide_graph.keys())[0]: list(slide_graph.values())[0]
            for slide_graph in results
        }
        joblib.dump(
            all_graphs, f"{output_dir}/{slide_type}_graphs.pkl")
    else:
        all_graphs = joblib.load(graphs_path)

    ######################################################################
    # ---- Fraction of cell type clusters (simultaneous clustering) ---- #
    ######################################################################

    # Spatially Hierarchical Constrained Clustering with all quantification of all cell types
    slide_clusters= Parallel(n_jobs=NUM_CORES)(delayed(clustering.schc_all)(predictions, all_graphs[id], id) for id in slide_submitter_ids)
    # Combine the tiles labeled with their cluster id for all slides
    tiles_all_schc = pd.concat(slide_clusters, axis=0)

    # Assign a cell type label based on the mean of all cluster means across all slides
    all_slide_clusters_characterized = clustering.characterize_clusters(tiles_all_schc)

    # Count the number of clusters per cell type for each slide
    num_clust_by_slide = features.n_clusters_per_cell_type(all_slide_clusters_characterized, cell_types=cell_types)

    ######################################################################################
    # ---- Fraction of highly abundant cell types (individual cell type clustering) ---- #
    ######################################################################################

    # Spatially Hierarchical Constrained Clustering with all quantification of all cell types for each individual cell type
    slide_indiv_clusters= Parallel(n_jobs=NUM_CORES)(delayed(clustering.schc_individual)(predictions, all_graphs[id], id) for id in slide_submitter_ids)
    all_slide_indiv_clusters = pd.concat(slide_indiv_clusters, axis=0)

    # Add metadata
    all_slide_indiv_clusters = pd.merge(predictions, all_slide_indiv_clusters, on="tile_ID")

    # Add abundance label 'high' or 'low' based on cluster means
    slide_indiv_clusters_labeled = clustering.label_cell_type_map_clusters(all_slide_indiv_clusters)

    # Count the fraction of 'high' clusters
    frac_high = features.n_high_clusters(slide_indiv_clusters_labeled)

    ##################################################################
    # ---- Compute proximity features (simultaneous clustering) ---- #
    ##################################################################

    ## Computing proximity for clusters derived with all cell types simultaneously
    clusters_all_schc_long = all_slide_clusters_characterized.melt(id_vars=["slide_submitter_id", "cluster_label"], value_name="is_assigned", var_name="cell_type")
    # remove all cell types that are not assigned to the cluster
    clusters_all_schc_long = clusters_all_schc_long[clusters_all_schc_long["is_assigned"]]
    clusters_all_schc_long = clusters_all_schc_long.drop(columns="is_assigned")

    results_schc_all= Parallel(n_jobs=NUM_CORES)(delayed(features.compute_proximity_clusters_pairs)(
        tiles=tiles_all_schc, slide_submitter_id=id, n_clusters=n_clusters, cell_types=cell_types, max_dist=max_dist, max_n_tiles_threshold=max_n_tiles_threshold, tile_size=tile_size, overlap=overlap, method="all") for id in slide_submitter_ids)
    prox_all_schc = pd.concat(results_schc_all)

    # Label clusters (a number) with the assigned cell types
    prox_all_schc = pd.merge(prox_all_schc, clusters_all_schc_long, left_on=["slide_submitter_id", "cluster1"], right_on=["slide_submitter_id", "cluster_label"])
    prox_all_schc = prox_all_schc.rename(columns={"cell_type": "cluster1_label"})
    prox_all_schc = prox_all_schc.drop(columns=["cluster_label"])

    prox_all_schc = pd.merge(prox_all_schc, clusters_all_schc_long, left_on=["slide_submitter_id", "cluster2"], right_on=["slide_submitter_id", "cluster_label"])
    prox_all_schc = prox_all_schc.rename(columns={"cell_type": "cluster2_label"})

    # Order doesn't matter: x <->
    prox_all_schc["pair"] = [f"{sorted([i, j])[0]}-{sorted([i, j])[1]}" for i, j in prox_all_schc[["cluster1_label", "cluster2_label"]].to_numpy()]
    prox_all_schc = prox_all_schc[((prox_all_schc.cluster1 == prox_all_schc.cluster2) & (prox_all_schc.cluster2_label != prox_all_schc.cluster1_label)) | (prox_all_schc.cluster1 != prox_all_schc.cluster2)]

    # slides = prox_all_schc[["MFP", "slide_submitter_id"]].drop_duplicates().to_numpy()
    slide_submitter_ids = list(set(prox_all_schc.slide_submitter_id))

    # Post Processing
    results_schc_all= Parallel(n_jobs=NUM_CORES)(delayed(features.post_processing_proximity)(prox_df=prox_all_schc, slide_submitter_id=id, method="all") for id in slide_submitter_ids)
    all_prox_df = pd.concat(results_schc_all)
    # Remove rows with a proximity of NaN
    all_prox_df = all_prox_df.dropna(axis=0)

    ##########################################################################
    # ---- Compute proximity features (individual cell type clustering) ---- #
    ##########################################################################

    ## Computing proximity for clusters derived for each cell type individually
    # Between clusters
    slide_submitter_ids = list(set(predictions.slide_submitter_id))
    results_schc_indiv= Parallel(n_jobs=NUM_CORES)(delayed(features.compute_proximity_clusters_pairs)(all_slide_indiv_clusters, slide_submitter_id=id, method="individual_between",n_clusters=n_clusters, cell_types=cell_types, max_dist=max_dist, max_n_tiles_threshold=max_n_tiles_threshold, tile_size=tile_size, overlap=overlap) for id in slide_submitter_ids)
    prox_indiv_schc = pd.concat(results_schc_indiv)

    # Formatting
    prox_indiv_schc = pd.merge(prox_indiv_schc,slide_indiv_clusters_labeled, left_on=["slide_submitter_id", "cluster1_label", "cluster1"], right_on=["slide_submitter_id", "cell_type_map", "cluster_label"])
    prox_indiv_schc = prox_indiv_schc.drop(columns=["cell_type_map", "cluster_label"])
    prox_indiv_schc = prox_indiv_schc.rename(columns={"is_high": "cluster1_is_high"})
    prox_indiv_schc = pd.merge(prox_indiv_schc,slide_indiv_clusters_labeled, left_on=["slide_submitter_id", "cluster2_label", "cluster2"], right_on=["slide_submitter_id", "cell_type_map", "cluster_label"])
    prox_indiv_schc = prox_indiv_schc.rename(columns={"is_high": "cluster2_is_high"})
    prox_indiv_schc = prox_indiv_schc.drop(columns=["cell_type_map", "cluster_label"])

    # Order matters
    prox_indiv_schc["ordered_pair"] = [f"{i}-{j}" for i, j in prox_indiv_schc[["cluster1_label", "cluster2_label"]].to_numpy()]
    prox_indiv_schc["comparison"] = [f"cluster1={i}-cluster2={j}" for i, j in prox_indiv_schc[["cluster1_is_high", "cluster2_is_high"]].to_numpy()]

    # Post-processing
    slide_submitter_ids = list(set(predictions.slide_submitter_id))
    results_schc_indiv= pd.concat(Parallel(n_jobs=NUM_CORES)(delayed(features.post_processing_proximity)(prox_df=prox_indiv_schc, slide_submitter_id=id, method="individual_between") for id in slide_submitter_ids))

    # Within clusters
    slide_submitter_ids = list(set(predictions.slide_submitter_id))
    results_schc_indiv_within= Parallel(n_jobs=NUM_CORES)(delayed(features.compute_proximity_clusters_pairs)(all_slide_indiv_clusters, slide_submitter_id=id, method="individual_within",n_clusters=n_clusters, cell_types=cell_types, max_dist=max_dist, max_n_tiles_threshold=max_n_tiles_threshold, tile_size=tile_size, overlap=overlap,) for id in slide_submitter_ids)
    prox_indiv_schc_within = pd.concat(results_schc_indiv_within)

    prox_indiv_schc_within = pd.merge(prox_indiv_schc_within,slide_indiv_clusters_labeled, left_on=["slide_submitter_id", "cell_type", "cluster1"], right_on=["slide_submitter_id", "cell_type_map", "cluster_label"])
    prox_indiv_schc_within = prox_indiv_schc_within.drop(columns=[ "cluster_label"])
    prox_indiv_schc_within = prox_indiv_schc_within.rename(columns={"is_high": "cluster1_is_high", "cell_type_map":"cell_type_map1"})
    prox_indiv_schc_within = pd.merge(prox_indiv_schc_within,slide_indiv_clusters_labeled, left_on=["slide_submitter_id", "cell_type", "cluster2"], right_on=["slide_submitter_id", "cell_type_map", "cluster_label"])
    prox_indiv_schc_within = prox_indiv_schc_within.rename(columns={"is_high": "cluster2_is_high", "cell_type_map": "cell_type_map2"})
    prox_indiv_schc_within = prox_indiv_schc_within.drop(columns=["cluster_label"])

    # Order doesn't matter (only same cell type combinations)
    prox_indiv_schc_within["pair"] = [f"{i}-{j}" for i, j in prox_indiv_schc_within[["cell_type_map1", "cell_type_map2"]].to_numpy()]
    prox_indiv_schc_within["comparison"] = [f"cluster1={sorted([i,j])[0]}-cluster2={sorted([i,j])[1]}" for i, j in prox_indiv_schc_within[["cluster1_is_high", "cluster2_is_high"]].to_numpy()]

    # Post-processing
    slide_submitter_ids = list(set(prox_indiv_schc_within.slide_submitter_id))
    results_schc_indiv_within= pd.concat(Parallel(n_jobs=NUM_CORES)(delayed(features.post_processing_proximity)(prox_df=prox_indiv_schc_within, slide_submitter_id=id,method="individual_within") for id in slide_submitter_ids))

    # Concatenate within and between computed proximity values
    prox_indiv_schc_combined = pd.concat([results_schc_indiv_within, results_schc_indiv])

    # Remove rows with a proximity of NaN
    prox_indiv_schc_combined = prox_indiv_schc_combined.dropna(axis=0)

    ####################################
    # ---- Compute shape features ---- #
    ####################################

    # Compute shape features based on clustering with all cell types simultaneously
    slide_submitter_ids = list(set(predictions.slide_submitter_id))
    all_slide_clusters_characterized = all_slide_clusters_characterized.rename(columns=dict(zip(cell_types, [f"is_{cell_type}_cluster" for cell_type in cell_types])))
    tiles_all_schc = pd.merge(tiles_all_schc, all_slide_clusters_characterized, on=["slide_submitter_id", "cluster_label"])

    res = pd.concat(Parallel(n_jobs=NUM_CORES)(delayed(features.compute_shape_features)(tiles=tiles_all_schc, slide_submitter_id=id, tile_size=tile_size, overlap=overlap, cell_types=cell_types) for id in slide_submitter_ids))
    res = res.drop(axis=1, labels=["cluster_label"])
    shape_feature_means = res.groupby(["slide_submitter_id", "cell_type"]).mean().reset_index()

    ##############################################
    # ---- Formatting all computed features ---- #
    ##############################################

    frac_high_sub = frac_high[frac_high["is_high"]].copy()
    frac_high_sub = frac_high_sub.drop(columns=["is_high", "n_clusters", "n_total_clusters"])

    frac_high_wide = frac_high_sub.pivot(index=["slide_submitter_id"], columns=["cell_type_map"])["fraction"]
    new_cols=[('fraction {0} clusters labeled high'.format(col)) for col in frac_high_wide.columns]
    frac_high_wide.columns = new_cols
    frac_high_wide = frac_high_wide.sort_index(axis="columns").reset_index()

    num_clust_by_slide_sub = num_clust_by_slide.copy()
    num_clust_by_slide_sub = num_clust_by_slide_sub.drop(columns=["is_assigned", "n_clusters"])

    num_clust_slide_wide = num_clust_by_slide_sub.pivot(index=[ "slide_submitter_id"], columns=["cell_type"])["fraction"]
    new_cols=[('fraction {0} clusters'.format(col)) for col in num_clust_slide_wide.columns]
    num_clust_slide_wide.columns = new_cols
    num_clust_slide_wide = num_clust_slide_wide.sort_index(axis="columns").reset_index()

    all_prox_df_wide = all_prox_df.pivot(index=["slide_submitter_id"], columns=["pair"])["proximity"]
    new_cols = [f'prox CC {col.replace("_", " ")} clusters' for col in all_prox_df_wide.columns]
    all_prox_df_wide.columns = new_cols
    all_prox_df_wide = all_prox_df_wide.reset_index()

    prox_indiv_schc_combined.comparison = prox_indiv_schc_combined.comparison.replace(dict(zip(['cluster1=True-cluster2=True', 'cluster1=True-cluster2=False',
        'cluster1=False-cluster2=True', 'cluster1=False-cluster2=False'], ["high-high", "high-low", "low-high", "low-low"])))
    prox_indiv_schc_combined["pair (comparison)"] = [f"{pair} ({comp})" for pair, comp in prox_indiv_schc_combined[["pair", "comparison"]].to_numpy()]
    prox_indiv_schc_combined = prox_indiv_schc_combined.drop(axis=1, labels=["pair", "comparison"])
    prox_indiv_schc_combined_wide = prox_indiv_schc_combined.pivot(index=[ "slide_submitter_id"], columns=["pair (comparison)"])["proximity"]
    new_cols = [f'prox CC {col.replace("_", " ")}' for col in prox_indiv_schc_combined_wide.columns]
    prox_indiv_schc_combined_wide.columns = new_cols
    prox_indiv_schc_combined_wide = prox_indiv_schc_combined_wide.reset_index()

    shape_feature_means_wide = shape_feature_means.pivot(index=["slide_submitter_id"], columns="cell_type")[["solidity", "roundness"]]
    new_cols = [f'prox CC {col.replace("_", " ")}' for col in prox_indiv_schc_combined_wide.columns]
    shape_feature_means_wide.columns = [f"{i.capitalize()} {j}"  for i, j in shape_feature_means_wide.columns]
    shape_feature_means_wide = shape_feature_means_wide.reset_index()

    # Store features
    all_features = pd.merge(frac_high_wide, num_clust_slide_wide, on=["slide_submitter_id"])
    all_features = pd.merge(all_features, all_prox_df_wide)
    all_features = pd.merge(all_features, prox_indiv_schc_combined_wide)
    all_features = pd.merge(all_features, shape_feature_means_wide)

    tiles_all_schc = tiles_all_schc.drop(axis=1, columns=cell_types) # drop the predicted probabilities
    all_slide_indiv_clusters = all_slide_indiv_clusters.drop(axis=1, columns=cell_types)# drop the predicted probabilities

    ################################
    # ---- Store all features ---- #
    ################################

    # tiles_all_schc (DataFrame): dataframe containing the metadata columns and the cluster_label (int)
    # all_slide_clusters_characterized (DataFrame): dataframe containing the slide_submitter_id, and the the columns for the cell types filled with booleans (True if the cluster is assigned with that cell type)
    # all_slide_indiv_clusters (DataFrame): dataframe containing the metadata columns and columns with to which cell type cluster the tile belongs to
    # slide_indiv_clusters_labeled (DataFrame): dataframe containing the slide_submitter_id, cell_type_map, cluster_label (int), and is_high (abundance)
    # all_prox_df (DataFrame): dataframe containing slide_submitter_id, pair, proximity
    # prox_indiv_schc_combined (DataFrame): dataframe containing slide_submitter_id, comparison (high/low abundance label), pair (cell type pair) and proximity
    # shape_features_mean (DataFrame): dataframe containing slide_submitter_id, cell_type, slide_submitter_id, solidity, roundness
    tiles_all_schc.to_csv(f"{output_dir}/{slide_type}_all_schc_tiles.csv", sep="\t", index=False)
    all_slide_clusters_characterized.to_csv(f"{output_dir}/{slide_type}_all_schc_clusters_labeled.csv", sep="\t", index=False)
    all_slide_indiv_clusters.to_csv(f"{output_dir}/{slide_type}_indiv_schc_tiles.csv", sep="\t", index=False)
    slide_indiv_clusters_labeled.to_csv(f"{output_dir}/{slide_type}_indiv_schc_clusters_labeled.csv", sep="\t", index=False)
    all_prox_df.to_csv(f"{output_dir}/{slide_type}_features_clust_all_schc_prox.csv", sep="\t", index=False)
    prox_indiv_schc_combined.to_csv(f"{output_dir}/{slide_type}_features_clust_indiv_schc_prox.csv", sep="\t", index=False)
    shape_feature_means.to_csv(f"{output_dir}/{slide_type}_features_clust_shapes.csv", sep="\t", index=False)
    all_features.to_csv(f"{output_dir}/{slide_type}_clustering_features.csv", sep="\t", index=False)


def post_processing(output_dir, slide_type="FF", metadata_path="", is_TCGA=False, merge_var="slide_submitter_id", sheet_name=None):
    """
    Combine network and clustering features into a single file. If metadata_path is not None, add the metadata as well, based on variable slide_submitter_id

    Args:
        output_dir (str): directory containing the graph and clustering features
        slide_type (str): slide type to identify correct files for merging, either "FF" or "FFPE" (default="FF")
        metadata_path (str): path to file containing metadata
        is_TCGA (bool): whether data is from TCGA
        merge_var (str): variable on which to merge (default: slide_submitter_id)

    """
    output_dir = f"{output_dir}/4_spatial_features"
    if not path.exists(output_dir):
        os.mkdir(output_dir)
    all_features_graph = pd.read_csv(f"{output_dir}/features/{slide_type}_all_graph_features.csv", sep="\t")
    all_features_clustering =  pd.read_csv(f"{output_dir}/features/{slide_type}_clustering_features.csv", sep="\t")

    all_features_combined = pd.merge(all_features_graph, all_features_clustering)

    # Add additional identifiers for TCGA
    if is_TCGA:
        all_features_combined["TCGA_patient_ID"] = all_features_combined.slide_submitter_id.str[0:12]
        all_features_combined["TCGA_sample_ID"] = all_features_combined.slide_submitter_id.str[0:15]
        all_features_combined["sample_submitter_id"] = all_features_combined.slide_submitter_id.str[0:16]

    if path.isfile(metadata_path):
        file_extension = metadata_path.split(".")[-1]
        if (file_extension.startswith("xls")):
            if sheet_name is None:
                metadata = pd.read_excel(metadata_path)
        elif (file_extension == "txt") or (file_extension == "csv"):
            metadata = pd.read_csv(metadata_path, sep="\t")
        all_features_combined = pd.merge(all_features_combined, metadata, on=merge_var, how="left")
    all_features_combined.to_csv(f"{output_dir}/all_features_combined.csv", sep="\t", index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Derive spatial features")
    parser.add_argument("--workflow_mode", type=int, help="Choose which steps to execute: all = 1, graph-based = 2, clustering-based = 3, combining features = 4 (default: 1)", default=1)
    parser.add_argument("--tile_quantification_path", type=str, help="Path to csv file with tile-level quantification (predictions)", required=True)
    parser.add_argument("--output_dir",type=str, help="Path to output folder to store generated files", required=True)
    parser.add_argument("--slide_type", type=str,help="Type of slides 'FFPE' or 'FF' used for naming generated files (default: 'FF')", default="FF")
    parser.add_argument("--cell_types_path", type=str,help="Path to file with list of cell types (default: CAFs, endothelial_cells, T_cells, tumor_purity)", default="")
    parser.add_argument("--graphs_path", type=str,help="Path to pkl with generated graphs in case this was done before (OPTIONAL) if not specified, graphs will be generated", default=None)

    parser.add_argument("--cutoff_path_length", type=int, help="Max path length for proximity based on graphs", default=2, required=False)
    parser.add_argument("--shapiro_alpha", type=float, help="Choose significance level alpha (default: 0.05)", default=0.05, required=False)
    parser.add_argument("--abundance_threshold", type=float, help="Threshold for assigning cell types based on predicted probability (default: 0.5)", default=0.5, required=False)

    parser.add_argument("--n_clusters", type=int, help="Number of clusters for SCHC (default: 8)", required=False, default=8)
    parser.add_argument("--max_dist", type=int, help="Maximum distance between clusters", required=False, default=None)
    parser.add_argument("--max_n_tiles_threshold", type=int, help="Number of tiles for computing max. distance between two points in two different clusters", default=2, required=False)
    parser.add_argument("--tile_size", type=int, help="Size of tile (default: 512)", default=512, required=False)
    parser.add_argument("--overlap", type=int, help="Overlap of tiles (default: 50)", default=50, required=False)

    parser.add_argument("--metadata_path", type=str,help="Path to tab-separated file with metadata", default="")
    parser.add_argument("--is_TCGA", type=bool, help="dataset is from TCGA (default: True)", default=True, required=False)
    parser.add_argument("--merge_var", type=str,help="Variable to merge metadata and computed features on", default=None)
    parser.add_argument("--sheet_name", type=str,help="Name of sheet for merging in case a path to xls(x) file is given for metadata_path", default=None)

    args=parser.parse_args()

    # Common variables
    workflow_mode = args.workflow_mode
    tile_quantification_path = args.tile_quantification_path
    slide_type = args.slide_type
    cell_types_path = args.cell_types_path
    graphs_path = args.graphs_path
    output_dir = args.output_dir

    # Variables for network/graph based features
    cutoff_path_length = args.cutoff_path_length
    shapiro_alpha = args.shapiro_alpha
    abundance_threshold = args.abundance_threshold

    # Variables clustering
    n_clusters = args.n_clusters
    max_dist=args.max_dist
    max_n_tiles_threshold=args.max_n_tiles_threshold
    tile_size =args.tile_size
    overlap = args.overlap

    # Variables post-processing
    metadata_path = args.metadata_path
    merge_var = args.merge_var
    sheet_name=args.sheet_name

    if path.isfile(cell_types_path):
        cell_types = pd.read_csv(cell_types_path, header=None).to_numpy().flatten()
    else:
        cell_types = DEFAULT_CELL_TYPES

    if (workflow_mode in [1, 2, 3]) & (graphs_path is None):
        predictions = pd.read_csv(tile_quantification_path, sep="\t")
        slide_submitter_ids = list(set(predictions.slide_submitter_id))
        results = Parallel(n_jobs=NUM_CORES)(
            delayed(graphs.construct_graph)(predictions=predictions, slide_submitter_id=id)
            for id in slide_submitter_ids
        )
        # Extract/format graphs
        all_graphs = {
            list(slide_graph.keys())[0]: list(slide_graph.values())[0]
            for slide_graph in results
        }
        joblib.dump(
            all_graphs, f"{output_dir}/{slide_type}_graphs.pkl")

        graphs_path = f"{output_dir}/{slide_type}_graphs.pkl"

    if workflow_mode == 1:
        print("Workflow mode: all steps")

        print("Compute network features...")
        compute_network_features(
        tile_quantification_path=tile_quantification_path,
        output_dir=output_dir,
        slide_type=slide_type,
        cell_types=cell_types,
        graphs_path=graphs_path, cutoff_path_length=cutoff_path_length, shapiro_alpha=shapiro_alpha, abundance_threshold=abundance_threshold)

        print("Compute clustering features...")
        compute_clustering_features(
            tile_quantification_path=tile_quantification_path,
            output_dir=output_dir,
            slide_type=slide_type,
            cell_types=cell_types,
            graphs_path=graphs_path, n_clusters=n_clusters, max_dist=max_dist, max_n_tiles_threshold=max_n_tiles_threshold, tile_size=tile_size, overlap=overlap)

        print("Post-processing: combining all features")
        post_processing(output_dir=output_dir, slide_type=slide_type, metadata_path=metadata_path, is_TCGA=False, merge_var=merge_var, sheet_name=sheet_name)

        print("Finished with all steps.")
    elif workflow_mode == 2:
        print("Compute network features...")
        compute_network_features(
              tile_quantification_path=tile_quantification_path,
        output_dir=output_dir,
        slide_type=slide_type,
        cell_types=cell_types,
        graphs_path=graphs_path, cutoff_path_length=cutoff_path_length, shapiro_alpha=shapiro_alpha, abundance_threshold=abundance_threshold)
        print("Finished.")
    elif workflow_mode == 3:
        print("Compute clustering features...")
        compute_clustering_features(
          tile_quantification_path=tile_quantification_path,
            output_dir=output_dir,
            slide_type=slide_type,
            cell_types=cell_types,
            graphs_path=graphs_path, n_clusters=n_clusters, max_dist=max_dist, max_n_tiles_threshold=max_n_tiles_threshold, tile_size=tile_size, overlap=overlap)
        print("Finished.")

    elif workflow_mode == 4:
        print("Post-processing: combining all features")

        post_processing(output_dir=output_dir, slide_type=slide_type, metadata_path=metadata_path, is_TCGA=False, merge_var=merge_var, sheet_name=sheet_name)
        print("Finished.")
    else:
        raise Exception("Invalid workflow mode, please choose one of the following (int): all = 1, graph-based = 2, clustering-based = 3, combining features = 4 (default: 1)")
