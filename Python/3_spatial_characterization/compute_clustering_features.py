import multiprocessing
import sys
import joblib
import pandas as pd
from joblib import Parallel, delayed
import argparse

sys.path.append(f"{os.path.dirname(os.getcwd())}/Python/libs")

import features.clustering as clustering # trunk-ignore(flake8/E402)
import features.features as features # trunk-ignore(flake8/E402)
import features.graphs as graphs # trunk-ignore(flake8/E402)

NUM_CORES = multiprocessing.cpu_count()

def compute_clustering_features(tile_quantification_path, output_dir, slide_type="FF", cell_types=None, graphs_path=None):

    if cell_types is None:
        cell_types = ["CAFs", "T_cells", "endothelial_cells", "tumor_purity"]

    predictions = pd.read_csv(tile_quantification_path, sep="\t", index_col=0)]

    #####################################
    # ---- Constructing the graphs ---- #
    #####################################

    if graphs_path is None:
        results = Parallel(n_jobs=NUM_CORES)(
            delayed(graphs.construct_graph)(predictions=predictions, slide_submitter_id=slide_submitter_id)
            for _, slide_submitter_id in slides.to_numpy()
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
    slide_clusters= Parallel(n_jobs=NUM_CORES)(delayed(clustering.schc_all)(predictions, all_graphs[slide_submitter_id], slide_submitter_id) for subtype, slide_submitter_id in slides.to_numpy())
    # Combine the tiles labeled with their cluster id for all slides
    tiles_all_schc = pd.concat(slide_clusters, axis=0)

    # Assign a cell type label based on the mean of all cluster means across all slides
    all_slide_clusters_characterized = clustering.characterize_clusters(tiles_all_schc)

    # Count the number of clusters per cell type for each slide
    num_clust_by_slide = features.n_clusters_per_cell_type(all_slide_clusters_characterized)

    ######################################################################################
    # ---- Fraction of highly abundant cell types (individual cell type clustering) ---- #
    ######################################################################################

    # Spatially Hierarchical Constrained Clustering with all quantification of all cell types for each individual cell type
    slide_indiv_clusters= Parallel(n_jobs=NUM_CORES)(delayed(clustering.schc_individual)(predictions, all_graphs[slide_submitter_id], slide_submitter_id) for subtype, slide_submitter_id in slides.to_numpy())
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
    clusters_all_schc_long = all_slide_clusters_characterized.melt(id_vars=["MFP", "slide_submitter_id", "cluster_label"], value_name="is_assigned", var_name="cell_type")
    # remove all cell types that are not assigned to the cluster
    clusters_all_schc_long = clusters_all_schc_long[clusters_all_schc_long["is_assigned"]]
    clusters_all_schc_long = clusters_all_schc_long.drop(columns="is_assigned")

    results_schc_all= Parallel(n_jobs=NUM_CORES)(delayed(features.compute_proximity_clusters_pairs)(tiles_all_schc, slide_submitter_id, method="all") for _, slide_submitter_id in slides.to_numpy())
    prox_all_schc = pd.concat(results_schc_all)

    # Label clusters (a number) with the assigned cell types
    prox_all_schc = pd.merge(prox_all_schc, clusters_all_schc_long, left_on=["slide_submitter_id", "cluster1"], right_on=["slide_submitter_id", "cluster_label"])
    prox_all_schc = prox_all_schc.rename(columns={"cell_type": "cluster1_label"})
    prox_all_schc = prox_all_schc.drop(columns=["cluster_label", "MFP"])

    prox_all_schc = pd.merge(prox_all_schc, clusters_all_schc_long, left_on=["slide_submitter_id", "cluster2"], right_on=["slide_submitter_id", "cluster_label"])
    prox_all_schc = prox_all_schc.rename(columns={"cell_type": "cluster2_label"})
    # prox_all_schc = prox_all_schc.drop(columns=["cluster_label", "cluster1", "cluster2" ])

    # Order doesn't matter: x <->
    prox_all_schc["pair"] = [f"{sorted([i, j])[0]}-{sorted([i, j])[1]}" for i, j in prox_all_schc[["cluster1_label", "cluster2_label"]].to_numpy()]
    prox_all_schc = prox_all_schc[((prox_all_schc.cluster1 == prox_all_schc.cluster2) & (prox_all_schc.cluster2_label != prox_all_schc.cluster1_label)) | (prox_all_schc.cluster1 != prox_all_schc.cluster2)]
    # prox_all_schc.to_csv(f"{output_dir}/{slide_type}_features_clust_all_schc_prox.txt", sep="\t")

    slides = prox_all_schc[["MFP", "slide_submitter_id"]].drop_duplicates().to_numpy()

    # Post Processing
    results_schc_all= Parallel(n_jobs=NUM_CORES)(delayed(features.post_processing_proximity)(prox_df=prox_all_schc, slide_submitter_id=slide_submitter_id,subtype=subtype, method="all") for subtype, slide_submitter_id in slides)
    all_prox_df = pd.concat(results_schc_all)
    # Remove rows with a proximity of NaN
    all_prox_df = all_prox_df.dropna(axis=0)

    ##########################################################################
    # ---- Compute proximity features (individual cell type clustering) ---- #
    ##########################################################################

    ## Computing proximity for clusters derived for each cell type individually
    # Between clusters
    slides = (
        predictions[["MFP", "slide_submitter_id"]].drop_duplicates().reset_index(drop=True))
    results_schc_indiv= Parallel(n_jobs=NUM_CORES)(delayed(features.compute_proximity_clusters_pairs)(all_slide_indiv_clusters, slide_submitter_id, method="individual_between") for _, slide_submitter_id in slides.to_numpy())
    prox_indiv_schc = pd.concat(results_schc_indiv)

    # Formatting
    prox_indiv_schc = pd.merge(prox_indiv_schc,slide_indiv_clusters_labeled, left_on=["slide_submitter_id", "cluster1_label", "cluster1"], right_on=["slide_submitter_id", "cell_type_map", "cluster_label"])
    prox_indiv_schc = prox_indiv_schc.drop(columns=["cell_type_map", "MFP", "cluster_label"])
    prox_indiv_schc = prox_indiv_schc.rename(columns={"is_high": "cluster1_is_high"})
    prox_indiv_schc = pd.merge(prox_indiv_schc,slide_indiv_clusters_labeled, left_on=["slide_submitter_id", "cluster2_label", "cluster2"], right_on=["slide_submitter_id", "cell_type_map", "cluster_label"])
    prox_indiv_schc = prox_indiv_schc.rename(columns={"is_high": "cluster2_is_high"})
    prox_indiv_schc = prox_indiv_schc.drop(columns=["cell_type_map", "cluster_label"])

    # Order matters
    prox_indiv_schc["ordered_pair"] = [f"{i}-{j}" for i, j in prox_indiv_schc[["cluster1_label", "cluster2_label"]].to_numpy()]
    prox_indiv_schc["comparison"] = [f"cluster1={i}-cluster2={j}" for i, j in prox_indiv_schc[["cluster1_is_high", "cluster2_is_high"]].to_numpy()]

    # Post-processing
    slides = prox_indiv_schc[["MFP", "slide_submitter_id"]].drop_duplicates().to_numpy()
    results_schc_indiv= pd.concat(Parallel(n_jobs=NUM_CORES)(delayed(features.post_processing_proximity)(prox_df=prox_indiv_schc, slide_submitter_id=slide_submitter_id,subtype=subtype, method="individual_between") for subtype, slide_submitter_id in slides))

    # Within clusters
    slides = (
        predictions[["MFP", "slide_submitter_id"]].drop_duplicates().reset_index(drop=True))
    results_schc_indiv_within= Parallel(n_jobs=NUM_CORES)(delayed(features.compute_proximity_clusters_pairs)(all_slide_indiv_clusters, slide_submitter_id, method="individual_within") for _, slide_submitter_id in slides.to_numpy())
    prox_indiv_schc_within = pd.concat(results_schc_indiv_within)

    prox_indiv_schc_within = pd.merge(prox_indiv_schc_within,slide_indiv_clusters_labeled, left_on=["slide_submitter_id", "cell_type", "cluster1"], right_on=["slide_submitter_id", "cell_type_map", "cluster_label"])
    prox_indiv_schc_within = prox_indiv_schc_within.drop(columns=["MFP", "cluster_label"])
    prox_indiv_schc_within = prox_indiv_schc_within.rename(columns={"is_high": "cluster1_is_high", "cell_type_map":"cell_type_map1"})
    prox_indiv_schc_within = pd.merge(prox_indiv_schc_within,slide_indiv_clusters_labeled, left_on=["slide_submitter_id", "cell_type", "cluster2"], right_on=["slide_submitter_id", "cell_type_map", "cluster_label"])
    prox_indiv_schc_within = prox_indiv_schc_within.rename(columns={"is_high": "cluster2_is_high", "cell_type_map": "cell_type_map2"})
    prox_indiv_schc_within = prox_indiv_schc_within.drop(columns=["cluster_label"])

    # Order doesn't matter (only same cell type combinations)
    prox_indiv_schc_within["pair"] = [f"{i}-{j}" for i, j in prox_indiv_schc_within[["cell_type_map1", "cell_type_map2"]].to_numpy()]
    prox_indiv_schc_within["comparison"] = [f"cluster1={sorted([i,j])[0]}-cluster2={sorted([i,j])[1]}" for i, j in prox_indiv_schc_within[["cluster1_is_high", "cluster2_is_high"]].to_numpy()]

    # prox_indiv_schc_within.to_csv(f"{output_dir}/{slide_type}_features_clust_indiv_schc_prox_within.txt", sep="\t")
    slides = prox_indiv_schc_within[["slide_submitter_id", "MFP"]].drop_duplicates().to_numpy()
    results_schc_indiv_within= pd.concat(Parallel(n_jobs=NUM_CORES)(delayed(features.post_processing_proximity)(prox_df=prox_indiv_schc_within, slide_submitter_id=slide_submitter_id,subtype=subtype, method="individual_within") for slide_submitter_id, subtype in slides))

    # Concatenate within and between computed proximity values
    prox_indiv_schc_combined = pd.concat([results_schc_indiv_within, results_schc_indiv])

    # Remove rows with a proximity of NaN
    prox_indiv_schc_combined = prox_indiv_schc_combined.dropna(axis=0)

    ####################################
    # ---- Compute shape features ---- #
    ####################################

    # Compute shape features based on clustering with all cell types simultaneously
    slides = (
        predictions[["MFP", "slide_submitter_id"]].drop_duplicates().reset_index(drop=True))

    all_slide_clusters_characterized = all_slide_clusters_characterized.rename(columns=dict(zip(cell_types, [f"is_{cell_type}_cluster" for cell_type in cell_types])))
    tiles_all_schc = pd.merge(tiles_all_schc, all_slide_clusters_characterized, on=["slide_submitter_id", "MFP", "cluster_label"])
    res = pd.concat(Parallel(n_jobs=NUM_CORES)(delayed(features.compute_shape_features)(tiles=tiles_all_schc, slide_submitter_id=slide_submitter_id,subtype=subtype) for subtype, slide_submitter_id in slides.to_numpy()))
    res = res.drop(axis=1, labels=["cluster_label"])
    shape_feature_means = res.groupby(["slide_submitter_id", "cell_type"]).mean().reset_index()

    ##############################################
    # ---- Formatting all computed features ---- #
    ##############################################

    frac_high_sub = frac_high[frac_high["is_high"]].copy()
    frac_high_sub = frac_high_sub.drop(columns=["is_high", "n_clusters", "n_total_clusters"])

    frac_high_wide = frac_high_sub.pivot(index=["MFP", "slide_submitter_id"], columns=["cell_type_map"])["fraction"]
    new_cols=[('fraction {0} clusters labeled high'.format(col)) for col in frac_high_wide.columns]
    frac_high_wide.columns = new_cols
    frac_high_wide = frac_high_wide.sort_index(axis="columns").reset_index()

    num_clust_by_slide_sub = num_clust_by_slide.copy()
    num_clust_by_slide_sub = num_clust_by_slide_sub.drop(columns=["is_assigned", "n_clusters"])

    num_clust_slide_wide = num_clust_by_slide_sub.pivot(index=["MFP", "slide_submitter_id"], columns=["cell_type"])["fraction"]
    new_cols=[('fraction {0} clusters'.format(col)) for col in num_clust_slide_wide.columns]
    num_clust_slide_wide.columns = new_cols
    num_clust_slide_wide = num_clust_slide_wide.sort_index(axis="columns").reset_index()

    all_prox_df_wide = all_prox_df.pivot(index=["MFP", "slide_submitter_id"], columns=["pair"])["proximity"]
    new_cols = [f'prox CC {col.replace("_", " ")} clusters' for col in all_prox_df_wide.columns]
    all_prox_df_wide.columns = new_cols
    all_prox_df_wide = all_prox_df_wide.reset_index()

    prox_indiv_schc_combined.comparison = prox_indiv_schc_combined.comparison.replace(dict(zip(['cluster1=True-cluster2=True', 'cluster1=True-cluster2=False',
        'cluster1=False-cluster2=True', 'cluster1=False-cluster2=False'], ["high-high", "high-low", "low-high", "low-low"])))
    prox_indiv_schc_combined["pair (comparison)"] = [f"{pair} ({comp})" for pair, comp in prox_indiv_schc_combined[["pair", "comparison"]].to_numpy()]
    prox_indiv_schc_combined = prox_indiv_schc_combined.drop(axis=1, labels=["pair", "comparison"])
    prox_indiv_schc_combined_wide = prox_indiv_schc_combined.pivot(index=["MFP", "slide_submitter_id"], columns=["pair (comparison)"])["proximity"]
    new_cols = [f'prox CC {col.replace("_", " ")}' for col in prox_indiv_schc_combined_wide.columns]
    prox_indiv_schc_combined_wide.columns = new_cols
    prox_indiv_schc_combined_wide = prox_indiv_schc_combined_wide.reset_index()

    shape_feature_means_wide = shape_feature_means.pivot(index=["slide_submitter_id"], columns="cell_type")[["solidity", "roundness"]]
    new_cols = [f'prox CC {col.replace("_", " ")}' for col in prox_indiv_schc_combined_wide.columns]
    shape_feature_means_wide.columns = [f"{i.capitalize()} {j}"  for i, j in shape_feature_means_wide.columns]
    shape_feature_means_wide = shape_feature_means_wide.reset_index()

    # Store features
    all_features = pd.merge(frac_high_wide, num_clust_slide_wide, on=["MFP", "slide_submitter_id"])
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
    frac_high.to_csv(f"{output_dir}/{slide_type}_num_clusters_per_cell_type_indiv_clustering.csv", sep="\t", index=False)
    num_clust_by_slide.to_csv(f"{output_dir}/{slide_type}_num_clusters_per_cell_type_all_clustering.csv", sep="\t", index=False)
    tiles_all_schc.to_csv(f"{output_dir}/{slide_type}_all_schc_tiles.csv", sep="\t", index=False)
    all_slide_clusters_characterized.to_csv(f"{output_dir}/{slide_type}_all_schc_clusters_labeled.csv", sep="\t", index=False)
    all_slide_indiv_clusters.to_csv(f"{output_dir}/{slide_type}_indiv_schc_tiles.csv", sep="\t", index=False)
    slide_indiv_clusters_labeled.to_csv(f"{output_dir}/{slide_type}_indiv_schc_clusters_labeled.csv", sep="\t", index=False)
    all_prox_df.to_csv(f"{output_dir}/{slide_type}_features_clust_all_schc_prox.csv", sep="\t", index=False)
    prox_indiv_schc_combined.to_csv(f"{output_dir}/{slide_type}_features_clust_indiv_schc_prox.csv", sep="\t", index=False)
    shape_feature_means.to_csv(f"{output_dir}/{slide_type}_features_clust_shapes.csv", sep="\t", index=False)
    all_features.to_csv(f"{output_dir}/{slide_type}_clustering_features.csv", sep="\t", index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute clustering features")
    parser.add_argument("--tile_quantification_path", help="Path to csv file with tile-level quantification (predictions)")
    parser.add_argument("--output_dir", help="Path to output folder to store generated files")
    parser.add_argument("--slide_type", help="Type of slides 'FFPE' or 'FF' used for naming generated files (by default='FF')", default="FF")
    parser.add_argument("--cell_types", help="List of cell types", default=None) # TODO: adapt to external file with the cell types, easier for parsing
    parser.add_argument("--graphs_path", help="Path to pkl with generated graphs in case this was done before (OPTIONAL) if not specified, graphs will be generated", default=None)

    args=parser.parse_args()

    compute_clustering_features(
        tile_quantification_path=args.tile_quantification_path,
        output_dir=args.output_dir,
        slide_type=args.slide_type,
        cell_types=args.cell_types,
        graphs_path=args.graphs_path)
