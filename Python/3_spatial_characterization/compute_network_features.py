import multiprocessing
import sys
import argparse
import joblib
from joblib import Parallel, delayed
import pandas as pd


import features.features as features # trunk-ignore(flake8/E402)
import features.graphs as graphs # trunk-ignore(flake8/E402)
import features.utils as utils # trunk-ignore(flake8/E402)
def compute_network_features(tile_quantification_path, output_dir, slide_type="FF", cell_types=None, graphs_path=None):
    NUM_CORES = multiprocessing.cpu_count()

    if cell_types is None:
        cell_types = ["CAFs", "T_cells", "endothelial_cells", "tumor_purity"]

    predictions = pd.read_csv(tile_quantification_path, sep="\t", index_col=0)

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

    #######################################################
    # ---- Compute connectedness and co-localization ---- #
    #######################################################

    all_largest_cc_sizes = []
    all_dual_nodes_frac = []
    for _, slide_submitter_id in slides.to_numpy():
        slide_data = utils.get_slide_data(predictions, slide_submitter_id)
        node_cell_types = utils.assign_cell_types(slide_data)
        lcc = features.determine_lcc(
            graph=all_graphs[slide_submitter_id], cell_type_assignments=node_cell_types
        )
        lcc["slide_submitter_id"] = slide_submitter_id
        all_largest_cc_sizes.append(lcc)

        dual_nodes_frac = features.compute_dual_node_fractions(node_cell_types)
        dual_nodes_frac["slide_submitter_id"] = slide_submitter_id
        all_dual_nodes_frac.append(dual_nodes_frac)

    all_largest_cc_sizes = pd.concat(all_largest_cc_sizes, axis=0)
    all_dual_nodes_frac = pd.concat(all_dual_nodes_frac, axis=0)

    #######################################################
    # ---- Compute N shortest paths with max. length ---- #
    #######################################################

    results = Parallel(n_jobs=NUM_CORES)(
        delayed(features.compute_n_shortest_paths_max_length)(
            predictions=predictions, slide_submitter_id=slide_submitter_id, graph=all_graphs[slide_submitter_id]
        )
        for _, slide_submitter_id in slides.to_numpy()
    )
    # Formatting and count the number of shortest paths with max length
    all_shortest_paths_thresholded = pd.concat(results, axis=0)
    all_shortest_paths_thresholded["n_paths"] = 1
    proximity_graphs = (
        all_shortest_paths_thresholded.groupby(
            ["slide_submitter_id", "source", "target"]
        )
        .sum()
        .reset_index()
    )
    # Post-processing
    proximity_graphs["pair"] = [f"{source}-{target}" for source, target in proximity_graphs[["source", "target"]].to_numpy()]
    proximity_graphs = proximity_graphs.drop(columns=["path_length"])

    ###############################################
    # ---- Compute ES based on ND difference ---- #
    ###############################################

    nd_results= Parallel(n_jobs=NUM_CORES)(delayed(features.node_degree_wrapper)(all_graphs[slide_submitter_id], predictions, slide_submitter_id) for _, slide_submitter_id in slides.to_numpy())

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
    # shapiro_tests = Parallel(n_jobs=NUM_CORES)(delayed(utils.test_normality)(sims_nd=all_sims_nd, slide_submitter_id=slide_submitter_id) for slide_submitter_id in all_sims_nd.slide_submitter_id.unique())
    # all_shapiro_tests = pd.concat(shapiro_tests, axis=0)
    # print(f"Number of samples from normal distribution { len(all_shapiro_tests) -  all_shapiro_tests.is_not_normal.sum()}/{len(all_shapiro_tests)}")

    # Computing Cohen's d effect size and perform t-test
    effect_sizes = Parallel(n_jobs=NUM_CORES)(delayed(features.compute_effect_size)(all_mean_nd_df, all_sims_nd, slide_submitter_id) for slide_submitter_id in all_sims_nd.slide_submitter_id.unique())
    all_effect_sizes = pd.concat(effect_sizes, axis=0)
    all_effect_sizes["pair"] = [f"{c}-{n}" for c, n in all_effect_sizes[["center", "neighbor"]].to_numpy()]

    ########################
    # ---- Formatting ---- #
    ########################
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

    ################################
    # ---- Store all features ---- #
    ################################

    # all_effect_sizes (DataFrame): dataframe containing the slide_submitter_id, center, neighbor, effect_size (Cohen's d), Tstat, pval, and the pair (string of center and neighbor)
    # all_sims_nd (DataFrame): dataframe containing slide_submitter_id, center, neighbor, simulation_nr and degree (node degree)
    # all_mean_nd_df (DataFrame): dataframe containing slide_submitter_id, center, neighbor, mean_sim (mean node degree across the N simulations), mean_obs
    # all_largest_cc_sizes (DataFrame): dataframe containing slide_submitter_id, cell type and type_spec_frac (fraction of LCC w.r.t. all tiles for cell type)
    # shortest_paths_slide (DataFrame): dataframe containing slide_submitter_id, source, target, pair and n_paths (number of shortest paths for a pair)
    # all_dual_nodes_frac (DataFrame): dataframe containing slide_submitter_id, pair, counts (absolute) and frac

    all_effect_sizes.to_csv(
        f"{output_dir}/{slide_type}_features_ND_ES.csv", sep="\t", index=False)
    all_sims_nd.to_csv(
        f"{output_dir}/{slide_type}_features_ND_sims.csv", sep="\t", index=False)
    all_mean_nd_df.to_csv(
        f"{output_dir}/{slide_type}_features_ND.csv", sep="\t", index=False)
    joblib.dump(example_simulations,
        f"{output_dir}/{slide_type}_features_ND_sim_assignments.pkl")

    all_largest_cc_sizes_wide.to_csv(f"{output_dir}/{slide_type}_features_lcc_fraction.csv", sep="\t", index=False)
    proximity_graphs.to_csv(f"{output_dir}/{slide_type}_features_shortest_paths_thresholded.csv", sep="\t", index=False)
    all_dual_nodes_frac.to_csv(f"{output_dir}/{slide_type}_features_coloc_fraction.csv", sep="\t", index=False)

    all_features.to_csv(f"{output_dir}/{slide_type}_all_graph_features.csv", sep="\t", index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute network features")
    parser.add_argument("--tile_quantification_path", help="Path to csv file with tile-level quantification (predictions)")
    parser.add_argument("--output_dir", help="Path to output folder to store generated files")
    parser.add_argument("--slide_type", help="Type of slides 'FFPE' or 'FF' used for naming generated files (by default='FF')", default="FF")
    parser.add_argument("--cell_types", help="List of cell types", default=None)
    parser.add_argument("--graphs_path", help="Path to pkl with generated graphs in case this was done before (OPTIONAL) if not specified, graphs will be generated", default=None)

    args=parser.parse_args()
    compute_network_features(
        tile_quantification_path=args.tile_quantification_path,
        output_dir=args.output_dir,
        slide_type=args.slide_type,
        cell_types=args.cell_types,
        graphs_path=args.graphs_path)
