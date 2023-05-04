import joblib
import numpy as np
import pandas as pd
import scipy.stats as stats
from scipy.stats import pearsonr, spearmanr
import os
import logging

def custom_spearmanr(x, y, in_gridsearch=True):
    """
    Compute spearman correlation (of absolute values) of main diagonal

    Args:
        x (array): vector containing features
        y (array): vector containing features

    Returns:
        (float): mean value of diagonal of spearman correlation matrix

        (array): array containing the diagonal of spearman correlation matrix
    """
    if in_gridsearch:
        return np.mean(abs(np.diag(stats.spearmanr(x, y, axis=0)[0], x.shape[1])))
    else:
        return np.diag(stats.spearmanr(x, y, axis=0)[0], x.shape[1])


def compute_aggregated_scores(
    y_true, y_pred, target_vars, level_ids, aggregation_level="slide_submitter_id"
):
    """
    Compute the aggregated scores on a defined level e.g. slide level

    Args:
        y_true (array): array containing the true values for the tasks
        y_pred (array): array containing the predicted values by the algorithm for the tasks
        target_vars (list): list of quantification method names
        level_ids (list): list of IDs, e.g slide_submitter_ids if aggregation_level="slide_submitter_id""
        aggregation_level (str)

    Returns:
        (list): nested list of 1) aggregated quantification of individual tasks for true abundances, 2) aggregated quantification of individual tasks for predicted abundances
    """
    y_pred_df = pd.DataFrame(y_pred, columns=target_vars)
    y_true_df = pd.DataFrame(y_true, columns=target_vars)

    # add IDs
    y_pred_df[aggregation_level] = level_ids.to_numpy()
    y_true_df[aggregation_level] = level_ids.to_numpy()

    # Compute score by level
    y_pred_agg = y_pred_df.groupby([aggregation_level]).mean()
    y_true_agg = y_true_df.groupby([aggregation_level]).mean()
    return [y_true_agg, y_pred_agg]


def compute_tile_predictions(
    cell_type,
    models_dir,
    var_names,
    n_outerfolds=5,
    prediction_mode='test',
    X=pd.DataFrame(),
    metadata=pd.DataFrame(),
    slide_type="FF"
):
    """
    Predict the scores using the trained multi-task lasso models
    1. prediction_mode=tcga_train_validation: predict for all tiles (validation+train), average cross cell-type signatures and average across folds
    2. prediction_mode=tcga_validation: predict only the tiles in the test sets
    3. prediction_mode=test: predict in external dataset

    Args:
        cell_type (str): one of the following CAFs, endothelial_cells, tumor_purity or T_cells
        models_dir (str): path pointing to the folder where the files associated with the trained models are stored
        var_names (dict): dictionary with keys CAFs, T_cells, tumor_purity, endothelial_cells, IDs, and tile IDs (created in 1_processing_transcriptomics.ipynb)

    Returns:
        all_scores (DataFrame): dataframe containing the mean predicted scores for the individual tasks and the mean combined score for the given cell type
    """


    # Check model directory
    if models_dir.find(cell_type) != -1: # models_dir = folder for specific cell type
        full_model_path = models_dir
    else: # models_dir = parent directory
        subfolders = os.listdir(models_dir)
        # Find subfolder for given cell type
        full_model_path = f"{models_dir}/{[s for s in subfolders if cell_type in s][0]}"

    # Load models and scalers
    outer_models = joblib.load(f"{full_model_path}/outer_models.pkl")
    x_scalers = joblib.load(f"{full_model_path}/x_train_scaler.pkl")

    # Check prediction mode
    if  any([prediction_mode == item for item in ['tcga_train_validation', 'tcga_validation']]):
        outer_fold_splits = joblib.load(f"{full_model_path}/cv_outer_splits.pkl")
        total_tile_selection = joblib.load(f"{full_model_path}/total_tile_selection.pkl")[["tile_ID", "slide_submitter_id"]]

    # Initialize variables
    combi_scores = []
    test_pred = dict.fromkeys(range(n_outerfolds))
    index_pred = dict.fromkeys(range(n_outerfolds))

    for i in range(n_outerfolds):
        # extract best model and scalers
        outer_model = outer_models[i]
        if any([prediction_mode == item for item in ['tcga_train_validation', 'test']]):
	        X_test = X
        elif prediction_mode == 'tcga_validation':
            # Only select all TILES from the slide submitter ids (slides) in the test set
            _, test_index = outer_fold_splits[i]
            test_slides = total_tile_selection.loc[test_index, "slide_submitter_id"].unique()
            if slide_type == "FF":
               tiles_subset = metadata[metadata["slide_submitter_id"].isin(test_slides)]
               X_test = X.loc[tiles_subset.tile_ID,:]
            elif slide_type == "FFPE":
               tiles_subset = metadata[metadata["slide_submitter_id"].isin(test_slides)]["tile_ID"]
               X_test = X.map_partitions(lambda x: x[x.index.isin(tiles_subset)])
               X_test = X_test.compute()

        # scale to training data
        x_scaler = x_scalers[i]
        X_z = x_scaler.transform(X_test)

        # Predict tasks for each tile and average their prediction across quantification methods
        test_pred[i] = outer_model.predict(X_z).mean(axis=1)
        index_pred[i] = X_test.index

    if any([prediction_mode == item for item in ['tcga_train_validation', 'test']]):
        combi_scores = np.concatenate((test_pred[0][:,np.newaxis],test_pred[1][:,np.newaxis],test_pred[2][:,np.newaxis], test_pred[3][:,np.newaxis], test_pred[4][:,np.newaxis]),axis=1)
	# Average across folds for final prediction
        combi_scores = combi_scores.mean(axis=1)
        combi_indexes = index_pred[1]

    elif prediction_mode == 'tcga_validation':
        combi_scores = np.concatenate((test_pred[0],test_pred[1],test_pred[2], test_pred[3], test_pred[4]),axis=0)
        combi_indexes = np.concatenate((index_pred[0],index_pred[1],index_pred[2], index_pred[3], index_pred[4]),axis=0)

    print("Finished with predictions for cell type: {}...".format(cell_type))

    combi_scores_df = pd.DataFrame(combi_scores)
    combi_scores_df.columns = [f"{cell_type} (combi)"]
    combi_scores_df.index = combi_indexes
    return combi_scores_df


def compute_correlations(df1, df2, corr_method="pearson"):
    """
    Compute correlations between two dataframes

    Args:
        df1 (DataFrame)
        df2 (DataFrame)
        corr_method (str): pearson or spearman

    Returns:
        (list): nested list of two dataframes containing the computed correlation coefficients and the corresponding p-values
    """
    if (corr_method != "pearson") | (corr_method != "spearman"):
        raise Exception("Choose a valid method: 'pearson' or 'spearman'")

    corr_coefs = pd.DataFrame(columns=df1.columns, index=df2.columns)
    pvalues = pd.DataFrame(columns=df1.columns, index=df2.columns)

    for r in df1.columns:
        for c in df2.columns:
            if corr_method == "spearman":
                pvalues[r][c] = spearmanr(df1[r], df2[c])[1]
                corr_coefs[r][c] = spearmanr(df1[r], df2[c])[0]

            elif corr_method == "pearson":
                pvalues[r][c] = pearsonr(df1[r], df2[c])[1]
                corr_coefs[r][c] = spearmanr(df1[r], df2[c])[0]

    return [corr_coefs.astype(float), pvalues.astype(float)]
