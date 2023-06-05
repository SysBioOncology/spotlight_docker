import pandas as pd
import scipy.stats as stats
import itertools
import sys
import os
from scipy.stats import pearsonr, spearmanr
import numpy as np

sys.path.append(f"{os.path.dirname(os.getcwd())}/Python/libs")

# Own modules
from model.constants import *

def get_slide_data(data, slide_submitter_id):
    """
    Select data for slide

    Args:
        data (DataFrame): dataframe containing at least the given slide_submitter_id
        slide_submitter_id (str)

    Returns:
        (DataFrame): dataframe containing only the data from the given slide id
    """
    slide_data = data[data["slide_submitter_id"] == slide_submitter_id]
    slide_data = slide_data.reset_index(drop=True)
    return slide_data


def assign_cell_types(
    slide_data,
    cell_types=None,
    threshold=0.5,
):
    """
    Assign nodes with cell types based on their value being greater than the threshold (by default=0.5)

    Args:
        slide_data (DataFrame): dataframe containing at least the cell type predictions for the tiles
        cell_types (list): list of cell types
        threshold (float): threshold to assign a cell type to a tile

    Returns
        (DataFrame): Tiles with assigned cell types len(cell_types) columns with True/False indicating whether the tile is assigned with the cell type and corresponding metadata.
    """
    if cell_types is None:
        cell_types = DEFAULT_CELL_TYPES

    has_cell_types = slide_data[cell_types] > threshold
    node_cell_types = slide_data.copy()
    node_cell_types[cell_types] = has_cell_types
    return node_cell_types


def test_normality(sims_nd, slide_submitter_id, alpha=0.05, cell_types=None):
    """
    Test for normality using the shapiro tests
    Args:
        sims_nd (DataFrame): dataframe containing the node degree for each individual simulated slide for all cell type pairs
        slide_submitter_id (str): e.g. TCGA-77-6842-01A-01-BS1
        alpha (float): significance level

    Returns:
        (DataFrame): containing slide_submitter_id, center, neighbor, resulting p-value and variable indicating the result of the test (True/False)
    """
    if cell_types is None:
        cell_types = DEFAULT_CELL_TYPES
    test_res = []
    cell_type_interactions = list(itertools.combinations(cell_types, 2)) + list(
        zip(cell_types, cell_types)
    )
    for center, neighbor in cell_type_interactions:
        mask = (
            (sims_nd.slide_submitter_id == slide_submitter_id)
            & (sims_nd.center == center)
            & (sims_nd.neighbor == neighbor)
        )
        if sum(mask) > 0:
            sim_degrees = sims_nd.loc[mask, "degree"]
            _, pval = stats.shapiro(sim_degrees)
            test_res.append([slide_submitter_id, center, neighbor, pval, pval < alpha])
    return pd.DataFrame(
        test_res,
        columns=["slide_submitter_id", "center", "neighbor", "pval", "is_not_normal"],
    )


def cohens_d(theoretical_mean, sample_mean, sample_sd):
    """
    Compute effect size (Cohen's d)

    Args:
        theoretical_mean (float): 'true' mean
        sample_mean (float): mean from the observations
        sample_sd (float): mean from the observations

    Returns:
        (float): effect size
    """
    return (sample_mean - theoretical_mean) / sample_sd


def compute_correlations(df1, df2, corr_method="pearson"):
    corr_coefs = pd.DataFrame(columns=df1.columns, index=df2.columns)
    pvalues = pd.DataFrame(columns=df1.columns, index=df2.columns)
    for r in df1.columns:
        for c in df2.columns:
            if r == c:
                sub = df1[r]
                sub[sub < -1e300] = np.nan
                sub = sub.dropna(how="any", axis=0).values
                if len(sub) >= 25:
                    pvalues[r][c] = pearsonr(sub, sub)[1]
                    corr_coefs[r][c] = pearsonr(sub, sub)[0]

            else:
                temp = pd.merge(df1[r], df2[c], left_index=True, right_index=True)
                temp[temp < -1e300] = np.nan
                temp = temp.dropna()
                if len(temp) >=25:
                    if corr_method == "spearman":
                        pvalues[r][c] = spearmanr(temp[r].values, temp[c].values)[1]
                        corr_coefs[r][c] = spearmanr(temp[r].values, temp[c].values)[0]

                    elif corr_method == "pearson":
                        pvalues[r][c] = pearsonr(temp[r].values, temp[c].values)[1]
                        corr_coefs[r][c] = pearsonr(temp[r].values, temp[c].values)[0]
    return corr_coefs, pvalues
