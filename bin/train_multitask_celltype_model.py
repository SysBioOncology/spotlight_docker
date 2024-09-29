#!/usr/bin/env python3
#  Module imports
import argparse
import os
import dask.dataframe as dd
import joblib
import numpy as np
import pandas as pd
from sklearn import linear_model, metrics
from sklearn.model_selection import GridSearchCV, GroupKFold
from sklearn.preprocessing import StandardScaler

# Custom imports
import model.evaluate as meval
import model.preprocessing as preprocessing
import model.utils as utils

from argparse import ArgumentParser as AP

from os.path import abspath
import time
from pathlib import Path
import multiprocessing


def get_args():
    # Script description
    description = """Train multi-task cell type model"""

    # Add parser
    parser = AP(
        description=description, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # Inputs
    inputs_group = parser.add_argument_group("Inputs")
    inputs_group.add_argument(
        "--tile_quantification_path",
        type=str,
        help="Path to csv file with tile-level quantification (predictions)",
        required=True,
    )
    inputs_group.add_argument(
        "--path_target_features",
        help="Path pointing to file containing the cell type abundances",
    )
    inputs_group.add_argument(
        "--path_bottleneck_features",
        help="Path pointing to file containing the histopathological features",
    )

    inputs_group.add_argument(
        "--var_names_path",
        help="Path pointing to pkl file with the tasks to use",
    )

    parser.add_argument(
        "--output_dir",
        type=str,
        help="Path to output folder to store generated files",
        default="",
    )
    parser.add_argument(
        "--slide_type",
        type=str,
        help="Type of slides 'FFPE' or 'FF' used for naming generated files (default: 'FF')",
        default="FF",
        choices=["FF", "FFPE"],
    )
    parser.add_argument("--prefix", type=str, help="Prefix for output file", default="")
    parser.add_argument(
        "--n_cores", type=int, help="Number of cores to use (parallelization)"
    )

    parser.add_argument("--category", help="Cell type", type=str, required=True)

    # GridSearch - grid
    gridsearch_group = parser.add_argument_group("GridSearch grid")
    gridsearch_group.add_argument(
        "--alpha_min",
        help="Min. value of hyperparameter alpha, converts to 10^alpha_min",
        type=int,
        default=-4,
    )
    gridsearch_group.add_argument(
        "--alpha_max",
        help="Max. value of hyperparameter alpha, converts to 10^alpha_max",
        type=int,
        default=-1,
    )
    gridsearch_group.add_argument(
        "--n_steps", help="Stepsize for grid [alpha_min, alpha_max]", default=40
    )

    # Cross validation
    crossvalid_group = parser.add_argument_group("Cross validation")
    crossvalid_group.add_argument(
        "--n_innerfolds", help="Number of inner loops", default=10
    )
    crossvalid_group.add_argument(
        "--n_outerfolds", help="Number of outer loops", default=5
    )

    # Setup splits
    splits_group = parser.add_argument_group("Setup splits for Cross-validation")
    splits_group.add_argument(
        "--n_tiles", help="Number of tiles to select per slide", default=50
    )
    splits_group.add_argument(
        "--split_level",
        help="Split level of slides for creating splits",
        default="sample_submitter_id",
    )
    parser.add_argument("--version", action="version", version="0.1.0")
    arg = parser.parse_args()
    arg.output_dir = abspath(arg.output_dir)

    if (arg.output_dir != "") & (not os.path.isdir(arg.output_dir)):
        arg.output_dir = Path(args.output_dir, "models", arg.category)
        os.mkdir(arg.output_dir)
    else:
        arg.output_dir = Path(args.output_dir, arg.category)
        os.mkdir(arg.output_dir)

    if arg.n_cores is None:
        arg.n_cores = multiprocessing.cpu_count()
    return arg


def nested_cv_multitask(
    bottleneck_features_path: str,
    category: str,
    alpha_min: int = -4,
    alpha_max: int = -1,
    n_steps: int = 40,
    n_outerfolds: int = 5,
    n_innerfolds: int = 10,
    n_tiles: int = 50,
    split_level: str = "sample_submitter_id",
    slide_type: str = "FF",
    n_jobs: int = -1,
    var_names_path: str = "task_selection_names.pkl",
    target_features_path: str = "ensembled_selected_tasks.csv",
):
    """Using Transfer Learning to quantify cell types on a tile-level
    Use a nested cross-validation strategy to train a multi-task lasso algorithm. Tuning and evaluation based on spearman correlation.

    Args:
        bottleneck_features_path (str): path to predicted bottleneck features from the PC-ChiP framework
        category (str): Cell type for which you want to train a model
        alpha_min (int, optional): Min. value of hyperparameter grid, converts to 10^alpha_min. Defaults to -4.
        alpha_max (int, optional): Max. value of hyperparameter grid, converts to 10^alpha_max. Defaults to -1.
        n_steps (int, optional): Stepsize for grid in logspace. Defaults to 40.
        n_outerfolds (int, optional): Number of outer loops to use. Defaults to 5.
        n_innerfolds (int, optional): Number of inner loops to use. Defaults to 10.
        n_tiles (int, optional): Number of tiles to select per slide. Defaults to 50.
        split_level (str, optional): Split level of slides for creating splits. Defaults to "sample_submitter_id".
        slide_type (str, optional): Type of slide. Defaults to "FF".
        n_jobs (int, optional): Number of jobs to run in parallel, by default using all cores. Defaults to -1.
        var_names_path (str, optional): Path to 'task_selection_names.pkl'. Defaults to "task_selection_names.pkl".
        target_features_path (str, optional): Path to computed transcriptomics features, i.e. the target features for the model to be trained. Defaults to "ensembled_selected_tasks.csv".

        Returns:
            _type_:
    """
    # Hyperparameter grid for tuning
    alphas = np.logspace(alpha_min, alpha_max, n_steps)
    scoring = meval.custom_spearmanr

    # Load data
    if Path(var_names_path).exists():
        var_names = joblib.load(var_names_path)
    # TODO place models.constants
    var_names["T_cells"] = [
        "Cytotoxic cells",
        "Effector cells",
        "CD8 T cells (quanTIseq)",
        "Immune score",
    ]

    # Load computed transcriptomics features
    target_features = pd.read_csv(
        target_features_path,
        sep="\t",
        index_col=0,
    )

    # Load bottleneck features
    if Path(bottleneck_features_path).suffix == "txt":
        # FF slides
        bottleneck_features = pd.read_csv(
            bottleneck_features_path,
            sep="\t",
            index_col=0,
        )
    else:
        # FFPE slides
        bottleneck_features = dd.read_parquet(bottleneck_features_path)

    target_vars = var_names[category]
    metadata_colnames = var_names["tile_IDs"] + var_names["IDs"] + var_names[category]

    # Preprocessing
    IDs = ["sample_submitter_id", "slide_submitter_id"]
    merged_data = preprocessing.clean_data(
        bottleneck_features, target_features.loc[:, target_vars + IDs], slide_type
    )
    total_tile_selection = utils.selecting_tiles(merged_data, n_tiles, slide_type)
    X, Y = utils.split_in_XY(
        total_tile_selection, metadata_colnames, var_names[category]
    )

    # TF learning
    # Create variables for storing
    model_learned = dict.fromkeys(range(n_outerfolds))
    x_train_scaler = dict.fromkeys(range(n_outerfolds))
    y_train_scaler = dict.fromkeys(range(n_outerfolds))

    multi_task_lasso = linear_model.MultiTaskLasso()

    # Setup nested cv
    sample_id = pd.factorize(total_tile_selection[split_level])[0]
    cv_outer = GroupKFold(n_splits=n_outerfolds)
    cv_inner = GroupKFold(n_splits=n_innerfolds)
    cv_outer_splits = list(cv_outer.split(X, Y, groups=sample_id))

    # Storing scores
    tiles_spearman_train = {}
    tiles_spearman_test = {}
    slides_spearman_train = {}
    slides_spearman_test = {}

    print("Feature matrix dimensions [tiles, features]:", X.shape)
    print("Response matrix dimensions:", Y.shape)

    # Run nested cross-validation
    for outerfold in range(n_outerfolds):
        print(f"Outerfold {outerfold}")
        train_index, test_index = cv_outer_splits[outerfold]
        x_train, x_test = X.iloc[train_index], X.iloc[test_index]
        y_train, y_test = Y.iloc[train_index], Y.iloc[test_index]

        # Standardizing predictors
        scaler_x = StandardScaler()
        scaler_x.fit(x_train)
        x_train_z = scaler_x.transform(x_train)
        x_test_z = scaler_x.transform(x_test)

        # Standardizing targets
        scaler_y = StandardScaler()
        scaler_y.fit(y_train)
        y_train_z = scaler_y.transform(y_train)
        y_test_z = scaler_y.transform(y_test)
        grid = GridSearchCV(
            estimator=multi_task_lasso,
            param_grid=[{"alpha": alphas}],
            cv=cv_inner,
            scoring=metrics.make_scorer(scoring),
            return_train_score=True,
            n_jobs=n_jobs,
        )
        grid.fit(x_train_z, y_train_z, groups=sample_id[train_index])

        # Evaluate on tile level (spearmanr)
        y_train_z_pred = grid.predict(x_train_z)
        y_test_z_pred = grid.predict(x_test_z)

        tiles_spearman_train[outerfold] = meval.custom_spearmanr(
            y_train_z, y_train_z_pred, in_gridsearch=False
        )
        tiles_spearman_test[outerfold] = meval.custom_spearmanr(
            y_test_z, y_test_z_pred, in_gridsearch=False
        )

        # Evaluate on Slide level
        # For aggregation for evaluation in gridsearch for hyper parameter choosing
        slide_IDs_train = total_tile_selection["slide_submitter_id"].iloc[train_index]
        slide_IDs_test = total_tile_selection["slide_submitter_id"].iloc[test_index]

        # 1. Aggregate on slide level (averaging)
        Y_train_true_agg, Y_train_pred_agg = meval.compute_aggregated_scores(
            y_train_z, y_train_z_pred, target_vars, slide_IDs_train
        )
        Y_test_true_agg, Y_test_pred_agg = meval.compute_aggregated_scores(
            y_test_z, y_test_z_pred, target_vars, slide_IDs_test
        )

        # 2. Compute spearman correlation between ground truth and predictions
        slides_spearman_train[outerfold] = meval.custom_spearmanr(
            Y_train_true_agg, Y_train_pred_agg, in_gridsearch=False
        )
        slides_spearman_test[outerfold] = meval.custom_spearmanr(
            Y_test_true_agg, Y_test_pred_agg, in_gridsearch=False
        )

        # Store scalers for future predictions/later use
        x_train_scaler[outerfold] = scaler_x
        y_train_scaler[outerfold] = scaler_y
        model_learned[outerfold] = grid

    return (
        cv_outer_splits,
        total_tile_selection,
        model_learned,
        x_train_scaler,
        y_train_scaler,
        slides_spearman_train,
        slides_spearman_test,
        tiles_spearman_train,
        tiles_spearman_test,
    )


def main(args):
    (
        cv_outer_splits,
        total_tile_selection,
        model_learned,
        x_train_scaler,
        y_train_scaler,
        slides_spearman_train,
        slides_spearman_test,
        tiles_spearman_train,
        tiles_spearman_test,
    ) = nested_cv_multitask(
        bottleneck_features_path=args.bottleneck_features_path,
        category=args.category,
        alpha_min=args.alpha_min,
        alpha_max=args.alpha_max,
        n_steps=args.n_steps,
        n_outerfolds=args.n_outerfolds,
        n_innerfolds=args.n_innerfolds,
        n_tiles=args.n_tiles,
        split_level=args.split_level,
        slide_type=args.slide_type,
        n_jobs=args.n_cores,
        var_names_path=args.var_names_path,
        target_features_path=args.target_features,
    )

    # Store for reproduction of outer scores
    joblib.dump(cv_outer_splits, Path(args.output_dir, "cv_outer_splits.pkl"))
    joblib.dump(total_tile_selection, Path(args.output_dir, "total_tile_selection.pkl"))
    joblib.dump(model_learned, Path(args.output_dir, "outer_models.pkl"))
    joblib.dump(x_train_scaler, Path(args.output_dir, "x_train_scaler.pkl"))
    joblib.dump(y_train_scaler, Path(args.output_dir, "y_train_scaler.pkl"))
    joblib.dump(
        slides_spearman_train, Path(args.output_dir, "outer_scores_slides_train.pkl")
    )
    joblib.dump(
        slides_spearman_test, Path(args.output_dir, "outer_scores_slides_test.pkl")
    )
    joblib.dump(
        tiles_spearman_train, Path(args.output_dir, "outer_scores_tiles_train.pkl")
    )
    joblib.dump(
        tiles_spearman_test, Path(args.output_dir, "outer_scores_tiles_test.pkl")
    )


if __name__ == "__main__":
    args = get_args()
    st = time.time()
    main(args)
    rt = time.time() - st
    print(f"Script finished in {rt // 60:.0f}m {rt % 60:.0f}s")
