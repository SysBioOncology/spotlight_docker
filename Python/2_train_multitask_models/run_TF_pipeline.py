# Module imports
import argparse
import os
import sys
import dask.dataframe as dd
import joblib
import git
import numpy as np
import pandas as pd
from sklearn import linear_model, metrics
from sklearn.model_selection import GridSearchCV, GroupKFold
from sklearn.preprocessing import StandardScaler
REPO_DIR = git.Repo(os.getcwd(), search_parent_directories=True).working_tree_dir
sys.path.append(f"{REPO_DIR}/Python/libs")

# Custom imports
import model.evaluate as meval
import model.preprocessing as preprocessing
import model.utils as utils

def nested_cv_multitask(
    output_dir,
    category,
    alpha_min,
    alpha_max,
    n_steps=40,
    n_outerfolds=5,
    n_innerfolds=10,
    n_tiles=50,
    split_level="sample_submitter_id",
    slide_type="FF"
):
    """
    Transfer Learning to quantify the cell types on a tile-level
    Use a nested cross-validation strategy to train a multi-task lasso algorithm. Tuning and evaluation based on spearman correlation.

    Args:
        output_dir (str): Path pointing to folder where models will be stored
        category (str): cell type
        alpha_min (int): Min. value of hyperparameter alpha
        alpha_max (int): Max. value of hyperparameter alpha
        n_steps (int): Stepsize for grid [alpha_min, alpha_max]
	slide_type (str): slide format (FF or FFPE)
        n_outerfolds (int): Number of outer loops
        n_innerfolds (int): Number of inner loops
        n_tiles (int): Number of tiles to select per slide
        split_level (str): Split level of slides for creating splits

    Returns:
        {output_dir}/: Pickle files containing the created splits, selected tiles, learned models, scalers, and evaluation of the slides and tiles using the spearman correlation for both train and test sets
    """
    # Hyperparameter grid for tuning
    alphas = np.logspace(int(alpha_min), int(alpha_max), int(n_steps))
    scoring = meval.custom_spearmanr
    N_JOBS = -1
    OUTPUT_PATH = f"{output_dir}/2_TF_training/models/{category}"

    print(slide_type)

    # Load data
    var_names = joblib.load("./task_selection_names.pkl")
    var_names['T_cells'] = ['Cytotoxic cells', 'Effector cells', 'CD8 T cells (quanTIseq)', 'Immune score']
    target_features = pd.read_csv(f"{output_dir}/2_TF_training/ensembled_selected_tasks.csv", sep="\t", index_col=0)
    if slide_type == "FF":
        bottleneck_features = pd.read_csv(f"{output_dir}/1_histopathological_features/features.txt", sep="\t", index_col=0)
    elif slide_type == "FFPE":
        bottleneck_features = dd.read_parquet(f"{output_dir}/1_histopathological_features/features.parquet")

    target_vars = var_names[category]
    metadata_colnames = var_names["tile_IDs"] + var_names["IDs"] + var_names[category]

    if os.path.exists(OUTPUT_PATH):
        print("Folder exists")
    else:
        os.makedirs(OUTPUT_PATH)

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
    ## Create variables for storing
    model_learned = dict.fromkeys(range(n_outerfolds))
    x_train_scaler = dict.fromkeys(range(n_outerfolds))
    y_train_scaler = dict.fromkeys(range(n_outerfolds))

    multi_task_lasso = linear_model.MultiTaskLasso()

    ## Setup nested cv
    sample_id = pd.factorize(total_tile_selection[split_level])[0]
    cv_outer = GroupKFold(n_splits=n_outerfolds)
    cv_inner = GroupKFold(n_splits=n_innerfolds)
    cv_outer_splits = list(cv_outer.split(X, Y, groups=sample_id))

    ## Storing scores
    tiles_spearman_train = {}
    tiles_spearman_test = {}
    slides_spearman_train = {}
    slides_spearman_test = {}

    print("Feature matrix dimensions [tiles, features]:", X.shape)
    print("Response matrix dimensions:", Y.shape)

    ## Run nested cross-validation
    for outerfold in range(n_outerfolds):
        print(f"Outerfold {outerfold}")
        train_index, test_index = cv_outer_splits[outerfold]
        x_train, x_test = X.iloc[train_index], X.iloc[test_index]
        y_train, y_test = Y.iloc[train_index], Y.iloc[test_index]

        ### Standardizing predictors
        scaler_x = StandardScaler()
        scaler_x.fit(x_train)
        x_train_z = scaler_x.transform(x_train)
        x_test_z = scaler_x.transform(x_test)

        ### Standardizing targets
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
            n_jobs=N_JOBS,
        )
        grid.fit(x_train_z, y_train_z, groups=sample_id[train_index])

        ### Evaluate on tile level (spearmanr)
        y_train_z_pred = grid.predict(x_train_z)
        y_test_z_pred = grid.predict(x_test_z)

        tiles_spearman_train[outerfold] = meval.custom_spearmanr(
            y_train_z, y_train_z_pred, in_gridsearch=False
        )
        tiles_spearman_test[outerfold] = meval.custom_spearmanr(
            y_test_z, y_test_z_pred, in_gridsearch=False
        )

        ### Evaluate on Slide level
        #### For aggregation for evaluation in gridsearch for hyper parameter choosing
        slide_IDs_train = total_tile_selection["slide_submitter_id"].iloc[train_index]
        slide_IDs_test = total_tile_selection["slide_submitter_id"].iloc[test_index]

        ##### 1. Aggregate on slide level (averaging)
        Y_train_true_agg, Y_train_pred_agg = meval.compute_aggregated_scores(
            y_train_z, y_train_z_pred, target_vars, slide_IDs_train
        )
        Y_test_true_agg, Y_test_pred_agg = meval.compute_aggregated_scores(
            y_test_z, y_test_z_pred, target_vars, slide_IDs_test
        )

        ###### 2. Compute spearman correlation between ground truth and predictions
        slides_spearman_train[outerfold] = meval.custom_spearmanr(
            Y_train_true_agg, Y_train_pred_agg, in_gridsearch=False
        )
        slides_spearman_test[outerfold] = meval.custom_spearmanr(
            Y_test_true_agg, Y_test_pred_agg, in_gridsearch=False
        )

        ###### Store scalers for future predictions/later use
        x_train_scaler[outerfold] = scaler_x
        y_train_scaler[outerfold] = scaler_y
        model_learned[outerfold] = grid

    ## Store for reproduction of outer scores
    joblib.dump(cv_outer_splits, f"{OUTPUT_PATH}/cv_outer_splits.pkl")
    joblib.dump(total_tile_selection, f"{OUTPUT_PATH}/total_tile_selection.pkl")
    joblib.dump(model_learned, f"{OUTPUT_PATH}/outer_models.pkl")
    joblib.dump(x_train_scaler, f"{OUTPUT_PATH}/x_train_scaler.pkl")
    joblib.dump(y_train_scaler, f"{OUTPUT_PATH}/y_train_scaler.pkl")
    joblib.dump(slides_spearman_train, f"{OUTPUT_PATH}/outer_scores_slides_train.pkl")
    joblib.dump(slides_spearman_test, f"{OUTPUT_PATH}/outer_scores_slides_test.pkl")
    joblib.dump(tiles_spearman_train, f"{OUTPUT_PATH}/outer_scores_tiles_train.pkl")
    joblib.dump(tiles_spearman_test, f"{OUTPUT_PATH}/outer_scores_tiles_test.pkl")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--path_target_features", help="Path pointing to file containing the cell type abundances",
    )
    parser.add_argument(
        "--output_dir", help="Path pointing to folder where models will be stored"
    )
    parser.add_argument(
        "--path_bottleneck_features", help="Path pointing to file containing the histopathological features",
    )
    parser.add_argument("--category", help="Cell type")
    parser.add_argument("--alpha_min", help="Min. value of hyperparameter alpha")
    parser.add_argument("--alpha_max", help="Max. value of hyperparameter alpha")
    parser.add_argument("--n_steps", help="Stepsize for grid [alpha_min, alpha_max]", default=40)
    parser.add_argument("--n_innerfolds", help="Number of inner loops", default=10)
    parser.add_argument("--n_outerfolds", help="Number of outer loops", default=5)
    parser.add_argument("--n_tiles", help="Number of tiles to select per slide", default=50)
    parser.add_argument("--split_level", help="Split level of slides for creating splits",
        default="sample_submitter_id",
    )
    parser.add_argument("--slide_type", help="Type of tissue slide (FF or FFPE)]")

    args = parser.parse_args()
    nested_cv_multitask(
        output_dir=args.output_dir,
        category=args.category,
        alpha_min=args.alpha_min,
        alpha_max=args.alpha_max,
        n_steps=args.n_steps,
        n_outerfolds=args.n_outerfolds,
        n_innerfolds=args.n_innerfolds,
        n_tiles=args.n_tiles,
        split_level=args.split_level,
	slide_type=args.slide_type
    )
