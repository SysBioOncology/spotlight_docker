#!/usr/bin/env python3

import os
import pandas as pd
import dask.dataframe as dd
import argparse
import joblib
import scipy.stats as stats
from pathlib import Path
from model.constants import DEFAULT_CELL_TYPES
from model.evaluate import compute_tile_predictions
import time
from argparse import ArgumentParser as AP

def get_args():
    # Script description
    description = """Tile-level cell type quantification"""

    # Add parser
    parser = AP(description=description,
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--models_dir", type=str,
                        help="Path to models directory", required=True)
    parser.add_argument("--output_dir", type=str,
                        help="Path to output directory", required=False, default = "")
    parser.add_argument("--histopatho_features_dir", type=str,
                        help="Path to histopathological features file", required=False, default = "")
    parser.add_argument("--var_names_path", type=str,
                        help="Path to variable names pkl file", required=True)
    parser.add_argument("--features_input", type = str, default = None)
    parser.add_argument("--prediction_mode", type=str,
                        help="Choose prediction mode 'performance' or 'all' (default='all')", default="all", required=False)
    parser.add_argument("--n_outerfolds", type=int, default=5,
                        help="Number of outer folds (default=5)", required=False)
    parser.add_argument("--cell_types_path", type=str, default="",
                        help="List of cell types by default=['T_cells','CAFs',  'tumor_purity','endothelial_cells']", required=False)
    parser.add_argument(
        "--slide_type", help="Type of tissue slide (FF or FFPE)", type=str, required=True)

    arg = parser.parse_args()

    if (arg.features_input is None):
        if arg.slide_type == "FF":
            arg.features_input = Path(arg.histopatho_features_dir, "features.txt")

        elif arg.slide_type == "FFPE":
            arg.features_input = Path(arg.histopatho_features_dir, "features_format_parquet")

    if (not Path(arg.features_input).exists()):
            raise Exception("Invalid argument, please check `features_input` or `histopatho_features_dir`")

    if ((arg.output_dir != "") & (not os.path.isdir(arg.output_dir))):
        # Create an empty folder for TF records if folder doesn't exist
        os.mkdir(arg.output_dir)
    return arg


def tile_level_quantification(features_input, models_dir, var_names_path, histopatho_features_dir="", prediction_mode="all", n_outerfolds=5, cell_types_path="", slide_type="FF"):
    """
    Quantify the cell type abundances for the different tiles. Creates three files:
    (1) z-scores and
    (2) probability scores
    for the individual tasks and the combined score (incl. metadata)

    (3) z-scores for only the combined scores (incl. metadata)
        Args:
            prediction_mode: choose 1) "performance" [predict only test folds] or 2) "all" [predict all tiles]
            n_outerfolds (int): number of outer loops
            cell_types (list)
            models_dir (str): path pointing to a folder containing the subfolders (automatically created previously) of the trained models
            output_dir (str): path pointing to a folder where the predictions can be stored
            var_names_path (str): path pointing to the file containing the various variable names
            histopatho_features (str): path pointing to the file with the extracted histopathological features

    """
    # Read data
    cell_types = DEFAULT_CELL_TYPES
    if os.path.isfile(cell_types_path):
        cell_types = pd.read_csv(
            cell_types_path, header=None).to_numpy().flatten()
    print(cell_types)

    var_names = joblib.load(var_names_path)
    print(var_names)

    if slide_type == "FF":
        histopatho_features = pd.read_csv(
            features_input, sep="\t", index_col=0)
    elif slide_type == "FFPE":
        histopatho_features = dd.read_parquet(features_input)

    print(histopatho_features.head())

    # Compute predictions based on bottleneck features
    tile_predictions = pd.DataFrame()
    bottleneck_features = histopatho_features.loc[:, [
        str(i) for i in range(1536)]]
    bottleneck_features.index = histopatho_features.tile_ID
    var_names['IDs'] = 'sample_submitter_id'
    var_names['tile_IDs'] = ['Coord_X', 'Coord_Y', 'tile_ID']
    var_names['tile_IDs'].append(var_names['IDs'])
    metadata = histopatho_features.loc[:, var_names["tile_IDs"]]
    if slide_type == "FFPE":
        metadata = metadata.compute()

    print("Computing tile predictions for each cell type...")

    ##############################################################################
    # If predicting on all FFPE slides, we do this by chunks:
    # if any([prediction_mode == item for item in ['tcga_train_validation', 'test']]):
    #     if slide_type == "FFPE":
    #        tmp=metadata.slide_submitter_id.unique().tolist()
    #        ffpe_slides = dict.fromkeys(range(n_outerfolds))
    #        for ii in range(n_outerfolds):
    #            ffpe_slides[ii] = tmp[(ii*75):75*(ii+1)]
    #            if ii == 4:
    #                ffpe_slides[ii] = tmp[(ii*75):len(tmp)]
    #
    # for key, chunk in ffpe_slides.items():
    #    print('Chunk:', key)
    #    tiles_subset = metadata[metadata["slide_submitter_id"].isin(chunk)]["tile_ID"]
    #    X = bottleneck_features.map_partitions(lambda x: x[x.index.isin(tiles_subset)])
    #    X = X.compute()
    #    predictions = pd.DataFrame()
    #    for cell_type in cell_types:
    #        cell_type_tile_predictions = compute_tile_predictions(
    #            cell_type=cell_type, models_dir=models_dir, n_outerfolds=n_outerfolds,
    #            prediction_mode=prediction_mode, X=X, metadata=metadata,var_names=var_names, slide_type=slide_type
    #        )
    #        predictions = pd.concat([predictions, cell_type_tile_predictions], axis=1)
    #
    #    print(predictions.shape)
    #    tile_predictions = pd.concat([tile_predictions, predictions], axis=0)
    #
    ##############################################################################

    for cell_type in cell_types:
        cell_type_tile_predictions = compute_tile_predictions(
            cell_type=cell_type, models_dir=models_dir, n_outerfolds=n_outerfolds,
            prediction_mode=prediction_mode, X=bottleneck_features, metadata=metadata, var_names=var_names, slide_type="FF"
        )
        tile_predictions = pd.concat(
            [tile_predictions, cell_type_tile_predictions], axis=1)

    print(tile_predictions.shape)

    # Remove slides with nan values
    tile_predictions = tile_predictions.dropna()

    #  Order tile_predictions according to metadata
    metadata = metadata[metadata.tile_ID.isin(tile_predictions.index)]
    metadata.index = metadata.tile_ID
    tile_predictions = tile_predictions.loc[metadata.index, :]

    # Convert predictions to probabilities using cdf.
    feature_names = tile_predictions.columns
    pred_proba = pd.DataFrame(data=stats.norm.cdf(
        tile_predictions), columns=feature_names, index=tile_predictions.index)
    pred_proba = pd.concat([pred_proba, metadata], axis=1)

    # Remove suffix '(combi)'
    pred_proba.columns = [col.replace(" (combi)", "")
                          for col in pred_proba.columns]
    tile_predictions.columns = [col.replace(
        " (combi)", "") for col in tile_predictions.columns]

    #  Change columnsynta" "sample" for "slide"
    tile_predictions = tile_predictions.rename(
        columns={'sample_submitter_id': 'slide_submitter_id'})
    pred_proba = pred_proba.rename(
        columns={'sample_submitter_id': 'slide_submitter_id'})
    return (tile_predictions, pred_proba)


def main(args):
    tile_predictions, pred_proba = tile_level_quantification(
        features_input=args.features_input,
        models_dir=args.models_dir,
        histopatho_features_dir=args.histopatho_features_dir,
        prediction_mode=args.prediction_mode,
        n_outerfolds=args.n_outerfolds,
        cell_types_path=args.cell_types_path,
        var_names_path=args.var_names_path,
        slide_type=args.slide_type)


    tile_predictions.to_csv(
    Path(args.output_dir, f"{args.prediction_mode}_tile_predictions_zscores.csv"), sep="\t", index=False)
    pred_proba.to_csv(
        Path(args.output_dir, f"{args.prediction_mode}_tile_predictions_proba.csv"), sep="\t", index=False)
    print("Finished tile predictions...")

if __name__ == "__main__":
    args = get_args()
    st = time.time()
    main(args)
    rt = time.time() - st
    print(f"Script finished in {rt // 60:.0f}m {rt % 60:.0f}s")
