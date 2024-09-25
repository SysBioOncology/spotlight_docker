#!/usr/bin/env python

# Module imports
import argparse
from argparse import ArgumentParser as AP
import os
import dask.dataframe as dd
import pandas as pd

#  Custom imports
import DL.utils as utils
import numpy as np
from os.path import abspath
import time
from pathlib import Path

def get_args():
    # Script description
    description = """Post-processing predictions"""

    # Add parser
    parser = AP(description=description,
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output_dir", help="Set output folder", default = ".")
    parser.add_argument("--create_parquet_subdir", help = "Whether to create a subdirectory called 'predictions_format_parquet' if slide_type == 'FFPE', default=False", default = False)
    parser.add_argument(
        "--slide_type", help="Type of tissue slide (FF or FFPE) (default='FF')", type = str, default = "FF")
    parser.add_argument(
        "--path_codebook", help="codebook.txt file", required=True, type=str)
    parser.add_argument(
        "--path_tissue_classes", help="Tissue_classes.csv file", required=True, type=str)
    parser.add_argument("--cancer_type", help = "Cancer type", required = True, type =str)
    parser.add_argument("--pred_train_file", help = "", type = str, default = None)
    arg = parser.parse_args()

    if arg.pred_train_file is None:
        arg.pred_train_file = Path(arg.output_dir, "pred_train.txt")

    if (arg.create_parquet_subdir):
        arg.output_dir = abspath(Path(arg.output_dir, "predictions_format_parquet"))

    if not os.path.isdir(arg.output_dir):
        os.mkdir(arg.output_dir)

    return arg



def handle_ff_slides(pred_train_file, codebook, tissue_classes, cancer_type):
    predictions_raw = pd.read_csv(pred_train_file, sep="\t", header=None)
    # Extract tile name incl. coordinates from path
    tile_names = [utils.get_tile_name(tile_path)
                    for tile_path in predictions_raw[0]]
    # Create output dataframe for post-processed data
    predictions = pd.DataFrame(tile_names, columns=["tile_ID"])
    # Get predicted probabilities for all 42 classes + rename columns
    pred_probabilities = predictions_raw.iloc[:, 2:]
    pred_probabilities.columns = codebook["class_id"]
    # Get predicted and true class ids
    predictions["pred_class_id"] = pred_probabilities.idxmax(
        axis="columns")
    predictions["true_class_id"] = 41
    # Get corresponding max probabilities to the predicted class
    predictions["pred_probability"] = pred_probabilities.max(axis=1)
    # Replace class id with class name
    predictions["true_class_name"] = predictions["true_class_id"].copy()
    predictions["pred_class_name"] = predictions["pred_class_id"].copy()
    found_class_ids = set(predictions["true_class_id"]).union(
        set(predictions["pred_class_id"]))
    for class_id in found_class_ids:
        predictions["true_class_name"].replace(
            class_id, codebook["class_name"][class_id], inplace=True
        )
        predictions["pred_class_name"].replace(
            class_id, codebook["class_name"][class_id], inplace=True
        )

    # Define whether prediction was right
    predictions["is_correct_pred"] = (
        predictions["true_class_id"] == predictions["pred_class_id"])
    predictions["is_correct_pred"] = predictions["is_correct_pred"].replace(
        False, "F")
    predictions.is_correct_pred = predictions.is_correct_pred.astype(str)
    # Get tumor and tissue ID
    temp = pd.DataFrame(
        {"tumor_type": predictions["true_class_name"].str[:-2]})
    temp = pd.merge(temp, tissue_classes, on="tumor_type", how="left")
    # Set of IDs for normal and tumor (because of using multiple classes)
    IDs_tumor = list(set(temp["ID_tumor"]))
    if list(set(temp.tumor_type.tolist()))[0] == cancer_type:
        # Probability for predicting tumor and normal label (regardless of tumor [tissue] type)
        predictions["tumor_label_prob"] = np.nan
        predictions["normal_label_prob"] = np.nan
        for ID_tumor in IDs_tumor:
            vals = pred_probabilities.loc[temp["ID_tumor"]
                                            == ID_tumor, ID_tumor]
            predictions.loc[temp["ID_tumor"] ==
                            ID_tumor, "tumor_label_prob"] = vals

        predictions["is_correct_pred_label"] = np.nan
    else:
        IDs_normal = list(set(temp["ID_normal"]))
        # Probability for predicting tumor and normal label (regardless of tumor [tissue] type)
        predictions["tumor_label_prob"] = np.nan
        predictions["normal_label_prob"] = np.nan
        for ID_tumor in IDs_tumor:
            vals = pred_probabilities.loc[temp["ID_tumor"]
                                            == ID_tumor, ID_tumor]
            predictions.loc[temp["ID_tumor"] ==
                            ID_tumor, "tumor_label_prob"] = vals

        for ID_normal in IDs_normal:
            vals = pred_probabilities.loc[temp["ID_normal"]
                                            == ID_normal, ID_normal]
            predictions.loc[temp["ID_normal"] ==
                            ID_normal, "normal_label_prob"] = vals

        # Check if the correct label (tumor/normal) is predicted
        temp_probs = predictions[["tumor_label_prob", "normal_label_prob"]]
        is_normal_label_prob = (
            temp_probs["normal_label_prob"] > temp_probs["tumor_label_prob"]
        )
        is_tumor_label_prob = (
            temp_probs["normal_label_prob"] < temp_probs["tumor_label_prob"]
        )
        is_normal_label = predictions["true_class_name"].str.find(
            "_N") != -1
        is_tumor_label = predictions["true_class_name"].str.find(
            "_T") != -1

        is_normal = is_normal_label & is_normal_label_prob
        is_tumor = is_tumor_label & is_tumor_label_prob

        predictions["is_correct_pred_label"] = is_normal | is_tumor
        predictions["is_correct_pred_label"].replace(
            True, "T", inplace=True)
        predictions["is_correct_pred_label"].replace(
            False, "F", inplace=True)
    return(predictions)

def handle_ffpe_slides(pred_train_file, codebook, tissue_classes, cancer_type):
        predictions_raw = dd.read_csv(pred_train_file, sep="\t", header=None)
        predictions_raw['tile_ID'] = predictions_raw.iloc[:, 0]
        predictions_raw.tile_ID = predictions_raw.tile_ID.map(
            lambda x: x.split("/")[-1])
        predictions_raw['tile_ID'] = predictions_raw['tile_ID'].str.replace(
            ".jpg'", "")
        predictions = predictions_raw.map_partitions(
            lambda df: df.drop(columns=[0, 1]))
        new_names = list(map(lambda x: str(x), codebook["class_id"]))
        new_names.append('tile_ID')
        predictions.columns = new_names
        predictions = predictions.map_partitions(lambda x: x.assign(
            pred_class_id=x.iloc[:, 0:41].idxmax(axis="columns")))
        predictions["true_class_id"] = 41
        predictions = predictions.map_partitions(lambda x: x.assign(
            pred_probability=x.iloc[:, 0:41].max(axis="columns")))
        predictions["true_class_name"] = predictions["true_class_id"].copy()
        predictions["pred_class_name"] = predictions["pred_class_id"].copy()
        predictions.pred_class_id = predictions.pred_class_id.astype(int)
        res = dict(zip(codebook.class_id, codebook.class_name))
        predictions = predictions.map_partitions(lambda x: x.assign(
            pred_class_name=x.loc[:, 'pred_class_id'].replace(res)))
        predictions = predictions.map_partitions(lambda x: x.assign(
            true_class_name=x.loc[:, 'true_class_id'].replace(res)))
        predictions["is_correct_pred"] = (
            predictions["true_class_id"] == predictions["pred_class_id"])
        predictions["is_correct_pred"] = predictions["is_correct_pred"].replace(
            False, "F")
        predictions.is_correct_pred = predictions.is_correct_pred.astype(str)
        temp = predictions.map_partitions(lambda x: x.assign(
            tumor_type=x["true_class_name"].str[:-2]))
        temp = temp.map_partitions(lambda x: pd.merge(
            x, tissue_classes, on="tumor_type", how="left"))
        if (temp['tumor_type'].compute() == cancer_type).any():
            # Probability for predicting tumor and normal label (regardless of tumor [tissue] type)
            predictions["tumor_label_prob"] = np.nan
            predictions["normal_label_prob"] = np.nan
            predictions = predictions.map_partitions(
                lambda x: x.assign(tumor_label_prob=x.loc[:, '41']))
            predictions["is_correct_pred_label"] = np.nan
        else:
            # TO DO
            predictions["tumor_label_prob"] = np.nan
            predictions["normal_label_prob"] = np.nan
            # predictions = predictions.map_partitions(lambda x: x.assign(tumor_label_prob=x.loc[:, '41']))
            # predictions = predictions.map_partitions(lambda x: x.assign(tumor_label_prob=x.loc[:, '41']))
        return predictions

def post_process_predictions(pred_train_file, slide_type, path_codebook, path_tissue_classes, cancer_type):
    """
    Format predicted tissue classes and derive tumor purity from pred.train.txt file generated by myslim/bottleneck_predict.py and
    The pred.train.txt file contains the tile ID, the true class id and the 42 predicted probabilities for the 42 tissue classes.

    Args:
        output_dir (str): path pointing to folder for storing all created files by script

    Returns:
        {output_dir}/predictions.txt containing the following columns
        - tile_ID,
        - pred_class_id and true_class_id: class ids defined in codebook.txt)
        - pred_class_name and true_class_name: class names e.g. LUAD_T, defined in codebook.txt)
        - pred_probability: corresponding probability
        - is_correct_pred (boolean): correctly predicted tissue class label
        - tumor_label_prob and normal_label_prob: probability for predicting tumor and normal label (regardless of tumor or tissue type)
        - is_correct_pred_label (boolean): correctly predicted 'tumor' or 'normal' tissue regardless of tumor or tissue type
        In the rows the tiles.
    """

    # Initialize
    codebook = pd.read_csv(path_codebook, delim_whitespace=True, header=None)
    codebook.columns = ["class_name", "class_id"]
    tissue_classes = pd.read_csv(path_tissue_classes, sep="\t")

    # Read predictions
    if slide_type == "FF":
        return(handle_ff_slides(pred_train_file=pred_train_file, codebook=codebook, tissue_classes=tissue_classes, cancer_type = cancer_type))
    #  Save features to .csv file
    elif slide_type == "FFPE":
        return(handle_ffpe_slides(pred_train_file=pred_train_file, codebook=codebook, tissue_classes= tissue_classes, cancer_type=cancer_type))
    else:
        raise Exception("Invalid `slide_type`, please choose 'FF' or 'FFPE' ")

def main(args):
    predictions = post_process_predictions(pred_train_file = args.pred_train_file, slide_type=args.slide_type, path_codebook=args.path_codebook,
                             path_tissue_classes=args.path_tissue_classes, cancer_type=args.cancer_type)
    if (args.slide_type == "FF"):
            predictions.to_csv(Path(args.output_dir, "predictions.txt"), sep="\t")
    elif (args.slide_type == "FFPE"):
        # Save features using parquet
        def name_function(x): return f"predictions-{x}.parquet"
        predictions.to_parquet(
            path=args.output_dir, compression='gzip', name_function=name_function)
    print("Finished post-processing of predictions...")


if __name__ == "__main__":
    args = get_args()
    st = time.time()
    main(args)
    rt = time.time() - st
    print(f"Script finished in {rt // 60:.0f}m {rt % 60:.0f}s")
