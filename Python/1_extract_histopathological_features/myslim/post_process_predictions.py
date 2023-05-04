# Module imports
import argparse
import os
import sys
import dask.dataframe as dd
import pandas as pd

# WINDOWS
# import git
#REPO_DIR= git.Repo('.', search_parent_directories=True).working_tree_dir
#sys.path.append(f"{REPO_DIR}\\Python\\libs")

# Point to folder with custom imports
sys.path.append(f"{os.path.dirname(os.getcwd())}/libs")

# Custom imports
import DL.utils as utils
import numpy as np


def post_process_predictions(output_dir, slide_type):
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
    path_codebook = f"{os.path.dirname(os.getcwd())}/1_extract_histopathological_features/codebook.txt"
    path_tissue_classes = f"{os.path.dirname(os.getcwd())}/1_extract_histopathological_features/tissue_classes.csv"
    codebook = pd.read_csv(path_codebook, delim_whitespace=True, header=None)
    codebook.columns = ["class_name", "class_id"]
    tissue_classes = pd.read_csv(path_tissue_classes, sep="\t")

    # Read predictions
    if slide_type == "FF":
        predictions_raw = pd.read_csv(output_dir + "/pred.train.txt", sep="\t", header=None)
        ## Extract tile name incl. coordinates from path
        tile_names = [utils.get_tile_name(tile_path) for tile_path in predictions_raw[0]]
        ## Create output dataframe for post-processed data
        predictions = pd.DataFrame(tile_names, columns=["tile_ID"])
        ## Get predicted probabilities for all 42 classes + rename columns
        pred_probabilities = predictions_raw.iloc[:, 2:]
        pred_probabilities.columns = codebook["class_id"]
        ## Get predicted and true class ids
        predictions["pred_class_id"] = pred_probabilities.idxmax(axis="columns")
        predictions["true_class_id"] = 41
        ## Get corresponding max probabilities to the predicted class
        predictions["pred_probability"] = pred_probabilities.max(axis=1)
        # Replace class id with class name
        predictions["true_class_name"] = predictions["true_class_id"].copy()
        predictions["pred_class_name"] = predictions["pred_class_id"].copy()
        found_class_ids = set(predictions["true_class_id"]).union(set(predictions["pred_class_id"]))
        for class_id in found_class_ids:
                predictions["true_class_name"].replace(
                    class_id, codebook["class_name"][class_id], inplace=True
                )
                predictions["pred_class_name"].replace(
                    class_id, codebook["class_name"][class_id], inplace=True
                )

        ## Define whether prediction was right
        predictions["is_correct_pred"] = (predictions["true_class_id"] == predictions["pred_class_id"])
        predictions["is_correct_pred"] = predictions["is_correct_pred"].replace(False, "F")
        predictions.is_correct_pred=predictions.is_correct_pred.astype(str)
        ## Get tumor and tissue ID
        temp = dd.DataFrame({"tumor_type": predictions["true_class_name"].str[:-2]})
        temp = pd.merge(temp, tissue_classes, on="tumor_type", how="left")
        ## Set of IDs for normal and tumor (because of using multiple classes)
        IDs_tumor = list(set(temp["ID_tumor"]))
        if list(set(temp.tumor_type.tolist()))[0] == 'SKCM':
            # Probability for predicting tumor and normal label (regardless of tumor [tissue] type)
            predictions["tumor_label_prob"] = np.nan
            predictions["normal_label_prob"] = np.nan
            for ID_tumor in IDs_tumor:
                vals = pred_probabilities.loc[temp["ID_tumor"] == ID_tumor, ID_tumor]
                predictions.loc[temp["ID_tumor"] == ID_tumor, "tumor_label_prob"] = vals

            predictions["is_correct_pred_label"] = np.nan
        else:
            IDs_normal = list(set(temp["ID_normal"]))
            # Probability for predicting tumor and normal label (regardless of tumor [tissue] type)
            predictions["tumor_label_prob"] = np.nan
            predictions["normal_label_prob"] = np.nan
            for ID_tumor in IDs_tumor:
                vals = pred_probabilities.loc[temp["ID_tumor"] == ID_tumor, ID_tumor]
                predictions.loc[temp["ID_tumor"] == ID_tumor, "tumor_label_prob"] = vals

            for ID_normal in IDs_normal:
                vals = pred_probabilities.loc[temp["ID_normal"] == ID_normal, ID_normal]
                predictions.loc[temp["ID_normal"] == ID_normal, "normal_label_prob"] = vals

            # Check if the correct label (tumor/normal) is predicted
            temp_probs = predictions[["tumor_label_prob", "normal_label_prob"]]
            is_normal_label_prob = (
                temp_probs["normal_label_prob"] > temp_probs["tumor_label_prob"]
            )
            is_tumor_label_prob = (
                temp_probs["normal_label_prob"] < temp_probs["tumor_label_prob"]
            )
            is_normal_label = predictions["true_class_name"].str.find("_N") != -1
            is_tumor_label = predictions["true_class_name"].str.find("_T") != -1

            is_normal = is_normal_label & is_normal_label_prob
            is_tumor = is_tumor_label & is_tumor_label_prob

            predictions["is_correct_pred_label"] = is_normal | is_tumor
            predictions["is_correct_pred_label"].replace(True, "T", inplace=True)
            predictions["is_correct_pred_label"].replace(False, "F", inplace=True)

        # Save features to .csv file
        predictions.to_csv(output_dir + "/predictions.txt", sep="\t")

    elif slide_type== "FFPE":
        predictions_raw = dd.read_csv(output_dir + "/pred.train.txt", sep = "\t", header=None)
        predictions_raw['tile_ID'] = predictions_raw.iloc[:,0].str[104:180]
        predictions_raw['tile_ID'] = predictions_raw['tile_ID'].str.replace(".jpg'","")
        predictions = predictions.map_partitions(lambda df: df.drop(columns=[0,1]))
        new_names=list(map(lambda x: str(x), codebook["class_id"]))
        new_names.append('tile_ID')
        predictions.columns = new_names
        predictions = predictions.map_partitions(lambda x: x.assign(pred_class_id=x.iloc[:,0:41].idxmax(axis="columns")))
        predictions["true_class_id"] = 41
        predictions = predictions.map_partitions(lambda x: x.assign(pred_probability=x.iloc[:,0:41].max(axis="columns")))
        predictions["true_class_name"] = predictions["true_class_id"].copy()
        predictions["pred_class_name"] = predictions["pred_class_id"].copy()
        predictions.pred_class_id=predictions.pred_class_id.astype(int)
        res = dict(zip(codebook.class_id, codebook.class_name))
        predictions = predictions.map_partitions(lambda x: x.assign(pred_class_name=x.loc[:,'pred_class_id'].replace(res)))
        predictions = predictions.map_partitions(lambda x: x.assign(true_class_name=x.loc[:,'true_class_id'].replace(res)))
        predictions["is_correct_pred"] = (predictions["true_class_id"] == predictions["pred_class_id"])
        predictions["is_correct_pred"] = predictions["is_correct_pred"].replace(False, "F")
        predictions.is_correct_pred=predictions.is_correct_pred.astype(str)
        temp = predictions.map_partitions(lambda x: x.assign(tumor_type=x["true_class_name"].str[:-2]))
        temp = temp.map_partitions(lambda x: pd.merge(x, tissue_classes, on="tumor_type", how="left"))
        if (temp['tumor_type'].compute() == 'SKCM').any():
            # Probability for predicting tumor and normal label (regardless of tumor [tissue] type)
            predictions["tumor_label_prob"] = np.nan
            predictions["normal_label_prob"] = np.nan
            predictions = predictions.map_partitions(lambda x: x.assign(tumor_label_prob=x.loc[:, '41']))
            predictions["is_correct_pred_label"] = np.nan
        else:
            # TO DO
            predictions["tumor_label_prob"] = np.nan
            predictions["normal_label_prob"] = np.nan
            #predictions = predictions.map_partitions(lambda x: x.assign(tumor_label_prob=x.loc[:, '41']))
            #predictions = predictions.map_partitions(lambda x: x.assign(tumor_label_prob=x.loc[:, '41']))

        # Save features using parquet
        name_function=lambda x: f"predictions-{x}.parquet"
        OUTPUT_PATH = f"{output_dir}/predictions_format_parquet"
        if os.path.exists(OUTPUT_PATH):
            print("Folder exists")
        else:
            os.makedirs(OUTPUT_PATH)

        predictions.to_parquet(path=OUTPUT_PATH, compression='gzip', name_function=name_function)

    print("Formatted all predicted tissue labels")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output_dir", help="Set output folder")
    parser.add_argument("--slide_type", help="Type of tissue slide (FF or FFPE)]")
    args = parser.parse_args()
    post_process_predictions(output_dir=args.output_dir, slide_type=args.slide_type)


#  $cur_dir/codebook.txt $cur_dir/tissue_classes.csv $output_dir
