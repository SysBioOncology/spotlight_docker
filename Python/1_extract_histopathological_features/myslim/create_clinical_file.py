import argparse
import os
import os.path
import numpy as np
import pandas as pd
import sys


def create_TCGA_clinical_file(
    class_names,
    clinical_files_dir,
    output_dir=None,
    tumor_purity_threshold=80,
    path_codebook=None
):
    """
    Create a clinical file based on the slide metadata downloaded from the GDC data portal
    1. Read the files and add classname and id based on CODEBOOK.txt
    2. Filter tumor purity
    3. Save file

    Args:
        class_names (str): single class name e.g. LUAD_T or path to file with class names
        clinical_files_dir (str): String with path to folder with subfolders pointing to the raw clinical files (slide.tsv)
        output_dir (str): Path to folder where the clinical file should be stored
        tumor_purity_threshold (int): default=80
        multi_class_path (str): path to file with class names to be merged into one clinical file

    Returns:
        {output_dir}/generated_clinical_file.txt" containing the slide_submitter_id, sample_submitter_id, image_file_name, percent_tumor_cells, class_name, class_id in columns and records (slides) in rows.

    """
    # ---- Setup parameters ---- #
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    if (os.path.isfile(class_names)):  # multi class names
        class_names = pd.read_csv(
            class_names, header=None).to_numpy().flatten()
    else:  # single class names
        class_name = class_names

    CODEBOOK = pd.read_csv(
        path_codebook,
        delim_whitespace=True,
        header=None, names=["class_name", "value"]
    )

    # ---- 1. Constructing a merged clinical file ---- #
    # Read clinical files
    # a) Single class
    if os.path.isfile(clinical_files_dir):
        clinical_file = pd.read_csv(clinical_files_dir, sep="\t")
        # only keep tissue (remove _T or _N) to check in filename
        clinical_file["class_name"] = class_name
        clinical_file["class_id"] = int(
            CODEBOOK.loc[CODEBOOK["class_name"] == class_name].values[0][1]
        )
        print(clinical_file)
        print(CODEBOOK)
    # b) Multiple classes
    elif os.path.isdir(clinical_files_dir) & (len(class_names) > 1):
        clinical_file_list = []
        # Combine all clinical raw files based on input
        for class_name in class_names:
            clinical_file_temp = pd.read_csv(
                f"{clinical_files_dir}/clinical_file_TCGA_{class_name[:-2]}.tsv",
                sep="\t",
            )
            # only keep tissue (remove _T or _N) to check in filename
            clinical_file_temp["class_name"] = class_name
            clinical_file_temp["class_id"] = int(
                CODEBOOK.loc[CODEBOOK["class_name"] == class_name].values[0][1]
            )
            clinical_file_list.append(clinical_file_temp)
        clinical_file = pd.concat(
            clinical_file_list, axis=0).reset_index(drop=True)

    # ---- 2) Filter: Availability of tumor purity (percent_tumor_cells) ---- #
    # Remove rows with missing tumor purity
    clinical_file["percent_tumor_cells"] = (
        clinical_file["percent_tumor_cells"]
        .replace("'--", np.nan, regex=True)
        .astype(float)
    )

    # Convert strings to numeric type
    clinical_file["percent_tumor_cells"] = pd.to_numeric(
        clinical_file["percent_tumor_cells"]
    )
    clinical_file = clinical_file.dropna(subset=["percent_tumor_cells"])
    clinical_file = clinical_file.where(
        clinical_file["percent_tumor_cells"] >= float(tumor_purity_threshold)
    )
    # ---- 3) Formatting and saving ---- #
    clinical_file["image_file_name"] = [
        f"{slide_submitter_id}.{str(slide_id).upper()}.svs"
        for slide_submitter_id, slide_id in clinical_file[
            ["slide_submitter_id", "slide_id"]
        ].to_numpy()
    ]

    clinical_file = clinical_file.dropna(how="all")
    clinical_file = clinical_file.drop_duplicates()
    clinical_file = clinical_file.drop_duplicates(subset="slide_submitter_id")
    clinical_file = clinical_file[
        [
            "slide_submitter_id",
            "sample_submitter_id",
            "image_file_name",
            "percent_tumor_cells",
            "class_name",
            "class_id",
        ]
    ]
    clinical_file = clinical_file.dropna(how="any", axis=0)
    clinical_file.to_csv(
        f"{output_dir}/generated_clinical_file.txt",
        index=False,
        sep="\t",
    )
    print("\nFinished creating a new clinical file")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--class_names",
        help="Either (a) single classname or (b) Path to file with classnames according to codebook.txt (e.g. LUAD_T)", required=True
    )
    parser.add_argument(
        "--clinical_files_dir",
        help="Path to folders containing subfolders for multiple tumor types.", required=True
    )
    parser.add_argument(
        "--tumor_purity_threshold",
        help="Integer for filtering tumor purity assessed by pathologists",
        default=80, required=False
    )
    parser.add_argument(
        "--output_dir", help="Path to folder for saving all created files", default=None, required=False
    )
    parser.add_argument(
        "--path_codebook", help="Path to codebook", default=None, required=False
    )
    args = parser.parse_args()

    create_TCGA_clinical_file(
        class_names=args.class_names,
        tumor_purity_threshold=args.tumor_purity_threshold,
        clinical_files_dir=args.clinical_files_dir,
        output_dir=args.output_dir,
        path_codebook=args.path_codebook
    )
