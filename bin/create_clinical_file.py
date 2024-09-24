#!/usr/bin/env python3
import argparse
from argparse import ArgumentParser as AP
import os
import os.path
import numpy as np
import pandas as pd

from os.path import abspath
import time
from pathlib import Path


def get_args():
    # Script description
    description = """Creating a clinical file for TCGA dataset(s)"""

    # Add parser
    parser = AP(description=description,
                formatter_class=argparse.RawDescriptionHelpFormatter)

    # Sections
    parser.add_argument(
        "--clinical_files_input",
        help="Path to either a folder for multiple cancer types or single txt file.", required=False,
        default=None
    )
    parser.add_argument("--out_file", type=str, required=False, default="generated_clinical_file.txt",
                        help="Output filename with .txt extension (default='generated_clinical_file.txt')")
    parser.add_argument(
        "--path_codebook",
        help="Path to codebook",
        default=None, required=False, type=str
    )
    parser.add_argument(
        "--output_dir", help="Path to folder for saving all created files", default="", required=False, type=str
    )
    parser.add_argument(
        "--class_name",
        help="Single classname or (b) Path to file with classnames according to codebook.txt (e.g. LUAD_T)", default=None,
        type=str
    )
    parser.add_argument("--class_names_path",
                             type=str,
                             help="Path to file with classnames according to codebook.txt",
                             default=None)
    parser.add_argument(
        "--tumor_purity_threshold",
        help="Integer for filtering tumor purity assessed by pathologists",
        default=80, required=False, type=int
    )

    parser.add_argument("--version", action="version", version="0.1.0")
    arg = parser.parse_args()
    arg.output_dir = abspath(arg.output_dir)

    if ((arg.output_dir != "") & (not os.path.isdir(arg.output_dir))):
        os.mkdir(arg.output_dir)
    return arg

def handle_single_class(input, class_name, codebook):
        input = pd.read_csv(input, sep="\t")
        # only keep tissue (remove _T or _N) to check in filename
        input["class_name"] = class_name
        input["class_id"] = int(
            codebook.loc[codebook["class_name"]
                            == class_name].values[0][1]
        )
        return (input)

def handle_multi_class(input, class_names, codebook):
    clinical_file_list = []
    # Combine all clinical raw files based on input
    for class_name in class_names:
        clinical_file_temp = pd.read_csv(
            f"{input}/clinical_file_TCGA_{class_name[:-2]}.tsv",
            sep="\t",
        )
        # only keep tissue (remove _T or _N) to check in filename
        clinical_file_temp["class_name"] = class_name
        clinical_file_temp["class_id"] = int(
            codebook.loc[codebook["class_name"]
                            == class_name].values[0][1]
        )
        clinical_file_list.append(clinical_file_temp)
    clinical_file = pd.concat(
        clinical_file_list, axis=0).reset_index(drop=True)
    return (clinical_file)


def filter_tumor_purity(df, threshold = 80):
    # ---- 2) Filter: Availability of tumor purity (percent_tumor_cells) ---- #
    # Remove rows with missing tumor purity
    df["percent_tumor_cells"] = (
        df["percent_tumor_cells"]
        .replace("'--", np.nan, regex=True)
        .astype(float)
    )

    # Convert strings to numeric type
    df["percent_tumor_cells"] = pd.to_numeric(
        df["percent_tumor_cells"]
    )
    df = df.dropna(subset=["percent_tumor_cells"])
    df = df.where(
        df["percent_tumor_cells"] >= float(
            threshold)
    )
    return(df)

def is_valid_class_name_input (input, codebook):
    res = None
    if (input is not None):
        if (input in codebook["class_name"].values):
            res = "single"
        elif(Path(input)):
            res = "multi"
    return(res)


def create_TCGA_clinical_file(
    class_name,
    class_names_path,
    clinical_files_input,
    tumor_purity_threshold=80,
    path_codebook=None
):
    """
    Create a clinical file based on the slide metadata downloaded from the GDC data portal
    1. Read the files and add classname and id based on codebook_df.txt
    2. Filter tumor purity
    3. Save file

    Args:
        class_names (str): single class name e.g. LUAD_T or path to file with class names
        clinical_files_input (str): String with path to folder with subfolders pointing to the raw clinical files (slide.tsv)
        tumor_purity_threshold (int): default=80
        multi_class_path (str): path to file with class names to be merged into one clinical file

    Returns:
        {output_dir}/generated_clinical_file.txt" containing the slide_submitter_id, sample_submitter_id, image_file_name, percent_tumor_cells, class_name, class_id in columns and records (slides) in rows.

    """
    codebook_df = pd.read_csv(
            path_codebook,
            delim_whitespace=True,
            header=None, names=["class_name", "value"]
        )
    init_check_single_class =is_valid_class_name_input(input = class_name, codebook= codebook_df)
    init_check_multi_class = is_valid_class_name_input(input = class_names_path, codebook = codebook_df)
    is_single_class = init_check_single_class == "single"
    is_multi_class = init_check_multi_class == "multi"

    passes_input_check = (is_single_class| is_multi_class) & (clinical_files_input is not None)

    if passes_input_check:
        if (is_multi_class):  # multi class names
            class_names = pd.read_csv(
                class_names_path, header=None).to_numpy().flatten()
            if os.path.isdir(clinical_files_input) & (len(class_names) > 1):
                clinical_file = handle_multi_class(input = clinical_files_input, class_names = class_names, codebook=codebook_df)
        elif (is_single_class):  # single class names
            # a) Single class
            if os.path.isfile(clinical_files_input):
                clinical_file = handle_single_class(input = clinical_files_input, class_name=class_name, codebook = codebook_df)

        clinical_file = filter_tumor_purity(df = clinical_file, threshold= tumor_purity_threshold)

        # ---- 3) Formatting ---- #
        clinical_file["image_file_name"] = [
                f"{slide_submitter_id}.{str(slide_id).upper()}.svs"
                for slide_submitter_id, slide_id in clinical_file[
                    ["slide_submitter_id", "slide_id"]
                ].to_numpy()
            ]

        clinical_file = clinical_file.dropna(how="all")
        clinical_file = clinical_file.drop_duplicates()
        clinical_file = clinical_file.drop_duplicates(
            subset="slide_submitter_id")
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
        return clinical_file


def main(args):
    # Generate clinical file
    clinical_file = create_TCGA_clinical_file(
        class_name=args.class_name,
        class_names_path = args.class_names_path,
        tumor_purity_threshold=args.tumor_purity_threshold,
        clinical_files_input=args.clinical_files_input,
        path_codebook=args.path_codebook,
    )
    # Save file
    clinical_file.to_csv(
        Path(args.output_dir, args.out_file),
        index=False,
        sep="\t",
    )


if __name__ == "__main__":
    args = get_args()
    st = time.time()
    main(args)
    rt = time.time() - st
    print(f"Script finished in {rt // 60:.0f}m {rt % 60:.0f}s")
