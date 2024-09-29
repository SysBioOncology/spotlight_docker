#!/usr/bin/env python3
# Module imports
import os
import argparse
import joblib
import numpy as np
import pandas as pd

import model.preprocessing as preprocessing
from model.constants import (
    TUMOR_PURITY,
    T_CELLS,
    ENDOTHELIAL_CELLS,
    CAFS,
    IDS,
    TILE_VARS,
)

from os.path import abspath
import time
from pathlib import Path
import multiprocessing
from argparse import ArgumentParser as AP


def get_args():
    # Script description
    description = """Preprocessing multi-task model target features"""

    # Add parser
    parser = AP(
        description=description, formatter_class=argparse.RawDescriptionHelpFormatter
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
    )
    parser.add_argument("--prefix", type=str, help="Prefix for output file", default="")
    parser.add_argument(
        "--n_cores", type=int, help="Number of cores to use (parallelization)"
    )

    parser.add_argument(
        "--cancer_type",
        help="Abbreviation of cancer type for naming of generated files",
    )
    parser.add_argument(
        "--clinical_file_path", help="Full path to clinical file", default=None
    )

    parser.add_argument("--tpm_path", help="Path to tpm file", type=str, required=True)

    parser.add_argument(
        "--output_dir", help="Path to folder for generated file", default=""
    )

    parser.add_argument("--version", action="version", version="0.1.0")
    arg = parser.parse_args()
    arg.output_dir = abspath(arg.output_dir)

    if (arg.output_dir != "") & (not os.path.isdir(arg.output_dir)):
        os.mkdir(arg.output_dir)

    if arg.n_cores is None:
        arg.n_cores = multiprocessing.cpu_count()
    return arg


def processing_transcriptomics(
    clinical_file_path: str,
    tpm_path: str,
    thorsson_signatures_path: str = "Thorsson_Scores_160_Signatures.tsv",
    estimate_signatures_path: str = "Yoshihara_ESTIMATE_XXX_RNAseqV2.txt",
    absolute_tumor_purity_path: str = "ABSOLUTE_tumor_purity.txt",
    gibbons_signatures_path: str = "Gibbons_supp1.xlsx",
    mcp_counter_path: str = "mcp_counter.csv",
    quantiseq_path: str = "quantiseq.csv",
    xcell_path: str = "xcell.csv",
    epic_path: str = "epic.csv",
):
    """Compute and combine cell type abundances from different quantification methods necessary for TF learning

    Args:
        cancer_type (str): abbreviation of cancer type
        clinical_file_path (str): path to clinical file
        tpm_path (str): path pointing to the tpm file (expression matrix)
        slide_type (str, optional): Type of slide, either FF or FFPE. Defaults to "FF".
        path_data (_type_, optional): _description_. Defaults to None.
        thorsson_signatures_path (str, optional): Path to file with Thorsson signatureds. Defaults to "Thorsson_Scores_160_Signatures.tsv".
        estimate_signatures_path (str, optional): Path to file with ESTIMATE signatures. Defaults to "Yoshihara_ESTIMATE_XXX_RNAseqV2.txt".
        absolute_tumor_purity_path (str, optional): Path to file with ABSOLUTE tumor purity. Defaults to "ABSOLUTE_tumor_purity.txt".
        gibbons_signatures_path (str, optional): Path to Gibbons signatures. Defaults to "Gibbons_supp1.xlsx".
        mcp_counter_path (str, optional): Path to results of MCP counter. Defaults to "mcp_counter.csv".
        quantiseq_path (str, optional): Path to results of Quantiseq. Defaults to "quantiseq.csv".
        xcell_path (str, optional): Path to results of Xcell. Defaults to "xcell.csv".
        epic_path (str, optional): path to results of EPIC. Defaults to "epic.csv".

    Returns:
        ./task_selection_names.pkl: pickle file containing variable names.
            {output_dir}/TCGA_{cancer_type}_ensembled_selected_tasks.csv" containing the following cell type quantification methods:
                tumor_purity = [
                    'tumor purity (ABSOLUTE)',
                    'tumor purity (estimate)',
                    'tumor purity (EPIC)'
                ]
                T_cells = [
                    'CD8 T cells (Thorsson)',
                    'Cytotoxic cells',
                    'Effector cells',
                    'CD8 T cells (quanTIseq)',
                    'TIL score',
                    'Immune score',
                ]
                endothelial_cells = [
                    'Endothelial cells (xCell)',
                    'Endothelial cells (EPIC)',
                    'Endothelium', ]
                CAFs = [
                    'Stromal score',
                    'CAFs (MCP counter)',
                    'CAFs (EPIC)',
                    'CAFs (Bagaev)',
                ]
    """

    var_dict = {
        "CAFs": CAFS,
        "T_cells": T_CELLS,
        "tumor_purity": TUMOR_PURITY,
        "endothelial_cells": ENDOTHELIAL_CELLS,
        "IDs": IDS,
        "tile_IDs": TILE_VARS,
    }
    joblib.dump(var_dict, "./task_selection_names.pkl")
    clinical_file = pd.read_csv(clinical_file_path, sep="\t")

    var_IDs = ["sample_submitter_id", "slide_submitter_id"]
    all_slide_features = clinical_file.loc[:, var_IDs]

    # Published Data
    Thorsson = pd.read_csv(thorsson_signatures_path, sep="\t")
    estimate = pd.read_csv(
        estimate_signatures_path,
        sep="\t",
    )
    tcga_absolute = pd.read_csv(absolute_tumor_purity_path, sep="\t")
    gibbons = pd.read_excel(
        gibbons_signatures_path,
        skiprows=2,
        sheet_name="DataFileS1 - immune features",
    )

    # Computed Data: Immunedeconv
    mcp_counter = pd.read_csv(mcp_counter_path, index_col=0, sep=",")
    quantiseq = pd.read_csv(quantiseq_path, index_col=0, sep=",")
    xCell = pd.read_csv(xcell_path, index_col=0, sep=",", header=[0])
    EPIC = pd.read_csv(epic_path, index_col=0, sep=",")

    # Re(compute) Fges scores with TPM
    Fges_computed = preprocessing.compute_gene_signature_scores(tpm_path)
    Fges_computed = Fges_computed.loc[:, ["Effector_cells", "Endothelium", "CAF"]]
    Fges_computed.columns = ["Effector cells", "Endothelium", "CAFs (Bagaev)"]

    Fges_computed = Fges_computed.reset_index()
    Fges_computed = Fges_computed.rename(columns={"index": "TCGA_sample"})

    # From immunedeconv
    quantiseq = preprocessing.process_immunedeconv(quantiseq, "quanTIseq")
    EPIC = preprocessing.process_immunedeconv(EPIC, "EPIC")
    mcp_counter = preprocessing.process_immunedeconv(mcp_counter, "MCP")
    xCell = preprocessing.process_immunedeconv(xCell, "xCell")

    # Merge cell fractions
    cellfrac = pd.merge(xCell, quantiseq, on=["TCGA_sample"])
    cellfrac = pd.merge(cellfrac, mcp_counter, on=["TCGA_sample"])
    cellfrac = pd.merge(cellfrac, EPIC, on=["TCGA_sample"])

    # Merge cell fractions
    cellfrac = pd.merge(xCell, quantiseq, on=["TCGA_sample"])
    cellfrac = pd.merge(cellfrac, mcp_counter, on=["TCGA_sample"])
    cellfrac = pd.merge(cellfrac, EPIC, on=["TCGA_sample"])

    # estimate data
    estimate = estimate.rename(columns={"ID": "TCGA_sample"})
    estimate = estimate.set_index("TCGA_sample")
    estimate.columns = ["Stromal score", "Immune score", "ESTIMATE score"]

    # According the tumor purity formula provided in the paper
    estimate["tumor purity (ESTIMATE)"] = np.cos(
        0.6049872018 + 0.0001467884 * estimate["ESTIMATE score"]
    )
    estimate = estimate.drop(columns=["ESTIMATE score"])

    # Thorsson data
    Thorsson = Thorsson.drop(columns="Source")
    Thorsson = Thorsson.set_index("SetName").T
    Thorsson = Thorsson.rename_axis(None, axis=1)
    Thorsson.index.name = "TCGA_aliquot"
    Thorsson = Thorsson.loc[:, ["LIexpression_score", "CD8_PCA_16704732"]]
    Thorsson.columns = ["TIL score", "CD8 T cells (Thorsson)"]

    # TCGA PanCanAtlas
    tcga_absolute = tcga_absolute.rename(
        columns={"purity": "tumor purity (ABSOLUTE)", "sample": "TCGA_aliquot"}
    )
    tcga_absolute = tcga_absolute.set_index("TCGA_aliquot")
    tcga_absolute = pd.DataFrame(tcga_absolute.loc[:, "tumor purity (ABSOLUTE)"])

    gibbons = gibbons.rename(columns={"Unnamed: 1": "id"})
    gibbons["slide_submitter_id"] = gibbons["id"].str[0:23]
    gibbons["Cytotoxic cells"] = gibbons["Cytotoxic cells"].astype(float)
    gibbons = gibbons.set_index("slide_submitter_id")

    all_slide_features["TCGA_sample"] = clinical_file["slide_submitter_id"].str[0:15]

    # add IDs
    Thorsson["TCGA_sample"] = Thorsson.index.str[0:15]
    tcga_absolute["TCGA_sample"] = tcga_absolute.index.str[0:15]
    gibbons["TCGA_sample"] = gibbons.index.str[0:15]

    tcga_absolute_merged = pd.merge(
        all_slide_features,
        tcga_absolute,
        on=[
            "TCGA_sample",
        ],
        how="left",
    )
    Thorsson_merged = pd.merge(
        all_slide_features,
        Thorsson,
        on=[
            "TCGA_sample",
        ],
        how="left",
    )
    gibbons_merged = pd.merge(
        all_slide_features, gibbons, on=["TCGA_sample"], how="left"
    )

    cellfrac_merged = pd.merge(
        all_slide_features, cellfrac, on=["TCGA_sample"], how="left"
    )
    estimate_merged = pd.merge(
        all_slide_features, estimate, on=["TCGA_sample"], how="left"
    )
    Fges_computed_merged = pd.merge(
        all_slide_features, Fges_computed, on=["TCGA_sample"], how="left"
    )

    # Combine in one dataframe
    all_merged = pd.merge(all_slide_features, tcga_absolute_merged, how="left")
    all_merged = pd.merge(all_merged, Thorsson_merged, how="left")
    all_merged = pd.merge(all_merged, gibbons_merged, how="left")
    all_merged = pd.merge(all_merged, estimate_merged, how="left")
    all_merged = pd.merge(all_merged, cellfrac_merged, how="left")
    all_merged = pd.merge(all_merged, Fges_computed_merged, how="left")

    # ---- Transform features to get a normal distribution (immunedeconv) ---- #
    featuresnames_transform = [
        "CAFs (MCP counter)",
        "CAFs (EPIC)",
    ]
    feature_data = all_merged.loc[:, CAFS].astype(float)
    data_log2_transformed = feature_data.copy()
    data_log2_transformed[featuresnames_transform] = np.log2(
        feature_data[featuresnames_transform] * 100 + 0.001
    )
    CAFs_transformed = data_log2_transformed

    featuresnames_transform = [
        "Endothelial cells (xCell)",
        "Endothelial cells (EPIC)",
    ]
    feature_data = all_merged.loc[:, ENDOTHELIAL_CELLS].astype(float)
    data_log2_transformed = feature_data.copy()
    data_log2_transformed[featuresnames_transform] = np.log2(
        feature_data[featuresnames_transform] * 100 + 0.001
    )
    endothelial_cells_transformed = data_log2_transformed

    feature_data = all_merged.loc[:, T_CELLS].astype(float)
    featuresnames_transform = ["CD8 T cells (quanTIseq)"]
    data_log2_transformed = feature_data.copy()
    data_log2_transformed[featuresnames_transform] = np.log2(
        feature_data[featuresnames_transform] * 100 + 0.001
    )
    T_cells_transformed = data_log2_transformed

    feature_data = all_merged.loc[:, TUMOR_PURITY].astype(float)
    featuresnames_transform = ["tumor purity (EPIC)"]
    data_log2_transformed = feature_data.copy()
    data_log2_transformed[featuresnames_transform] = np.log2(
        feature_data[featuresnames_transform] * 100 + 0.001
    )
    tumor_cells_transformed = data_log2_transformed

    # Store processed data
    IDs = ["slide_submitter_id", "sample_submitter_id", "TCGA_sample"]
    metadata = all_merged[IDs]
    merged = pd.concat(
        [
            metadata,
            CAFs_transformed,
            endothelial_cells_transformed,
            T_cells_transformed,
            tumor_cells_transformed,
        ],
        axis=1,
    )
    merged = merged.fillna(np.nan)

    # Remove slides if there are no values at all
    merged = merged.dropna(
        axis=0, subset=T_CELLS + CAFS + ENDOTHELIAL_CELLS + TUMOR_PURITY, how="all"
    )
    return (merged, var_dict)


def main(args):
    tasks, var_dict = processing_transcriptomics(
        clinical_file_path=args.clinical_file_path,
        tpm_path=args.tpm_path,
        thorsson_signatures_path=args.thorsson_signatures_path,
        estimate_signatures_path=args.estimate_signatures_path,
        absolute_tumor_purity_path=args.absolute_tumor_purity_path,
        gibbons_signatures_path=args.gibbons_signatures_path,
        mcp_counter_path=args.mcp_counter_path,
        quantiseq_path=args.quantiseq_path,
        xcell_path=args.xcell_path,
        epic_path=args.epic_path,
    )
    tasks.to_csv(Path(args.output_dir, "ensembled_selected_tasks.csv"), sep="\t")
    joblib.dump(var_dict, Path(args.output_dir, "task_selection_names.pkl"))


if __name__ == "__main__":
    args = get_args()
    st = time.time()
    main(args)
    rt = time.time() - st
    print(f"Script finished in {rt // 60:.0f}m {rt % 60:.0f}s")
