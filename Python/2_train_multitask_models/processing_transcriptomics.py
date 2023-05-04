# Module imports
import os
import sys
import argparse
import joblib

import numpy as np
import pandas as pd
import git
REPO_DIR= git.Repo('.', search_parent_directories=True).working_tree_dir
sys.path.append(f"{REPO_DIR}/Python/libs")
import model.preprocessing as preprocessing
from model.constants import TUMOR_PURITY, T_CELLS, ENDOTHELIAL_CELLS, CAFS, IDS, TILE_VARS

# cancer_type="SKCM"
# slide_type="FFPE"
# # clinical_file = pd.read_csv("../../data/SKCM/slide.tsv", sep="\t")
# clinical_file = pd.read_csv("../../data/FFPE_generated_clinical_file.txt", sep="\t")

# # Set paths: 1) folder with all published data, 2) folder with all computed data, 3) folder for storing the ensembled tasks
# path_published_data = "../../data/published"
# path_computed_features= "../../data/"
# output_dir = "../../data/"

def processing_transcriptomics(cancer_type, slide_type, clinical_file_path, tpm_path, output_dir,  path_data=None,):
    """ Compute and combine cell type abundances from different quantification methods necessary for TF learning
    Args:
        cancer_type (str): abbreviation of cancer_type
        slide_type (str): type of slide either 'FFPE' or 'FF' for naming and necessary for merging data
        clinical_file_path (str): clinical_file_path: path to clinical_file
        tpm_path (str): path pointing to the tpm file
        output_dir (str): path pointing to a folder where the dataframe containing all features should be stored, stored as .txt file


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
    full_output_dir = f"{output_dir}/2_TF_training"
    if not os.path.exists(full_output_dir):
        os.makedirs(full_output_dir)

    var_dict = {
        "CAFs": CAFS,
        "T_cells": T_CELLS,
        "tumor_purity": TUMOR_PURITY,
        "endothelial_cells": ENDOTHELIAL_CELLS,
        "IDs":IDS,
        "tile_IDs": TILE_VARS
    }
    joblib.dump(var_dict, "./task_selection_names.pkl")
    clinical_file = pd.read_csv(clinical_file_path, sep="\t")

    var_IDs = ['sample_submitter_id','slide_submitter_id']
    all_slide_features = clinical_file.loc[:,var_IDs]

    # Published Data
    Thorsson =  pd.read_csv(f"{REPO_DIR}/data/published/Thorsson_Scores_160_Signatures.tsv",  sep="\t")
    estimate = pd.read_csv(f"{REPO_DIR}/data/published/Yoshihara_ESTIMATE_{cancer_type}_RNAseqV2.txt", sep="\t")
    tcga_absolute = pd.read_csv(f"{REPO_DIR}/data/published/TCGA_ABSOLUTE_tumor_purity.txt", sep="\t")
    gibbons = pd.read_excel(f"{REPO_DIR}/data/published/Gibbons_supp1.xlsx", skiprows=2, sheet_name="DataFileS1 - immune features")

    # Computed Data: Immunedeconv
    mcp_counter = pd.read_csv(f"{output_dir}/immunedeconv/mcp_counter.csv", index_col=0, sep=",")
    quantiseq = pd.read_csv(f"{output_dir}/immunedeconv/quantiseq.csv", index_col=0, sep=",")
    xCell = pd.read_csv(f"{output_dir}/immunedeconv/xcell.csv", index_col=0, sep=",", header=[0])
    EPIC = pd.read_csv(f"{output_dir}/immunedeconv/epic.csv", index_col=0, sep=",")

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
        0.6049872018 + .0001467884 * estimate["ESTIMATE score"])
    estimate = estimate.drop(columns=["ESTIMATE score"])

    # Thorsson data
    Thorsson = Thorsson.drop(columns="Source")
    Thorsson = Thorsson.set_index("SetName").T
    Thorsson = Thorsson.rename_axis(None, axis=1)
    Thorsson.index.name="TCGA_aliquot"
    Thorsson = Thorsson.loc[:, ["LIexpression_score", "CD8_PCA_16704732"]]
    Thorsson.columns = ["TIL score", "CD8 T cells (Thorsson)"]

    # TCGA PanCanAtlas
    tcga_absolute = tcga_absolute.rename(columns = {"purity": "tumor purity (ABSOLUTE)", "sample": "TCGA_aliquot"})
    tcga_absolute = tcga_absolute.set_index("TCGA_aliquot")
    tcga_absolute = pd.DataFrame(tcga_absolute.loc[:, "tumor purity (ABSOLUTE)"])

    gibbons = gibbons.rename(columns={'Unnamed: 1': "id"})
    gibbons["slide_submitter_id"] =  gibbons["id"].str[0:23]
    gibbons["Cytotoxic cells"] = gibbons["Cytotoxic cells"].astype(float)
    gibbons = gibbons.set_index("slide_submitter_id")

    all_slide_features["TCGA_sample"] = clinical_file["slide_submitter_id"].str[0:15]

    # add IDs
    Thorsson["TCGA_sample"] = Thorsson.index.str[0:15]
    tcga_absolute["TCGA_sample"] = tcga_absolute.index.str[0:15]
    gibbons["TCGA_sample"] = gibbons.index.str[0:15]

    tcga_absolute_merged = pd.merge(all_slide_features, tcga_absolute, on=["TCGA_sample", ], how="left")
    Thorsson_merged = pd.merge(all_slide_features, Thorsson, on=["TCGA_sample",], how="left")
    gibbons_merged = pd.merge(all_slide_features, gibbons, on=["TCGA_sample"], how="left")

    cellfrac_merged = pd.merge(all_slide_features, cellfrac, on=["TCGA_sample"], how="left")
    estimate_merged = pd.merge(all_slide_features, estimate, on=["TCGA_sample" ], how="left")
    Fges_computed_merged = pd.merge(all_slide_features, Fges_computed, on=["TCGA_sample"], how="left")

    # Combine in one dataframe
    all_merged = pd.merge(all_slide_features, tcga_absolute_merged, how="left")
    all_merged = pd.merge(all_merged, Thorsson_merged ,how="left")
    all_merged = pd.merge(all_merged, gibbons_merged,how="left")
    all_merged = pd.merge(all_merged, estimate_merged,  how="left")
    all_merged = pd.merge(all_merged, cellfrac_merged,  how="left")
    all_merged = pd.merge(all_merged, Fges_computed_merged,  how="left")

    # ---- Transform features to get a normal distribution (immunedeconv) ---- #
    featuresnames_transform = ["CAFs (MCP counter)",
        'CAFs (EPIC)',]
    feature_data = all_merged.loc[:, CAFS].astype(float)
    data_log2_transformed = feature_data.copy()
    data_log2_transformed[featuresnames_transform] = np.log2(feature_data[featuresnames_transform] * 100 + 0.001)
    CAFs_transformed = data_log2_transformed

    featuresnames_transform = ["Endothelial cells (xCell)",
        "Endothelial cells (EPIC)",]
    feature_data = all_merged.loc[:, ENDOTHELIAL_CELLS].astype(float)
    data_log2_transformed = feature_data.copy()
    data_log2_transformed[featuresnames_transform] = np.log2(feature_data[featuresnames_transform] * 100 + 0.001)
    endothelial_cells_transformed = data_log2_transformed

    feature_data = all_merged.loc[:, T_CELLS].astype(float)
    featuresnames_transform = ['CD8 T cells (quanTIseq)']
    data_log2_transformed =  feature_data.copy()
    data_log2_transformed[featuresnames_transform] = np.log2(feature_data[featuresnames_transform] * 100 + 0.001)
    T_cells_transformed = data_log2_transformed

    feature_data = all_merged.loc[:, TUMOR_PURITY].astype(float)
    featuresnames_transform = ["tumor purity (EPIC)"]
    data_log2_transformed =  feature_data.copy()
    data_log2_transformed[featuresnames_transform] = np.log2(feature_data[featuresnames_transform] * 100 + 0.001)
    tumor_cells_transformed = data_log2_transformed

    # Store processed data
    IDs = ['slide_submitter_id', 'sample_submitter_id', "TCGA_sample"]
    metadata = all_merged[IDs]
    merged = pd.concat([
        metadata,
    CAFs_transformed, endothelial_cells_transformed, T_cells_transformed, tumor_cells_transformed], axis=1)
    merged = merged.fillna(np.nan)

    # Remove slides if there are no values at all
    merged = merged.dropna(axis=0, subset=T_CELLS + CAFS + ENDOTHELIAL_CELLS + TUMOR_PURITY, how="all")
    merged.to_csv(f"{full_output_dir}/ensembled_selected_tasks.csv", sep="\t")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process transcriptomics data for use in TF learning')
    parser.add_argument(
        "--cancer_type",
        help="Abbreviation of cancer type for naming of generated files",
    )
    parser.add_argument(
        "--clinical_file_path",
        help="Full path to clinical file", default=None
    )
    parser.add_argument(
        "--slide_type",
        help="Type of pathology slides either 'FF' (fresh frozen) or 'FFPE' (Formalin-Fixed Paraffin-Embedded) by default 'FF'",
        type=str, required=True
    )
    parser.add_argument(
        "--tpm_path", help="Path to tpm file", type=str, required=True
    )

    parser.add_argument(
            "--output_dir", help="Path to folder for generated file")
    args = parser.parse_args()

    # old_stdout = sys.stdout
    # log_file = open(f"{REPO_DIR}/logs/processing_transcriptomics.log", "w")
    # sys.stdout = log_file

    processing_transcriptomics(
        cancer_type=args.cancer_type,
        slide_type=args.slide_type,
        tpm_path=args.tpm_path,
        path_data=args.path_data,
        clinical_file_path=args.clinical_file_path,
        output_dir=args.output_dir,
    )

    # sys.stdout = old_stdout
    # log_file.close()
