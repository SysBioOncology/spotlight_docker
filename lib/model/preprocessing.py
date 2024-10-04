import numpy as np
import pandas as pd
from MFP.portraits.utils import median_scale, read_gene_sets, ssgsea_formula


def process_immunedeconv(data, method, return_all_features=False):
    def format(data, method):
        """
        Formatting data from immunedeconv

        Args:
            data: dataframe with rows = features, columns = samples
            method (str): EPIC, MCP, quanTIseq, xCell
        Returns:
         data: dataframe with rows = samples and columns = features
        """
        if (method == "EPIC") | (method == "quanTIseq"):
            data = data.T
            data = data.set_axis(data.index.str.replace("\\.", "-", regex=True), axis=0)

        elif (method == "xCell") | (method == "MCP"):
            data = data.T
            data = data.set_axis(data.index.str.replace("\\.", "-", regex=True), axis=0)
            data = data.rename_axis(None, axis=1)
        data.index.name = "TCGA_sample"
        return data

    def task_selection(data, method):
        """
        Args:
            data: dataframe with rows = samples and columns = features

        Returns:
            data: dataframe with selected tasks defined by investigation of correlations using initial (large) set of tasks
        """
        if method == "quanTIseq":
            # Renaming
            data = data.rename(columns={"T cell CD8+": "CD8 T cells (quanTIseq)"})
            return pd.DataFrame(data.loc[:, "CD8 T cells (quanTIseq)"])
        elif method == "EPIC":
            data = data.rename(
                columns={
                    "CAFs": "CAFs (EPIC)",
                    "Endothelial": "Endothelial cells (EPIC)",
                    "otherCells": "tumor purity (EPIC)",
                }
            )
            return data.loc[
                :, ["CAFs (EPIC)", "Endothelial cells (EPIC)", "tumor purity (EPIC)"]
            ]
        elif method == "MCP":
            data = data.rename(
                columns={"Cancer associated fibroblast": "CAFs (MCP counter)"}
            )
            return pd.DataFrame(data.loc[:, "CAFs (MCP counter)"])
        elif method == "xCell":
            data = data.rename(
                columns={"Endothelial cell": "Endothelial cells (xCell)"}
            )
            return pd.DataFrame(data.loc[:, "Endothelial cells (xCell)"])

    if return_all_features:
        return add_suffix_to_featurename(format(data, method), method)

    else:
        return task_selection(format(data, method), method)


def add_suffix_to_featurename(data, method):
    """
    Add a suffix to the name

    Args:
        data (DataFrame)
        method (str): suffix to add
    """
    data.columns = ["{} ({})".format(name, method) for name in data.columns]
    return data


def compute_gene_signature_scores(
    tpm_path: str,
    gmt_signatures_path: str,
):
    """
    Re(computing) gene signature scores for Fges Bagaev

    Args:
        tpm_path (str): path pointing to the tpm file .txt
        gmt_signatures_path (str): l path pointing to gmt file containing the gene signatures

    Returns:
        dataframe containing the computed signatures in columns and rows representing the RNAseq samples
    """
    # Read signatures
    gmt = read_gene_sets(gmt_signatures_path)  # GMT format like in MSIGdb

    # Read expressions
    tpm = pd.read_csv(tpm_path, sep="\t", index_col=0, engine="pyarrow")
    tpm_log_transformed = np.log2(tpm + 1)

    # Calc signature scores
    signature_scores = ssgsea_formula(tpm_log_transformed, gmt)

    # Scale signatures
    return median_scale(signature_scores)


def clean_data(bottleneck_features, target_features, slide_type: str):
    """
    Merge histopathological features with the tasks (quantification methods) for TCGA data

    Args:
        bottleneck_features (DataFrame): dataframe containing tiles in rows and the extracted histopathological features (+ ID variables) in the columns
        target_features (DataFrame): dataframe containing the samples (patients) in the rows and the quantification methods of the cell type in the columns.


    """
    if slide_type == "FF":
        all_features = pd.merge(
            bottleneck_features,
            target_features,
            how="inner",
            on=["sample_submitter_id", "slide_submitter_id"],
        )
        all_features = all_features.drop_duplicates()
        all_features = all_features.dropna(axis=0, how="any")

    elif slide_type == "FFPE":
        all_features = bottleneck_features.map_partitions(
            lambda x: pd.merge(
                x,
                target_features,
                how="inner",
                on=["sample_submitter_id", "slide_submitter_id"],
            )
        )
        all_features = all_features.map_partitions(lambda x: x.drop_duplicates())
        all_features = all_features.map_partitions(
            lambda x: x.dropna(axis=0, how="any")
        )

    return all_features
