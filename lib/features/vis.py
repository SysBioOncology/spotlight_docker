import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import sys
import git
from statsmodels.stats.multitest import multipletests
import itertools
import sys
import pandas as pd
from scipy.stats import mannwhitneyu
from statannotations.Annotator import Annotator


REPO_DIR= git.Repo('.', search_parent_directories=True).working_tree_dir
sys.path.append(f"{REPO_DIR}/Python/libs")

# Own modules
from model.constants import *

# TODO: test this function
def proba_maps(predictions, output_dir,  slide_submitter_id=None, vmin=0.20, vmax=0.80, cell_types=None, cmap_cell_types=None, fontsize=12):
    """ Create cell type probability maps """
    if cell_types is None:
        cell_types = DEFAULT_CELL_TYPES
    if cmap_cell_types is None:
        cmap_cell_types = {"T_cells": "Greens", "CAFs": "Blues", "tumor_purity":"Reds", "endothelial_cells":"Purples"}


    if len(np.unique(predictions.slide_submitter_id)) > 1:
        if slide_submitter_id is None:
            raise Exception("If slide_submitter_id is not specified then predictions should contain only data for one slide")
        else:
            predictions = predictions[predictions.slide_submitter_id == slide_submitter_id]

    fig, axes = plt.subplots(1, len(cell_types), sharex=False, sharey=False, figsize=(20, 5))

    # Plot cell type maps
    for j, cell_type in enumerate(cell_types):
        temp = predictions.pivot("Coord_Y", "Coord_X" ,cell_type)
        temp = temp.sort_index(ascending=False)
        sns.heatmap(temp,cmap=cmap_cell_types[cell_type], vmin=vmin, vmax=vmax,ax=axes[j], cbar=False)
        axes[j].set_xticklabels("")
        axes[j].set_yticklabels("")
        axes[j].set_xlabel("")
        axes[j].set_ylabel("")
        axes[j].get_xaxis().set_ticks([])
        axes[j].get_yaxis().set_ticks([])

        for spine in ['top', 'right', "left", "bottom"]:
            axes[j].spines[spine].set_visible(False)

    # Add custom colorbar
    norm = plt.Normalize(vmin, vmax)
    for j, cell_type in enumerate(cell_types):
        axes[j].set_ylabel(cell_type.replace("_", " ").replace("purity", "cells"), size=fontsize)
        sm = plt.cm.ScalarMappable(cmap=cmap_cell_types[cell_type], norm=norm, )
        sm.set_array([])
        axins = inset_axes(axes[j],
                        width="100%",
                        height="3%",
                        loc='lower center',
                        borderpad=-2
                    )
        fig.colorbar(sm, cax=axins, orientation="horizontal",label="probability",
        ticks=np.linspace(vmin, vmax, 4),
        spacing='proportional',)

    if slide_submitter_id is None:
        slide_submitter_id = np.unique(predictions["slide_submitter_id"])[0]

    plt.savefig(f"{output_dir}/{slide_submitter_id}_cell_type_maps.pdf", bbox_inches="tight", dpi=600)


def create_boxplots(feature_data, feature="value", level="variable", orient="h", xlabel="", ylabel="", alpha=0.05, figsize=(6, 4), labelpad=10, fontsize=12, output_dir=None, plot_name=None, compare_subtypes=True):
    """ Create boxplots with multiple testing correction and statistical annotations

    """

    MFP_order = ["IE", "IE/F", "F", "D"]

    # Setup for visualization
    labelpad = 10
    fontsize=12
    subtype_palette =dict(zip(["F", "IE/F", "IE", "D"], sns.color_palette("pastel")[:4]))

    # TODO adapt later for the use given cell types (might even not be necessary)
    feature_data = feature_data.replace({
    "_":" ",
    "purity": "cells",
    "endothelial": "endo",
    }, regex=True)


    if orient == "v":
        x = level
        y = feature
    elif orient== "h":
        x = feature
        y = level

    common_params = {
        'data': feature_data,
        "x": (orient == "v") * level or (orient == "h") * feature,
        "y": (orient == "v") * feature or (orient == "h") * level,
        "order": sorted(feature_data[level].unique()),
        "hue": (compare_subtypes * "MFP") or None,
        "hue_order": (compare_subtypes * MFP_order) or None,
        }

    box_plot_params = common_params.copy()
    box_plot_params.update({
        "flierprops":{"marker":"o"},
        "fliersize": 1,
        "showfliers": False,
        "orient": orient,
        "palette": subtype_palette if compare_subtypes else {"T cells":"forestgreen", "CAFs": "royalblue", "endo cells": "lightskyblue", "tumor cells": "red"} ,
    })

    strip_params = common_params.copy()
    strip_params.update({
        "palette": "dark:black",
        "size": 1,
        "dodge": True
    })


    # Create new plot
    fig, ax = plt.subplots(figsize=figsize)

    # Add distribution
    g = sns.boxplot(ax=ax, **box_plot_params)
    handles, labels = ax.get_legend_handles_labels()

    # Add individual datapoints
    sns.stripplot(ax=ax, **strip_params)

    if compare_subtypes:
        all_pairs = []
        H_values = []
        p_values = []

        for x in sorted(feature_data[level].unique()):
            temp = list(itertools.product([x], MFP_order))
            subtype_combinations = [list(comparison) for comparison in list(itertools.combinations(temp, 2))]
            # Investigate difference between subtypes per cell type or pair (Adjust per cell type or pair)
            for (pair1, subtype1), (pair2, subtype2) in subtype_combinations:
                H, p = mannwhitneyu(feature_data.loc[(feature_data[level] == pair1) & (feature_data.MFP == subtype1), feature], feature_data.loc[(feature_data[level] == pair2) & (feature_data.MFP == subtype2), feature])
                H_values.append(H)
                p_values.append(p)
            all_pairs += [list(comparison) for comparison in list(itertools.combinations(temp, 2))]

        # Statistical testing
        p_values = multipletests(p_values, alpha=alpha, method="fdr_bh")[1]
        sign_pvals = pd.Series(p_values)[list(np.array(p_values) < alpha)].to_list()
        sign_pairs = pd.Series(all_pairs)[list(np.array(p_values) < alpha)].to_list()

        if len(sign_pairs) > 0:
            # Add annotations: pvalues
            annot = Annotator(g, sign_pairs, verbose=1, **box_plot_params)
            annot.new_plot(ax=g, pairs=sign_pairs,
                        **box_plot_params)
            (annot
            .configure(test=None, test_short_name="mannwhitneyu", verbose=0)
            .set_pvalues(pvalues=sign_pvals)
            .annotate())

        # Configure axis/legends
        plt.legend(handles, labels, loc="upper left", bbox_to_anchor=(1,1), frameon=False, title="Subtype")
    else:
        pairs = list(itertools.combinations(feature_data.cell_type.unique(),2))
        annot = Annotator(ax, pairs, data=feature_data, x=level, y=feature,  order=sorted(feature_data[level].unique()))
        annot.configure(test='Mann-Whitney', text_format='star', loc='inside', verbose=2, comparisons_correction="BH")
        annot.apply_test()
        ax, _ = annot.annotate()


    ax.set_xlabel(xlabel, size=fontsize, labelpad=labelpad)
    ax.set_ylabel(ylabel, size=fontsize)

    if orient == "v":
        ax.set_xticklabels(sorted(feature_data[level].unique()), size=fontsize)
    elif orient == "h":
        ax.set_yticklabels(sorted(feature_data[level].unique()), size=fontsize)

    if output_dir is not None:
        if plot_name is None:
            plt.savefig(f"{output_dir}/boxplots.pdf", bbox_inches="tight", dpi=600)
        else:
            plt.savefig(f"{output_dir}/boxplots_{plot_name}.pdf", bbox_inches="tight", dpi=600)

    plt.show()
