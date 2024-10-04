import numpy as np
import pandas as pd


def selecting_tiles(all_features, n_tiles, slide_type):
    """
    Selecting a random number of tiles for each slide in the dataset used for training of the multi-task models

    Args:
        all_features (DataFrame): dataframe containing the 1,536 features extracted per tile and at least the slide_submitter_id.
        n_tiles (int): number of tiles to select per slide

    Returns:
        total_tile_selection (DataFrame): dataframe containing the randomly selected n_tiles tiles per slide.
    """
    column_names = list(all_features.columns)
    all_features["dummy"] = 1  # for counting number of tiles
    total_tile_selection = pd.DataFrame(columns=column_names)  # Empty dataframe
    selected_tiles = pd.Series(dtype="str")  # Empty series

    if slide_type == "FF":
        tiles_per_slide = all_features.groupby("slide_submitter_id")["dummy"].count()
        slide_IDs = tiles_per_slide[tiles_per_slide >= n_tiles].index.values
        all_features = all_features.drop(columns="dummy")
        # Retrieve all tiles
        all_tiles = all_features.tile_ID
        all_tiles_list = all_tiles.to_list()
        tmp = [i[0:23] for i in all_tiles_list]
        all_tiles_list = np.array(all_tiles_list)

        for slide_ID in slide_IDs:
            print("slide_ID", slide_ID, "\n")
            # 1. Select tiles for corresponding slide
            sel = np.array([i == slide_ID for i in tmp])
            tiles_subset = all_tiles_list[sel]
            # 2. Get 50 random selected tiels
            tiles_subset = pd.Series(tiles_subset).sample(n=50)
            # 3. Add to series
            selected_tiles = pd.concat([selected_tiles, tiles_subset])

        # Get features for all slides corresponding tiles
        total_tile_selection = all_features[all_features.tile_ID.isin(selected_tiles)]

    elif slide_type == "FFPE":
        tiles_per_slide = all_features.map_partitions(
            lambda x: x.groupby("slide_submitter_id").count()
        )
        tiles_per_slide = tiles_per_slide.compute().loc[:, "dummy"]
        tiles_per_slide = (
            tiles_per_slide.to_frame().groupby("slide_submitter_id")["dummy"].sum()
        )
        slide_IDs = tiles_per_slide[tiles_per_slide >= n_tiles].index.values
        all_features = all_features.drop(columns="dummy")
        # Retrieve all tiles
        all_tiles = all_features.tile_ID.compute()
        all_tiles_list = all_tiles.to_list()
        tmp = [i[0:23] for i in all_tiles_list]
        all_tiles_list = np.array(all_tiles_list)

        for slide_ID in slide_IDs:
            print("slide_ID", slide_ID, "\n")
            # 1. Select tiles for corresponding slide
            sel = np.array([i == slide_ID for i in tmp])
            tiles_subset = all_tiles_list[sel]
            # 2. Get 50 random selected tiels
            tiles_subset = pd.Series(tiles_subset).sample(n=50)
            # 3. Add to series
            selected_tiles = pd.concat([selected_tiles, tiles_subset])

        # Get features for all slides corresponding tiles
        features_slide = all_features.map_partitions(
            lambda x: x.loc[x["tile_ID"].isin(selected_tiles)]
        )
        total_tile_selection = features_slide.compute()

    total_tile_selection = total_tile_selection.reset_index(drop=True)
    return total_tile_selection


def split_in_XY(tile_selection, metadata_colnames, category):
    """
    Separate the histopathological features from the metadata (IDs, coordinates etc.) used for training of the multi-task models

    Args:
        tile_selection (DataFrame): dataframe containing the randomly selected n_tiles per slide
        metadata_colnames (list): list of variables from task_seleciton_names.pkl file
        category (str): cell type, one of T_cells, tumor_purity, CAFs or endothelial_cells

    Returns:
        X (DataFrame): dataframe containing the histopathological features for the selected slides
        Y (DataFrame): dataframe containing the cell type quantification for the specified cell type
    """
    X = tile_selection.loc[:, ~tile_selection.columns.isin(metadata_colnames)]
    Y = tile_selection.loc[:, category]
    return [X, Y]
