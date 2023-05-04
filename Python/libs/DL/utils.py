def get_tile_name(tile_path):
    """
    Return tile names including x and y coordinates

    Args:
        tile_path (str): string to be formatted, which is a path pointing to the tile image e.g. b'C:/Users/20171880/PC_disk/LUNG_TEST/tiles/TCGA-43-A475-01A-01-TSA_10165_15709.jpg'

    Returns:
        tile_name (str): tile name including position of tile, e.g. TCGA-43-A475-01A-01-TSA_10165_15709

    """
    filename = tile_path.split("/")[-1]
    tile_name = filename.split(".")[0]
    return tile_name


def get_slide_submitter_id(tile_filename):
    """ "
    Returns slide_submitter_id from tile filename.

    Args:
        tile_filename (str): e.g. TCGA-05-4250-01A-01-BS1_10165_10165.jpg

    Returns:
        slide_name (str): e.g. TCGA-05-4250-01A-01-BS1
    """
    tilename = tile_filename.strip(".jpg")
    slidename = tilename.split("_")[0]
    return slidename
