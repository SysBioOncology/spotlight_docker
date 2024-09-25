import multiprocessing
import sys

NUM_CORES = multiprocessing.cpu_count()
METADATA_COLS = ['tile_ID',  'slide_submitter_id', 'Section',
                    'Coord_X', 'Coord_Y', 'TCGA_patient_ID', ]
DEFAULT_CELL_TYPES = ["CAFs", "T_cells", "endothelial_cells", "tumor_purity"]
DEFAULT_SLIDE_TYPE = "FF"
ISMACOS = sys.platform == "darwin"


# Final feature selection
TUMOR_PURITY = [
    'tumor purity (ABSOLUTE)',
    'tumor purity (ESTIMATE)',
    'tumor purity (EPIC)'
]

T_CELLS = [
    # 'CD8 T cells (Thorsson)',
    'Cytotoxic cells',
    'Effector cells',
    'CD8 T cells (quanTIseq)',
    # 'TIL score',
    'Immune score',
]

ENDOTHELIAL_CELLS = [
    'Endothelial cells (xCell)',
    'Endothelial cells (EPIC)',
    'Endothelium', ]

CAFS = [
    'Stromal score',
    'CAFs (MCP counter)',
    'CAFs (EPIC)',
    'CAFs (Bagaev)',
]

IDS = ['slide_submitter_id', 'sample_submitter_id']
TILE_VARS = ['Section', 'Coord_X', 'Coord_Y', "tile_ID"]
