//
// Subworkflow with functionality specific to the nf-core/spotlight pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GENERATE_GRAPHS } from '../../modules/local/generategraphs.nf'
include {COMPUTE_CONNECTEDNESS} from '../../modules/local/computeconnectedness.nf'
include { COMPUTE_COLOCALIZATION} from '../../modules/local/computecolocalization.nf'
include { COMPUTE_NODE_DEGREE_WITH_ES} from '../../modules/local/computenodedegreewithes.nf'
include {COMPUTE_N_SHORTEST_PATHS_WITH_MAX_LENGTH} from '../../modules/local/computenshortestwithmaxlength.nf'
include {CLUSTERING_SCHC_SIMULTANEOUS} from '../../modules/local/clusteringschcsimultaneous.nf'
include { CLUSTERING_SCHC_INDIVIDUAL} from '../../modules/local/clusteringschcindividual.nf'
include { COMPUTE_NCLUSTERS } from '../../modules/local/computenclusters.nf'
include {COMPUTE_FRAC_HIGH } from '../../modules/local/computefrachigh.nf'
include {COMPUTE_PROXIMITY_FROM_SIMULTANEOUS_SCHC} from '../../modules/local/computeproximityfromsimultaneousschc.nf'
include {COMPUTE_PROXIMITY_FROM_INDIV_SCHC_BETWEEN} from '../../modules/local/computeproximityfromindivschcbetween.nf'
include {COMPUTE_PROXIMITY_FROM_INDIV_SCHC_WITHIN} from '../../modules/local/computeproximityfromindivschcwithin.nf'
include { COMPUTE_PROXIMITY_FROM_INDIV_SCHC_COMBINE} from '../../modules/local/computeproximityfromindivschccombine.nf'
include {COMBINE_CLUSTERING_FEATURES} from '../../modules/local/combineclusteringfeatures.nf'
include {COMBINE_NETWORK_FEATURES} from '../../modules/local/combinenetworkfeatures.nf'
include { COMBINE_ALL_SPATIAL_FEATURES} from '../../modules/local/combineallspatialfeatures.nf'


workflow DERIVE_SPATIAL_FEATURES {
    take:
        tile_level_cell_type_quantification
        out_prefix
        slide_type
        abundance_threshold
        cell_types
        shapiro_alpha
        cutoff_path_length
        n_clusters
        max_dist
        max_n_tiles_threshold
        tile_size
        overlap
        metadata_path
        is_tcga
        merge_var
        sheet_name

    main:

    GENERATE_GRAPHS(
        tile_quantification_path = tile_level_cell_type_quantification,
        out_prefix = out_prefix,
        slide_type = slide_type
    )

    // COMPUTE_SPATIAL_FEATURES(
    //     tile_quantification_path = tile_level_cell_type_quantification
    //     )

    COMPUTE_CONNECTEDNESS(
        tile_quantification_path = tile_level_cell_type_quantification,
        cell_types = cell_types,
        graphs_path = GENERATE_GRAPHS.out.pkl,
        abundance_threshold = abundance_threshold,
        slide_type = slide_type,
        out_prefix = out_prefix
    )

    COMPUTE_COLOCALIZATION(
        tile_quantification_path = tile_level_cell_type_quantification,
        cell_types = cell_types,
        graphs_path = GENERATE_GRAPHS.out.pkl,
        abundance_threshold = abundance_threshold,
        slide_type = slide_type,
        out_prefix = out_prefix
    )

    COMPUTE_NODE_DEGREE_WITH_ES(
        tile_quantification_path = tile_level_cell_type_quantification,
        cell_types = cell_types,
        graphs_path = GENERATE_GRAPHS.out.pkl,
        shapiro_alpha = shapiro_alpha,
        slide_type = slide_type,
        out_prefix = out_prefix
    )

    COMPUTE_N_SHORTEST_PATHS_WITH_MAX_LENGTH(
        tile_quantification_path = tile_level_cell_type_quantification,
        cell_types = cell_types,
        graphs_path = GENERATE_GRAPHS.out.pkl,
        cutoff_path_length = cutoff_path_length,
        slide_type = slide_type,
        out_prefix = out_prefix

    )

    CLUSTERING_SCHC_SIMULTANEOUS(
        tile_quantification_path = tile_level_cell_type_quantification,
        cell_types = cell_types,
        graphs_path = GENERATE_GRAPHS.out.pkl,
        slide_type = slide_type,
        out_prefix = out_prefix

    )

    CLUSTERING_SCHC_INDIVIDUAL(
        tile_quantification_path = tile_level_cell_type_quantification,
        cell_types = cell_types,
        graphs_path = GENERATE_GRAPHS.out.pkl,
        slide_type = slide_type,
        out_prefix = out_prefix

    )

    COMPUTE_NCLUSTERS(
        CLUSTERING_SCHC_SIMULTANEOUS.out.csv,
        cell_types = cell_types,
        slide_type = slide_type,
        out_prefix = out_prefix
    )

    COMPUTE_FRAC_HIGH(
        CLUSTERING_SCHC_INDIVIDUAL.out.csv,
        slide_type = slide_type,
        out_prefix = out_prefix
    )

    COMPUTE_PROXIMITY_FROM_SIMULTANEOUS_SCHC(
        CLUSTERING_SCHC_SIMULTANEOUS.out.csv,
        cell_types = cell_types,
        n_clusters,
        max_dist,
        max_n_tiles_threshold,
        tile_size,
        overlap,
        slide_type = slide_type,
        out_prefix = out_prefix
    )


    COMPUTE_PROXIMITY_FROM_INDIV_SCHC_WITHIN(
        CLUSTERING_SCHC_INDIVIDUAL.out.csv,
        cell_types = cell_types,
        n_clusters,
        max_dist,
        max_n_tiles_threshold,
        tile_size,
        overlap,
        slide_type = slide_type,
        out_prefix = out_prefix
    )


    COMPUTE_PROXIMITY_FROM_INDIV_SCHC_BETWEEN(
        CLUSTERING_SCHC_INDIVIDUAL.out.csv,
        cell_types = cell_types,
        n_clusters,
        max_dist,
        max_n_tiles_threshold,
        tile_size,
        overlap,
        slide_type = slide_type,
        out_prefix = out_prefix
    )

    COMPUTE_PROXIMITY_FROM_INDIV_SCHC_COMBINE(
        prox_between = COMPUTE_PROXIMITY_FROM_INDIV_SCHC_BETWEEN.out.csv,
        prox_within = COMPUTE_PROXIMITY_FROM_INDIV_SCHC_WITHIN.out.csv,
        slide_type = slide_type,
        out_prefix = out_prefix

    )

    COMBINE_CLUSTERING_FEATURES(
        frac_high_wide = COMPUTE_FRAC_HIGH.out.csv,
        num_clust_slide_wide = COMPUTE_NCLUSTERS.out.csv,
        all_prox_df_wide  = COMPUTE_PROXIMITY_FROM_SIMULTANEOUS_SCHC.out.csv,
        prox_indiv_schc_combined_wide = COMPUTE_PROXIMITY_FROM_INDIV_SCHC_COMBINE.out.csv,
        slide_type = slide_type,
        out_prefix = out_prefix
    )

    COMBINE_NETWORK_FEATURES(
    all_largest_cc_sizes_wide = COMPUTE_CONNECTEDNESS.out.csv,
    shortest_paths_wide = COMPUTE_N_SHORTEST_PATHS_WITH_MAX_LENGTH.out.csv,
    colocalization_wide = COMPUTE_COLOCALIZATION.out.csv,
    slide_type = slide_type,
    out_prefix = out_prefix
    )

    COMBINE_ALL_SPATIAL_FEATURES(
        graph_features = COMBINE_NETWORK_FEATURES.out.csv,
        clustering_features = COMBINE_CLUSTERING_FEATURES.out.csv,
        metadata_path = metadata_path,
        is_tcga = is_tcga,
        merge_var = merge_var,
        sheet_name = sheet_name,
        slide_type = slide_type,
        out_prefix = out_prefix
    )

    // COMPUTE_SHAPE_FEATURES(
    //     CLUSTERING_SCHC_SIMULTANEOUS.out.csv,
    //     cell_types = cell_types,
    //     tile_size = tile_size,
    //     overlap = overlap,
    //     slide_type = slide_type,
    //     out_prefix = out_prefix
    // )
    // COMPUTE_NETWORK_FEATURES() {
    //     tile_level_cell_type_quantification = tile_level_cell_type_quantification,
    //     slide_type = slide_type,
    //     graphs_path = GENERATE_GRAPHS.pkl,
    //     abundance_threshold = abundance_threshold,
    //     shapiro_alpha = shapiro_alpha,
    //     cutoff_path_length = cutoff_path_length,
    // }



}
