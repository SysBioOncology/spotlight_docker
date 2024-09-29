//
// Subworkflow with functionality specific to the SysBioOncology/spotlight_docker pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include {   CREATE_CLINICAL_FILE  } from '../../modules/local/createclinicalfile.nf'
include {   CREATE_LIST_AVAIL_SLIDES } from '../../modules/local/createlistavailableslides.nf'
include {   TILING_SINGLE_SLIDE } from '../../modules/local/tilingsingleslide.nf'
include {   FORMAT_TILE_DATA_STRUCTURE  } from '../../modules/local/formattiledatastructure.nf'
include {   PREPROCESSING_SLIDES } from '../../modules/local/preprocessingslides.nf'
include {   BOTTLENECK_PREDICT  } from '../../modules/local/bottleneckpredict.nf'
include {   POST_PROCESSING_FEATURES    } from '../../modules/local/postprocessingfeatures.nf'
include {   POST_PROCESSING_PREDICTIONS } from '../../modules/local/postprocessingpredictions.nf'



workflow EXTRACT_HISTOPATHO_FEATURES {
    take:
        clinical_files_input
        path_codebook
        class_name
        out_file
        tumor_purity_threshold
        is_tcga
        image_dir
        gradient_mag_filter
        n_shards
        bot_out
        pred_out
        model_name
        checkpoint_path
        slide_type
        path_tissue_classes

    main:

    CREATE_CLINICAL_FILE(
        clinical_files_input = clinical_files_input,
        class_name = class_name,
        out_file = out_file,
        path_codebook = path_codebook,
        tumor_purity_threshold = tumor_purity_threshold,
        is_tcga = is_tcga,
        image_dir = image_dir
    )
    CREATE_LIST_AVAIL_SLIDES(
        clinical_file_path = CREATE_CLINICAL_FILE.out.txt,
        image_dir = image_dir
    )

    avail_img_to_process = CREATE_LIST_AVAIL_SLIDES.out.csv \
                                    | splitCsv(header:true) \
                                    | map { row -> tuple(row.slide_id, row.slide_filename, file("${image_dir}/${row.slide_filename}")) }
    TILING_SINGLE_SLIDE (
        avail_img_to_process,
        gradient_mag_filter = gradient_mag_filter,
    )

    FORMAT_TILE_DATA_STRUCTURE(
        all_tiles               = TILING_SINGLE_SLIDE.out.jpg.collect(),
        clinical_file_path      = CREATE_CLINICAL_FILE.out.txt,
        image_dir               = image_dir,
        is_tcga                 = is_tcga
    )

    PREPROCESSING_SLIDES(
        file_info_train         = FORMAT_TILE_DATA_STRUCTURE.out.txt,
        n_shards                = n_shards
    )

    BOTTLENECK_PREDICT(
        bot_out                 = bot_out,
        pred_out                = pred_out,
        tf_records              = PREPROCESSING_SLIDES.out.tfrecords.collect(),
        model_name              = model_name,
        checkpoint_path         = checkpoint_path
    )

    POST_PROCESSING_FEATURES(
        bot_train_file = BOTTLENECK_PREDICT.out.bot_txt,
        slide_type = slide_type,
        is_tcga = is_tcga
    )

    POST_PROCESSING_PREDICTIONS(
        path_codebook = path_codebook,
        path_tissue_classes = path_tissue_classes,
        pred_train_file = BOTTLENECK_PREDICT.out.pred_txt,
        cancer_type = class_name,
        slide_type = slide_type
    )

    emit:
    features = POST_PROCESSING_FEATURES.out.txt_parquet
    predictions = POST_PROCESSING_PREDICTIONS.out.txt_parquet
}
