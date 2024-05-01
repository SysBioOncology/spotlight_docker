import argparse
import os
from myslim.post_process_features import post_process_features
from myslim.post_process_predictions import post_process_predictions


def execute_postprocessing(output_dir, slide_type, path_codebook, path_tissue_classes, data_source):
    """
    1. Format extracted histopathological features
    2. Format predictions of the 42 classes

    Args:
        output_dir (str): path pointing to folder for storing all created files by script

    Returns:
        {output_dir}/features.txt
        {output_dir}/predictions.txt
    """
    post_process_features(output_dir=output_dir,
                          slide_type=slide_type, data_source=data_source)
    post_process_predictions(output_dir=output_dir, slide_type=slide_type,
                             path_codebook=path_codebook, path_tissue_classes=path_tissue_classes)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output_dir", help="Set output folder", type=str)
    parser.add_argument(
        "--slide_type", help="Type of tissue slide (FF or FFPE)", required=True, type=str)
    parser.add_argument(
        "--path_codebook", help="codebook.txt file", required=True, type=str)
    parser.add_argument(
        "--path_tissue_classes", help="Tissue_classes.csv file", required=True, type=str)
    parser.add_argument(
        "--data_source", help="Data source, default='TCGA'")
    args = parser.parse_args()

    if os.path.exists(args.output_dir):
        print("Output folder exists")
    else:
        os.makedirs(args.output_dir)

    execute_postprocessing(
        output_dir=args.output_dir,
        slide_type=args.slide_type,
        path_codebook=args.path_codebook,
        path_tissue_classes=args.path_tissue_classes,
        data_source=args.data_source
    )
