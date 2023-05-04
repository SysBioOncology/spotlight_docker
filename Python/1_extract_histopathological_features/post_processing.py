import argparse
import os
from myslim.post_process_features import post_process_features
from myslim.post_process_predictions import post_process_predictions


def execute_postprocessing(output_dir, slide_type):
    """
    1. Format extracted histopathological features
    2. Format predictions of the 42 classes

    Args:
        output_dir (str): path pointing to folder for storing all created files by script

    Returns:
        {output_dir}/features.txt
        {output_dir}/predictions.txt
    """
    post_process_features(output_dir=output_dir, slide_type=slide_type)
    post_process_predictions(output_dir=output_dir, slide_type=slide_type)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output_dir", help="Set output folder", type=str)
    parser.add_argument("--slide_type", help="Type of tissue slide (FF or FFPE)", required=True, type=str)
    args = parser.parse_args()

    full_output_dir = f"{args.output_dir}/1_histopathological_features"
    if not os.path.exists(args.output_dir):
        output_dir = full_output_dir

    execute_postprocessing(
        output_dir=full_output_dir,
	    slide_type=args.slide_type
    )
