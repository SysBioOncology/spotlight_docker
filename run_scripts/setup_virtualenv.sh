#!/bin/bash

# ----------------------------------- #
# ---- Setup virtual environment ---- #
# ----------------------------------- #

# Define repo_path
#repo_path="/Users/joankant/Library/CloudStorage/OneDrive-TUEindhoven/spotlight"
repo_path="/home/olapuent/Desktop/spatial_localization/repo_manuscript/spotlight/"
cd $repo_path/run_scripts

# python3 -m pip install virtualenv (virtualenv needs to be installed)
if [ ! -d "spotlight_env" ]
then
	python3 -m virtualenv spotlight_env
fi
# for bash
. spotlight_env/bin/activate
#source $repo_dir/run_scripts/spotlight_env/bin/activate
spotlight_env/bin/pip install -r $repo_path/requirements.txt --default-timeout=500
