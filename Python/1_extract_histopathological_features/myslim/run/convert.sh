#!/bin/bash
# Create a virtual environment for tensorflow 
virtualenv tensorflow

cd tensorflow/

virtualenv tfgpu

cd ..

source ./tensorflow/tfgpu/bin/activate

#TODO Check later if this is really happening, because I think it just installs the latest version.
pip install tensorflow==1.4   
# Error: module 'tensorflow' has no attribute 'placeholder'. Om die op te lossen, heb ik er ==1.4 achter gezet (was 2.4.1)
# omdat 'placeholder' alleen in versie 1 gebruikt kan worden.

CurDir="$( cd "$(dirname "$0")" ; pwd -P )"

outputDir=$1
num_shards=$2

python $CurDir/../datasets/convert.py $outputDir $num_shards

echo "Done with conversion"


