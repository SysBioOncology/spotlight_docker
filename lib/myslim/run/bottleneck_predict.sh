#!/bin/bash
# ACTION Activate virtual environment called tensorflow (WINDOWS)
# virtualenv tensorflow

# cd tensorflow/

# virtualenv tfgpu

# cd ..

# source ./tensorflow/tfgpu/bin/activate

#source /tensorflow/bin/activate
CurDir="$( cd "$(dirname "$0")" ; pwd -P )"

CHECKPOINT_PATH=$1
num_classes=$2		
DATASET_DIR=$3
bot_out=$4
pred_out=$5
model_name=$6

# echo $CHECKPOINT_PATH 
# echo $num_classes 
# echo $DATASET_DIR 
# echo $pred_out $bot_out

python $CurDir/../bottleneck_predict.py \
    --num_classes=$num_classes \
    --bot_out=$bot_out \
    --pred_out=$pred_out \
    --model_name=$model_name \
    --checkpoint_path=$CHECKPOINT_PATH \
    --filedir=$DATASET_DIR \
    --eval_image_size=299 \

echo "Extracted all histopathological features"