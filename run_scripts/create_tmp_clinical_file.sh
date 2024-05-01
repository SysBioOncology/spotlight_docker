#!/bin/bash
clinical_file="$data_example"
ls $slides_dir | tee list_images.txt
awk -v a=81 -v b="SKCM_T" -v c=41 'FNR==NR{print; next}{split($1, tmp, "."); print tmp[1], tmp[1], $1, a, b, c}' $clinical_file list_images.txt > $clinical_file/final_clinical_file.txt

