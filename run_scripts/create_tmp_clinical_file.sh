#!/bin/bash

awk -v var="$slides" -v a=81 -v b="SKCM_T" -v c=41 -v OFS='\t' 'NR==1{print};NR>1{split(var, tmp, "."); print tmp[1], tmp[1], var, a, b, c}' tmp_clinical_file.txt
