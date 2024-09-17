#! /bin/bash

export PATH=/home/project/11003054/e1124313/env/GRM/bin/:$PATH
# conda activate GRM 没报错

while IFS= read -r cell_type; do
    mkdir -p /home/users/nus/e1124313/scratch/eqtl/0720_input/${cell_type}/GRM
    if [ -f "/home/users/nus/e1124313/scratch/eqtl/0720_input/${cell_type}/chr/ATGC_${cell_type}_Asian_sle.fam" ]; then
        echo "File exists."
        /home/users/nus/e1124313/project/e1124313/env/GRM/bin/createSparseGRM.R \
        --plinkFile=/home/users/nus/e1124313/scratch/eqtl/0720_input/${cell_type}/chr/ATGC_${cell_type}_Asian_sle \
        --outputPrefix=/home/users/nus/e1124313/scratch/eqtl/0720_input/${cell_type}/GRM/ATGC_${cell_type}_Asian_sle     \
        --numRandomMarkerforSparseKin=2000      \
        --relatednessCutoff=0.125 \
        --nThreads=2
    else
        echo "File does not exist."
    fi

done < "/home/users/nus/e1124313/scratch/eqtl/input/cell_types.csv"