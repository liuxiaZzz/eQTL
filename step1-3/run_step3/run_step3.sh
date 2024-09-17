#! /bin/bash

cell_type=$1
gene_id=$2
chrom=$3
start=$4
end=$5

mkdir -p /home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step3/$cell_type
mkdir -p /home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step3/$cell_type/output
mkdir -p /home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step3/$cell_type/oe

step2prefix=/home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step2/${cell_type}/output/cis_${cell_type}_$gene_id

export PATH=/home/project/11003054/e1124313/env/RSAIGE/bin/:$PATH

Rscript /home/users/nus/e1124313/scratch/qtl/extdata/step3_gene_pvalue_qtl.R       \
--assocFile=$step2prefix \
--geneName=$gene_id \
--genePval_outputFile=${step2prefix}_pval >> /home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step3/${cell_type}/oe/step3_${cell_type}_${gene_id}.log
