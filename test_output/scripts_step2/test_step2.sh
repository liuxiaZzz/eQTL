#! /bin/bash

# cell_type="B_mem"
# gene_id="ENSG00000150093"
# chrom="10"
# start="32887273"
# end="33005792"

cell_type=$1
gene_id=$2
chrom=$3
start=$4
end=$5

# 强制转换为整数
start=${start%.*}
end=${end%.*}
start=$((start - 1000000))
end=$((end + 1000000))

mkdir -p /home/users/nus/e1124313/scratch/eqtl/test_output/step2/
mkdir -p /home/users/nus/e1124313/scratch/eqtl/test_output/step2/${cell_type}
mkdir -p /home/users/nus/e1124313/scratch/eqtl/test_output/step2/${cell_type}/oe
mkdir -p /home/users/nus/e1124313/scratch/eqtl/test_output/step2/${cell_type}/output
mkdir -p /home/users/nus/e1124313/scratch/eqtl/test_output/step2/${cell_type}/region

step2prefix=/home/users/nus/e1124313/scratch/eqtl/test_output/step2/${cell_type}/output/cis_${cell_type}_$gene_id
step1prefix=/home/users/nus/e1124313/scratch/eqtl/test_output/step1/${cell_type}/output/${cell_type}_$gene_id
region="${chrom}\t${start}\t${end}"
regionFile=/home/users/nus/e1124313/scratch/eqtl/test_output/step2/${cell_type}/region/${cell_type}_${gene_id}_range.txt
echo -e $region > $regionFile

export PATH=/home/project/11003054/e1124313/env/RSAIGE/bin/:$PATH

Rscript /home/users/nus/e1124313/scratch/qtl/extdata/step2_tests_qtl.R      \
--bedFile=/home/users/nus/e1124313/scratch/eqtl/0920_input/${cell_type}/chr/ATGC_${cell_type}_Asian_sle_chr${chrom}.bed      \
--bimFile=/home/users/nus/e1124313/scratch/eqtl/0920_input/${cell_type}/chr/ATGC_${cell_type}_Asian_sle_chr${chrom}.bim      \
--famFile=/home/users/nus/e1124313/scratch/eqtl/0920_input/${cell_type}/chr/ATGC_${cell_type}_Asian_sle_chr${chrom}.fam      \
--SAIGEOutputFile=${step2prefix}     \
--chrom=$chrom       \
--minMAF=0.05 \
--minMAC=20 \
--LOCO=FALSE    \
--GMMATmodelFile=${step1prefix}.rda     \
--SPAcutoff=2 \
--varianceRatioFile=${step1prefix}.varianceRatio.txt    \
--rangestoIncludeFile=${regionFile}     \
--markers_per_chunk=10000 > /home/users/nus/e1124313/scratch/eqtl/test_output/step2/${cell_type}/oe/step2_${cell_type}_${gene_id}.log
