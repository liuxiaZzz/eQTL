#! /bin/bash

# this script is to run manually

# cell_type="B_mem"
# gene_id="ENSG00000150093"

cell_type=$1
gene_id=$2

mkdir -p /home/users/nus/e1124313/scratch/eqtl/test_output/step1
mkdir -p /home/users/nus/e1124313/scratch/eqtl/test_output/step1/${cell_type}
mkdir -p /home/users/nus/e1124313/scratch/eqtl/test_output/step1/${cell_type}/oe
mkdir -p /home/users/nus/e1124313/scratch/eqtl/test_output/step1/${cell_type}/output

export PATH=/home/project/11003054/e1124313/env/RSAIGE/bin/:$PATH

### use GRM
Rscript /home/users/nus/e1124313/scratch/qtl/extdata/step1_fitNULLGLMM_qtl.R \
--useSparseGRMtoFitNULL=TRUE \
--sparseGRMSampleIDFile=/home/users/nus/e1124313/scratch/eqtl/0720_input/${cell_type}/GRM/ATGC_${cell_type}_Asian_sle_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
--sparseGRMFile=/home/users/nus/e1124313/scratch/eqtl/0720_input/${cell_type}/GRM/ATGC_${cell_type}_Asian_sle_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx \
--phenoFile=/home/users/nus/e1124313/scratch/eqtl/0920_input/${cell_type}/${cell_type}.txt \
--phenoCol=$gene_id \
--covarColList=Age,Sex,geno_PC1,geno_PC2,geno_PC3,geno_PC4,geno_PC5,geno_PC6 \
--sampleCovarColList=Age,Sex,geno_PC1,geno_PC2,geno_PC3,geno_PC4,geno_PC5,geno_PC6 \
--sampleIDColinphenoFile=Sample_ID \
--traitType=count \
--outputPrefix=/home/users/nus/e1124313/scratch/eqtl/test_output/step1/${cell_type}/output/${cell_type}_$gene_id \
--skipVarianceRatioEstimation=FALSE \
--isRemoveZerosinPheno=FALSE \
--isCovariateOffset=FALSE \
--isCovariateTransform=TRUE \
--skipModelFitting=FALSE \
--tol=0.00001 \
--plinkFile=/home/users/nus/e1124313/scratch/eqtl/0920_input/${cell_type}/chr/ATGC_${cell_type}_Asian_sle \
--nThreads=1 \
--IsOverwriteVarianceRatioFile=TRUE > /home/users/nus/e1124313/scratch/eqtl/test_output/step1/${cell_type}/oe/step1_${cell_type}_${gene_id}.log

# ### don't use GRM
# Rscript /home/users/nus/e1124313/scratch/qtl/extdata/step1_fitNULLGLMM_qtl.R \
# --useSparseGRMtoFitNULL=TRUE \
# --sparseGRMSampleIDFile=/home/users/nus/e1124313/scratch/eqtl/0720_input/${cell_type}/GRM/ATGC_${cell_type}_Asian_sle_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
# --sparseGRMFile=/home/users/nus/e1124313/scratch/eqtl/0720_input/${cell_type}/GRM/ATGC_${cell_type}_Asian_sle_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx \
# --phenoFile=/home/users/nus/e1124313/scratch/eqtl/0720_input/${cell_type}/${cell_type}.txt \
# --phenoCol=$gene_id \
# --covarColList=Age,Sex,geno_PC1,geno_PC2,geno_PC3,geno_PC4,geno_PC5,geno_PC6 \
# --sampleCovarColList=Age,Sex,geno_PC1,geno_PC2,geno_PC3,geno_PC4,geno_PC5,geno_PC6 \
# --sampleIDColinphenoFile=Sample_ID \
# --traitType=count \
# --outputPrefix=/home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step1/${cell_type}/output/${cell_type}_$gene_id \
# --skipVarianceRatioEstimation=FALSE \
# --isRemoveZerosinPheno=FALSE \
# --isCovariateOffset=FALSE \
# --isCovariateTransform=TRUE \
# --skipModelFitting=FALSE \
# --tol=0.00001 \
# --plinkFile=/home/users/nus/e1124313/scratch/eqtl/0720_input/${cell_type}/chr/ATGC_${cell_type}_Asian_sle \
# --nThreads=1 \
# --IsOverwriteVarianceRatioFile=TRUE > /home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step1/${cell_type}/oe/step1_${cell_type}_${gene_id}.log