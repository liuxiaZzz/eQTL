#! /bin/bash

cell_type=$1
gene_id=$2

mkdir -p /home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step1/${cell_type}
mkdir -p /home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step1/${cell_type}/oe
mkdir -p /home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step1/${cell_type}/output

export PATH=/home/project/11003054/e1124313/env/RSAIGE/bin/:$PATH

Rscript /home/users/nus/e1124313/scratch/qtl/extdata/step1_fitNULLGLMM_qtl.R \
--useSparseGRMtoFitNULL=FALSE \
--useGRMtoFitNULL=FALSE \
--phenoFile=/home/users/nus/e1124313/scratch/eqtl/input/${cell_type}/${cell_type}.txt \
--phenoCol=$gene_id \
--covarColList=Age,Sex \
--sampleCovarColList=Age,Sex \
--sampleIDColinphenoFile=Sample_ID \
--traitType=count \
--outputPrefix=/home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step1/${cell_type}/output/${cell_type}_$gene_id \
--skipVarianceRatioEstimation=FALSE \
--isRemoveZerosinPheno=FALSE \
--isCovariateOffset=FALSE \
--isCovariateTransform=TRUE \
--skipModelFitting=FALSE \
--tol=0.00001 \
--plinkFile=/home/users/nus/e1124313/scratch/eqtl/input/${cell_type}/chr/ATGC_${cell_type}_Asian_sle \
--nThreads=1 \
--IsOverwriteVarianceRatioFile=TRUE >> /home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step1/${cell_type}/oe/step1_${cell_type}_${gene_id}.log
