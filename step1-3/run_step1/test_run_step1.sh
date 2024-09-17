#! /bin/bash

cell_type='myeloid'
gene_id='ENSG00000000419'

mkdir -p /home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step1/${cell_type}
mkdir -p /home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step1/${cell_type}/oe
mkdir -p /home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step1/${cell_type}/output

export PATH=/home/project/11003054/e1124313/env/RSAIGE/bin/:$PATH

# plink --bfile /home/users/nus/e1124313/scratch/eqtl/input/myeloid/chr/ATGC_myeloid_Asian_sle --make-rel square --out /home/users/nus/e1124313/scratch/eqtl/input/myeloid/chr/ATGC_myeloid_Asian_sle

# gcta64 --bfile /home/users/nus/e1124313/scratch/eqtl/raw_plink/Asian_sle --make-grm --maf 0.01 --out /home/users/nus/e1124313/scratch/eqtl/raw_plink/GRM/Asian_sle.grm

# Rscript /home/users/nus/e1124313/scratch/qtl/extdata/createSparseGRM.R       \
# --plinkFile=/home/users/nus/e1124313/scratch/eqtl/raw_plink/Asian_sle \
# --outputPrefix=/home/users/nus/e1124313/scratch/eqtl/raw_plink/GRM/Asian_sle     \
# --numRandomMarkerforSparseKin=2000      \
# --relatednessCutoff=0.125 \
# --nThreads=1

# conda activate GRM 没报错
mkdir -p /home/users/nus/e1124313/scratch/eqtl/input/${cell_type}/GRM
/home/users/nus/e1124313/project/e1124313/env/GRM/bin/createSparseGRM.R \
--plinkFile=/home/users/nus/e1124313/scratch/eqtl/input/myeloid/chr/ATGC_myeloid_Asian_sle \
--outputPrefix=/home/users/nus/e1124313/scratch/eqtl/raw_plink/GRM/ATGC_myeloid_Asian_sle     \
--numRandomMarkerforSparseKin=2000      \
--relatednessCutoff=0.125 \
--nThreads=2

plink --bfile /home/users/nus/e1124313/scratch/eqtl/input/myeloid/chr/ATGC_myeloid_Asian_sle --pca 20 --out /home/users/nus/e1124313/scratch/eqtl/input/myeloid/ATGC_myeloid_Asian_sle

Rscript /home/users/nus/e1124313/scratch/qtl/extdata/step1_fitNULLGLMM_qtl.R \
--useSparseGRMtoFitNULL=TRUE \
--sparseGRMSampleIDFile=/home/users/nus/e1124313/scratch/eqtl/raw_plink/GRM/Asian_sle_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
--sparseGRMFile=/home/users/nus/e1124313/scratch/eqtl/raw_plink/GRM/Asian_sle_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx \
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
--IsOverwriteVarianceRatioFile=TRUE > /home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step1/${cell_type}/oe/step1_${cell_type}_${gene_id}.log