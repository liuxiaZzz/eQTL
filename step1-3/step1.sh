#!/bin/bash
echo $(hostname) >> $HOME/hostname.txt

crontab -l | grep -v 'step1.sh' | crontab -

export PATH=/opt/pbs/bin/:$PATH

# 细胞类型文件路径
CELL_TYPES_FILE="/home/users/nus/e1124313/scratch/eqtl/input/cell_types.csv"
BASIC_PATH="/home/users/nus/e1124313/scratch/eqtl/input/"
WORKDIR="/home/users/nus/e1124313/scratch/eqtl/"
STATE_FILE="/home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step1/submitted_jobs.txt"
# 如果状态文件不存在，创建一个空文件
if [ ! -f "$STATE_FILE" ]; then
  touch "$STATE_FILE"
fi

# 获取当前正在运行的任务数量
current_jobs=$(qstat | wc -l)
current_jobs=$((current_jobs - 2))
max_jobs=98

# 读取细胞类型列表并循环处理
while IFS= read -r cell_type; do
  GENE_LIST="${BASIC_PATH}${cell_type}/${cell_type}_var_names.txt"
  mkdir -p /home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step1/$cell_type
  mkdir -p /home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step1/$cell_type/output
  mkdir -p /home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step1/$cell_type/oe
  # 读取基因列表并循环处理
  while IFS= read -r gene_id; do
    job_id="step1_${cell_type}_${gene_id}"
    # 检查当前正在运行的任务数量是否小于最大任务数量
    if [ "$current_jobs" -lt "$max_jobs" ]; then
      # 检查任务是否已经提交
      if ! grep -q "$job_id" "$STATE_FILE"; then
        # 生成PBS脚本内容
        pbs_script="#!/bin/bash
#PBS -N step1_${cell_type}_$gene_id
#PBS -o /home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step1/${cell_type}/oe/step1_${cell_type}_${gene_id}.out
#PBS -e /home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step1/${cell_type}/oe/step1_${cell_type}_${gene_id}.err
#PBS -j oe
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=2:mem=32gb

export PATH=/home/project/11003054/e1124313/env/RSAIGE/bin/:$PATH

cd $WORKDIR

echo \"cell_type: ${cell_type}\gene_id: ${gene_id}\"

Rscript /home/users/nus/e1124313/scratch/qtl/extdata/step1_fitNULLGLMM_qtl.R \
        --useSparseGRMtoFitNULL=FALSE  \
        --useGRMtoFitNULL=FALSE \
        --phenoFile=/home/users/nus/e1124313/scratch/eqtl/input/${cell_type}/${cell_type}.txt \
        --phenoCol=$gene_id \
        --covarColList=Age \
        --sampleCovarColList=Age \
        --sampleIDColinphenoFile=Sample_ID \
        --traitType=count \
        --outputPrefix=/home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step1/${cell_type}/output/${cell_type}_$gene_id \
        --skipVarianceRatioEstimation=FALSE \
        --isRemoveZerosinPheno=FALSE \
        --isCovariateOffset=FALSE \
        --isCovariateTransform=TRUE \
        --skipModelFitting=FALSE \
        --tol=0.00001 \
        --plinkFile=/home/users/nus/e1124313/scratch/eqtl/input/${cell_type}/chr/${cell_type}_id_updated_ok_chr2 \
        --nThreads=24  \
        --IsOverwriteVarianceRatioFile=TRUE
"

        # 保存PBS脚本到文件
        pbs_file="/home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step1/${cell_type}/step1_${cell_type}_${gene_id}.pbs"
        echo "$pbs_script" > $pbs_file

        # 提交PBS脚本
        qsub $pbs_file
        
        # 增加当前运行的任务数量
        current_jobs=$((current_jobs + 1))

        rm $pbs_file

        # 将任务ID保存到状态文件
        echo "$job_id" >> "$STATE_FILE"
      fi
    else
      # 如果达到最大任务数量，退出循环
      break
    fi 
  done < "$GENE_LIST"
done < "$CELL_TYPES_FILE"


## 0528报错找不到基因，和之前一样，不知道为什么大的矩阵好像就不行，取500个基因就可以跑
## test1：换chr2测试， 结果失败
## test2：我怀疑是存在重复基因，但之前看过其实没有，adata.var_names_make_unique()，
# 于是重新跑了一遍preparation.py，结果发现Asian.txt从1GB多变成了3GB，说明在输出txt时有概率输出残缺文件，原因不明

## 打算不使用adata_sampled，或许它自己就不完整，改用成功过的

## 0529 重新用adata_sampled.h5ad生成0528_combined_data.txt，结果和combined_data.txt一样，说明adata_sampled.h5ad没有问题
# 那么问题应该出在preparation.py，我准备去掉进程池，改用普通循环， 普通循环也失败了
# 0529晚上，更改源代码发现没效果，开始重装环境，同时看到一个github，或许将txt分成多个小文件可行，后续准备试试
# 0529晚上，问题出在Raw barcode中间的空格，导致读取时出错

# 0603又出现新bug，直接将obs_list改为手写的特定列


cell_type='B_mem'
gene_id='ENSG00000150093'
Rscript /home/users/nus/e1124313/scratch/qtl/extdata/step1_fitNULLGLMM_qtl.R \
--useSparseGRMtoFitNULL=TRUE  \
--useGRMtoFitNULL=FALSE \
--phenoFile=/home/users/nus/e1124313/scratch/eqtl/0720_input/${cell_type}/${cell_type}.txt \
--phenoCol=$gene_id \
--covarColList=Age,Sex,geno_PC1,geno_PC2,geno_PC3,geno_PC4,geno_PC5,geno_PC6 \
--sampleCovarColList=Age,Sex,geno_PC1,geno_PC2,geno_PC3,geno_PC4,geno_PC5,geno_PC6 \
--sampleIDColinphenoFile=Sample_ID \
--traitType=count \
--outputPrefix=/home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step1/${cell_type}/output/${cell_type}_$gene_id \
--skipVarianceRatioEstimation=FALSE \
--isRemoveZerosinPheno=FALSE \
--isCovariateOffset=FALSE \
--isCovariateTransform=TRUE \
--skipModelFitting=FALSE \
--tol=0.00001 \
--plinkFile=/home/users/nus/e1124313/scratch/eqtl/0720_input/${cell_type}/chr/ATGC_${cell_type}_Asian_sle \
--nThreads=1 \
--IsOverwriteVarianceRatioFile=TRUE > /home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step1/${cell_type}/oe/step1_${cell_type}_${gene_id}.log


