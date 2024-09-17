#!/bin/bash
export PATH=/opt/pbs/bin/:$PATH

# 细胞类型文件路径
CELL_TYPES_FILE="/home/users/nus/e1124313/scratch/eqtl/input/cell_types.csv"
BASIC_PATH="/home/users/nus/e1124313/scratch/eqtl/input/"
WORKDIR="/home/users/nus/e1124313/scratch/eqtl/"
STATE_FILE="/home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step3/submitted_jobs.txt"
# 如果状态文件不存在，创建一个空文件
if [ ! -f "$STATE_FILE" ]; then
  touch "$STATE_FILE"
fi

# 获取当前正在运行的任务数量
current_jobs=$(qstat | wc -l)
# current_jobs=$((current_jobs - 2))
max_jobs=98

# 读取细胞类型列表并循环处理
while IFS= read -r cell_type; do
  GENE_LIST="${BASIC_PATH}${cell_type}/${cell_type}_gene_locations.txt"
  mkdir -p /home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step3/$cell_type
  mkdir -p /home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step3/$cell_type/output
  mkdir -p /home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step3/$cell_type/oe
  # 读取基因列表并循环处理
  while IFS=$'\t' read -r gene_id chrom start end; do
    step2prefix=/home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step2/${cell_type}/output/cis_${cell_type}_$gene_id

    job_id="cis_step3_${cell_type}_${gene_id}"
    # 检查当前正在运行的任务数量是否小于最大任务数量
    if [ "$current_jobs" -lt "$max_jobs" ]; then
      # 检查任务是否已经提交
      if ! grep -q "$job_id" "$STATE_FILE"; then
        # 生成PBS脚本内容
        pbs_script="#!/bin/bash
#PBS -N cis_step3_${cell_type}_$gene_id
#PBS -o /home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step3/${cell_type}/oe/step3_${cell_type}_${gene_id}.out
#PBS -e /home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step3/${cell_type}/oe/step3_${cell_type}_${gene_id}.err
#PBS -j oe
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=2:mem=128gb

export PATH=/home/project/11003054/e1124313/env/RSAIGE/bin/:$PATH

cd $WORKDIR

Rscript /home/users/nus/e1124313/scratch/qtl/extdata/step3_gene_pvalue_qtl.R \
        --assocFile=$step2prefix        \
        --geneName=$gene_id       \
        --genePval_outputFile=/home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step3/$cell_type/output/cis_${cell_type}_${gene_id}_pval
"

        # 保存PBS脚本到文件
        pbs_file="/home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step3/${cell_type}/step3_${cell_type}_${gene_id}.pbs"
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
  done < <(tail -n +2 "$GENE_LIST")  # 使用 process substitution 代替管道
done < "$CELL_TYPES_FILE"


cell_type='B_mem'
gene_id='ENSG00000150093'
chrom='10'
start='31887273'
end='34005792'

mkdir -p /home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step3/$cell_type
mkdir -p /home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step3/$cell_type/output
mkdir -p /home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step3/$cell_type/oe
step2prefix=/home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step2/${cell_type}/output/cis_${cell_type}_$gene_id
step3prefix=/home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step3/${cell_type}/output/cis_${cell_type}_$gene_id

Rscript /home/users/nus/e1124313/scratch/qtl/extdata/step3_gene_pvalue_qtl.R       \
--assocFile=$step2prefix \
--geneName=$gene_id \
--genePval_outputFile=${step3prefix}_pval > /home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step3/${cell_type}/oe/step3_${cell_type}_${gene_id}.log


