#!/bin/bash
export PATH=/opt/pbs/bin/:$PATH

# 细胞类型文件路径
CELL_TYPES_FILE="/home/users/nus/e1124313/scratch/eqtl/input/cell_types.csv"
BASIC_PATH="/home/users/nus/e1124313/scratch/eqtl/input/"
WORKDIR="/home/users/nus/e1124313/scratch/eqtl/"
STATE_FILE="/home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step2/submitted_jobs.txt"
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
  mkdir -p /home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step2/$cell_type
  mkdir -p /home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step2/$cell_type/output
  mkdir -p /home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step2/$cell_type/oe
  mkdir -p /home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step2/$cell_type/region
  # 读取基因列表并循环处理
  while IFS=$'\t' read -r gene_id chrom start end; do
    step2prefix=/home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step2/${cell_type}/output/cis_${cell_type}_$gene_id
    step1prefix=/home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step1/${cell_type}/output/${cell_type}_$gene_id
    region="${chrom}\t${start}\t${end}"
    regionFile=/home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step2/${cell_type}/region/${cell_type}_${gene_id}_range.txt
    echo -e $region > $regionFile

    job_id="cis_step2_${cell_type}_${gene_id}"
    # 检查当前正在运行的任务数量是否小于最大任务数量
    if [ "$current_jobs" -lt "$max_jobs" ]; then
      # 检查任务是否已经提交
      if ! grep -q "$job_id" "$STATE_FILE"; then
        # 生成PBS脚本内容
        pbs_script="#!/bin/bash
#PBS -N cis_step2_${cell_type}_$gene_id
#PBS -o /home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step2/${cell_type}/oe/step2_${cell_type}_${gene_id}.out
#PBS -e /home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step2/${cell_type}/oe/step2_${cell_type}_${gene_id}.err
#PBS -j oe
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=2:mem=24gb

export PATH=/home/project/11003054/e1124313/env/RSAIGE/bin/:$PATH

cd $WORKDIR

Rscript /home/users/nus/e1124313/scratch/qtl/extdata/step2_tests_qtl.R       \
        --bedFile=/home/users/nus/e1124313/scratch/eqtl/input/${cell_type}/chr/ATGC_${cell_type}_Asian_sle_chr${chrom}.bed      \
        --bimFile=/home/users/nus/e1124313/scratch/eqtl/input/${cell_type}/chr/ATGC_${cell_type}_Asian_sle_chr${chrom}.bim      \
        --famFile=/home/users/nus/e1124313/scratch/eqtl/input/${cell_type}/chr/ATGC_${cell_type}_Asian_sle_chr${chrom}.fam      \
        --SAIGEOutputFile=${step2prefix}     \
        --chrom=$chrom       \
        --minMAF=0 \
        --minMAC=20 \
        --LOCO=FALSE    \
        --GMMATmodelFile=${step1prefix}.rda     \
        --SPAcutoff=2 \
        --varianceRatioFile=${step1prefix}.varianceRatio.txt    \
        --rangestoIncludeFile=${regionFile}     \
        --markers_per_chunk=10000
"

        # 保存PBS脚本到文件
        pbs_file="/home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step2/${cell_type}/step2_${cell_type}_${gene_id}.pbs"
        echo "$pbs_script" > $pbs_file

        # 提交PBS脚本
        qsub $pbs_file
        
        # 增加当前运行的任务数量
        current_jobs=$((current_jobs + 1))
        # echo $current_jobs

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
# 强制转换为整数
start=${start%.*}
end=${end%.*}
start=$((start - 1000000))
end=$((end + 1000000))

step2prefix=/home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step2/${cell_type}/output/cis_${cell_type}_$gene_id
step1prefix=/home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step1/${cell_type}/output/${cell_type}_$gene_id
region="${chrom}\t${start}\t${end}"
regionFile=/home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step2/${cell_type}/region/${cell_type}_${gene_id}_range.txt

Rscript /home/users/nus/e1124313/scratch/qtl/extdata/step2_tests_qtl.R       \
--bedFile=/home/users/nus/e1124313/scratch/eqtl/input/${cell_type}/chr/ATGC_${cell_type}_Asian_sle_chr${chrom}.bed      \
--bimFile=/home/users/nus/e1124313/scratch/eqtl/input/${cell_type}/chr/ATGC_${cell_type}_Asian_sle_chr${chrom}.bim      \
--famFile=/home/users/nus/e1124313/scratch/eqtl/input/${cell_type}/chr/ATGC_${cell_type}_Asian_sle_chr${chrom}.fam      \
--SAIGEOutputFile=${step2prefix}     \
--chrom=$chrom       \
--minMAF=0.05 \
--minMAC=20 \
--LOCO=FALSE    \
--GMMATmodelFile=${step1prefix}.rda     \
--SPAcutoff=2 \
--varianceRatioFile=${step1prefix}.varianceRatio.txt    \
--rangestoIncludeFile=${regionFile}     \
--markers_per_chunk=10000
