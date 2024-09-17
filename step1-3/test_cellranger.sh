#!/bin/bash
# export PATH=/opt/pbs/bin/:$PATH

export PATH=/usr/bin:/bin:/usr/sbin:/sbin:/home/users/nus/e1124313/scratch/eqtl/step1-3


# 文件夹路径
PBS_FOLDER="/home/users/nus/e1124313/project/ly/share/cellranger"
STATE_FILE="/home/users/nus/e1124313/scratch/eqtl/step1-3/submitted_jobs.txt"

# 如果状态文件不存在，创建一个空文件
if [ ! -f "$STATE_FILE" ]; then
  touch "$STATE_FILE"
fi

# 获取当前正在运行的任务数量
current_jobs=$(qstat | awk 'NR>2' | wc -l)
# current_jobs=$((current_jobs - 2))
max_jobs=90

# 遍历PBS脚本文件夹并提交每个PBS脚本
for pbs_file in "$PBS_FOLDER"/*.pbs; do
  job_id=$(basename "$pbs_file" .pbs)

  # 检查当前正在运行的任务数量是否小于最大任务数量
  if [ "$current_jobs" -lt "$max_jobs" ]; then
    # 检查任务是否已经提交
    if ! grep -q "$job_id" "$STATE_FILE"; then
      # 提交PBS脚本
      qsub "$pbs_file"
      
      # 增加当前运行的任务数量
      current_jobs=$((current_jobs + 1))

      # 将任务ID保存到状态文件
      echo "$job_id" >> "$STATE_FILE"
    fi
  else
    # 如果达到最大任务数量，退出循环
    break
  fi
done
