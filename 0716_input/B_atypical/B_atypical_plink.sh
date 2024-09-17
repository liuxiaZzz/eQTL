#!/bin/bash
mkdir -p /home/users/nus/e1124313/scratch/eqtl/0716_input/B_atypical/chr/
#cp /home/users/nus/e1124313/scratch/eqtl/raw_plink/split/chr* /home/users/nus/e1124313/scratch/eqtl/0716_input/B_atypical/chr/
plink --bfile /home/users/nus/e1124313/scratch/eqtl/raw_plink/Asian_sle --keep /home/users/nus/e1124313/scratch/eqtl/0716_input/B_atypical/B_atypical_keep_samples.txt --make-bed --out /home/users/nus/e1124313/scratch/eqtl/0716_input/B_atypical/chr/ATGC_B_atypical_Asian_sle
plink --bfile /home/users/nus/e1124313/scratch/eqtl/0716_input/B_atypical/chr/ATGC_B_atypical_Asian_sle --pca 20 --out /home/users/nus/e1124313/scratch/eqtl/0716_input/B_atypical/chr/ATGC_B_atypical_Asian_sle
# 遍历所有22个染色体
for i in {1..22}; do
    # 构建plink命令
    #plink --bfile /home/users/nus/e1124313/scratch/eqtl/0716_input/B_atypical/chr/chr$i --update-ids /home/users/nus/e1124313/scratch/eqtl/0716_input/B_atypical/B_atypical_updated_id.txt --make-bed --out /home/users/nus/e1124313/scratch/eqtl/0716_input/B_atypical/chr/B_atypical_id_updated_chr$i #临时使用
    plink --bfile /home/users/nus/e1124313/scratch/eqtl/raw_plink/update_id_split/Asian_sle_chr$i --keep /home/users/nus/e1124313/scratch/eqtl/0716_input/B_atypical/B_atypical_keep_samples.txt --make-bed --out /home/users/nus/e1124313/scratch/eqtl/0716_input/B_atypical/chr/ATGC_B_atypical_Asian_sle_chr$i 
done
