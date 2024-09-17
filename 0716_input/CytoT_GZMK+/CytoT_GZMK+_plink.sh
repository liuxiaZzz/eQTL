#!/bin/bash
mkdir -p /home/users/nus/e1124313/scratch/eqtl/0716_input/CytoT_GZMK+/chr/
#cp /home/users/nus/e1124313/scratch/eqtl/raw_plink/split/chr* /home/users/nus/e1124313/scratch/eqtl/0716_input/CytoT_GZMK+/chr/
plink --bfile /home/users/nus/e1124313/scratch/eqtl/raw_plink/Asian_sle --keep /home/users/nus/e1124313/scratch/eqtl/0716_input/CytoT_GZMK+/CytoT_GZMK+_keep_samples.txt --make-bed --out /home/users/nus/e1124313/scratch/eqtl/0716_input/CytoT_GZMK+/chr/ATGC_CytoT_GZMK+_Asian_sle
plink --bfile /home/users/nus/e1124313/scratch/eqtl/0716_input/CytoT_GZMK+/chr/ATGC_CytoT_GZMK+_Asian_sle --pca 20 --out /home/users/nus/e1124313/scratch/eqtl/0716_input/CytoT_GZMK+/chr/ATGC_CytoT_GZMK+_Asian_sle
# 遍历所有22个染色体
for i in {1..22}; do
    # 构建plink命令
    #plink --bfile /home/users/nus/e1124313/scratch/eqtl/0716_input/CytoT_GZMK+/chr/chr$i --update-ids /home/users/nus/e1124313/scratch/eqtl/0716_input/CytoT_GZMK+/CytoT_GZMK+_updated_id.txt --make-bed --out /home/users/nus/e1124313/scratch/eqtl/0716_input/CytoT_GZMK+/chr/CytoT_GZMK+_id_updated_chr$i #临时使用
    plink --bfile /home/users/nus/e1124313/scratch/eqtl/raw_plink/update_id_split/Asian_sle_chr$i --keep /home/users/nus/e1124313/scratch/eqtl/0716_input/CytoT_GZMK+/CytoT_GZMK+_keep_samples.txt --make-bed --out /home/users/nus/e1124313/scratch/eqtl/0716_input/CytoT_GZMK+/chr/ATGC_CytoT_GZMK+_Asian_sle_chr$i 
done
