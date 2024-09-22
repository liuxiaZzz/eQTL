#!/bin/bash
    mkdir -p /home/users/nus/e1124313/scratch/eqtl/0920_input/T4_naive/chr/
    #cp /home/users/nus/e1124313/scratch/eqtl/raw_plink/split/chr* /home/users/nus/e1124313/scratch/eqtl/0920_input/T4_naive/chr/
    plink --bfile /home/users/nus/e1124313/scratch/eqtl/raw_plink/Asian_sle --keep /home/users/nus/e1124313/scratch/eqtl/0920_input/T4_naive/T4_naive_keep_samples.txt --make-bed --out /home/users/nus/e1124313/scratch/eqtl/0920_input/T4_naive/chr/ATGC_T4_naive_Asian_sle
    plink --bfile /home/users/nus/e1124313/scratch/eqtl/0920_input/T4_naive/chr/ATGC_T4_naive_Asian_sle --pca 20 --out /home/users/nus/e1124313/scratch/eqtl/0920_input/T4_naive/chr/ATGC_T4_naive_Asian_sle
    # 遍历所有22个染色体
    for i in {1..22}; do
        # 构建plink命令
        #plink --bfile /home/users/nus/e1124313/scratch/eqtl/0920_input/T4_naive/chr/chr$i --update-ids /home/users/nus/e1124313/scratch/eqtl/0920_input/T4_naive/T4_naive_updated_id.txt --make-bed --out /home/users/nus/e1124313/scratch/eqtl/0920_input/T4_naive/chr/T4_naive_id_updated_chr$i #临时使用
        plink --bfile /home/users/nus/e1124313/scratch/eqtl/raw_plink/update_id_split/Asian_sle_chr$i --keep /home/users/nus/e1124313/scratch/eqtl/0920_input/T4_naive/T4_naive_keep_samples.txt --make-bed --out /home/users/nus/e1124313/scratch/eqtl/0920_input/T4_naive/chr/ATGC_T4_naive_Asian_sle_chr$i 
    done
    