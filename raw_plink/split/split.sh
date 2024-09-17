#!/bin/bash

#在这之前使用了change_bim.py

# # 定义输入文件的基本名称（不带扩展名）
# base_file="/home/users/nus/e1124313/scratch/eqtl/raw_plink/chr_all_ok"

# # 循环处理每条染色体
# for chr in {1..22}
# do
#     # 使用 plink 工具生成单独染色体文件
#     plink --bfile ${base_file} --chr ${chr} --make-bed --out /home/users/nus/e1124313/scratch/eqtl/raw_plink/split/chr${chr}

#     echo "Chromosome ${chr} files generated: chr${chr}.bed, chr${chr}.bim, chr${chr}.fam"
# done

# 定义输入文件的基本名称（不带扩展名）
# base_file="/home/users/nus/e1124313/scratch/eqtl/raw_plink/chr_all"

# # 循环处理每条染色体
# for chr in {1..22}
# do
#     # 使用 plink 工具生成单独染色体文件
#     plink --bfile ${base_file} --chr ${chr} --make-bed --out /home/users/nus/e1124313/scratch/eqtl/raw_plink/raw_split/chr${chr}

#     echo "Chromosome ${chr} files generated: chr${chr}.bed, chr${chr}.bim, chr${chr}.fam"
# done

# for chr in {1..22}
# do
#     # 使用 plink 工具生成单独染色体文件
#     plink --bfile /home/users/nus/e1124313/scratch/eqtl/raw_plink/update_id_split/keeped_chr_all --chr ${chr} --make-bed --out /home/users/nus/e1124313/scratch/eqtl/raw_plink/update_id_split/keeped_chr${chr}

#     echo "Chromosome ${chr} files generated: chr${chr}.bed, chr${chr}.bim, chr${chr}.fam"
# done

for chr in {1..22}
do
    # 使用 plink 工具生成单独染色体文件
    plink --bfile /home/users/nus/e1124313/scratch/eqtl/raw_plink/Asian_sle --chr ${chr} --make-bed --out /home/users/nus/e1124313/scratch/eqtl/raw_plink/update_id_split/Asian_sle_chr${chr}

    echo "Chromosome ${chr} files generated: chr${chr}.bed, chr${chr}.bim, chr${chr}.fam"
done
