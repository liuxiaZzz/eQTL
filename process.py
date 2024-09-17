'''
:@Author: liuxia
:@Date: 5/22/2024, 11:26:36 AM
:@LastEditors: liuxia
:@LastEditTime: 5/22/2024, 11:26:36 AM
:Description: 
'''
# import re
# import csv

# # 打开文件并读取所有行
# with open('/home/users/nus/e1124313/scratch/eqtl/GSE196829_series_matrix.txt', 'r') as f:
#     lines = f.readlines()

# # 提取第41，52，65行
# line41 = lines[40].strip().replace('"', '').split('\t')
# line52 = lines[51].strip().replace('"', '').split('\t')
# line65 = lines[64].strip().replace('"', '').split('\t')

# # 先选取/分割后的最后一个元素，再选取“_Grn"前面的元素
# processed_elements = [element.split('/')[-1].split('_Grn')[0] for element in line65]

# # 将结果写入到一个新的csv文件中
# with open('/home/users/nus/e1124313/scratch/eqtl/output.csv', 'w', newline='') as f:
#     writer = csv.writer(f)
#     # 分别写入三行
#     writer.writerow(line41)
#     writer.writerow(line52)
#     writer.writerow(processed_elements)

#############
import os
import csv

# 文件路径
input_file_path = '/home/users/nus/e1124313/scratch/eqtl/GSE196829_series_matrix.txt'
output_file_path = '/home/users/nus/e1124313/scratch/eqtl/output.csv'

# 要提取的行索引
lines_to_extract = [40, 51]  # 这些是基于0索引的行号
special_lines_to_extract = [64]

with open(input_file_path, 'r') as f:
    lines = f.readlines()

# 提取行
extracted_lines = [lines[i].strip().replace('"', '').split('\t') for i in lines_to_extract]

processed_special_lines = []
for i in special_lines_to_extract:
    special_line = lines[i].strip().replace('"', '').split('\t')
    processed_elements = [element.split('/')[-1].split('_Grn')[0] for element in special_line]
    processed_special_lines.append(processed_elements)

with open(output_file_path, 'w', newline='') as f:
    writer = csv.writer(f)
    # 写入普通行
    for line in extracted_lines:
        writer.writerow(line)
    # 写入处理过的特殊行
    for processed_elements in processed_special_lines:
        writer.writerow(processed_elements)

#############
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

adata = sc.read("/home/users/nus/e1124313/project/ly/3/data/combined_3.h5ad")
# 使用 value_counts 统计每个 sample_ID 的细胞数量
cell_counts_per_sample = adata.obs['sample_ID'].value_counts()
# 绘制直方图
plt.figure(figsize=(10, 6))
cell_counts_per_sample.plot(kind='bar')
plt.xlabel('Sample ID')
plt.ylabel('Number of Cells')
plt.title('Number of Cells per Sample ID')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig('/home/users/nus/e1124313/scratch/eqtl/cell_counts_per_sample.pdf', dpi=300, bbox_inches='tight')

############# 取30万细胞，更改plink文件
import numpy as np
np.random.seed(42)
sample_size = 300000
sampled_indices = np.random.choice(adata.n_obs, size=sample_size, replace=False)

adata_sampled = adata[sampled_indices, :]
adata_sampled.obs.rename(columns={'Sample ID': 'Sample_ID', 'Batch ID': 'Batch_ID'}, inplace=True)
adata_sampled.write('/home/users/nus/e1124313/scratch/eqtl/adata_sampled.h5ad')

fam_df = pd.read_csv('/home/users/nus/e1124313/scratch/eqtl/chr2_ok.fam', sep=' ', header=None)
num_samples = adata_sampled.obs['Sample_ID'].nunique()
subset_fam_df = fam_df.iloc[:num_samples]
# 提取旧的 ID 和新的 ID
old_fids = subset_fam_df[0]
old_iids = subset_fam_df[1]
new_iids = adata_sampled.obs['Sample_ID'].unique()  # 这将作为新的 Sample ID
# 创建 update_ids.txt 文件内容
update_ids = pd.DataFrame({
    "old_fid": old_fids,
    "old_iid": old_iids,
    "new_fid": new_iids,  
    "new_iid": new_iids  # 新的 Individual ID
})
update_ids_file = '/home/users/nus/e1124313/scratch/eqtl/update_ids.txt'
update_ids.to_csv(update_ids_file, sep="\t", index=False, header=False)

keep_ids = pd.DataFrame({
    "keep_ID": new_iids,
    "keep_fID": new_iids
})
keep_ids_file = '/home/users/nus/e1124313/scratch/eqtl/keep.txt'
keep_ids.to_csv(keep_ids_file, sep="\t", index=False, header=False)

plink --bfile /home/users/nus/e1124313/scratch/eqtl/chr2_ok --update-ids /home/users/nus/e1124313/scratch/eqtl/update_ids.txt --make-bed --out /home/users/nus/e1124313/scratch/eqtl/id_updated_chr2_ok
plink --bfile /home/users/nus/e1124313/scratch/eqtl/id_updated_chr2_ok --keep keep.txt --make-bed --out /home/users/nus/e1124313/scratch/eqtl/id_updated_chr2_ok
#还需要

############### 生成表达矩阵
adata_sampled = sc.read('/home/users/nus/e1124313/scratch/eqtl/adata_sampled.h5ad')
sample_ids = adata_sampled.obs['Sample_ID']
covariates = adata_sampled.obs[['Batch_ID', 'Age', 'Sex', 'Ethnicity_1', 'Ethnicity_2', 'Status', 'Disease_status', 'Assay_method', 'CT_1', 'CT_2']]
gene_data = pd.DataFrame(adata_sampled.X.toarray(), columns=adata_sampled.var_names)
# gene_data = pd.DataFrame(adata_sampled.X.toarray()[:,:500], columns=adata_sampled.var_names[:500])
# gene_data = pd.DataFrame(adata_sampled.X.toarray()[:,-500:], columns=adata_sampled.var_names[-500:])
gene_data.index = sample_ids.index
# gene_data_head = gene_data.head(5)

# 根据基因拆分成多个文件
# for i in len(gene_data.columns):
#     combined_df = pd.concat([sample_ids, covariates, gene_data.iloc[:, i]], axis=1)
#     combined_df.rename(columns={'Sample ID': 'IND_ID'}, inplace=True)
#     output_file = f'/home/users/nus/e1124313/scratch/eqtl/input/combined_data_gene_{i}.txt'
#     combined_df.to_csv(output_file, sep='\t', index=False)

combined_df = pd.concat([sample_ids, covariates, gene_data], axis=1)
# combined_df.rename(columns={'Sample ID': 'IND_ID'}, inplace=True) #####一定要做，反正不能有空格
output_file = '/home/users/nus/e1124313/scratch/eqtl/input/combined_data.txt'
combined_df.to_csv(output_file, sep='\t', index=False)

df = pd.read_csv('/home/users/nus/e1124313/scratch/eqtl/input/combined_data.txt', sep='\t')
############### Rscrpit Step1
## 这个也成功了
Rscript /home/users/nus/e1124313/project/e1124313/qtl/extdata/step1_fitNULLGLMM_qtl.R \
        --useSparseGRMtoFitNULL=FALSE  \
        --useGRMtoFitNULL=FALSE \
        --phenoFile=/home/users/nus/e1124313/scratch/eqtl/input/combined_data.txt	\
        --phenoCol=ENSG00000273492       \
        --covarColList=Age    \
        --sampleCovarColList=Age      \
        --sampleIDColinphenoFile=Sample_ID \
        --traitType=count \
        --outputPrefix=/home/users/nus/e1124313/scratch/eqtl/output/combined_data_ENSG00000273492 \
        --skipVarianceRatioEstimation=FALSE  \
        --isRemoveZerosinPheno=FALSE \
        --isCovariateOffset=FALSE  \
        --isCovariateTransform=TRUE  \
        --skipModelFitting=FALSE  \
        --tol=0.00001   \
        --plinkFile=/home/users/nus/e1124313/scratch/eqtl/id_updated_chr2_ok       \
        --nThreads=24  \
        --IsOverwriteVarianceRatioFile=TRUE
       
# 这个成功了
Rscript /home/users/nus/e1124313/project/e1124313/qtl/extdata/step1_fitNULLGLMM_qtl.R \
        --useSparseGRMtoFitNULL=FALSE  \
        --useGRMtoFitNULL=FALSE \
        --phenoFile=/home/users/nus/e1124313/scratch/eqtl/input/test_combined_data.txt	\
        --phenoCol=ENSG00000000003       \
        --covarColList=Age    \
        --sampleCovarColList=Age      \
        --sampleIDColinphenoFile=Sample_ID \
        --traitType=count \
        --outputPrefix=/home/users/nus/e1124313/scratch/eqtl/output/combined_data_gene1 \
        --skipVarianceRatioEstimation=FALSE  \
        --isRemoveZerosinPheno=FALSE \
        --isCovariateOffset=FALSE  \
        --isCovariateTransform=TRUE  \
        --skipModelFitting=FALSE  \
        --tol=0.00001   \
        --plinkFile=/home/users/nus/e1124313/scratch/eqtl/id_updated_chr2_ok       \
        --nThreads=24  \
        --IsOverwriteVarianceRatioFile=TRUE

################ 这个也成功了，我认为就是文件完整不完整的问题
Rscript /home/users/nus/e1124313/project/e1124313/qtl/extdata/step1_fitNULLGLMM_qtl.R \
        --useSparseGRMtoFitNULL=FALSE  \
        --useGRMtoFitNULL=FALSE \
        --phenoFile=/home/users/nus/e1124313/scratch/eqtl/input/combined_data.txt	\
        --phenoCol=ENSG00000000003       \
        --covarColList=Age    \
        --sampleCovarColList=Age      \
        --sampleIDColinphenoFile=Sample_ID \
        --traitType=count \
        --outputPrefix=/home/users/nus/e1124313/scratch/eqtl/output/combined_data_ENSG00000000003 \
        --skipVarianceRatioEstimation=FALSE  \
        --isRemoveZerosinPheno=FALSE \
        --isCovariateOffset=FALSE  \
        --isCovariateTransform=TRUE  \
        --skipModelFitting=FALSE  \
        --tol=0.00001   \
        --plinkFile=/home/users/nus/e1124313/scratch/eqtl/id_updated_chr2_ok       \
        --nThreads=24  \
        --IsOverwriteVarianceRatioFile=TRUE

### 换combined_data就成功，换Asian.txt就失败，Asian.txt=2.84GB
Rscript /home/users/nus/e1124313/project/e1124313/qtl/extdata/step1_fitNULLGLMM_qtl.R \
        --useSparseGRMtoFitNULL=FALSE  \
        --useGRMtoFitNULL=FALSE \
        --phenoFile=/home/users/nus/e1124313/scratch/eqtl/input/combined_data.txt	\
        --phenoCol=ENSG00000000003       \
        --covarColList=Age    \
        --sampleCovarColList=Age      \
        --sampleIDColinphenoFile=Sample_ID \
        --traitType=count \
        --outputPrefix=/home/users/nus/e1124313/scratch/eqtl/output/Asian_ENSG00000000003 \
        --skipVarianceRatioEstimation=FALSE  \
        --isRemoveZerosinPheno=FALSE \
        --isCovariateOffset=FALSE  \
        --isCovariateTransform=TRUE  \
        --skipModelFitting=FALSE  \
        --tol=0.00001   \
        --plinkFile=/home/users/nus/e1124313/scratch/eqtl/input/Asian/chr/Asian_id_updated_ok_chr2       \
        --nThreads=24  \
        --IsOverwriteVarianceRatioFile=TRUE

### 换0528_combined_data也成功，说明0528_combined_data和combined_data是一样的
Rscript /home/users/nus/e1124313/project/e1124313/qtl/extdata/step1_fitNULLGLMM_qtl.R \
        --useSparseGRMtoFitNULL=FALSE  \
        --useGRMtoFitNULL=FALSE \
        --phenoFile=/home/users/nus/e1124313/scratch/eqtl/input/0528_combined_data.txt  \
        --phenoCol=ENSG00000000003       \
        --covarColList=Age    \
        --sampleCovarColList=Age      \
        --sampleIDColinphenoFile=Sample_ID \
        --traitType=count \
        --outputPrefix=/home/users/nus/e1124313/scratch/eqtl/output/combined_data_ENSG00000000003 \
        --skipVarianceRatioEstimation=FALSE  \
        --isRemoveZerosinPheno=FALSE \
        --isCovariateOffset=FALSE  \
        --isCovariateTransform=TRUE  \
        --skipModelFitting=FALSE  \
        --tol=0.00001   \
        --plinkFile=/home/users/nus/e1124313/scratch/eqtl/id_updated_chr2_ok       \
        --nThreads=24  \
        --IsOverwriteVarianceRatioFile=TRUE

### 注意！0529下午重新装了一次包，所以上面的路径是老路径
### 换Asian.txt还是失败
Rscript /home/users/nus/e1124313/scratch/qtl/extdata/step1_fitNULLGLMM_qtl.R \
        --useSparseGRMtoFitNULL=FALSE  \
        --useGRMtoFitNULL=FALSE \
        --phenoFile=/home/users/nus/e1124313/scratch/eqtl/input/Asian/Asian.txt  \
        --phenoCol=ENSG00000000419       \
        --covarColList=Age    \
        --sampleCovarColList=Age      \
        --sampleIDColinphenoFile=Sample_ID \
        --traitType=count \
        --outputPrefix=/home/users/nus/e1124313/scratch/eqtl/output/Asian_ENSG00000000419 \
        --skipVarianceRatioEstimation=FALSE  \
        --isRemoveZerosinPheno=FALSE \
        --isCovariateOffset=FALSE  \
        --isCovariateTransform=TRUE  \
        --skipModelFitting=FALSE  \
        --tol=0.00001   \
        --plinkFile=/home/users/nus/e1124313/scratch/eqtl/input/Asian/chr/Asian_id_updated_ok_chr1       \
        --nThreads=24  \
        --IsOverwriteVarianceRatioFile=TRUE

system.file(package = "SAIGEQTL")

### 现在试一下不把allele变成1 2可以不
plink --bfile /home/users/nus/e1124313/scratch/eqtl/raw_plink/chr_all --update-ids /home/users/nus/e1124313/scratch/eqtl/rubbish/update_ids.txt --make-bed --out /home/users/nus/e1124313/scratch/eqtl/raw_plink/update_id_split/id_updated_chr_all
plink --bfile /home/users/nus/e1124313/scratch/eqtl/raw_plink/update_id_split/id_updated_chr_all --keep /home/users/nus/e1124313/scratch/eqtl/rubbish/keep.txt --make-bed --out /home/users/nus/e1124313/scratch/eqtl/raw_plink/update_id_split/keeped_chr_all
for chr in {1..22}
do
    # 使用 plink 工具生成单独染色体文件
    plink --bfile /home/users/nus/e1124313/scratch/eqtl/raw_plink/update_id_split/keeped_chr_all --chr ${chr} --make-bed --out /home/users/nus/e1124313/scratch/eqtl/raw_plink/update_id_split/keeped_chr${chr}

    echo "Chromosome ${chr} files generated: chr${chr}.bed, chr${chr}.bim, chr${chr}.fam"
done

#### 测试不把allele换成12，成功了
Rscript /home/users/nus/e1124313/scratch/qtl/extdata/step1_fitNULLGLMM_qtl.R \
        --useSparseGRMtoFitNULL=FALSE  \
        --useGRMtoFitNULL=FALSE \
        --phenoFile=/home/users/nus/e1124313/scratch/eqtl/input/0528_combined_data.txt  \
        --phenoCol=ENSG00000000003       \
        --covarColList=Age    \
        --sampleCovarColList=Age      \
        --sampleIDColinphenoFile=Sample_ID \
        --traitType=count \
        --outputPrefix=/home/users/nus/e1124313/scratch/eqtl/output/combined_data_ENSG00000000003 \
        --skipVarianceRatioEstimation=FALSE  \
        --isRemoveZerosinPheno=FALSE \
        --isCovariateOffset=FALSE  \
        --isCovariateTransform=TRUE  \
        --skipModelFitting=FALSE  \
        --tol=0.00001   \
        --plinkFile=/home/users/nus/e1124313/scratch/eqtl/raw_plink/update_id_split/keeped_chr1       \
        --nThreads=24  \
        --IsOverwriteVarianceRatioFile=TRUE

#### 测试手动选取obs_names
Rscript /home/users/nus/e1124313/scratch/qtl/extdata/step1_fitNULLGLMM_qtl.R \
        --useSparseGRMtoFitNULL=FALSE  \
        --useGRMtoFitNULL=FALSE \
        --phenoFile=/home/users/nus/e1124313/scratch/eqtl/input/Asian/Asian.txt \
        --phenoCol=ENSG00000000419 \
        --covarColList=Age \
        --sampleCovarColList=Age \
        --sampleIDColinphenoFile=Sample_ID \
        --traitType=count \
        --outputPrefix=/home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step1/Asian/output/Asian_ENSG00000000419 \
        --skipVarianceRatioEstimation=FALSE \
        --isRemoveZerosinPheno=FALSE \
        --isCovariateOffset=FALSE \
        --isCovariateTransform=TRUE \
        --skipModelFitting=FALSE \
        --tol=0.00001 \
        --plinkFile=/home/users/nus/e1124313/scratch/eqtl/input/Asian/chr/Asian_id_updated_ok_chr2 \
        --nThreads=24  \
        --IsOverwriteVarianceRatioFile=TRUE

### 测试是否跟step1选取的染色体有关,结果是有一点变化
Rscript /home/users/nus/e1124313/scratch/qtl/extdata/step1_fitNULLGLMM_qtl.R \
        --useSparseGRMtoFitNULL=FALSE  \
        --useGRMtoFitNULL=FALSE \
        --phenoFile=/home/users/nus/e1124313/scratch/eqtl/input/Asian/Asian.txt \
        --phenoCol=ENSG00000000419 \
        --covarColList=Age \
        --sampleCovarColList=Age \
        --sampleIDColinphenoFile=Sample_ID \
        --traitType=count \
        --outputPrefix=/home/users/nus/e1124313/scratch/eqtl/output/step1_chr20_Asian_ENSG00000000419 \
        --skipVarianceRatioEstimation=FALSE \
        --isRemoveZerosinPheno=FALSE \
        --isCovariateOffset=FALSE \
        --isCovariateTransform=TRUE \
        --skipModelFitting=FALSE \
        --tol=0.00001 \
        --plinkFile=/home/users/nus/e1124313/scratch/eqtl/input/Asian/chr/Asian_id_updated_ok_chr20 \
        --nThreads=24  \
        --IsOverwriteVarianceRatioFile=TRUE

Rscript /home/users/nus/e1124313/scratch/qtl/extdata/step2_tests_qtl.R       \
        --bedFile=/home/users/nus/e1124313/scratch/eqtl/input/Asian/chr/Asian_id_updated_ok_chr20.bed      \
        --bimFile=/home/users/nus/e1124313/scratch/eqtl/input/Asian/chr/Asian_id_updated_ok_chr20.bim      \
        --famFile=/home/users/nus/e1124313/scratch/eqtl/input/Asian/chr/Asian_id_updated_ok_chr20.fam      \
        --SAIGEOutputFile=/home/users/nus/e1124313/scratch/eqtl/output/step2_chr20_Asian_ENSG00000000419     \
        --chrom=20       \
        --minMAF=0 \
        --minMAC=20 \
        --LOCO=FALSE    \
        --GMMATmodelFile=/home/users/nus/e1124313/scratch/eqtl/output/step1_chr20_Asian_ENSG00000000419.rda     \
        --SPAcutoff=2 \
        --varianceRatioFile=/home/users/nus/e1124313/scratch/eqtl/output/step1_chr20_Asian_ENSG00000000419.varianceRatio.txt    \
        --rangestoIncludeFile=/home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step2/Asian/region/Asian_ENSG00000000419_range.txt     \
        --markers_per_chunk=10000

### 测试conditional analysis
Rscript /home/users/nus/e1124313/scratch/qtl/extdata/step2_tests_qtl.R       \
        --bedFile=/home/users/nus/e1124313/scratch/eqtl/input/Asian/chr/ATGC_Asian_id_updated_ok_chr20.bed      \
        --bimFile=/home/users/nus/e1124313/scratch/eqtl/input/Asian/chr/ATGC_Asian_id_updated_ok_chr20.bim      \
        --famFile=/home/users/nus/e1124313/scratch/eqtl/input/Asian/chr/ATGC_Asian_id_updated_ok_chr20.fam      \
        --SAIGEOutputFile=/home/users/nus/e1124313/scratch/eqtl/output/conditional_step2_chr20_Asian_ENSG00000000419     \
        --chrom=20       \
        --minMAF=0 \
        --minMAC=20 \
        --LOCO=FALSE    \
        --GMMATmodelFile=/home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step1/Asian/output/Asian_ENSG00000000419.rda    \
        --SPAcutoff=2 \
        --varianceRatioFile=/home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step1/Asian/output/Asian_ENSG00000000419.varianceRatio.txt    \
        --rangestoIncludeFile=/home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step2/Asian/region/Asian_ENSG00000000419_range.txt     \
        --markers_per_chunk=10000 \
        --condition=20:50936367:A:G,20:50937815:A:T


filtered_adata = adata[(adata.obs['Ethnicity_1'] == 'Asian') & (adata.obs['Status'] == 'SLE')]
filtered_adata.obs['CT_2'] = filtered_adata.obs['CT_2'].cat.add_categories('myeloid').fillna('myeloid')
filtered_adata.write('/home/users/nus/e1124313/scratch/eqtl/asian_sle.h5ad')
qsub -I -l select=1:ncpus=12:mem=1024g -l walltime=24:00:00 -P 11003054 -q normal

import pandas as pd

df = pd.read_csv('/home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step2/B_mem/output/cis_B_mem_ENSG00000150093', sep='\t')
column_data = df['p.value']
count_less_than_005 = (column_data < 0.05).sum()
total_count = len(column_data)
percentage = (count_less_than_005 / total_count) * 100

print(f"小于 0.05 的值占比: {percentage:.2f}%")


