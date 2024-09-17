'''
:@Author: liuxia
:@Date: 5/23/2024, 4:00:32 PM
:@LastEditors: liuxia
:@LastEditTime: 5/23/2024, 4:00:32 PM
:Description: 最新更新的功能：添加了归一化和regress，并且使用plink计算了genotype的PCA，输出到保存到txt和h5ad中，后续用于计算PEER
下一步应该是计算GRM和PEER
'''
# conda activate go
# module load hdf5/1.12.1-nvhpc
# module load openmpi/4.1.5-gcc11
######### read adata

# import rpy2.rinterface_lib.callbacks
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri                       #先加载这些再加载后面的scanpy那些就不会报错，不知道为什么
from rpy2.robjects.conversion import localconverter
# from rpy2.robjects.packages import importr
import anndata2ri
pandas2ri.activate()
anndata2ri.activate()

import scanpy as sc
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
import numpy as np
import os
import subprocess
import gffutils
from pyensembl import EnsemblRelease
import scipy.sparse as sp

print("Start reading data...")
adata = sc.read("/home/users/nus/e1124313/scratch/eqtl/asian_sle.h5ad") #替换成完整数据路径
adata.obs.rename(columns={'Sample ID': 'Sample_ID', 'Batch ID': 'Batch_ID'}, inplace=True)
adata.var_names_make_unique() #测试能否解决找不到基因的问题
adata.obs['Sex'] = adata.obs['Sex'].map({'FEMALE': 0, 'MALE': 1})


############### 临时，用后删
# adata.obs['cell_type'] = adata.obs['Ethnicity_1']
adata.obs.rename(columns={'Raw Barcode': 'Raw_Barcode'}, inplace=True)
SAMPLE_LIST_doc = pd.read_csv("/home/users/nus/e1124313/scratch/eqtl/raw_plink/Asian_sle.fam", sep=' ', header=None)
SAMPLE_LIST = SAMPLE_LIST_doc[0].tolist()
adata = adata[adata.obs['Sample_ID'].isin(SAMPLE_LIST)]
###############

# cell_types = adata.obs['cell_type'].unique()
cell_types = adata.obs['CT_2'].unique()
# Save cell_types as a CSV file
cell_types_df = pd.DataFrame(cell_types, columns=['cell_type'])
cell_types_df.to_csv('/home/users/nus/e1124313/scratch/eqtl/input/cell_types.csv', index=False, header=False)

subsets = {ctype: adata[adata.obs['CT_2'] == ctype].copy() for ctype in cell_types}
subsets.keys()

gtf_path = "/home/users/nus/e1124313/scratch/eqtl/Homo_sapiens.GRCh38.112.chr.gtf"
db = gffutils.create_db(gtf_path, dbfn=':memory:', force=True, keep_order=True, disable_infer_genes=True, disable_infer_transcripts=True)

def process_input(adata_subset, min_genes=200, min_cells=3, min_counts=5000, path='/home/users/nus/e1124313/scratch/eqtl/0720_input'):
    # if adata_subset.obs['cell_type'].startswith('B'):
    #     min_genes = 300
    #     min_cells = 5
    #     min_counts = 10000

    adata_subset.raw = adata_subset
    celltype = adata_subset.obs['CT_2'].values[0]
    output_path = f"{path}/{celltype}"
    os.makedirs(output_path, exist_ok=True)
    # 基础的质控过滤
    sc.pp.filter_cells(adata_subset, min_genes=min_genes)
    sc.pp.filter_genes(adata_subset, min_cells=min_cells)
    
    ensembl = EnsemblRelease(111) # 初始化Ensembl数据库
    genes = ensembl.genes() # 使用合适的版本号
    mt_genes = [gene for gene in genes if gene.contig == 'MT'] # 获取所有线粒体基因
    mt_genes_in_adata = [gene.gene_id for gene in mt_genes if gene.gene_id in adata_subset.var_names] # 筛选在AnnData中存在的线粒体基因
    adata_subset.obs['mt_expr_sum'] = adata_subset[:, mt_genes_in_adata].X.sum(axis=1) # 计算线粒体基因的总表达
    adata_subset.obs['total_expr_sum'] = adata_subset.X.sum(axis=1) # 计算总表达
    adata_subset.obs['pct_counts_mt'] = (adata_subset.obs['mt_expr_sum'] / adata_subset.obs['total_expr_sum']) * 100 # 计算线粒体基因表达占比
    print(f"{celltype} QC finishied!")

    #过滤
    sample_counts = adata_subset.obs['Sample_ID'].value_counts()
    valid_samples = sample_counts[sample_counts >= 100].index
    adata_subset = adata_subset[adata_subset.obs['Sample_ID'].isin(valid_samples)]
    if adata_subset.obs['Sample_ID'].nunique() > 6:
        # 归一化和计算变异基因
        adata_subset.layers["counts"] = adata_subset.X.copy()

        # use scanpy
        # sc.pp.normalize_total(adata_subset, target_sum=1e4)
        # sc.pp.log1p(adata_subset)
        # sc.pp.regress_out(adata_subset, ['pct_counts_mt', 'Batch_ID'])
        # sc.pp.highly_variable_genes(adata_subset, min_mean=0.0125, max_mean=3, min_disp=0.5)

        # use SCTransform
        ro.r('library("tidyverse")')
        ro.r('library("here")')
        ro.r('library("Seurat")')
        ro.r('library(SeuratDisk)')

        # 转换所有分类列为字符串类型
        for column in adata_subset.obs.select_dtypes(include=['category']).columns:
            adata_subset.obs[column] = adata_subset.obs[column].astype(str)
        
        ro.r.assign("adata_r", adata_subset)
        ro.r('seurat_object <- as.Seurat(adata_r, data="X")')
        ro.r('seurat_object <- SCTransform(object = seurat_object, assay="originalexp", vars.to.regress = c("pct_counts_mt", "Batch_ID"), conserve.memory = TRUE, return.only.var.genes = FALSE)')
        # 从Seurat对象中获取调整后的数据
        ro.r('data_matrix <- GetAssayData(seurat_object, slot = "data", assay = "SCT")')
        # 转换数据回Python
        data_matrix = ro.r('data_matrix')
        data_matrix = data_matrix.T
        data_matrix = sp.csr_matrix(data_matrix, dtype=np.float32)

        # 将矩阵转换为CSR格式，如果需要的话
        seurat_genes = np.array(ro.r('rownames(seurat_object[["SCT"]]@data)'))
        adata_subset = adata_subset[: ,adata_subset.var_names.isin(seurat_genes)] #这里去掉了一些基因
        adata_subset.layers['SCT'] = data_matrix

        #保存筛选过的人群
        sample_ids = adata_subset.obs['Sample_ID']
        keep_samples = pd.concat([sample_ids, sample_ids], axis=1)
        keep_samples.to_csv(f"{output_path}/{celltype}_keep_samples.txt", sep='\t', index=False, header=False)
        print(f'{celltype} keep_samples saved!')

        ################################################################################# 临时对齐使用
        # fam_df = pd.read_csv(f'/home/users/nus/e1124313/scratch/eqtl/raw_plink/chr_all_ok.fam', sep=' ', header=None)
        # num_samples = adata_subset.obs['Sample_ID'].nunique()
        # subset_fam_df = fam_df.iloc[:num_samples]
        # old_fids = subset_fam_df[0]
        # old_iids = subset_fam_df[1]
        # new_iids = adata_subset.obs['Sample_ID'].unique()  # 这将作为新的 Sample ID
        # # 创建 update_ids.txt 文件内容
        # update_ids = pd.DataFrame({
        #     "old_fid": old_fids,
        #     "old_iid": old_iids,
        #     "new_fid": new_iids,  
        #     "new_iid": new_iids  # 新的 Individual ID
        # })
        # update_ids_file = f'{output_path}/{celltype}_updated_id.txt' 
        # update_ids.to_csv(update_ids_file, sep="\t", index=False, header=False)
        # print(f'{celltype} update_ids saved!')
        ##########################################################################################

        # 保存var_names,并且和gtf文件进行比对,保存对应的染色体位置
        var_names = adata_subset.var_names
        gene_list = pd.DataFrame(var_names, columns=['gene_id'])
        # gene_list.to_csv(f"{output_path}/{celltype}_var_names.txt", sep='\t', index=False, header=False)
        locations = []
        for gene_id in gene_list['gene_id']:
            try:
                gene = db[gene_id]
                locations.append([gene_id, gene.chrom, gene.start, gene.end])
            except gffutils.exceptions.FeatureNotFoundError:
                locations.append([gene_id, 'NA', 'NA', 'NA'])

        location_df = pd.DataFrame(locations, columns=['gene_id', 'chrom', 'start', 'end'])
        
        location_df = location_df[~location_df['chrom'].isin(['X', 'Y'])]  # Remove rows with chrom X or Y
        location_df.to_csv(f"{output_path}/{celltype}_gene_locations.txt", sep='\t', index=False)
        # 保存需要测试的基因名，去除了XY染色体上的基因
        gene_list = location_df['gene_id']
        gene_list.to_csv(f"{output_path}/{celltype}_var_names.txt", sep='\t', index=False, header=False)
        print(f'{celltype} var_names saved!')
        print(f'{celltype} gene_locations saved!')
            
        # 生成plink脚本并执行
        plink_script = f"""#!/bin/bash
    mkdir -p {output_path}/chr/
    #cp /home/users/nus/e1124313/scratch/eqtl/raw_plink/split/chr* {output_path}/chr/
    plink --bfile /home/users/nus/e1124313/scratch/eqtl/raw_plink/Asian_sle --keep {output_path}/{celltype}_keep_samples.txt --make-bed --out {output_path}/chr/ATGC_{celltype}_Asian_sle
    plink --bfile {output_path}/chr/ATGC_{celltype}_Asian_sle --pca 20 --out {output_path}/chr/ATGC_{celltype}_Asian_sle
    # 遍历所有22个染色体
    for i in {{1..22}}; do
        # 构建plink命令
        #plink --bfile {output_path}/chr/chr$i --update-ids {output_path}/{celltype}_updated_id.txt --make-bed --out {output_path}/chr/{celltype}_id_updated_chr$i #临时使用
        plink --bfile /home/users/nus/e1124313/scratch/eqtl/raw_plink/update_id_split/Asian_sle_chr$i --keep {output_path}/{celltype}_keep_samples.txt --make-bed --out {output_path}/chr/ATGC_{celltype}_Asian_sle_chr$i 
    done
    """
        with open(f"{output_path}/{celltype}_plink.sh", 'w') as f:
            f.write(plink_script)
        subprocess.run(['chmod', '+x', output_path + f'/{celltype}_plink.sh'])
        subprocess.run([output_path + f'/{celltype}_plink.sh'], capture_output=True, text=True)
        print(f'{celltype} plink finished!')

        # 保存PCA数据
        eigenvec_path = f"{output_path}/chr/ATGC_{celltype}_Asian_sle.eigenvec"
        pca_df = pd.read_csv(eigenvec_path, sep=' ', header=None)
        
        pca_df.columns = ['FID'] + ['Sample_ID'] + [f'geno_PC{i+1}' for i in range(pca_df.shape[1] - 2)]
        pca_df_reduced = pca_df[['Sample_ID', 'geno_PC1', 'geno_PC2', 'geno_PC3', 'geno_PC4', 'geno_PC5', 'geno_PC6']]
        adata_subset.obs = adata_subset.obs.merge(pca_df_reduced, on='Sample_ID', how='left')

        # 保存adata
        adata_subset.write(f'{output_path}/{celltype}.h5ad') ## 第一次保存
        print(f'{celltype} adata_subset saved!')

        # 保存输入矩阵
        adata_subset = sc.read(f'{output_path}/{celltype}.h5ad') # 这一步很重要，因为前面的保存会将obs中的string转换为category，这样后面生成矩阵才不会为空
        obs_list = adata_subset.obs.columns.tolist()
        covariates = adata_subset.obs[['Batch_ID', 'Age', 'Sex', 'Ethnicity_1', 'Ethnicity_2', 'Status', 'Disease_status', 'Assay_method', 'CT_1', 'CT_2', 'geno_PC1', 'geno_PC2', 'geno_PC3', 'geno_PC4', 'geno_PC5', 'geno_PC6']]
        gene_data = pd.DataFrame(adata_subset.layers['SCT'].toarray(), columns=adata_subset.var_names) # 如果做了normalization，这里的数据是array，不需要toarray()
        sample_ids = adata_subset.obs['Sample_ID'] # 这里得重新定义一下，不然会和covariates的index不一致，不知道前面哪里影响了
        gene_data.index = sample_ids.index
        combined_df = pd.concat([sample_ids, covariates, gene_data], axis=1)
        combined_df.to_csv(f"{output_path}/{celltype}.txt", sep='\t', index=False)
        print(f'{celltype} input matrix saved!')
    
    else:
        print(f"{celltype} has less than 6 samples, skipped!")

    # 返回处理后的数据
    return adata_subset


def main():
    # 使用进程池来并行处理每个细胞类型的子集
    with ProcessPoolExecutor(16) as executor:
        # 创建一个future到adata子集的字典
        futures = {executor.submit(process_input, subset): ctype for ctype, subset in subsets.items()}
        # 为每个执行完成的future收集结果
    results = {futures[future]: future.result() for future in futures}

    return results

inputs = main()

# def main():
#     results = {}
#     for ctype, subset in subsets.items():
#         result = process_input(subset)
#         results[ctype] = result

#     return results

# inputs = main()
