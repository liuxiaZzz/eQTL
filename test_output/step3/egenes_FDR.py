'''
:@Author: liuxia
:@Date: 9/23/2024, 1:15:41 PM
:@LastEditors: liuxia
:@LastEditTime: 9/23/2024, 1:15:41 PM
:Description: 
'''
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests
import os

import concurrent.futures
from tqdm import tqdm

cell_type_list = ['B_mem']#,'B_naive','CytoT_GZMH+','CytoT_GZMK+','myeloid','NK_dim','Progen','T4_em','T4_reg','T8_naive']

def compute_gene_p_value_FDR(cell_type):
    # Define the directory containing the files
    # path = "/home/users/nus/e1124313/scratch/eqtl/test_output"
    path = "/home/users/nus/e1124313/scratch/eqtl/step1-3/cis"
    directory = f"{path}/step2/{cell_type}/output"

    gene_p_values = []
    gene_list = []
    top_snp_pval_list = []
    top_snp_marker_list = []

    pval_files = [filename for filename in os.listdir(directory) if filename.endswith("_pval")] #进度条

    for filename in tqdm(pval_files, desc=f"Processing {cell_type}", unit="file"):
        # Read the file
        df = pd.read_csv(os.path.join(directory, filename), sep='\t')
        # Extract the p-values
        gene_p_values.extend(df['ACAT_p'].values)
        gene_list.extend(df['gene'].values)
        top_snp_pval_list.extend(df['top_pval'].values)
        top_snp_marker_list.extend(df['top_MarkerID'].values)

    # Compute the FDR
    egene_bool, gene_FDR, _, _ = multipletests(gene_p_values, method='fdr_bh', alpha=0.05)
    egene_list = [gene_list[i] for i in range(len(gene_list)) if egene_bool[i]]
    top_snp_pval_selected = [top_snp_pval_list[i] for i in range(len(gene_list)) if egene_bool[i]]
    top_snp_marker_selected = [top_snp_marker_list[i] for i in range(len(gene_list)) if egene_bool[i]]


    output_file = f"{path}/step3/{cell_type}_egene_list.csv"
    pd.DataFrame({'gene': egene_list,
                  'top_pval': top_snp_pval_selected,
                  'top_MarkerID': top_snp_marker_selected}).to_csv(output_file, index=False)
    print(f"Saved {len(egene_list)} egenes for {cell_type} to {output_file}")
    


# Use concurrent.futures to run the process_cell_type function in parallel
with concurrent.futures.ProcessPoolExecutor() as executor:
    executor.map(compute_gene_p_value_FDR, cell_type_list)