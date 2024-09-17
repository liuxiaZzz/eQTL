'''
:@Author: liuxia
:@Date: 8/6/2024, 5:35:23 PM
:@LastEditors: liuxia
:@LastEditTime: 8/6/2024, 5:35:23 PM
:Description: 
'''
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests
import os

# cell_type = 'T4_em'
import concurrent.futures
import logging

cell_type_list = ['B_mem','B_naive','CytoT_GZMH+','CytoT_GZMK+','myeloid','NK_dim','Progen','T4_em','T4_reg','T8_naive']

def process_cell_type(cell_type):
    # Define the directory containing the files
    directory = f"/home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step2/{cell_type}/output"

    # Create an empty DataFrame to store results
    results = {}

    # Iterate through all the files in the directory
    for filename in os.listdir(directory):
        # Check if the file does not have ".index" in its name
        if ".index" not in filename:
            print(f"Processing {filename}")
            # Read the data
            file_path = os.path.join(directory, filename)
            try:
                df = pd.read_csv(file_path, sep="\t")
                num_p_0 = (df['p.value'] == 0).sum()
                results[filename] = num_p_0
            except Exception as e:
                logging.error(f"Error processing {filename}: {e}")
                continue

    # Save the results to a new file if needed
    results_df = pd.DataFrame(results.items(), columns=['filename', 'num_p_0'])
    results_df.to_csv(f"/home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step2/{cell_type}/test_0_results.csv", sep='\t', index=False)

# Use concurrent.futures to run the process_cell_type function in parallel
with concurrent.futures.ProcessPoolExecutor() as executor:
    executor.map(process_cell_type, cell_type_list)


import scanpy as sc
adata = sc.read("/home/users/nus/e1124313/scratch/eqtl/input/B_mem/B_mem.h5ad")

expression = adata[:, 'ENSG00000112715'].X
for sample in adata.obs['Sample_ID'].unique():
    sum_exrp = expression[adata.obs['Sample_ID'] == sample].sum()
    print(sample, sum_exrp)
