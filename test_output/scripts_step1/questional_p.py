import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests
import os

# cell_type = 'T4_em'
import concurrent.futures
import logging

cell_type_list = ['B_mem']#,'B_naive','CytoT_GZMH+','CytoT_GZMK+','myeloid','NK_dim','Progen','T4_em','T4_reg','T8_naive']

import scanpy as sc
for cell_type in cell_type_list:
    adata = sc.read(f"/home/users/nus/e1124313/scratch/eqtl/0920_input/{cell_type}/{cell_type}.h5ad")
    # Read the test_0_results.csv file
    test_0_results = pd.read_csv(f"/home/users/nus/e1124313/scratch/eqtl/test_output/step2/{cell_type}/test_0_results.csv", sep='\t')

    # Get the gene names with non-zero num_p_0 values
    non_zero_genes = test_0_results[test_0_results['num_p_0'] != 0]['filename'].str.split('_').str[-1]

    sum_expression = {}
    for gene in non_zero_genes:
        phenotype = adata[:, gene].layers['SCT']
        sum_expression[gene] = {}
        for sample in adata.obs['Sample_ID'].unique():
            total_expression = phenotype[adata.obs['Sample_ID'] == sample].sum()
            num_cells = (adata.obs['Sample_ID'] == sample).sum()
            # sum_expression[gene][sample] = phenotype[adata.obs['Sample_ID'] == sample].sum()
            sum_expression[gene][sample] = num_cells

        # Convert the sum_expression dictionary to a DataFrame
    
    sum_expression_df = pd.DataFrame(sum_expression)
    sum_expression_df.columns
    sum_expression_df.to_csv(f"/home/users/nus/e1124313/scratch/eqtl/test_output/step2/{cell_type}/questional_sum_SCT_2.csv", sep='\t')