'''
:@Author: liuxia
:@Date: 8/6/2024, 3:02:10 PM
:@LastEditors: liuxia
:@LastEditTime: 8/6/2024, 3:02:10 PM
:Description: 
'''
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests
import os

# cell_type = 'T4_em'
import concurrent.futures

cell_type_list = ['B_mem','B_naive','CytoT_GZMH+','CytoT_GZMK+','myeloid','NK_dim','Progen','T4_em','T4_reg','T8_naive']

def process_cell_type(cell_type):
    # Define the directory containing the files
    directory = f"/home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step2/{cell_type}/output"

    # Create an empty DataFrame to store results
    results_df = pd.DataFrame()

    # Iterate through all the files in the directory
    for filename in os.listdir(directory):
        # Check if the file does not have ".index" in its name
        if ".index" not in filename:
            print(f"Processing {filename}")
            # Read the data
            file_path = os.path.join(directory, filename)
            df = pd.read_csv(file_path, sep="\t")
            p_values = df["p.value"].values
            if len(p_values) != 0:
                fdr_results = multipletests(p_values, method="fdr_bh", alpha=0.05)
                filtered_df = df[fdr_results[0]]
                filtered_df = filtered_df[['CHR', 'POS', 'MarkerID', 'Allele1', 'Allele2']]
                filtered_df['Adjusted_pvalue'] = fdr_results[1][fdr_results[0]]
            
            # If there are any hits, extract the gene name from the filename
            if not filtered_df.empty:
                gene_name = filename.split('_')[-1]  # Adjust based on how the gene name is included in the filename
                filtered_df['Gene'] = gene_name  # Add gene name to the DataFrame
                
                # Append the filtered data to the results DataFrame
                results_df = pd.concat([results_df, filtered_df], ignore_index=True)

    # Save the results to a new file if needed
    results_df.to_csv(f"/home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step2/{cell_type}/FDR_results.csv", sep='\t', index=False)

    # Show the first few results
    results_df.head()

# Use concurrent.futures to run the process_cell_type function in parallel
with concurrent.futures.ProcessPoolExecutor() as executor:
    executor.map(process_cell_type, cell_type_list)