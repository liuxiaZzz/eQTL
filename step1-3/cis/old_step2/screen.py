import pandas as pd
import os

# cell_type = 'T4_em'
cell_type_list = ['B_mem','B_naive','CytoT_GZMH+','CytoT_GZMK+','myeloid','NK_dim','Progen','T4_em','T4_reg','T8_naive']
for cell_type in cell_type_list:

    # Define the directory containing the files
    directory = f"/home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step2/{cell_type}/output"

    # Create an empty DataFrame to store results
    results_df = pd.DataFrame()

    # Iterate through all the files in the directory
    for filename in os.listdir(directory):
        # Check if the file does not have ".index" in its name
        if ".index" not in filename:
            # Read the data
            file_path = os.path.join(directory, filename)
            df = pd.read_csv(file_path, delim_whitespace=True)
            
            # Filter the rows where p.value is <= 1e-4
            # filtered_df = df[(df["p.value"] <= 1e-4) & (df["p.value"] > 0)]
            filtered_df = df[df["p.value"] > 0]
            
            # If there are any hits, extract the gene name from the filename
            if not filtered_df.empty:
                gene_name = filename.split('_')[2]  # Adjust based on how the gene name is included in the filename
                filtered_df['Gene'] = gene_name  # Add gene name to the DataFrame
                
                # Append the filtered data to the results DataFrame
                results_df = pd.concat([results_df, filtered_df], ignore_index=True)

    # Save the results to a new file if needed
    results_df.to_csv(f"/home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step2/{cell_type}/all_results.csv", sep='\t', index=False)

    # Show the first few results
    results_df.head()

