from multiprocessing import Pool
import subprocess
import pandas as pd
import os

cell_type = os.getenv('CELL_TYPE') 
tasks = []
df = pd.read_csv(f'/home/users/nus/e1124313/scratch/eqtl/0720_input/{cell_type}/{cell_type}_gene_locations.txt', sep='\t')

# num_rows = len(df)
# chunk_size = num_rows // 3

# df_chunks = [df[i:i+chunk_size] for i in range(0, num_rows, chunk_size)]

for index, row in df.iterrows():
    gene_id = row['gene_id']
    tasks.append((cell_type, gene_id))

def run(args):
    cell_type, gene_id = args
    subprocess.call(['/home/users/nus/e1124313/scratch/eqtl/step1-3/run_step1/run_step1.sh', cell_type, gene_id])

if __name__ == "__main__":
     with Pool(128) as pool:
        pool.map(run, tasks)
