from multiprocessing import Pool
import subprocess
import pandas as pd
import os

cell_type = os.getenv('CELL_TYPE') 
tasks = []
df = pd.read_csv(f'/home/users/nus/e1124313/scratch/eqtl/0720_input/{cell_type}/{cell_type}_gene_locations.txt', sep='\t')

num_rows = len(df)
chunk_size = num_rows // 3

df_chunks = [df[i:i+chunk_size] for i in range(0, num_rows, chunk_size)]

if len(df_chunks) > 3:
    df_chunks[-2] = pd.concat([df_chunks[-2], df_chunks[-1]])
    df_chunks.pop()

for index, row in df_chunks[1].iterrows():
    gene_id = row['gene_id']
    tasks.append((cell_type, gene_id))

def run(args):
    cell_type, gene_id = args
    subprocess.call(['/home/users/nus/e1124313/scratch/eqtl/test_output/scripts_step1/test_step1.sh', cell_type, gene_id])

if __name__ == "__main__":
     with Pool(128) as pool:
        pool.map(run, tasks)
