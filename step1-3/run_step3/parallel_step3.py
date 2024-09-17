from multiprocessing import Pool
import subprocess
import pandas as pd
import os

cell_type = os.getenv('CELL_TYPE') 
tasks = []
df = pd.read_csv(f'/home/users/nus/e1124313/scratch/eqtl/input/{cell_type}/{cell_type}_gene_locations.txt', sep='\t')

for index, row in df.iterrows():
    gene_id = row['gene_id']
    chrom = row['chrom']
    start = row['start']
    end = row['end']
    tasks.append((cell_type, gene_id, chrom, start, end))

def run(args):
    cell_type, gene_id, chrom, start, end = args
    subprocess.call(['/home/users/nus/e1124313/scratch/eqtl/step1-3/run_step3/run_step3.sh', cell_type, gene_id, str(chrom), str(start), str(end)])

if __name__ == "__main__":
     with Pool(128) as pool:
        pool.map(run, tasks)
