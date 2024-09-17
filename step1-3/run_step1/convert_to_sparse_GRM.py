'''
:@Author: liuxia
:@Date: 7/15/2024, 12:18:25 PM
:@LastEditors: liuxia
:@LastEditTime: 7/15/2024, 12:18:25 PM
:Description: 失败了，还是先不用这个
'''
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.io import mmwrite

def read_grm_bin(prefix):
    # Construct file names
    bin_file_name = f"{prefix}.grm.bin"
    n_file_name = f"{prefix}.grm.N.bin"
    id_file_name = f"{prefix}.grm.id"

    # Read sample IDs
    ids = pd.read_csv(id_file_name, header=None, delim_whitespace=True, usecols=[1])

    # Number of individuals
    num_individuals = len(ids)
    
    # Read GRM bin file
    grm = np.fromfile(bin_file_name, dtype=np.float32)
    n_values = np.fromfile(n_file_name, dtype=np.float32)

    # Calculate positions in the symmetric matrix
    triu_indices = np.triu_indices(num_individuals)

    # Create sparse matrix
    sparse_grm = coo_matrix((grm, triu_indices), shape=(num_individuals, num_individuals))

    # Symmetrize the matrix
    diagonal_matrix = coo_matrix((sparse_grm.diagonal(), (np.arange(num_individuals), np.arange(num_individuals))), shape=(num_individuals, num_individuals))
    sparse_grm = sparse_grm + sparse_grm.T - diagonal_matrix

    return ids, sparse_grm, n_values

def save_sparse_grm(ids, sparse_grm, output_prefix):
    # Save IDs
    ids.to_csv(f"{output_prefix}_sampleIDs.txt", header=False, index=False, sep='\t')

    # Save sparse matrix
    mmwrite(f"{output_prefix}_sparseGRM.mtx", sparse_grm, symmetry='symmetric')

# Usage
prefix = '/home/users/nus/e1124313/scratch/eqtl/raw_plink/GRM/Asian_sle.grm'
output_prefix = '/home/users/nus/e1124313/scratch/eqtl/raw_plink/GRM/Asian_sle.grm'
ids, sparse_grm, n_values = read_grm_bin(prefix)
save_sparse_grm(ids, sparse_grm, output_prefix)


