'''
:@Author: liuxia
:@Date: 7/18/2024, 4:08:50 PM
:@LastEditors: liuxia
:@LastEditTime: 7/18/2024, 4:08:50 PM
:Description: 
'''
import peer
import pandas as pd


expr = pd.read_csv('/home/users/nus/e1124313/scratch/eqtl/0716_input/myeloid/myeloid.txt', index_col=False, sep='\t')
expr2 = expr.drop(columns=['Sample_ID', 'Batch_ID', 'Age', 'Sex', 'Ethnicity_1', 'Ethnicity_2', 'Status', 'Disease_status', 'Assay_method', 'CT_1', 'CT_2', 'geno_PC1', 'geno_PC2', 'geno_PC3', 'geno_PC4', 'geno_PC5', 'geno_PC6'])
covs = expr[['Age', 'Sex', 'geno_PC1', 'geno_PC2', 'geno_PC3', 'geno_PC4', 'geno_PC5', 'geno_PC6']]

model = peer.PEER()
model.setPhenoMean(expr2)
model.setCovariates(covs)

model.setNk(10)
model.getNK()
model.update()







