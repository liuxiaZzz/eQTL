import pandas as pd
from pandas_plink import read_plink
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# 本脚本用于绘制基因型和表达量之间的关系

celltype = 'myeloid'
df = pd.read_csv(f"/home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step2/{celltype}/all_results.csv", sep='\t')
filtered_df = df[(df['p.value'] > 0) & (df['p.value'] < 1e-4)]

adata = sc.read(f"/home/users/nus/e1124313/scratch/eqtl/input/{celltype}/{celltype}.h5ad")

# 'bim' 包含SNP信息
# 'fam' 包含样本信息
# 'bed' 是一个NumPy数组，包含基因型数据

for index, row in filtered_df.iterrows():
    plt.figure()
    
    gene = row['Gene']
    snp_id = row['MarkerID']
    chr = row['CHR']

    (bim, fam, bed) = read_plink(f'/home/users/nus/e1124313/scratch/eqtl/input/{celltype}/chr/ATGC_{celltype}_Asian_sle_chr{chr}')
    
    # 获取SNP索引
    snp_index = bim[bim.snp == snp_id].index.values[0]
    genotypes = bed[snp_index, :].compute().T

    # 创建基因型数据框
    genotype_df = pd.DataFrame({
        'Sample_ID': fam.iid,
        'Genotype': genotypes
    })

    # 创建表达数据框
    expression_df = pd.DataFrame({
        'Sample_ID': adata.obs['Sample_ID'],
        'Expression': adata[:, gene].X.toarray().ravel()
    })

    # 合并数据
    merged_df = pd.merge(genotype_df, expression_df, on='Sample_ID')


    # 计算每个基因型组的平均表达
    average_expression_per_sample = merged_df.groupby('Sample_ID')['Expression'].mean().reset_index()

    merged_df = merged_df.drop('Expression', axis=1)
    merged_df = merged_df.drop_duplicates()
    merged_df = pd.merge(merged_df, average_expression_per_sample, on='Sample_ID')

    # merged_df.to_csv(f'/home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step2/{celltype}/snp_gex/{snp_id}.csv', index=False)

    # 假设数据已经加载到DataFrame中，并且命名为df
    # df['Genotype'] 包含基因型数据 "AA", "AB", "BB"
    # df['Expression'] 包含基因表达水平数据

    # # 使用seaborn的lmplot绘制散点图和回归线
    # ax = sns.scatterplot(x='Genotype', y='Expression', data=merged_df, alpha=0.5, color='blue', label='Individual Data')

    # # 计算每个基因型的平均表达量，并作为点图绘制
    # sns.pointplot(x='Genotype', y='Expression', data=merged_df, join=True, color='red', markers='o', ci=None, ax=ax, label='Mean Expression')
    # sns.lmplot(x='Genotype', y='Expression', data=merged_df, height=5, aspect=1,
    #            fit_reg=True, markers='o', scatter_kws={'s': 50, 'alpha': 0.6})
    # 画小提琴图
    sns.violinplot(x='Genotype', y='Expression', data=merged_df, inner=None, edgecolor='lightgreen', fill=False)

    # 画盒形图
    sns.boxplot(x='Genotype', y='Expression', data=merged_df, width=0.2, showfliers=False, boxprops={'facecolor': 'none'}, whiskerprops={'linewidth':0})

    # 添加线性回归线
    model = np.polyfit([0, 1, 2], merged_df.groupby('Genotype')['Expression'].mean(), 1)
    plt.plot([0, 1, 2], np.polyval(model, [0, 1, 2]), color="red", linewidth=2)

    plt.title('Gene Expression Levels by SNP Genotype')
    plt.xlabel('SNP X Genotype')
    plt.ylabel(f'{gene} Expression levels')
    # plt.xticks(rotation=45)  # 基因型标签可能需要旋转以便更好地显示
    plt.savefig(f'/home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step2/{celltype}/snp_gex/{snp_id}.png', dpi=300, bbox_inches='tight')
    
    plt.close()