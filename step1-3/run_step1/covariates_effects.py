'''
:@Author: liuxia
:@Date: 7/15/2024, 7:15:54 PM
:@LastEditors: liuxia
:@LastEditTime: 7/15/2024, 7:15:54 PM
:Description: 
'''
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import seaborn as sns
import scipy.stats as stats

data_dir = '/home/users/nus/e1124313/scratch/eqtl/input'
celltype_dir = os.path.join(data_dir, 'cell_types.csv')
celltypes = pd.read_csv(celltype_dir, header=None)
celltypes = celltypes[0].tolist()
for celltype in celltypes:
    adata = sc.read(os.path.join(data_dir, f'{celltype}', f'{celltype}.h5ad'))
    adata.obs['Age'] = adata.obs['Age'].astype(int)
    df = pd.DataFrame({
        'age': adata.obs.groupby('Sample_ID')['Age'].first(),
        'sex': adata.obs.groupby('Sample_ID')['Sex'].first(), 
    })
    print(df['age'])
    bins = [0, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 999]
    labels = ['<20', '20-24', '25-29', '30-34', '35-39', '40-44', '45-49', '50-54', '55-59', '60-64', '65-69', '70-74', '75-79', '80-84', '85-89', '90+']
    df['age_group'] = pd.cut(df['age'], bins=bins, labels=labels, right=False)

    grouped = df.groupby(['age_group', 'sex']).size().unstack(fill_value=0)
    fig, ax = plt.subplots(figsize=(10, 8))
    width = 0.35  # 条形的宽度
    x = np.arange(len(grouped.index))  # x轴的标签位置

    for i, sex in enumerate(grouped.columns):
        ax.bar(x - width/2 + i*width, grouped[sex], width, label=sex)

    # 设置x轴标签和图例
    ax.set_xlabel('Age Group')
    ax.set_ylabel('Number of Participants')
    ax.set_title('Age and Sex Distribution')
    ax.set_xticks(x)
    ax.set_xticklabels(grouped.index)
    ax.legend(title='Sex')

    # 显示图形
    plt.savefig(os.path.join(data_dir, 'age_sex_distribution', f'{celltype}.png'), dpi=300, bbox_inches='tight')
    plt.close()


    df_2 = pd.DataFrame({
        'Age': adata.obs.groupby('Sample_ID')['Age'].first(),
        'cell_counts': adata.obs['Sample_ID'].value_counts(),
        'Sex': adata.obs.groupby('Sample_ID')['Sex'].first(),
    })
    # 创建图形
    plt.figure(figsize=(6, 4))
    x=df_2['Age']
    y=df_2['cell_counts']
    plt.scatter(x, y, color='purple', alpha=0.5)  # 散点图
    
    # 计算拟合线
    if len(x) > 1:
        coefficients = np.polyfit(x, y, 1)  # 1表示线性拟合
        poly = np.poly1d(coefficients)
        plt.plot(np.unique(x), poly(np.unique(x)), color='black')  # 拟合线
    
    # 图形定制
    plt.title(f'{celltype} Cell Counts vs. Age')
    plt.xlabel('Age')
    plt.ylabel('Cell counts')
    plt.grid(True)
    
    # 保存图像
    plt.savefig(os.path.join(data_dir, 'age_sex_distribution', f'{celltype}_age_counts_scatter.png'), dpi=300, bbox_inches='tight')
    plt.close()  # 关闭图形以释放内存

    # 绘制箱线图
    plt.figure(figsize=(6, 4))
    sns.boxplot(data=df_2, x='Sex', y='cell_counts', palette='Set2')
    # 添加标题和轴标签
    plt.title(f'{celltype} Counts by Sex')
    plt.xlabel('Sex')
    plt.ylabel('Cell Counts')

    # 计算男性和女性之间的t-test p值
    male_counts = df_2[df_2['Sex'] == 1]['cell_counts']
    female_counts = df_2[df_2['Sex'] == 0]['cell_counts']
    t_stat, p_val = stats.ttest_ind(male_counts, female_counts)

    # 添加p值
    if not df_2['cell_counts'].empty:
        plt.text(0.5, max(df_2['cell_counts']) * 0.95, f'p={p_val:.2e}', horizontalalignment='center', color='black')

    # 保存和显示图形
    plt.savefig(os.path.join(data_dir, 'age_sex_distribution', f'{celltype}_sex_counts_boxplot.png'), dpi=300, bbox_inches='tight')
    plt.close()





