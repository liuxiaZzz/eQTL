import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# 本脚本用于绘制曼哈顿图

# df = pd.read_csv("/home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step2/myeloid/screen_results.csv", sep='\t')
# df['-log10(p-value)'] = -np.log10(df['p.value'])
# fig, ax = plt.subplots(figsize=(10, 6))
# df['CHR'] = df['CHR'].astype(str)  # Ensure chromosome column is string for categorical plotting
# df = df.sort_values(['CHR', 'POS'])

# df['ind'] = range(len(df))
# df_grouped = df.groupby('CHR')

# for i, (name, group) in enumerate(df_grouped):
#     ax.scatter(group['ind'], group['-log10(p-value)'], s=10, label=name if i % 2 == 0 else "", color='blue' if i % 2 == 0 else 'red')

# ax.set_xlabel('Index')
# ax.set_ylabel('-log10(p-value)')
# ax.legend(title='Chromosome')
# plt.savefig("/home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step2/screen_results.pdf", bbox_inches='tight')

############################################################

cell_type = 'T4_em'

df = pd.read_csv(f"/home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step2/{cell_type}/all_results.csv", sep='\t')

# -log_10(pvalue)
df['minuslog10pvalue'] = -np.log10(df['p.value'])
df.chromosome = df['CHR'].astype('category')
df.chromosome = df.chromosome.cat.set_categories(['ch-%i' % i for i in range(22)], ordered=True)
df = df.sort_values('CHR')

# How to plot gene vs. -log10(pvalue) and colour it by chromosome?
df['ind'] = range(len(df))
df_grouped = df.groupby(('CHR'))

# manhattan plot
fig = plt.figure(figsize=(10,4),dpi=100) 
ax = fig.add_subplot(111)

colors = ["#30A9DE","#EFDC05","#E53A40","#090707"]
x_labels = []
x_labels_pos = []
for num, (name, group) in enumerate(df_grouped):
    group.plot(kind='scatter', x='ind', y='minuslog10pvalue',color=colors[num % len(colors)], ax=ax)
    x_labels.append(name)
    x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0])/2))
# add grid
ax.grid(axis="y",linestyle="--",linewidth=.5,color="gray")
ax.tick_params(direction='in',labelsize=13)
ax.set_xticks(x_labels_pos)
ax.set_xticklabels(x_labels)
ax.set_title(f'Manhattan Plot for {cell_type}')
ax.set_xlim([0, len(df)])
# ax.set_ylim([0, 10])
# x axis label
ax.set_xlabel('Chromosome',size=14)
plt.savefig(f"/home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step2/{cell_type}/{cell_type}_all_results.pdf",dpi=300,bbox_inches='tight',facecolor='white')