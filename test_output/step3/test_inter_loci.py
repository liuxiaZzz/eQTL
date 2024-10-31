import pandas as pd

# 读取第一个文件的第一列（gene 列）
file1 = "/home/users/nus/e1124313/scratch/eqtl/test_output/step3/B_mem_egene_list.csv"  # 替换为实际文件路径
file2 = "/home/users/nus/e1124313/scratch/eqtl/step1-3/cis/step3/B_mem_egene_list.csv"  # 替换为实际文件路径

# 读取文件并提取第一列
genes_file1 = pd.read_csv(file1, sep=',')['gene']
genes_file2 = pd.read_csv(file2, sep=',')['gene']

# 将两个基因列表转换为集合
set_file1 = set(genes_file1)
set_file2 = set(genes_file2)

# 计算交集
intersection_genes = set_file1.intersection(set_file2)

# 打印交集结果
print("交集中的基因：")
for gene in intersection_genes:
    print(gene)

# 如果需要将结果保存到文件，可以使用以下代码
output_file = "/home/users/nus/e1124313/scratch/eqtl/test_output/step3/intersection_genes.txt"
with open(output_file, 'w') as f:
    for gene in intersection_genes:
        f.write(f"{gene}\n")

print(f"交集基因已保存到 {output_file}")
