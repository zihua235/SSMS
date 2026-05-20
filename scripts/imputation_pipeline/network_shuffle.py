import pandas as pd
import os
import random

input_file = "adj_matrix/adj_Italy2_spieceasi.txt"
output_dir = "adj_matrix/adj_Italy2_spieceasi_label_shuffled"
num_shuffles = 100
base_seed = 68

os.makedirs(output_dir, exist_ok=True)

# 读取邻接矩阵
adj = pd.read_csv(input_file, sep="\t", index_col=0)

# 检查：行列标签是否一致
if not set(adj.index) == set(adj.columns):
    raise ValueError("unmatched index")

for i in range(num_shuffles):
    random.seed(base_seed + i)

    shuffled_labels = list(adj.index)
    random.shuffle(shuffled_labels)

    # 打乱标签，但不改变结构
    shuffled_matrix = adj.copy()
    shuffled_matrix.index = shuffled_labels
    shuffled_matrix.columns = shuffled_labels

    # 保存打乱标签的矩阵
    output_path = os.path.join(output_dir, f"{i}.txt")
    shuffled_matrix.to_csv(output_path, sep="\t")

    print(f"Saved label-shuffled matrix {i} to {output_path}")
