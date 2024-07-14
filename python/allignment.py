import pandas as pd
import numpy as np
import os

# load data
current_dir = os.getcwd()

org_EC = pd.read_csv(os.path.join(current_dir, 'merged_level4ec_relab_unstratified.tsv'),
                     sep='\t', index_col=0)
EC_aj = pd.read_csv(os.path.join(current_dir, 'SparCC', 'BackhedF_2015_0.03_0.05_network.adj.txt'),
                    sep='\t', index_col=0)

# number of overlap EC number
set1 = set(org_EC.index)
set2 = set(EC_aj.index)
intersection = set1.intersection(set2) # EC number that exit in both org_EC and EC_aj
set3 = set1 - intersection # EC number that only exit in org_EC
set4 = set2 - intersection # EC number that only exit in EC_aj

print("length of org_EC", len(set1))
print("length of EC_aj", len(set2))
print("Number of overlapping elements:", len(intersection))
print("length of unique EC number in org_EC", len(set3))
print("length of unique EC number in EC_aj", len(set4))

# select overlapped EC
org_EC_overlap = org_EC.loc[np.array(intersection)]

# modify both org_EC and EC_aj
EC_aj_overlap = EC_aj.loc[np.array(intersection)]
EC_aj_overlap = EC_aj_overlap[np.array(intersection)]
print("dim of EC_aj_overlap", EC_aj_overlap.shape)

EC_overlap = org_EC.loc[np.array(intersection)]

print("length of EC_overlap", len(EC_overlap.index))
# predict EC
def predict_features(adj_matrix, feature_matrix):
    n_features, n_samples = feature_matrix.shape
    predicted_features = np.zeros((n_features, n_samples))

    for feature_index in range(n_features):
        neighbors = np.nonzero(adj_matrix[feature_index])[0]  # 找到该特征的邻居索引
        for sample_index in range(n_samples):
            neighbor_features = feature_matrix[neighbors, sample_index]# 获取邻居的特征值
            if feature_matrix[feature_index, sample_index] > 0: # 不改变原有的值
                predicted_features[feature_index, sample_index] = feature_matrix[feature_index, sample_index]
            elif len(neighbor_features) > 0:
                predicted_features[feature_index, sample_index] = np.mean(neighbor_features)  # 计算邻居特征的平均值

    return predicted_features


# predition
EC_aj_overlap_array = EC_aj_overlap.values
EC_overlap_array = EC_overlap.values

predicted_features = predict_features(EC_aj_overlap_array, EC_overlap_array)
print("number of predicted Features", len(predicted_features))

# merge non-overlapped EC number
org_EC_non_overlap = org_EC.loc[np.array(set3)].values
predicted_features_final = np.concatenate((predicted_features, org_EC_non_overlap), axis=0)
col_names = org_EC.columns
index = EC_overlap.index.tolist() + org_EC.loc[np.array(set3)].index.tolist()
predicted_features_final = pd.DataFrame(predicted_features_final, columns=col_names, index=index)
print("dim of output", predicted_features_final.shape)
# save data
predicted_features_final.to_csv('/Users/zihuahuang/Documents/SSMS/BackhedF_2015_predict_level4ec_SparCC_0.03_0.05.csv',
                                index=True, header=True)
