import numpy as np

def predict_features(adj_matrix, feature_matrix):
    n_features, n_samples = feature_matrix.shape
    predicted_features = np.zeros((n_features, n_samples))

    for feature_index in range(n_features):
        neighbors = np.nonzero(adj_matrix[feature_index])[0]  # 找到该特征的邻居索引
        for sample_index in range(n_samples):
            neighbor_features = feature_matrix[neighbors, sample_index]  # 获取邻居的特征值
            if len(neighbor_features) > 0:
                predicted_features[feature_index, sample_index] = np.mean(neighbor_features)  # 计算邻居特征的平均值

    return predicted_features

# 示例数据
adj_matrix = np.array([[0, 1, 1],
                       [1, 0, 1],
                       [1, 1, 0]])

feature_matrix = np.array([[1, 2, 3],
                           [4, 5, 6],
                           [7, 8, 9]])

predicted_features = predict_features(adj_matrix, feature_matrix)
print("Predicted Features:")
print(predicted_features)