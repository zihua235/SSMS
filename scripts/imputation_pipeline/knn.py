from sklearn.neighbors import NearestNeighbors 
import numpy as np
import pandas as pd

def _align_column(df1:pd.DataFrame, df2:pd.DataFrame):
    """
    Align the row and column of two dataframes. Let both databases have the same index

    Args:
    - feature: pd.dataframe

    Returns:
    - pd.DataFrame, the overlapped index of two dataframes
    """
    df1_columns = df1.columns
    df2_columns = df2.columns

    intersection = df1_columns.intersection(df2_columns) # 
    
    df1_overlap = df1.loc[:, intersection]
    df2_overlap = df2.loc[:, intersection]
    df2_overlap = df2_overlap.reindex(index=df1_overlap.index, columns=df1_overlap.columns) # Ensure df2 has the same index as df1

    if not df1_overlap.index.equals(df2_overlap.index):
        raise ValueError("Error: two dataframes have different indexes or columns! Terminating program.")
    
    if not df1_overlap.columns.equals(df2_overlap.columns):
        raise ValueError("Error: two dataframes have different indexes or columns! Terminating program.")
    
    return df1_overlap, df2_overlap
    
def find_neighbors(*, samples, n_neighbors=5):
    """
    Find the nearest neighbors of each sample in the given data. Note that ach "row" in the data is considered as a sample.

    Parameters:
    n_neighbors: int
    algorithm: 'auto' as defoult, can also choose 'ball_tree', 'kd_tree', 'brute' 
    
    Returns:
    neighbors_indices: list of arrays. Each array contains the indices of the nearest neighbors of the corresponding sample
    """
    if samples.shape[0] >= samples.shape[1]:
        raise ValueError("The number of samples should be larger than the number of features")

    samples_filled = samples.loc[:,(samples > 0).all(axis=0)] # only consider the features without missing values and 0 values.
    knn = NearestNeighbors(n_neighbors=n_neighbors, algorithm='auto')
    knn.fit(samples_filled) #shape (n_samples, n_features)
    distances, indices = knn.kneighbors(samples_filled) 

    return distances, indices  

def recover_gene_from_neighbors(*, samples, n_neighbors=5):
    """
    recover gene content from the nearest neighbors. Each "row" in the data is considered as a sample.

    Parameters:
    n_neighbors: int
    algorithm: 'auto' as defoult, can also choose 'ball_tree', 'kd_tree', 'brute' 
    
    Returns:
    recovered_samples: pd.dataframe. weighted imputed samples
    """
    distances, indices = find_neighbors(samples=samples, n_neighbors=n_neighbors)
    
    n_features, n_samples = samples.shape
    recovered_samples = pd.DataFrame(np.nan, index=range(n_features), columns=range(n_samples))  # Initialize a zero matrix for predictions
    
    for i, (dist, ind) in enumerate(zip(distances, indices)):
        
        for gene_index in range(samples.shape[1]):            
            if samples.iloc[i, gene_index] > 0 or any(x == 0 for x in dist):
                recovered_samples.iloc[i, gene_index] = samples.iloc[i, gene_index]
            else:
                recovered_samples.iloc[i, gene_index] = np.average(samples.iloc[np.array(ind), gene_index], weights=1/dist)

    recovered_samples_transposed = recovered_samples.T    
    
    recovered_samples_transposed.index = samples.columns
    recovered_samples_transposed.columns = samples.index    

    return recovered_samples_transposed #shape (ec_number, sample)
