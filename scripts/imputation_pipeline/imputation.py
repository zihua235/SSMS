import pandas as pd
import numpy
import os
import sys
import numpy as np
import argparse
import logging
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LinearRegression
from knn import find_neighbors, recover_gene_from_neighbors
from joblib import Parallel, delayed
import multiprocessing
from tqdm import tqdm
logging.basicConfig(level=logging.DEBUG)

current_dir = os.path.dirname(os.path.abspath(__file__))
output_path = os.path.dirname(current_dir)


def _create_parser():
    """
    Create an argument parser for the script.
    """
    parser = argparse.ArgumentParser(description="Perform basic arithmetic operations")

    # 添加命令行参数
    parser.add_argument('--ssms', type=str, help="path to the shallow metagenomics file")
    parser.add_argument('--mgs', type=str, help="path to the metagenomics file")
    parser.add_argument('--adj', type=str, help="path to the adj matrix file")
    parser.add_argument('--model', choices=['RF', 'LR'], default='RF', help="model type, RF or LR, default is RF")  
    parser.add_argument('--output', type=str, help="output file name, csv format")
    return parser

def _LoadData(file_path:str):
    """
    Load data (including csv, tsv, txt, xlsx) from path, and detect if there's index and column

    Args:
    - feature: str, path to the feature file

    Returns:
    - pd.DataFrame, the feature matrix
    """

    file_ext = os.path.splitext(file_path)[-1].lower()
    
    if file_ext == ".csv":
        df = pd.read_csv(file_path, index_col=0, header=0)  
    elif file_ext == ".tsv":
        df = pd.read_csv(file_path, delimiter="\t", index_col=0, header=0) 
    elif file_ext == ".txt":
        df = pd.read_csv(file_path, delimiter="\t", index_col=0, header=0)  
    elif file_ext == ".xlsx":
        df = pd.read_excel(file_path, index_col=0, header=0)
    else:
        raise ValueError("unsupported file format {file_ext}, please provide .csv, .tsv, .txt, .xlsx file")

    return df

def _align_row_column(adj:pd.DataFrame, df:pd.DataFrame):
    """
    Align the row and column of two dataframes. Let both databases have the same index

    Args:
    - feature: pd.dataframe

    Returns:
    - pd.DataFrame, the overlapped index of two dataframes
    """
    adj_index = adj.index
    df_index = df.index

    intersection = adj_index.intersection(df_index) # sort first, EC number that exit in both org_EC and EC_aj
    
    adj_overlap = adj.loc[pd.Index(intersection), pd.Index(intersection)]
    df_overlap = df.loc[intersection]
    df_overlap = df_overlap.reindex(adj_overlap.index) # Ensure df2 has the same index as df1

    if not adj_overlap.index.equals(df_overlap.index):
        raise ValueError("Error: adj matrix and ssms/mgs have different indexes! Terminating program.")
    
    return adj_overlap, df_overlap

def create_model(model_type: str, correlated_ec: pd.DataFrame, predicted_ec: pd.DataFrame):
    if model_type == 'RF':
        model = RandomForestRegressor(n_estimators=500, random_state=666)
    elif model_type == 'LR':
        model = LinearRegression()
    else:
        raise ValueError(f"Unsupported model type: {model_type}")
    model.fit(correlated_ec, predicted_ec)
    return model

def _predict_features(adj_matrix: pd.DataFrame, ssms: pd.DataFrame, mgs: pd.DataFrame, model_type: str):

    adj_df, ssms_df = _align_row_column(adj_matrix, ssms)
    adj_df, mgs_df = _align_row_column(adj_df, mgs)
    logging.debug(f"CHECK if adj and ssms have same index: {adj_df.index.isin(ssms_df.index).all()}")
    logging.debug(f"CHECK if adj and mgs have same index: {adj_df.index.isin(mgs_df.index).all()}")


    n_features, n_samples = ssms_df.shape
    predicted_features = pd.DataFrame(np.nan, index=range(n_features), columns=range(n_samples))  # Initialize a zero matrix for predictions

    for feature in adj_df.index:
        logging.debug(f"feature: {feature}")
        # Find the neighbors of the current feature (those with a non-zero adjacency)
        selected_row = adj_df.loc[feature]
        neighbors = selected_row[selected_row > 0].index
        logging.debug(f"neighbors: {neighbors}")
        logging.debug(f"length of neighbors: {len(neighbors)}")
        if len(neighbors) == 0:  # If there are no neighbors, skip this feature
             model = None
        else: # Otherwise, create a model using the neighbor features
             logging.debug(f"CHECK if all neighbors are in mgs index: {neighbors.isin(mgs.index).all()}")
             #logging.debug(f"model_input: {mgs_df.loc[neighbors].T.shape}")
             #logging.debug(f"model_output: {mgs_df.loc[feature].T.shape}")
             model = create_model(model_type, mgs_df.loc[neighbors].T, mgs_df.loc[feature].T)  # Create RF model
            
        # Iterate over each sample in the ssms
        for sample in ssms_df.columns:
            logging.debug(f"sample: {sample}")
            if len(neighbors) == 0:   # there are no neighbors
                predicted_features.loc[feature, sample] = ssms_df.loc[feature, sample]
            elif ssms.loc[feature, sample] > 0:  # Keep the existing value if non-zero
                predicted_features.loc[feature, sample] = ssms_df.loc[feature, sample]
            elif model is None:  # If the value for the current feature and sample is missing (assumed to be zero here)
                predicted_features.loc[feature, sample] = ssms_df.loc[feature, sample]
            else:
                predicted_features.loc[feature, sample] = model.predict([ssms_df.loc[neighbors, sample]])[0]
                logging.debug(f"imputated value: {model.predict([ssms_df.loc[neighbors, sample]])[0]}")

    # Return the predicted feature matrix as a DataFrame
    return pd.DataFrame(predicted_features, index=ssms_df.index, columns=ssms_df.columns)


def _predict_features_quick(adj_matrix: pd.DataFrame, ssms: pd.DataFrame, mgs: pd.DataFrame, model_type: str, n_jobs: int = 64):
    # 对齐数据
    adj_df, ssms_df = _align_row_column(adj_matrix, ssms)
    adj_df, mgs_df = _align_row_column(adj_df, mgs)

    def predict_one_feature(feature):
        try:
            selected_row = adj_df.loc[feature]
            neighbors = selected_row[selected_row > 0].index

            pred_row = ssms_df.loc[feature].copy()

            if len(neighbors) == 0 or not neighbors.isin(mgs_df.index).all():
                return feature, pred_row

            model = create_model(model_type, mgs_df.loc[neighbors].T, mgs_df.loc[feature].T)

            for sample in ssms_df.columns:
                if ssms_df.loc[feature, sample] > 0:
                    continue
                sample_data = ssms_df.loc[neighbors, sample].to_frame().T
                pred_row[sample] = model.predict(sample_data)[0]


            return feature, pred_row
        except Exception as e:
            logging.warning(f"Failed to process feature {feature}: {e}")
            return feature, ssms_df.loc[feature]

    features = list(adj_df.index)

    n_jobs = min(64, multiprocessing.cpu_count())
    results = Parallel(n_jobs=n_jobs)(
    delayed(predict_one_feature)(feature) for feature in tqdm(features, desc="Predicting features"))

    predicted_df = pd.DataFrame({f: row for f, row in results}).T
    predicted_df.index.name = "EC"
    predicted_df.columns = ssms_df.columns

    return predicted_df



def _predict_features_with_recovry(adj_matrix: pd.DataFrame, ssms: pd.DataFrame, mgs: pd.DataFrame, model_type: str, n_neighbors:int):
    
    ssms_transposed = ssms.T
    recoverd_ssms = recover_gene_from_neighbors(samples=ssms_transposed, n_neighbors=n_neighbors)
    assert recoverd_ssms.shape == ssms.shape

    adj_df, recoverd_ssms_df = _align_row_column(adj_matrix, recoverd_ssms)
    adj_df, ssms_df = _align_row_column(adj_matrix, ssms)
    adj_df, mgs_df = _align_row_column(adj_df, mgs)
    logging.debug(f"CHECK if adj and ssms have same index: {adj_df.index.isin(ssms_df.index).all()}")
    logging.debug(f"CHECK if adj and mgs have same index: {adj_df.index.isin(mgs_df.index).all()}")


    n_features, n_samples = ssms_df.shape
    predicted_features = pd.DataFrame(np.nan, index=range(n_features), columns=range(n_samples))  # Initialize a zero matrix for predictions

    for feature in adj_df.index:
        logging.debug(f"feature: {feature}")
        # Find the neighbors of the current feature (those with a non-zero adjacency)
        selected_row = adj_df.loc[feature]
        neighbors = selected_row[selected_row > 0].index
        logging.debug(f"neighbors: {neighbors}")
        logging.debug(f"length of neighbors: {len(neighbors)}")
        if len(neighbors) == 0:  # If there are no neighbors, skip this feature
             model = None
        else: # Otherwise, create a model using the neighbor features
             logging.debug(f"CHECK if all neighbors are in mgs index: {neighbors.isin(mgs.index).all()}")
             #logging.debug(f"model_input: {mgs_df.loc[neighbors].T.shape}")
             #logging.debug(f"model_output: {mgs_df.loc[feature].T.shape}")
             model = create_model(model_type, mgs_df.loc[neighbors].T, mgs_df.loc[feature].T)  # Create RF model
            
        # Iterate over each sample in the ssms
        for sample in ssms_df.columns:
            logging.debug(f"sample: {sample}")
            if len(neighbors) == 0:   # there are no neighbors
                predicted_features.loc[feature, sample] = ssms_df.loc[feature, sample]
            elif ssms.loc[feature, sample] > 0:  # Keep the existing value if non-zero
                predicted_features.loc[feature, sample] = ssms_df.loc[feature, sample]
            elif model is None:  # If the value for the current feature and sample is missing (assumed to be zero here)
                predicted_features.loc[feature, sample] = ssms_df.loc[feature, sample]
            else:
                predicted_features.loc[feature, sample] = model.predict([recoverd_ssms_df.loc[neighbors, sample]])[0]
                logging.debug(f"imputated value: {model.predict([recoverd_ssms_df.loc[neighbors, sample]])[0]}")

    return pd.DataFrame(predicted_features, index=ssms_df.index, columns=ssms_df.columns)


def _predict_features_without_missing(adj_matrix: pd.DataFrame, ssms: pd.DataFrame, mgs: pd.DataFrame, model_type: str, n_neighbors:int):
    
    ssms_transposed = ssms.T
    recoverd_ssms = recover_gene_from_neighbors(samples=ssms_transposed, n_neighbors=n_neighbors)
    assert recoverd_ssms.shape == ssms.shape

    adj_df, recoverd_ssms_df = _align_row_column(adj_matrix, recoverd_ssms)
    adj_df, ssms_df = _align_row_column(adj_matrix, ssms)
    adj_df, mgs_df = _align_row_column(adj_df, mgs)
    logging.debug(f"CHECK if adj and ssms have same index: {adj_df.index.isin(ssms_df.index).all()}")
    logging.debug(f"CHECK if adj and mgs have same index: {adj_df.index.isin(mgs_df.index).all()}")


    n_features, n_samples = ssms_df.shape
    predicted_features = pd.DataFrame(np.nan, index=range(n_features), columns=range(n_samples))  # Initialize a zero matrix for predictions

    for feature in adj_df.index:
        logging.debug(f"feature: {feature}")
        # Find the neighbors of the current feature (those with a non-zero adjacency)
        selected_row = adj_df.loc[feature]
        neighbors = selected_row[selected_row > 0].index
        logging.debug(f"neighbors: {neighbors}")
        logging.debug(f"length of neighbors: {len(neighbors)}")
        if len(neighbors) == 0:  # If there are no neighbors, skip this feature
             model = None
        else: # Otherwise, create a model using the neighbor features
             logging.debug(f"CHECK if all neighbors are in mgs index: {neighbors.isin(mgs.index).all()}")
             #logging.debug(f"model_input: {mgs_df.loc[neighbors].T.shape}")
             #logging.debug(f"model_output: {mgs_df.loc[feature].T.shape}")
             model = create_model(model_type, mgs_df.loc[neighbors].T, mgs_df.loc[feature].T)  # Create RF model
            
        # Iterate over each sample in the ssms
        for sample in ssms_df.columns:
            logging.debug(f"sample: {sample}")
            if len(neighbors) == 0:   # there are no neighbors
                predicted_features.loc[feature, sample] = ssms_df.loc[feature, sample]
            elif ssms.loc[feature, sample] > 0:  # Keep the existing value if non-zero
                predicted_features.loc[feature, sample] = ssms_df.loc[feature, sample]
            elif model is None:  # If the value for the current feature and sample is missing (assumed to be zero here)
                predicted_features.loc[feature, sample] = ssms_df.loc[feature, sample]
            elif (recoverd_ssms_df.loc[neighbors, sample] > 0).all():
                predicted_features.loc[feature, sample] = model.predict([recoverd_ssms_df.loc[neighbors, sample]])[0]
                logging.debug(f"imputated value: {model.predict([recoverd_ssms_df.loc[neighbors, sample]])[0]}") 
            elif not (recoverd_ssms_df.loc[neighbors, sample] > 0).all():
                for neighbor in neighbors:
                    if recoverd_ssms_df.loc[neighbor, sample] == 0:
                        neighbors_of_neigbor_row = adj_df.loc[neighbor]
                        neighbors_of_neigbor = neighbors_of_neigbor_row[neighbors_of_neigbor_row  > 0].index
                        if len(neighbors_of_neigbor) == 0:  # If there are no neighbors, skip this feature
                            neighbor_model = None
                            recoverd_ssms_df.loc[neighbors, sample] = recoverd_ssms_df.loc[neighbors, sample]
                        else: # Otherwise, create a model using the neighbor features
                            neighbor_model = create_model(model_type, mgs_df.loc[neighbors_of_neigbor].T, mgs_df.loc[neighbor].T)  # Create RF model
                            predicted_neighbor = neighbor_model.predict([recoverd_ssms_df.loc[neighbors_of_neigbor, sample]])[0]
                            recoverd_ssms_df.loc[neighbor, sample] = predicted_neighbor
                    else:
                        pass
                predicted_features.loc[feature, sample] = model.predict([recoverd_ssms_df.loc[neighbors, sample]])[0]
                logging.debug(f"imputated value: {model.predict([recoverd_ssms_df.loc[neighbors, sample]])[0]}") 
            else:
                predicted_features.loc[feature, sample] = ssms_df.loc[feature, sample]
                
                
    return pd.DataFrame(predicted_features, index=ssms_df.index, columns=ssms_df.columns)


def _predict_features_with_recovry_quick(adj_matrix: pd.DataFrame,
                                            ssms: pd.DataFrame,
                                            mgs: pd.DataFrame,
                                            model_type: str,
                                            n_neighbors: int,
                                            n_jobs: int = 64):
    # Step 1: Recover missing values from neighbors
    ssms_transposed = ssms.T
    recoverd_ssms = recover_gene_from_neighbors(samples=ssms_transposed, n_neighbors=n_neighbors)
    assert recoverd_ssms.shape == ssms.shape

    # Step 2: Align all matrices
    adj_df, recoverd_ssms_df = _align_row_column(adj_matrix, recoverd_ssms)
    adj_df, ssms_df = _align_row_column(adj_matrix, ssms)
    adj_df, mgs_df = _align_row_column(adj_df, mgs)

    # Step 3: Prediction function for one feature
    def predict_one_feature(feature):
        try:
            selected_row = adj_df.loc[feature]
            neighbors = selected_row[selected_row > 0].index

            if len(neighbors) == 0 or not neighbors.isin(mgs_df.index).all():
                model = None
            else:
                model = create_model(model_type, mgs_df.loc[neighbors].T, mgs_df.loc[feature].T)

            # Predict for each sample
            pred_values = []
            for sample in ssms_df.columns:
                if len(neighbors) == 0:
                    pred_values.append(ssms_df.loc[feature, sample])
                elif ssms_df.loc[feature, sample] > 0:
                    pred_values.append(ssms_df.loc[feature, sample])
                elif model is None:
                    pred_values.append(ssms_df.loc[feature, sample])
                else:
                    sample_data = recoverd_ssms_df.loc[neighbors, sample]
                    pred = model.predict([sample_data])[0]
                    pred_values.append(pred)
            return feature, pred_values
        except Exception as e:
            logging.warning(f"Failed to process feature {feature}: {e}")
            return feature, list(ssms_df.loc[feature, :])  # fallback

    # Step 4: Parallel execution
    features = list(adj_df.index)
    n_jobs = min(n_jobs, multiprocessing.cpu_count())

    results = Parallel(n_jobs=n_jobs)(
        delayed(predict_one_feature)(feature) for feature in tqdm(features, desc="Predicting features")
    )

    # Step 5: Assemble results into DataFrame
    predicted_df = pd.DataFrame({f: row for f, row in results}).T
    predicted_df.index = ssms_df.index
    predicted_df.columns = ssms_df.columns

    return predicted_df

def _predict_features_without_missing_quick(adj_matrix: pd.DataFrame,
                                               ssms: pd.DataFrame,
                                               mgs: pd.DataFrame,
                                               model_type: str,
                                               n_neighbors: int,
                                               n_jobs: int = 64):
    # Step 1: 基于邻居恢复输入
    ssms_transposed = ssms.T
    recoverd_ssms = recover_gene_from_neighbors(samples=ssms_transposed, n_neighbors=n_neighbors)
    assert recoverd_ssms.shape == ssms.shape

    # Step 2: 对齐矩阵
    adj_df, recoverd_ssms_df = _align_row_column(adj_matrix, recoverd_ssms)
    adj_df, ssms_df = _align_row_column(adj_matrix, ssms)
    adj_df, mgs_df = _align_row_column(adj_df, mgs)

    # Step 3: 每个 feature 的预测逻辑封装为函数
    def predict_feature(feature):
        try:
            selected_row = adj_df.loc[feature]
            neighbors = selected_row[selected_row > 0].index
            pred_row = []

            if len(neighbors) == 0 or not neighbors.isin(mgs_df.index).all():
                model = None
            else:
                model = create_model(model_type, mgs_df.loc[neighbors].T, mgs_df.loc[feature].T)

            for sample in ssms_df.columns:
                if len(neighbors) == 0:
                    pred_row.append(ssms_df.loc[feature, sample])
                    continue

                if ssms_df.loc[feature, sample] > 0 or model is None:
                    pred_row.append(ssms_df.loc[feature, sample])
                    continue

                # 如果所有 neighbors 都有值，直接预测
                if (recoverd_ssms_df.loc[neighbors, sample] > 0).all():
                    pred = model.predict([recoverd_ssms_df.loc[neighbors, sample]])[0]
                    pred_row.append(pred)
                else:
                    # 否则对缺失 neighbor 进行回溯补全
                    for neighbor in neighbors:
                        if recoverd_ssms_df.loc[neighbor, sample] == 0:
                            n_row = adj_df.loc[neighbor]
                            n_neighbors2 = n_row[n_row > 0].index
                            if len(n_neighbors2) == 0:
                                continue
                            try:
                                n_model = create_model(model_type,
                                                       mgs_df.loc[n_neighbors2].T,
                                                       mgs_df.loc[neighbor].T)
                                pred_neighbor = n_model.predict([recoverd_ssms_df.loc[n_neighbors2, sample]])[0]
                                recoverd_ssms_df.loc[neighbor, sample] = pred_neighbor
                            except Exception as e:
                                logging.warning(f"Failed to impute neighbor {neighbor}: {e}")

                    # 回补完后再预测目标值
                    pred = model.predict([recoverd_ssms_df.loc[neighbors, sample]])[0]
                    pred_row.append(pred)

            return feature, pred_row

        except Exception as e:
            logging.warning(f"Failed to process feature {feature}: {e}")
            return feature, list(ssms_df.loc[feature, :])

    # Step 4: 并行处理每个 feature
    features = list(adj_df.index)
    n_jobs = min(n_jobs, multiprocessing.cpu_count())

    results = Parallel(n_jobs=n_jobs)(
        delayed(predict_feature)(feature) for feature in tqdm(features, desc="Predicting (with recovry)")
    )

    # Step 5: 重建矩阵
    predicted_df = pd.DataFrame({f: row for f, row in results}).T
    predicted_df.index = ssms_df.index
    predicted_df.columns = ssms_df.columns

    return predicted_df

def _calculate_missing_value_portion(adj_matrix: pd.DataFrame, ssms: pd.DataFrame, mgs: pd.DataFrame, model_type: str, n_neighbors:int):

    #adj_df, ssms_df = _align_row_column(adj_matrix, ssms)
    #adj_df, mgs_df = _align_row_column(adj_df, mgs)
    #logging.debug(f"CHECK if adj and ssms have same index: {adj_df.index.isin(ssms_df.index).all()}")
    #logging.debug(f"CHECK if adj and mgs have same index: {adj_df.index.isin(mgs_df.index).all()}")

    ssms_transposed = ssms.T
    recoverd_ssms = recover_gene_from_neighbors(samples=ssms_transposed, n_neighbors=n_neighbors)
    assert recoverd_ssms.shape == ssms.shape

    adj_df, recoverd_ssms_df = _align_row_column(adj_matrix, recoverd_ssms)
    adj_df, ssms_df = _align_row_column(adj_matrix, ssms)
    adj_df, mgs_df = _align_row_column(adj_df, mgs)
    logging.debug(f"CHECK if adj and ssms have same index: {adj_df.index.isin(ssms_df.index).all()}")
    logging.debug(f"CHECK if adj and mgs have same index: {adj_df.index.isin(mgs_df.index).all()}")

    n_features, n_samples = ssms_df.shape
    missing_features = []

    for feature in adj_df.index:
        logging.debug(f"feature: {feature}")
        # Find the neighbors of the current feature (those with a non-zero adjacency)
        selected_row = adj_df.loc[feature]
        neighbors = selected_row[selected_row > 0].index
        logging.debug(f"neighbors: {neighbors}")
        logging.debug(f"length of neighbors: {len(neighbors)}")
        if len(neighbors) == 0:  # If there are no neighbors, skip this feature
             model = None
        else: # Otherwise, create a model using the neighbor features
             logging.debug(f"CHECK if all neighbors are in mgs index: {neighbors.isin(mgs.index).all()}")
             #logging.debug(f"model_input: {mgs_df.loc[neighbors].T.shape}")
             #logging.debug(f"model_output: {mgs_df.loc[feature].T.shape}")
             #model = create_model(model_type, mgs_df.loc[neighbors].T, mgs_df.loc[feature].T)  # Create RF model
            
        # Iterate over each sample in the ssms
        for sample in ssms_df.columns:
            logging.debug(f"sample: {sample}")
            if len(neighbors) == 0:   # there are no neighbors
                pass
            elif ssms.loc[feature, sample] > 0:  # Keep the existing value if non-zero
                pass
            else:
                missing_portion = (recoverd_ssms_df.loc[neighbors, sample] == 0).sum() / len(neighbors)
                missing_features.append(missing_portion)

    # Return the predicted feature matrix as a DataFrame
    return pd.DataFrame(missing_features)

