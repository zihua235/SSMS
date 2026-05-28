import pandas as pd
import numpy
import os
import sys
import numpy as np
import argparse
import logging
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LinearRegression
from joblib import Parallel, delayed
import multiprocessing
from tqdm import tqdm
logging.basicConfig(level=logging.DEBUG)


current_dir = os.path.dirname(os.path.abspath(__file__))
output_path = os.path.dirname(current_dir)

def _create_parser():

    parser = argparse.ArgumentParser()

    parser.add_argument('--ssms', type=str, help="path to the unstratified relative level4 ec number profile of shallow metagenomics")
    parser.add_argument('--mgs', type=str, help="path to the unstratified relative level4 ec number profile of metagenomics")
    parser.add_argument('--adj', type=str, help="path to the adjacency matrix file")
    parser.add_argument('--model', choices=['RF', 'LR'], default='RF', help="model type, Random Forest (RF) or Linear Regression (LR), default is RF")  
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

def _align_index(adj:pd.DataFrame, df:pd.DataFrame):
    """
    Align the indexof two dataframes. Let both databases have the same index

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

def _predict_features_quick(adj_matrix: pd.DataFrame, ssms: pd.DataFrame, mgs: pd.DataFrame, model_type: str, n_jobs: int = 64):
    
    """
    for ssms, row is feature, column is sample.

    for mgs, row is feature, column is sample.
    """
     
    # Align the row of two dataframes. Let both databases have the same index
    adj_df, ssms_df = _align_index(adj_matrix, ssms)
    adj_df, mgs_df = _align_index(adj_df, mgs)

    def predict_one_feature(feature):
        try:
            selected_row = adj_df.loc[feature]
            neighbors = selected_row[selected_row > 0].index

            pred_row = ssms_df.loc[feature].copy()

            if len(neighbors) == 0 or not neighbors.isin(mgs_df.index).all():
                return feature, pred_row

            # input format of the model: X: neighbor features, shape (n_samples, n_neighbors), y: target feature, shape (n_samples,1)
            model = create_model(model_type, mgs_df.loc[neighbors].T, mgs_df.loc[feature].T)

            for sample in ssms_df.columns:
                # if the value > 0, keep the original value
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

