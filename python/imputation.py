import pandas as pd
import numpy
import os
import sys
import numpy as np
import argparse
import logging
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LinearRegression
logging.basicConfig(level=logging.DEBUG)
print("load module successfully")

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


def main():

    parser = _create_parser()
    args = parser.parse_args()

    # Load data
    adj_matrix = _LoadData(args.adj)
    print('shape of adj',adj_matrix.shape)
    ssms = _LoadData(args.ssms)
    print('shape of ssms',ssms.shape)
    mgs = _LoadData(args.mgs)
    print('shape of mgs',mgs.shape)
    model = args.model
    print('model',model)
    output = args.output
  
    # Predict features
    predicted_features = _predict_features(adj_matrix, ssms, mgs, model)
    print('shape of predicted_features',predicted_features.index)
    # Merge non-overlapped EC number
    ssms_non_overlap = ssms.loc[list(set(ssms.index) - set(predicted_features.index))]

    predicted_features_final = pd.concat([predicted_features, ssms_non_overlap], axis=0)

    predicted_features_final.to_csv(output, sep='\t')


if __name__ == "__main__":
    main()

