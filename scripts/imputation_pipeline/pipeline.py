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

from knn import find_neighbors, recover_gene_from_neighbors

from imputation import _create_parser, _LoadData, _align_row_column, _predict_features_quick, _predict_features_with_recovry, _predict_features_without_missing, _calculate_missing_value_portion,_predict_features_with_recovry_quick,_predict_features_without_missing_quick
import multiprocessing
import time
start = time.time()
print("needed time：", time.time() - start)
print(multiprocessing.cpu_count())
print(os.environ.get("SLURM_CPUS_PER_TASK"))

current_dir = os.path.dirname(os.path.abspath(__file__))
output_path = os.path.dirname(current_dir)

def _create_extended_parser():
    """
    Create an extended argument parser with an additional parameter.
    """
    parser = _create_parser()  # 调用原来的 _create_parser 函数
    
    # 添加额外的命令行参数
    parser.add_argument('--n_neighbors', type=int, default=5, help="number of neighbors to consider. Default is 5")
    
    return parser

def main():

    parser = _create_extended_parser()
    args = parser.parse_args()

    # Load data
    adj_matrix = _LoadData(args.adj)
    ssms = _LoadData(args.ssms)
    mgs = _LoadData(args.mgs)
    model = args.model
    output = args.output
    n_neighbors = args.n_neighbors

    # Predict features

    #predicted_features = _predict_features_with_recovry(adj_matrix, ssms, mgs, model, n_neighbors)
    #predicted_features = _predict_features_with_recovry_quick(adj_matrix, ssms, mgs, model, n_neighbors)
    #  
    #predicted_features = _predict_features_without_missing(adj_matrix, ssms, mgs, model, n_neighbors)
    predicted_features = _predict_features_without_missing_quick(adj_matrix, ssms, mgs, model, n_neighbors)

    #predicted_features = _calculate_missing_value_portion(adj_matrix, ssms, mgs, model, n_neighbors)
    #predicted_features = _predict_features_quick(adj_matrix, ssms, mgs, model, n_neighbors)

    # Merge non-overlapped EC number
    ssms_non_overlap = ssms.loc[list(set(ssms.index) - set(predicted_features.index))]

    predicted_features_final = pd.concat([predicted_features, ssms_non_overlap], axis=0)

    predicted_features_final.to_csv(output, sep='\t')
    #predicted_features.to_csv(output, sep='\t')

if __name__ == "__main__":
    main()
