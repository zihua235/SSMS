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
import multiprocessing
from joblib import Parallel, delayed
from tqdm import tqdm
import time

from imputation import _create_parser, _LoadData, _predict_features_quick

start = time.time()
print("needed time：", time.time() - start)
print(multiprocessing.cpu_count())
print(os.environ.get("SLURM_CPUS_PER_TASK"))

current_dir = os.path.dirname(os.path.abspath(__file__))
output_path = os.path.dirname(current_dir)


def main():

    parser = _create_parser()
    args = parser.parse_args()

    # Load data
    adj_matrix = _LoadData(args.adj)
    ssms = _LoadData(args.ssms)
    mgs = _LoadData(args.mgs)
    model = args.model
    output = args.output

    # Predict features
    predicted_features = _predict_features_quick(adj_matrix, ssms, mgs, model, n_jobs=multiprocessing.cpu_count())
    
    # Merge non-overlapped EC number
    ssms_non_overlap = ssms.loc[list(set(ssms.index) - set(predicted_features.index))]

    predicted_features_final = pd.concat([predicted_features, ssms_non_overlap], axis=0)

    predicted_features_final.to_csv(output, sep='\t')
    #predicted_features.to_csv(output, sep='\t')

if __name__ == "__main__":
    main()
