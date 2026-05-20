import pandas as pd
import numpy
import os
import sys
import numpy as np
import argparse
import logging
from sklearn.metrics.pairwise import cosine_similarity
#import matplotlib.pyplot as plt

from knn import _align_column, find_neighbors, recover_gene_from_neighbors

Italy1 = pd.read_csv('Italy1_500k_level4ec_overlap_completed.tsv', delimiter="\t", index_col=0, header=0) 
Italy2 = pd.read_csv('Italy2_500K_genefamily_level4ec.tsv', delimiter="\t", index_col=0, header=0) 
Austria = pd.read_csv('Austria_500k_level4ec_overlap_completed.tsv', delimiter="\t", index_col=0, header=0) 

Italy1_mgs = pd.read_csv('Italy1_level4ec.tsv', delimiter="\t", index_col=0, header=0)
Italy2_mgs = pd.read_csv('Italy2_gene_family_unstratified_level4ec.tsv', delimiter="\t", index_col=0, header=0)
Italy2_mgs = Italy2_mgs.rename(columns=lambda x: x.replace(".5", ".3"))
Austria_mgs = pd.read_csv('Austria_level4ec_overlap.tsv', delimiter="\t", index_col=0, header=0)

Italy1, Italy1_mgs = _align_column(Italy1, Italy1_mgs)
Italy2, Italy2_mgs = _align_column(Italy2, Italy2_mgs)
Austria, Austria_mgs = _align_column(Austria, Austria_mgs)



results = []
for data in ['Italy1', 'Italy2', 'Austria']: 
    
    #for n_neighbors in range(1, 51):
    for n_neighbors in [15]:

        if data == 'Italy1': 
            ssms_transposed = Italy1.T
            recoverd_ssms = recover_gene_from_neighbors(samples=ssms_transposed, n_neighbors=n_neighbors)
            
            ssms_cor = Italy1.corrwith(Italy1_mgs, method='spearman').mean() 
            ssms_rmse = np.sqrt(((Italy1 - Italy1_mgs) ** 2).mean().mean())

            recover_cor = recoverd_ssms.corrwith(Italy1_mgs, method='spearman').mean() 
            recover_rmse = np.sqrt(((recoverd_ssms - Italy1_mgs) ** 2).mean().mean())

            recover_genes = recoverd_ssms - Italy1
            nonzero_counts = (recover_genes != 0).sum(axis=1)
            nonzero_counts.to_csv('nonzero_counts_Italy1.csv', sep='\t')

        elif data == 'Italy2':
            ssms_transposed = Italy2.T
            recoverd_ssms = recover_gene_from_neighbors(samples=ssms_transposed, n_neighbors=n_neighbors)

            ssms_cor = Italy2.corrwith(Italy2_mgs, method='spearman').mean() 
            ssms_rmse = np.sqrt(((Italy2 - Italy2_mgs) ** 2).mean().mean())

            recover_cor = recoverd_ssms.corrwith(Italy2_mgs, method='spearman').mean() 
            recover_rmse = np.sqrt(((recoverd_ssms - Italy2_mgs) ** 2).mean().mean())

            recover_genes = recoverd_ssms - Italy2
            nonzero_counts = (recover_genes != 0).sum(axis=1)
            nonzero_counts.to_csv('nonzero_counts_Italy2.csv', sep='\t')

        elif data == 'Austria':
            ssms_transposed = Austria.T
            recoverd_ssms = recover_gene_from_neighbors(samples=ssms_transposed, n_neighbors=n_neighbors)
    
            ssms_cor = Austria.corrwith(Austria_mgs, method='spearman').mean() 
            ssms_rmse = np.sqrt(((Austria - Austria_mgs) ** 2).mean().mean())

            recover_cor = recoverd_ssms.corrwith(Austria_mgs, method='spearman').mean() 
            recover_rmse = np.sqrt(((recoverd_ssms - Austria_mgs) ** 2).mean().mean())

            recover_genes = recoverd_ssms - Austria
            nonzero_counts = (recover_genes != 0).sum(axis=1)
            nonzero_counts.to_csv('nonzero_counts_Austria.csv', sep='\t')

        results.append((data, n_neighbors, ssms_cor, recover_cor, ssms_rmse, recover_rmse))

results_df = pd.DataFrame(results, columns=['data', 'n_neighbors', 'ssms_cor', 'recover_cor', 'ssms_rmse', 'recover_rmse'])
#results_df.to_csv('knn_tuning_results.csv')

#recover_genes = recoverd_ssms - Italy1
#nonzero_counts = (recover_genes != 0).sum(axis=1)
#nonzero_counts.to_csv('nonzero_counts.csv', sep='\t')
