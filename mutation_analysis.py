# -*- coding: utf-8 -*-
"""
Created on Thu Jan 22 07:33:32 2015

@author: nitin
"""
from __future__ import division
import sys
import os
import math
import pandas as pd
import numpy as np
import scipy.stats as stats
from utils import *
import clean_data as clean
import pathway_analysis as pathway

def generate_mutation_scores(can_type):
    """ Compute and save pathway level mutation scores """
    results_dir = '../results/' + can_type
    makedirs(results_dir)
    
    # Load preprocessed mutation data    
    data_dir = '../data/processed/' + can_type 
    mut_path = data_dir + os.sep + can_type + '_mutation.txt'
    mut_df = pd.read_table(mut_path, index_col=0, header=0)

    # Compute pathway score for silent and non-silent mutations
    silent_s = mut_df.replace([-1, 0], [np.nan, np.nan]).count()
    nsilent_s = mut_df.replace([1, 0], [np.nan, np.nan]).count()

    opt_s = dict(zip(['nsilent', 'silent'], [nsilent_s, silent_s]))   
    
    for opt, s in opt_s.iteritems():
        score_df = pathway.score_pathways(mut_df, data_type='mut', opt=opt)
        score_df.to_csv(results_dir + os.sep + opt +'_pathway_score.txt', index=True, index_label='Pathway', sep='\t')
        print("Wrote " + can_type + " " + opt +'_pathway_score.txt %d x %d' %score_df.shape )

    
def analyze_pathway_mutations(can_type, refresh=0):
    """ Compute correlation between pathway wise mutation scores and non-silent
    mutation count for the given patient population.
    """
    results_dir = '../results/' + can_type
    result_file = 'pathway_mutation_score_vs_mutation_count.txt'
    if not os.path.exists(results_dir):
        os.makedirs(os.path.abspath(results_dir))
    
    # Use previously generated calculations, if exists
    prefix = 'rho_' 
    result_path = results_dir + os.sep + prefix + result_file
    if os.path.exists(result_path) and not refresh:
        pathway_corr_df =  clean.read_matrix_format_data(result_path, index_col=0, header=0)
        return pathway_corr_df
    
    data_dir = '../data/processed/' + can_type 
    mut_path = data_dir + os.sep + can_type + '_mutation.txt'
    
    path_df = pathway.preprocess_pathway_data()
    mut_df = clean.read_matrix_format_data(mut_path, index_col=0)

    # Now compute correlation of TOTAL mutations in a sample vs. mutation count 
    # in specific pathways
    silent_s = mut_df.replace([-1, 0], [np.nan, np.nan]).count()
    nsilent_s = mut_df.replace([1, 0], [np.nan, np.nan]).count()

    opt_s = dict(zip(['nsilent', 'silent'], [nsilent_s, silent_s]))   
    pathway_corr_df = pd.DataFrame(index=path_df.columns)
    for opt, s in opt_s.iteritems():
        score_path = results_dir + os.sep + opt +'_pathway_score.txt'
        score_df = pd.read_table(score_path, header=0, index_col=0) 

        for p in score_df.index:
            rho_pval = stats.spearmanr(score_df.loc[p], s)
            pathway_corr_df.loc[p, 'rho'+'_'+ opt] = rho_pval[0]
            pathway_corr_df.loc[p, 'pval'+'_'+ opt] = -1*math.log10(rho_pval[1])

    pathway_corr_df.insert(0, 'pathway', pathway_corr_df.index)
    save_df(pathway_corr_df, result_file, results_dir, prefix=prefix)
    
    return pathway_corr_df
    

if __name__ == '__main__':
    
    """ Gen pathway mutation correlation with non-silent mutation counts"""
    can_type = ['LUAD', 'LUSC']
    for can in can_type:
        df = generate_mutation_scores(can)
        #analyze_pathway_mutations(can, 0)

    
