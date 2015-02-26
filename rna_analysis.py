# -*- coding: utf-8 -*-
"""
Created on Thu Jan 22 07:33:32 2015

@author: nitin
"""
from __future__ import division
import sys
import os
import pandas as pd
import numpy as np
import scipy.stats as stats
from utils import error, makedirs
import pathway_analysis as pathway
from clean_data import CANCER_TYPES

def generate_rna_scores():
    """ Compute and save pathway level scores """
    for can_type in CANCER_TYPES:
        results_dir = '../results/' + can_type
        makedirs(results_dir)
        
        # Load preprocessed sampation data    
        data_dir = '../data/processed/' + can_type 
        sample_categories = ['PT', 'NB']
        
        path_df = pathway.preprocess_pathway_data()
        
        for samp in sample_categories:
            samp_path = data_dir + os.sep + can_type + '_' + samp + '_rnaseq.txt'
            if not os.path.exists(samp_path):
                print("processed sampation data is not found: %s" %can_type)
                continue
                
            samp_df = pd.read_table(samp_path, index_col=0, header=0)
            samp_df = (samp_df - samp_df.mean())/samp_df.std()
            print can_type, samp_df.shape
            score_df = pathway.score_pathways(path_df, samp_df, data_type='rna')
            score_df.dropna(how='all', inplace=True)
            score_df.to_csv(results_dir + os.sep + samp +'_rnaseq_pathway_score.txt', index=True, index_label='Pathway', sep='\t')
            print("Wrote " + can_type + " " + samp +'_rnaseq_pathway_score.txt %d x %d' %score_df.shape )

# Entry point, test code
if __name__ == '__main__':
    generate_rna_scores()
