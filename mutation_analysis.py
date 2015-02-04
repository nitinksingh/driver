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
    mut_categories = ['silent', 'nsilent']
    
    for mut in mut_categories:
        mut_path = data_dir + os.sep + can_type + '_' + mut + '_mutation.txt'
        if not os.path.exists(mut_path):
            print("processed mutation data is not found: %s" %can_type)
            return
            
        mut_df = pd.read_table(mut_path, index_col=0, header=0)
    
        score_df = pathway.score_pathways(mut_df, data_type='mut')
        score_df.to_csv(results_dir + os.sep + mut +'_mutation_pathway_score.txt', index=True, index_label='Pathway', sep='\t')
        print("Wrote " + can_type + " " + mut +'_mutation_pathway_score.txt %d x %d' %score_df.shape )

    

def get_pancancer_mutation_summary(agg_axis=0):
    """ Summarize mutation data across cancer types. 
    Parameters
    -----------------
    agg_axis: 0 Aggregate across genes for each sample.
    agg_axis: 1 Aggregate across samples for each GENE.    
    """
    cohort = "ACC BLCA BRCA CESC CHOL COAD COADREAD DLBC ESCA FPPP GBM GBMLGG HNSC KICH KIRC KIRP LAML LGG LIHC LUAD LUSC MESO OV PAAD PRAD READ SARC SKCM STAD TGCT THCA THYM UCEC UCS UVM"
    ns_summary = pd.DataFrame()
    s_summary = pd.DataFrame()
    for can in cohort.split():
        f1 = '../data/processed/' + can + os.sep + can + '_nsilent_mutation.txt'
        f2 = '../data/processed/' + can + os.sep + can + '_silent_mutation.txt'
        if not os.path.exists(f1) or not os.path.exists(f2):
            #print("Skipping %s " %can)
            continue

        nsdf = pd.read_table(f1, header=0, index_col=0)
        nsilent = nsdf.replace([0], [np.nan]).count(axis=agg_axis)

        sdf = pd.read_table(f2, header=0, index_col=0)
        silent = sdf.replace([0], [np.nan]).count(axis=agg_axis)

        if agg_axis == 1:
            silent = silent/sdf.shape[1]
            nsilent = nsilent/nsdf.shape[1]
        
        if ns_summary.empty:
            ns_summary[can] = nsilent
            s_summary[can] = silent
        else:
            ns_summary = ns_summary.join(pd.Series(nsilent, name=can), how='outer')
            s_summary = s_summary.join(pd.Series(silent, name=can), how='outer')

            
    return (ns_summary, s_summary)

if __name__ == '__main__':
    
    """ Gen pathway mutation correlation with non-silent mutation counts"""
    cohort = "ACC BLCA BRCA CESC CHOL COAD COADREAD DLBC ESCA FPPP GBM GBMLGG HNSC KICH KIRC KIRP LAML LGG LIHC LUAD LUSC MESO OV PAAD PRAD READ SARC SKCM STAD TGCT THCA THYM UCEC UCS UVM"
    can_type = cohort.split()
    
    for can in can_type:
        df = generate_mutation_scores(can)
        #analyze_pathway_mutations(can, 0)

    
