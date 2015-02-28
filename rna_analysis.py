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
        
        # Load preprocessed mrna data    
        data_dir = '../data/processed/' + can_type 
        sample_categories = ['PT', 'NB']
        
        path_df = pathway.preprocess_pathway_data()
        
        for samp in sample_categories:
            samp_path = data_dir + os.sep + can_type + '_' + samp + '_rnaseq.txt'
            if not os.path.exists(samp_path):
                print("processed mrna data is not found: %s" %can_type)
                continue
                
            samp_df = pd.read_table(samp_path, index_col=0, header=0)
            samp_df = (samp_df - samp_df.mean())/samp_df.std()
            print can_type, samp_df.shape
            score_df = pathway.score_pathways(path_df, samp_df, data_type='rna')
            score_df.dropna(how='all', inplace=True)
            score_df.to_csv(results_dir + os.sep + samp +'_rnaseq_pathway_score.txt', index=True, index_label='Pathway', sep='\t', float_format='%.2f')
            print("Wrote " + can_type + " " + samp +'_rnaseq_pathway_score.txt %d x %d' %score_df.shape )

def generate_rna_scores_centroid():
    """ Compute and save pathway score wrt to centroid of the normal samples
    """
    for can_type in CANCER_TYPES:
        results_dir = '../results/' + can_type
        makedirs(results_dir)
        
        # Load preprocessed data    
        data_dir = '../data/processed/' + can_type 
        
        pt_path = data_dir + os.sep + can_type + '_PT_rnaseq.txt'
        nb_path = data_dir + os.sep + can_type + '_NB_rnaseq.txt'
        if not (os.path.exists(pt_path) and os.path.exists(nb_path)):
            print("processed mrna data is not found: %s" %can_type)
            continue
        # Read and normalize mRNASeq RSEM values, replace NA with 0
        nb_df  = pd.read_table(nb_path, index_col=0, header=0)
        nb_df = (nb_df - nb_df.mean())/nb_df.std()
        nb_df.fillna(0, inplace=True)

        pt_df  = pd.read_table(pt_path, index_col=0, header=0)
        pt_df = (pt_df - pt_df.mean())/pt_df.std()
        pt_df.fillna(0, inplace=True)
        print can_type, pt_df.shape, nb_df.shape
            
        path_df = pathway.preprocess_pathway_data()

        score_df = pd.DataFrame(index=path_df.columns, columns=pt_df.columns)
        for p in path_df.columns:
            try:
                # Get the genes that are present in data_df (gexp/cna/mut)
                ppt_df = pt_df.loc[path_df[p].dropna()]
                pnb_df = nb_df.loc[path_df[p].dropna()]
            except KeyError:
                avail_genes = filter(lambda x: x in pt_df.index, path_df[p].dropna())
                ppt_df = pt_df.loc[avail_genes]
                pnb_df = nb_df.loc[avail_genes]
            # Compute Euclidean distance from the centroid of normal samples
            center = pnb_df.mean(axis=1)
            dist_df = ppt_df.sub(center, axis=0)
            score_s = np.sqrt(dist_df.applymap(lambda x: x**2).sum(axis=0))
            score_df.loc[p] = score_s/dist_df.shape[0]
        score_df.dropna(how='all', inplace=True)
        score_df.to_csv(results_dir + os.sep + 'PT_centroid_rnaseq_pathway_score.txt', index=True, index_label='Pathway', sep='\t', float_format='%.2f')
        print("Wrote " + can_type + " " + 'PT_centroid_rnaseq_pathway_score.txt %d x %d' %score_df.shape )
# Entry point, test code
if __name__ == '__main__':
    generate_rna_scores_centroid()
    generate_rna_scores()
