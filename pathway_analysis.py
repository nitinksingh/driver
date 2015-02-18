
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 18 16:27:20 2015

@author: Nitin Singh

Clean pathways from MSigDB
"""
from __future__ import division
import os
import pandas as pd
import numpy as np
from utils import *


def score_pathways(data_df, data_type, opt=''):
    """ Score each pathway against a given gene by samples matrix of 
    mutation/mRNASeq/CNA data. Gene should be index of the data frame i.e.
    there should not be any column with gene ids.
    """
    path_df = preprocess_pathway_data()
    score_df = pd.DataFrame(index=path_df.columns, columns=data_df.columns)
    for p in path_df.columns:
        try:
            # Get the genes that are present in data_df (gexp/cna/mut)
            avail_genes = filter(lambda x: x in data_df.index, path_df[p].dropna())
            p_df = data_df.loc[avail_genes]
            score_s = scorer(p_df, data_type, opt)
            score_df.loc[p] = score_s
        except KeyError:
            error()

    return score_df
        
def scorer(df, data_type, opt):
    """ Call individual scoring functions that are specific to data type"""
    if data_type == 'mut':
        sum_score = df.abs().sum(axis=0)
        ret_s = sum_score/df.shape[0]

        return ret_s
        
    if data_type == 'cna':
        error('Scoring Not implemented yet')

    if data_type == 'rna':
        error('Scoring Not implemented yet')
        
    return df
    
def preprocess_pathway_data():
    """Load GSEA MSigDB Broad Pathway DB"""
    input_dir = '../data/pathways'
    filename = input_dir + os.sep + 'kegg_biocarta_pid_positional.txt'
    pathways = {}
    with open(filename, 'r') as f:
        for line in f:
            p = line.strip().split('\t')
            pathways[p[0]] = p[2:]

    df = pd.DataFrame.from_dict(pathways, orient='index').transpose()
        
    return df

# Entry point, test code
if __name__ == '__main__':
    df = preprocess_pathway_data()