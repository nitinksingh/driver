# -*- coding: utf-8 -*-
"""
Created on Thu Jan 22 07:33:32 2015

@author: Nitin Singh
"""
from __future__ import division
import sys
import os
import pandas as pd
import numpy as np
import scipy.stats as stats
from utils import makedirs
import pathway_analysis as pathway
from clean_data import CANCER_TYPES

def generate_cna_scores(can_type):
    """ Compute and save pathway level cna scores """
    results_dir = '../results/' + can_type
    makedirs(results_dir)

    # Load preprocessed cna data
    data_dir = '../data/processed/' + can_type
    path_df = pathway.preprocess_pathway_data()
    cna_path = data_dir + os.sep + can_type + '_PT_cna.txt'

    try:
        cna_df = pd.read_table(cna_path, index_col=0, header=0)
    except IOError:
        print("processed cna data is not found: %s" %can_type)
        return

    score_df = pathway.score_pathways(path_df, cna_df, data_type='cna')
    score_df.dropna(how='all', inplace=True)
    score_df.to_csv(results_dir + os.sep + 'cna_pathway_score.txt', index=True, index_label='Pathway', sep='\t')
    print('Wrote ' + can_type + os.sep + 'cna_pathway_score.txt %d x %d' %score_df.shape )



def get_pancancer_cna_summary(agg_axis=0, refresh=False):
    """ Summarize cna data across cancer types.
    Parameters
    -----------------
    agg_axis: 0 or 1 Aggregate across Samples or Genes
    refresh: True/False Recompute cna frequencies from all cancer types
                or reload from previous computaion
    """
    outdir = '../data/pancancer'
    makedirs(outdir)
    if agg_axis:
        outf = outdir + os.sep + 'pancancer_cna_genes.txt'
    else:
        outf = outdir + os.sep + 'pancancer_cna_samples.txt'


    if os.path.exists(outf) and not refresh:
        summary = pd.read_table(outf, sep='\t', index_col=0, header=0)
        return summary

    summary = pd.DataFrame()
    for can in CANCER_TYPES:
        f1 = '../data/processed/' + can + os.sep + can + '_PT_cna.txt'
        if not os.path.exists(f1):
            #print("Skipping %s " %can)
            continue

        # Read CNAs and count extreme events ie -2 or +2, ignore +1, -1 CNAs
        df = pd.read_table(f1, header=0, index_col=0)
        df[df == -1] = 0
        df[df == 1] = 0
        cna_count = df.abs().sum(axis=agg_axis)/2


        # Compute cna frequencies (%) of gene across samples
        if agg_axis == 1:
            cna_count = cna_count*100/df.shape[1]

        if summary.empty:
            summary[can] = cna_count

        else:
            summary = summary.join(pd.Series(cna_count, name=can), how='outer')



    # Save for future usage
    summary.to_csv(outf, sep='\t', header=True, index=True)
    print("Wrote %s %d x %d matrix" %(os.path.basename(outf), len(summary.index), len(summary.columns)))


    return summary

if __name__ == '__main__':
    """ Gen pathway cna correlation with non-silent cna counts"""

    for can in CANCER_TYPES:
        generate_cna_scores(can)



