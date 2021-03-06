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
from utils import *
import pathway_analysis as pathway
from clean_data import CANCER_TYPES

def generate_mutation_scores():
    """ Compute and save pathway level mutation scores """
    for can_type in CANCER_TYPES:
        results_dir = '../results/' + can_type
        makedirs(results_dir)

        # Load preprocessed mutation data
        data_dir = '../data/processed/' + can_type
        mut_categories = ['silent', 'nsilent']

        path_df = pathway.preprocess_pathway_data()

        for mut in mut_categories:
            mut_path = data_dir + os.sep + can_type + '_' + mut + '_mutation.txt'
            if not os.path.exists(mut_path):
                print("processed mutation data is not found: %s" %can_type)
                continue

            mut_df = pd.read_table(mut_path, index_col=0, header=0)
            print can_type, mut_df.shape
            score_df = pathway.score_pathways(path_df, mut_df, data_type='mut')
            score_df.dropna(how='all', inplace=True)
            score_df = score_df.applymap(lambda x: '%.3f' %x)
            score_df.to_csv(results_dir + os.sep + mut +'_mutation_pathway_score.txt', index=True, index_label='Pathway', sep='\t', float_format='%.3f')
            print("Wrote " + can_type + " " + mut +'_mutation_pathway_score.txt %d x %d' %score_df.shape )



def get_pancancer_mutation_summary(agg_axis=0, refresh=False):
    """ Summarize mutation data across cancer types.
    Parameters
    -----------------
    agg_axis: 0 or 1 Aggregate across Samples or Genes
    refresh: True/False Recompute mutation frequencies from all cancer types
                or reload from previous computaion
    """
    outdir = '../data/pancancer'
    makedirs(outdir)
    if agg_axis:
        ns_outf = outdir + os.sep + 'non_silent_mutations_genes.txt'
        s_outf = outdir + os.sep + 'silent_mutations_genes.txt'
    else:
        ns_outf = outdir + os.sep + 'non_silent_mutations_samples.txt'
        s_outf = outdir + os.sep + 'silent_mutations_samples.txt'

    if os.path.exists(ns_outf) and not refresh:
        ns_summary = pd.read_table(ns_outf, sep='\t', index_col=0, header=0)
        s_summary = pd.read_table(s_outf, sep='\t', index_col=0, header=0)

        return (ns_summary, s_summary)


    ns_summary = pd.DataFrame()
    s_summary = pd.DataFrame()
    for can in CANCER_TYPES:
        f1 = '../data/processed/' + can + os.sep + can + '_nsilent_mutation.txt'
        f2 = '../data/processed/' + can + os.sep + can + '_silent_mutation.txt'
        if not os.path.exists(f1) or not os.path.exists(f2):
            #print("Skipping %s " %can)
            continue

        # Read non-silent mutations and compute mutation frequencies
        nsdf = pd.read_table(f1, header=0, index_col=0)
        nsilent = nsdf.sum(axis=agg_axis)

        sdf = pd.read_table(f2, header=0, index_col=0)
        silent = sdf.sum(axis=agg_axis)

        if ns_summary.empty:
            ns_summary[can] = nsilent
            s_summary[can] = silent
        else:
            ns_summary = ns_summary.join(pd.Series(nsilent, name=can), how='outer')
            s_summary = s_summary.join(pd.Series(silent, name=can), how='outer')


    # Save for future usage
    ns_summary.to_csv(ns_outf, sep='\t', header=True, index=True, float_format='%.3f')
    s_summary.to_csv(s_outf, sep='\t', header=True, index=True, float_format='%.3f')
    print("Wrote %s %d x %d matrix" %(os.path.basename(ns_outf), len(ns_summary.index), len(ns_summary.columns)))
    print("Wrote %s %d x %d matrix" %(os.path.basename(s_outf), len(s_summary.index), len(s_summary.columns)))

    return (ns_summary, s_summary)

if __name__ == '__main__':
    """ Gen pathway mutation correlation with non-silent mutation counts"""
    #generate_mutation_scores()
    get_pancancer_mutation_summary(0, True)
