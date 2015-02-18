# -*- coding: utf-8 -*-
"""
Created on Sun Feb  8 13:21:37 2015

@author: Nitin Singh
"""
from __future__ import division
import os
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as stats


from sklearn.cluster import KMeans, SpectralClustering
from sklearn.cluster.bicluster import SpectralCoclustering, SpectralBiclustering

import sklearn.metrics as metrics

def spectral_clustering(X, n_clusters, affinity='linear'):
    """ Spectral clustering of columns/samples based on rows/features.
    Parameters
    -------------
    X: A (m,n) shape matrix where m: samples, n: features
    
    affinity: Option to similarity Kernel. Valid values are
        [‘rbf’, ‘sigmoid’, ‘polynomial’, ‘poly’, ‘linear’, ‘cosine’]
        
    Return
    -------------
    cluster: m dim vector with cluster labels/ids for each sample.
    """

    sc = SpectralClustering(n_clusters=n_clusters, affinity=affinity)#, assign_labels='discretize')
    sc.fit_predict(X)
    
    return sc.labels_

def KMeans_clustering(X, n_clusters):
    """ Spectral clustering of columns/samples based on rows/features.
    Parameters
    -------------
    X: A (m,n) shape matrix where m: samples, n: features
    
    Return
    -------------
    cluster: m dim vector with cluster labels/ids for each sample.
    """
    km = KMeans(n_clusters=n_clusters, n_jobs=-1)
    km.fit_predict(X)
    
    return km.labels_

def spectral_bicoclustering(X, n_clusters):
    
    data = X
    model = SpectralCoclustering(n_clusters=n_clusters, random_state=0)
    model.fit(data)
    
    fit_data = data[np.argsort(model.row_labels_)]
    fit_data = fit_data[:, np.argsort(model.column_labels_)]

    plt.matshow(fit_data, cmap=plt.cm.Blues)
    plt.title("After biclustering; rearranged to show biclusters")

    plt.show() 

def get_data():
    can = 'LUAD'
    
    # Load clinical file drop patients with no followup or death time available
    clncl_fpath = '../data/processed/' + can + os.sep + can + '_clinical.csv'
    clinical_df = pd.read_table(clncl_fpath, index_col=0, header=0, sep='\t').transpose()
    
    death_followup = clinical_df[['patient.days_to_death', 'patient.days_to_last_followup']].dropna(axis=1, how='all')
    time = death_followup.fillna(0).astype(float).sum(axis=1)
    clinical_df = clinical_df.loc[death_followup.index]
    clinical_df.insert(0, 'time', time)
    clinical_df.insert(1, 'right_sensor', clinical_df['patient.vital_status'] != 'dead')
    
    # Now load pathway scores
    opt = 'silent'
    input_fpath = '../results/' + can + os.sep + opt +'_mutation_pathway_score.txt'
    sdf = pd.read_table(input_fpath, sep='\t', header=0, index_col=0)
    opt = 'nsilent'
    input_fpath = '../results/' + can + os.sep + opt +'_mutation_pathway_score.txt'
    nsdf = pd.read_table(input_fpath, sep='\t', header=0, index_col=0)
    
    res = pd.Series(index=sdf.index)
    for p in sdf.index:
        res.loc[p] = -1*math.log10(stats.ttest_ind(sdf.loc[p], nsdf.loc[p], equal_var=False)[1])
        
    res.sort(ascending=False)
    
    sig_nsdf = nsdf.loc[res.index].astype(float)
    
    return (sig_nsdf, clinical_df)

if __name__ == "__main__":
    (sig_nsdf, clinical_df) = get_data()
    X = sig_nsdf[:50].transpose().values
    
    spectral_bicoclustering(X, 2)