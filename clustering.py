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
from utils import error

from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from sklearn.cluster import KMeans, SpectralClustering
from sklearn.cluster.bicluster import SpectralCoclustering, SpectralBiclustering

import sklearn.metrics as metrics

# Adopted from to OxanaSachenkova/hclust-python at github
def hierarchical_clustering(data_array, labels=''):
    if type(labels) == str:
        labels = map(str, range(1, data_array.shape[0]+1))
    data_dist = pdist(data_array) # computing the distance
    data_link = linkage(data_dist) # computing the linkage

    # Compute and plot first dendrogram.

    fig = plt.figure(figsize=(22,15))
    # x ywidth height
    ax1 = fig.add_axes([0.05,0.1,0.2,0.6])
    Y = linkage(data_dist, method='single')
    Z1 = dendrogram(Y, orientation='right',labels=labels, leaf_font_size=1) # adding/removing the axes
    ax1.set_xticks([])

    # Compute and plot second dendrogram.
    ax2 = fig.add_axes([0.3,0.71,0.6,0.2])
    Z2 = dendrogram(Y)
    ax2.set_xticks([])
    ax2.set_yticks([])

    #Compute and plot the heatmap
    axmatrix = fig.add_axes([0.3,0.1,0.6,0.6])
    idx1 = Z1['leaves']
    idx2 = Z2['leaves']
    D = squareform(data_dist)
    D = D[idx1,:]
    D = D[:,idx2]
    im = axmatrix.matshow(D, aspect='auto', origin='lower', cmap=plt.cm.YlGnBu)
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])

    # Plot colorbar.
    axcolor = fig.add_axes([0.91,0.1,0.02,0.6])
    plt.colorbar(im, cax=axcolor)

    plt.show()

    # Gather cluster member count
    membership = pd.Series(fcluster(Y, 10, 'maxclust'), index=labels, name='groups')
    return membership




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

def get_data(can='LUAD', data_type='mut', nsdf_norm_factor=1):
    # Load clinical file drop patients with no followup or death time available
    clncl_fpath = '../data/processed/' + can + os.sep + can + '_clinical.csv'
    clinical_df = pd.read_table(clncl_fpath, index_col=0, header=0, sep='\t').transpose()

    death_followup = clinical_df[['patient.days_to_death', 'patient.days_to_last_followup']].dropna(axis=1, how='all')
    time = death_followup.fillna(0).astype(float).sum(axis=1)
    clinical_df = clinical_df.loc[death_followup.index]
    clinical_df.insert(0, 'time', time)
    clinical_df.insert(1, 'right_sensor', clinical_df['patient.vital_status'] != 'dead')

    # Now load pathway scores
    if data_type == 'mut':
        opt = 'silent'
        input_fpath = '../results/' + can + os.sep + opt +'_mutation_pathway_score.txt'
        sdf = pd.read_table(input_fpath, sep='\t', header=0, index_col=0)
        opt = 'nsilent'
        input_fpath = '../results/' + can + os.sep + opt +'_mutation_pathway_score.txt'
        nsdf = pd.read_table(input_fpath, sep='\t', header=0, index_col=0)
    elif data_type == 'rna':
        opt = 'NB'
        input_fpath = '../results/' + can + os.sep + opt +'_rnaseq_pathway_score.txt'
        sdf = pd.read_table(input_fpath, sep='\t', header=0, index_col=0)
        opt = 'PT'
        input_fpath = '../results/' + can + os.sep + opt +'_rnaseq_pathway_score.txt'
        nsdf = pd.read_table(input_fpath, sep='\t', header=0, index_col=0)
    elif data_type == 'rna_centroid':
        # Sort the pathways per mutation data significant
        mut_df, clinical_df = get_data(can, data_type='mut')
        input_fpath = '../results/' + can + os.sep + 'PT_centroid_rnaseq_pathway_score.txt'
        pt_df = pd.read_table(input_fpath, sep='\t', header=0, index_col=0)

        try:
            sorted_pt_df = pt_df.loc[mut_df.index]
        except KeyError:
            sorted_mut_path = filter(lambda x: x in pt_df.index, mut_df.index)
            sorted_pt_df = pt_df.loc[sorted_mut_path]
        return (sorted_pt_df, clinical_df)
    else:
        error('Unknown data_type')

    nsdf = nsdf/nsdf_norm_factor
    common = nsdf.index.intersection(sdf.index)
    res = pd.Series(index=common)
    for p in common:
        res.loc[p] = -1*math.log10(stats.ttest_ind(sdf.loc[p], nsdf.loc[p], equal_var=False)[1])

    res.sort(ascending=False)

    sig_nsdf = nsdf.loc[res.index].astype(float)

    return (sig_nsdf, clinical_df)

if __name__ == "__main__":
    (sig_nsdf, clinical_df) = get_data()
    X = sig_nsdf[:50].transpose().values
    hierarchical_clustering(X[:10,:5])
    #spectral_bicoclustering(X, 2)
