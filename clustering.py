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
import pathway_analysis as pathway

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


# Adopted from to OxanaSachenkova/hclust-python at github
def combined_heatmap(cancer_type, num_pathways=5, num_genes=10, groups=None):
    
    # Load preprocessed mutation data
    data_dir = '../data/processed/' + cancer_type
    path_df = pathway.preprocess_pathway_data()
    mut_path = data_dir + os.sep + cancer_type + '_nsilent_mutation.txt'
    mut_df = pd.read_table(mut_path, index_col=0, header=0)

    # Load mutation data aggregated at the level of pathways
    (sig_nsdf, clinical_df) = get_data(cancer_type, 'mut', nsdf_norm_factor=1)
    clinical_df = clinical_df.loc[clinical_df.time > 0]
    common_samples = sig_nsdf.columns.intersection(clinical_df.index)
    sig_nsdf = sig_nsdf[common_samples]
    clinical_df = clinical_df.loc[common_samples] # redaundant as of now
    df = sig_nsdf.iloc[:num_pathways]
    labels = sig_nsdf.columns

    # Do heirarchical clustering on samples
    data_array = df.transpose().values
    if type(labels) == str:
        labels = map(str, range(1, data_array.shape[0]+1))
    data_dist = pdist(data_array) # computing the distance
    Y = linkage(data_dist, method='single')

    x0 = 0.2; w=0.7; y00 = 0.93; y01 = 0.8; y = 25
    # Compute and plot first dendrogram.
    fig = plt.figure(figsize=(20,y))
    ax2 = fig.add_axes([x0,y00,w,0.06])
    Z2 = dendrogram(Y,labels=labels, leaf_font_size=8)#, color_list=color_list.values)
    idx2 = Z2['leaves']
    ax2.set_yticks([])
   
    # Cluster Ids
    ax3 = fig.add_axes([x0,y01+0.1005,w,0.002])
    if type(groups) is type(None):
        clust_mat = 10*np.ones((1, df.shape[1]))
    else:
        clust_mat = np.reshape(groups.iloc[idx2].values, (1,len(groups)))

    im = ax3.matshow(clust_mat, aspect='auto', origin='upper', cmap='Oranges')
    ax3.set_yticks([])
    ax3.set_xticks([])

    #Compute and plot the heatmap
    axmatrix = fig.add_axes([x0,y01,w,0.1])
    D = df.loc[:,df.columns[idx2]].values
    im = axmatrix.matshow(D, aspect='auto', origin='upper', cmap=plt.cm.YlGnBu)
  
    axmatrix.set_yticklabels([''] + list(df.index), rotation=0)
    axmatrix.set_xticks([])

    
    # Plot colorbar.
    axcolor = fig.add_axes([0.91,y01,0.02,0.1])
    plt.colorbar(im, cax=axcolor)
    
    cmaps = ['Oranges', 'binary', 'Spectral']
    h = y01/df.shape[0] - df.shape[0]*0.01
    for i in range(df.shape[0]):
        y0 = y01-h*(i+1)
    
        # Plot individual mutation
        ax3 = fig.add_axes([x0,y0,w, h] )
        path_genes = path_df[df.index[i]].dropna()
        path_genes_df = mut_df.loc[path_genes, df.columns[idx2]].fillna(0)
        path_genes_df = path_genes_df.loc[path_genes_df.sum(axis=1) != 0, :]
        sorted_df = path_genes_df.sum(axis=1).order(ascending=False)
        path_genes_df = path_genes_df.loc[sorted_df.index[:num_genes]]
#         print D.shape, path_genes_df.shape, df.index[i]
        D2 = path_genes_df.values
        im2 = ax3.matshow(D2, aspect='auto', origin='upper', cmap=cmaps[i%len(cmaps)])
        ax3.set_yticks(range(len(path_genes_df.index)))
        ax3.set_yticklabels(list(path_genes_df.index), fontsize=9)
        ax3.set_xticks([])
        # Plot colorbar.
        axcolor2 = fig.add_axes([0.91,y0,0.02,h])
        plt.colorbar(im2, cax=axcolor2)
        
#     print path_genes_df.head()
    plt.show()

    return idx2
    
if __name__ == "__main__":
    (sig_nsdf, clinical_df) = get_data()
    X = sig_nsdf[:50].transpose().values
    hierarchical_clustering(X[:10,:5])
    #spectral_bicoclustering(X, 2)
