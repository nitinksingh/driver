# -*- coding: utf-8 -*-
"""
Created on Sun Jan 18 16:27:17 2014

@author: Nitin Singh

Utility functions to read and load Firehose/GDAC pre-processed TCGA data
"""
from __future__ import division
import sys
import os
from utils import *
import pandas as pd
import tarfile


def read_matrix_format_data(filename, keycol=0, **read_table_args):
    """ This is a generic function that reads data from a file which is in matrix
    format. Returns pandas dataframe.
    """
    input_fpath = os.path.abspath(filename)
    # Read data into a data frame. First row is header, pass caller functions
    # arguments. This argument will typically tell which columns to select ie.
    # "usecols" and rows to skip "skiprows".
    if not read_table_args:
        read_table_args = {'header': 0}

    if 'header' not in read_table_args:
        read_table_args['header'] = 0

    try:
        df = pd.read_table(input_fpath, sep='\t', **read_table_args)
    except:
        print "Unable to read input file: ", filename
        error()

    # Drop if key column is missing i.e first column has NaN in the data frame
    df = df.dropna(subset = [df.columns[keycol]])
    
    
    return df


def save_df(df, filename, dest_dir='', prefix=''):
    # Save the desired dataset in matrix format   
    input_fpath = os.path.abspath(filename)
    if dest_dir:
        of_matrix = dest_dir + os.sep + prefix + os.path.basename(input_fpath) 
    else:
        of_matrix = make_filename(input_fpath, append='matrix')
    
    df.to_csv(of_matrix, index=False, sep='\t') 
    print("Wrote %s %dx %d matrix" %(os.path.basename(of_matrix), len(df.index), len(df.columns)-1))   

def preprocess_pathway_data():
    """Load GSEA MSigDB Broad Pathway DB"""
    input_dir = '/Users/nitin/research/driver/data/pathways'
    filename = input_dir + os.sep + 'kegg_biocarta_pid.txt'

    pathways = {}
    with open(filename, 'r') as f:
        for line in f:
            p = line.strip().split('\t')
            pathways[p[0]] = p[2:]

    df = pd.DataFrame.from_dict(pathways, orient='index').transpose()

    return df
    
def preprocess_cna_data(filename, dest_dir):
    """ Load CNA data, shorten long TCGA""" 
    df = read_matrix_format_data(filename)    
    # drop locus id and cytoband columns
    df.drop(df.columns[[1, 2]], axis=1, inplace=True)

    categorize_samples(df, filename, dest_dir)


def preprocess_rnaseq_data(filename, dest_dir):
    """ Fix Gene Id and Symbol mix-up column"""
    df = read_matrix_format_data(filename)
    gene_sym = [gid_sym.split('|')[0] for gid_sym in df['gene']]
    df['gene'] = gene_sym

    categorize_samples(df, filename, dest_dir)

def categorize_samples(df, filename, dest_dir):
    """ Check for duplicates, normal samples and then shorten """
    headers = df.columns[1:]
    samples = ['-'.join(x.split('-')[:3]) for x in headers] 

    sample_type = pd.Series([int(x.split('-')[3][:2]) for x in headers])
    summary = sample_type.value_counts().sort_index()
    print("type\tcount")
    pd_print_full(summary)

    primary_tumors = []; normals = []
    for s in headers:
        stype = int(s.split('-')[3][:2])
        if stype == 1:
            primary_tumors.append(s)

        if stype == 10:
            normals.append(s)

    if primary_tumors:
        p = ['-'.join(x.split('-')[:3]) for x in primary_tumors]
        if len(set(p)) != len(p):
            error("duplicate primary tumors?, total: %d, uniq: %d" %(len(p), len(set(p))))
        
        pt_df = df[[df.columns[0]] + primary_tumors]
        save_df(df, filename, dest_dir, 'PT-')


    if normals:
        p = ['-'.join(x.split('-')[:3]) for x in normals]
        if len(set(p)) != len(p):
            error("duplicate normals?, total: %d, uniq: %d" %(len(p), len(set(p))))
       
        pt_df = df[[df.columns[0]] + normals]
        save_df(df, filename, dest_dir, 'NB-')

    return

    all_samples = {}; non_tumor = {}; tumor = {}; duplicates = {}
    for x in headers:
        
        dict_append(all_samples, p, x)

        s = x.split('-')[3][:2]
        if int(s) not in range(10):
            dict_append(non_tumor, p, x) 
        else:
            dict_append(tumor, p, x)

    for k, v in tumor.iteritems():
        if len(v) > 1:
            print 'tumor', k, v
    for k, v in non_tumor.iteritems():
        print 'non_tumor', k, v
    print("tumor: %d, non-tumor: %d, total: %d" %(len(tumor), len(non_tumor), len(headers)))
    # Set operation keeps only unique objects
    if len(samples) != len(uniq_samples):
        print('Before %d, after: %d' %(len(samples), len(uniq_samples)))
        

    return samples

def preprocess_clinical_data(filename, dest_dir):
    df = read_matrix_format_data(filename, 0, header=15, skiprows=[0])
    save_df(df, filename, dest_dir)
    
    return df


def preprocess_gdac_data():
    """ Gather filenames for GDAC downloaded TCGA data """
    input_dir = '/Users/nitin/research/driver/data'
    output_dir = input_dir + os.sep + 'processed'
    CANCER_TYPES = ['LUAD', 'LUSC']
    GDAC_PREFIX = 'gdac.broadinstitute.org_'

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
        

    for c in CANCER_TYPES:
        can_dir = output_dir + os.sep + c
        if not os.path.exists(can_dir):
            os.mkdir(can_dir)

        # CNA 
        cna_tar = input_dir + os.sep + 'analyses__2014_10_17' + \
                    os.sep + c + os.sep + '20141017' + os.sep + GDAC_PREFIX + \
                    c + '-TP.CopyNumber_Gistic2.Level_4.2014101700.0.0.tar.gz'
        cna_file = 'all_thresholded.by_genes.txt'
                    
        
        print("%s CNA" %c)
        extracted_cna_file = tar_extract(cna_tar, cna_file, can_dir) 
        preprocess_cna_data(extracted_cna_file, can_dir)

        
                
        # RNA-Seq
        rnaseq_tar = input_dir + os.sep + 'stddata__2014_12_06' + os.sep + \
                    c + os.sep + '20141206' + os.sep + GDAC_PREFIX + c + \
                    '.mRNAseq_Preprocess.Level_3.2014120600.0.0.tar.gz'
        rnaseq_file = c + '.uncv2.mRNAseq_RSEM_normalized_log2.txt'

        print("%s RNA-Seq" %c)
        extracted_rnaseq_file = tar_extract(rnaseq_tar, rnaseq_file, can_dir) 
        preprocess_rnaseq_data(extracted_rnaseq_file, can_dir)

        continue 
        # Clinical
        clinical_tar = input_dir + os.sep + 'stddata__2014_12_06' + os.sep + \
                    c + os.sep + '20141206' + os.sep + GDAC_PREFIX + c + \
                    '.Merge_Clinical.Level_1.2014120600.0.0.tar.gz'
        clinical_file = c + '.clin.merged.txt'

        print("%s Clincal" %c)
        extracted_clinical = tar_extract(clinical_tar, clinical_file, can_dir) 
        df = preprocess_clinical_data(extracted_clinical, can_dir)
        


def tar_extract(tf, filename, dest_dir):
    extratced_file= tf.split('.tar.gz')[0] + os.sep + filename
    
    if os.path.exists(extratced_file):
        #print("The file %s is already extracted" %extratced_file)
        return extratced_file
        
    if os.path.exists(tf):
        tarfile.open(tf, 'r').extract(fullname, dest_dir)
    else:
        error('Tar %s not found' %(tf))
    
    return extratced_file


# Program entry point    
if __name__ == "__main__":
    df = preprocess_gdac_data()    
    #df = preprocess_pathway_data()
