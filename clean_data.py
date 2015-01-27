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
import glob


def read_matrix_format_data(filename, keycol=0, **read_table_args):
    """ This is a generic function that reads data from a file which is in matrix
    format. Returns pandas dataframe. Common arguments: index_col, names, 
    skip_rows, na_values, usecols, delimiter
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

    
def preprocess_cna_data(filename, dest_dir):
    """ Load data, drop chromosome and cytoband columns and check if there are 
    samples other than primary tumor type.
    """
    df = read_matrix_format_data(filename)    
    # drop locus id and cytoband columns
    df.drop(df.columns[[1, 2]], axis=1, inplace=True)
    df = df.rename(columns = {'Gene Symbol':'Gene_Symbol'})
    categorize_samples(df, filename, dest_dir, prefix='CNA_')


def preprocess_rnaseq_data(filename, dest_dir):
    """ Load data fix Gene Id and Symbol mix-up column, drop genes with '?' 
    symbols and check if there are samples other than primary tumor type.
    """
    df = read_matrix_format_data(filename)
    gene_sym = [gid_sym.split('|')[0] for gid_sym in df['gene']]
    df['gene'] = gene_sym
    print "rna seq", df.shape
    df = df[df.gene != '?']
    df = df.rename(columns = {'gene':'Gene_Symbol'})
    
    categorize_samples(df, filename, dest_dir, prefix='')

def categorize_samples(df, filename, dest_dir, prefix=''):
    """ Check for duplicates, normal samples and then shorten TCGA id """
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

        if stype in range(10, 15):
            normals.append(s)

    if primary_tumors:
        short_id = ['-'.join(x.split('-')[:3]) for x in primary_tumors]
        if len(set(short_id)) != len(short_id):
            error("duplicate primary tumors?, total: %d, uniq: %d" %(len(p), len(set(p))))
        
        pt_df = df[[df.columns[0]] + primary_tumors]
        pt_df.columns = [pt_df.columns[0]] + short_id
        save_df(pt_df, filename, dest_dir, prefix+'TP_')


    if normals:
        short_id = ['-'.join(x.split('-')[:3]) for x in normals]
        if len(set(short_id)) != len(short_id):
            error("duplicate normals?, total: %d, uniq: %d" %(len(p), len(set(p))))
       
        nb_df = df[[df.columns[0]] + normals]
        nb_df.columns = [nb_df.columns[0]] + short_id
        save_df(nb_df, filename, dest_dir, prefix+'NB_')

    return


def preprocess_clinical_data(filename, dest_dir):
    df = read_matrix_format_data(filename, 0, header=15, skiprows=[0])
    save_df(df, filename, dest_dir)
    
    return df

def preprocess_mutation_data(mut_dir, dest_dir, prefix=''):
    mut_dir += os.sep + '*.maf.txt'
    mut_files = glob.glob(mut_dir)

    samples = [os.path.basename(f).split('.')[0] for f in mut_files]
    sample_type = pd.Series([int(x.split('-')[3][:2]) for x in samples])
    summary = sample_type.value_counts().sort_index()
    print("TSS Code \t Sample Count")
    for i in summary.index:
        print("%d \t\t %d\n" %(i, summary.loc[i]))

    p = ['-'.join(x.split('-')[:3]) for x in samples]
    if len(set(p)) != len(p):
            error('Duplicate samples?,  total: %d, uniq: %d' %(len(p), len(set(p))))

    if any([0 if x ==1 else 0  for x in sample_type]):
        error('Non primary tumor samples?')

    joint_dict = {}
    # Some info variables
    unknown_count = 0; both_mut_count = 0;
    for f in mut_files:
        df = read_matrix_format_data(f, 0, usecols = [0, 8])
        # Substitute silent mutations: 1, non-silent: -1 (deleterious) 
        sub_dict = dict.fromkeys(df['Variant_Classification'].unique(), -1)
        sub_dict.update({'Silent':1})
        df.replace(sub_dict.keys(), sub_dict.values(), inplace=True)
        # If there are multple same type of mutations in the same gene
        df.drop_duplicates(inplace=True)
        df.set_index('Hugo_Symbol', drop=True, inplace=True)
        if 'Unknown' in df.index:
            unknown_count += df.loc["Unknown"].shape[0]
            df.drop("Unknown", axis=0, inplace=True)
        # Keep non-silent if a gene has both silent and non-silent mutations
        both_mut_genes = df.index.get_duplicates()
        df.loc[both_mut_genes] == -1
        df.drop_duplicates(inplace=True)
        both_mut_count += len(both_mut_genes)
        
        x = os.path.basename(f).split('.')[0] 
        df.columns = ['-'.join(x.split('-')[:3])]
        joint_dict.update(df.to_dict())        
        joint_df = pd.DataFrame.from_dict(joint_dict)
        
    print("dropped 'Unknown' genes count: %d" %unknown_count)
    print("Total incidence of both silent and non-silent mutations occuring in the same gene: %d" %both_mut_count)
    print("Saving in: %s" %dest_dir)
    joint_df.insert(0, 'Gene_Symbol', joint_df.index)
    joint_df.fillna(0, inplace=True)    
    save_df(joint_df, 'mutation.txt', dest_dir, prefix+'_')
    
    return joint_df
    
def preprocess_gdac_data():
    """ Gather filenames for GDAC downloaded TCGA data """
    input_dir = '../data'
    output_dir = input_dir + os.sep + 'processed'
    CANCER_TYPES = ['LUAD', 'LUSC']
    GDAC_PREFIX = 'gdac.broadinstitute.org_'

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
        

    for c in CANCER_TYPES:
        can_dir = output_dir + os.sep + c
        if not os.path.exists(can_dir):
            os.mkdir(can_dir)

        print('\n'+'*'*50)
        print(' '*20 +  c   + ' '*20)
        print('*'*50)
        # Mutation
        mut_dir = input_dir + os.sep + 'stddata__2014_12_06' + os.sep + \
                    c + os.sep + '20141206' + os.sep + GDAC_PREFIX + c + \
                    '.Mutation_Packager_Calls.Level_3.2014120600.0.0'
        print("\n\t MUTATION")
        print("-"*40)
        df = preprocess_mutation_data(mut_dir, can_dir, c)
        
        continue
        # CNA 
        cna_tar =   input_dir + os.sep + 'analyses__2014_10_17' + \
                    os.sep + c + os.sep + '20141017' + os.sep + GDAC_PREFIX + \
                    c + '-TP.CopyNumber_Gistic2.Level_4.2014101700.0.0.tar.gz'
        cna_file = 'all_thresholded.by_genes.txt'
                    
        
        print("\n\t CNA")
        print("-"*40)

        extracted_cna_file = tar_extract(cna_tar, cna_file, can_dir) 
        preprocess_cna_data(extracted_cna_file, can_dir)

        
                
        # RNA-Seq
        rnaseq_tar = input_dir + os.sep + 'stddata__2014_12_06' + os.sep + \
                    c + os.sep + '20141206' + os.sep + GDAC_PREFIX + c + \
                    '.mRNAseq_Preprocess.Level_3.2014120600.0.0.tar.gz'
        rnaseq_file = c + '.uncv2.mRNAseq_RSEM_normalized_log2.txt'

        print("\n\t mRNA-Seq")
        print("-"*40)

        extracted_rnaseq_file = tar_extract(rnaseq_tar, rnaseq_file, can_dir) 
        preprocess_rnaseq_data(extracted_rnaseq_file, can_dir)

        continue 
        # Clinical
        clinical_tar = input_dir + os.sep + 'stddata__2014_12_06' + os.sep + \
                    c + os.sep + '20141206' + os.sep + GDAC_PREFIX + c + \
                    '.Merge_Clinical.Level_1.2014120600.0.0.tar.gz'
        clinical_file = c + '.clin.merged.txt'

        print("\n\t Clinical")
        print("-"*40)

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
