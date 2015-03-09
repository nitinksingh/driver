# -*- coding: utf-8 -*-
"""
Created on Sun Jan 18 16:27:17 2014

@author: Nitin Singh

Utility functions to read and load Firehose/GDAC pre-processed TCGA data
"""
from __future__ import division
import sys
import os
from utils import error, save_df
import pandas as pd
import tarfile
import glob


COHORTS = "ACC BLCA BRCA CESC CHOL COAD COADREAD DLBC ESCA FPPP GBM GBMLGG HNSC KICH KIRC KIRP LAML LGG LIHC LUAD LUSC MESO OV PAAD PRAD READ SARC SKCM STAD TGCT THCA THYM UCEC UCS UVM"
CANCER_TYPES = COHORTS.split()
TSS_CODE = dict.fromkeys(CANCER_TYPES, 1)
TSS_CODE.update({'LAML': 3, 'SKCM': 6})
TSS_SYM = {1:'TP', 3:'TB', 6:'TM'}
GDAC_PREFIX = 'gdac.broadinstitute.org_'


def _summarize_sample_types(cols):
    sample_type = [int(c.split('-')[3][:2]) for c in cols]
    summary = pd.Series(sample_type)

    summary = summary.value_counts().sort_index()

    print("TSS Code \t Sample Count")
    for i in summary.index:
        print("%d \t\t %d\n" %(i, summary.loc[i]))



def preprocess_cna_data(filename, dest_dir, can):
    """ Load data, drop chromosome and cytoband columns and check if there are
    samples other than primary tumor type.
    """

    try:
        df = pd.read_table(filename, header=0)
    except IOError:
        print('****** Unable to open {} ******'.format(os.path.basename(filename)))
        return

    print "Raw CNA shape", df.shape
    # drop locus id and cytoband columns
    df.drop(df.columns[[1, 2]], axis=1, inplace=True)
    df = df.rename(columns = {'Gene Symbol':'Gene_Symbol'})
    pt_df, nb_df = categorize_samples(df, can)

    if type(pt_df) != int:
        save_df(pt_df, 'cna.txt', dest_dir, can+'_PT_')
    if type(nb_df) != int:
        save_df(nb_df, 'cna.txt', dest_dir, can+'_NB_')

def preprocess_rnaseq_data(filename, dest_dir, can):
    """ Load data fix Gene Id and Symbol mix-up column, drop genes with '?'
    symbols and check if there are samples other than primary tumor type.
    """

    try:
        df = pd.read_table(filename, header=0)
    except IOError:
        print('****** Unable to open {} ******'.format(os.path.basename(filename)))
        return

    df.columns = ['Gene_Symbol'] + list(df.columns[1:])

    gene_sym = [gid_sym.split('|')[0] for gid_sym in df['Gene_Symbol']]
    df['Gene_Symbol'] = gene_sym
    print "Raw rna seq shape", df.shape
    df = df[df['Gene_Symbol'] != '?']


    pt_df, nb_df = categorize_samples(df, can)
    if type(pt_df) != int:
        save_df(pt_df, 'rnaseq.txt', dest_dir, can+'_PT_')
    if type(nb_df) != int:
        save_df(nb_df, 'rnaseq.txt', dest_dir, can+'_NB_')


def categorize_samples(df, can):
    """ Check for duplicates, normal samples and then shorten TCGA id """
    headers = df.columns[1:]
    _summarize_sample_types(headers)
    primary_tumors = []; normals = []
    for s in headers:
        stype = int(s.split('-')[3][:2])
        if stype == TSS_CODE[can]:
            primary_tumors.append(s)

        if stype in range(10, 15):
            normals.append(s)

    pt_df = 0; nb_df = 0;
    if primary_tumors:
        short_id = ['-'.join(x.split('-')[1:3]) for x in primary_tumors]
        if len(set(short_id)) != len(short_id):
            error("duplicate primary tumors?, total: %d, uniq: %d" %(len(p), len(set(p))))

        pt_df = df[[df.columns[0]] + primary_tumors]
        pt_df.columns = [pt_df.columns[0]] + short_id



    if normals:
        short_id = ['-'.join(x.split('-')[1:3]) for x in normals]
        if len(set(short_id)) != len(short_id):
            error("duplicate normals?, total: %d, uniq: %d" %(len(p), len(set(p))))

        nb_df = df[[df.columns[0]] + normals]
        nb_df.columns = [nb_df.columns[0]] + short_id

    return (pt_df, nb_df)

total_samples = 0
def preprocess_clinical_data(filename, dest_dir):
    try:
        df = pd.read_table(filename, index_col=0, header=0)
    except IOError:
        print('****** Unable to open {} ******'.format(os.path.basename(filename)))
        return

    df.columns = df.loc['patient.bcr_patient_barcode']
    df.columns = ['-'.join(x.upper().split('-')[1:]) for x in df.columns]
    df.dropna(axis=(0, 1), how='all', inplace=True)
    df.insert(0, 'patient_info', df.index)

    params = ['patient.days_to_death','patient.days_to_last_followup']
    params += ['patient.days_to_birth', 'patient.vital_status', 'patient.gender']
    params += ['patient.primary_pathology.histological_type']
#    params += ['patient.prior_dx', 'patient.race', 'patient.radiation_therapy']
#    params += ['patient.primary_therapy_outcome_success']

    # Take the params that are available
    params = df.index.intersection(params)
    df = df.loc[params]

    # Check for auxilary information
    can = os.path.basename(dest_dir)
    if can in ['CESC', 'HNSC']:
        aux_file = os.path.dirname(filename) + os.sep + can + '.merged_only_auxiliary_clin_format.txt'
        aux_df = pd.read_table(aux_file, index_col=0, header=0)
        aux_s = aux_df.loc['patient.hpv_test_results.hpv_test_result.hpv_status']
        aux_s.name = 'HPV'
        aux_s.index = [ '-'.join(x.split('-')[1:]).upper() for x in aux_df.loc['patient.bcr_patient_barcode']]
        aux_s['patient_info'] = 'Auxiliary'
        df = df.append(aux_s)

    if can in ['COAD', 'READ', 'COADREAD', 'STAD', 'UCEC', 'UCS']:
        aux_file = os.path.dirname(filename) + os.sep + can + '.merged_only_auxiliary_clin_format.txt'
        aux_df = pd.read_table(aux_file, index_col=0, header=0)
        aux_s = aux_df.loc['patient.microsatellite_instability_test_results.microsatellite_instability_test_result.mononucleotide_and_dinucleotide_marker_panel_analysis_status']
        aux_s.name = 'MSI'
        aux_s.index = [ '-'.join(x.split('-')[1:]).upper() for x in aux_df.loc['patient.bcr_patient_barcode']]
        aux_s['patient_info'] = 'Auxiliary'
        df = df.append(aux_s)

    clinical_of = dest_dir + os.sep + os.path.basename(dest_dir) + '_clinical.csv'
    df.to_csv(clinical_of, index=False, sep='\t')
    print("Wrote %s %d x %d matrix" %(os.path.basename(clinical_of), len(df.index), len(df.columns)))

    global total_samples
    total_samples += len(df.columns)
    return df



def _group_mutation_count(xdf, can):
    """ Helper function to process mutation data"""
    grouped = xdf.groupby('patient')
    ret_df = pd.DataFrame()

    for p, i in grouped.groups.iteritems():
        x = xdf.ix[i]
        y = pd.DataFrame(x['Hugo_Symbol'].value_counts(), columns=[p])
        ret_df = ret_df.join(y, how='outer')

    ret_df.fillna(0, inplace=True)

    # If there are any sample that is of non-primary type, raise error.
    # Otherwise, shortent the patient ids in the column
    sample_type = [int(c.split('-')[3][:2]) for c in ret_df.columns]
    summary = pd.Series(sample_type)
    if summary[summary != TSS_CODE[can]].any():
        _summarize_sample_types(ret_df.columns)
        error('Non primary type sample found')
    else:
       sample_ids = ['-'.join(c.split('-')[1:3]) for c in ret_df.columns]
       ret_df.columns = sample_ids

    ret_df.insert(0, 'Hugo_Symbol', ret_df.index)
    ret_df.reset_index(drop=True, inplace=True)

    return ret_df

def process_mutation_data_mutsig(filename, dest_dir, can):
    try:
        df = pd.read_table(filename, header=0)
    except IOError:
        print('****** Unable to open {} ******'.format(os.path.basename(filename)))
        return

    df = df[[u'Hugo_Symbol', u'Variant_Classification', u'Variant_Type', u'Mutation_Status', u'patient']]
    sdf = df[(df.Mutation_Status == 'Somatic') & (df.Variant_Classification =='Silent')][[u'Hugo_Symbol', u'patient']]
    nsdf = df[(df.Mutation_Status == 'Somatic') & (df.Variant_Classification !='Silent')] [[u'Hugo_Symbol',u'patient']]

    nsdf = _group_mutation_count(nsdf, can)
    save_df(nsdf, 'mutation.txt', dest_dir, can+'_nsilent_')

    sdf = _group_mutation_count(sdf, can)
    save_df(sdf, 'mutation.txt', dest_dir, can+'_silent_')

def preprocess_gdac_data():
    """ Construct filenames for GDAC downloaded TCGA data and call preprocessing
    functions for Mutation/CNA/mRNASeq/Clinical data respectively """
    input_dir = '../data'
    output_dir = input_dir + os.sep + 'processed'


    if not os.path.exists(output_dir):
        os.mkdir(output_dir)


    for c in CANCER_TYPES:
        can_dir = output_dir + os.sep + c
        if not os.path.exists(can_dir):
            os.mkdir(can_dir)

        print('\n'+'*'*50)
        print(' '*20 +  c   + ' '*20)
        print('*'*50)

        # Clinical
        print("-"*40)
        print("\n\t Clinical")
        print("-"*40)

        clinical_dir = input_dir + os.sep + 'stddata__2014_12_06' + os.sep + \
                    c + os.sep + '20141206' + os.sep + GDAC_PREFIX + c + \
                    '.Merge_Clinical.Level_1.2014120600.0.0'
        clinical_file = clinical_dir + os.sep + c + '.merged_only_clinical_clin_format.txt'


        if os.path.exists(clinical_dir + '.tar.gz'):
	        tar_extract(clinical_dir + '.tar.gz', os.path.dirname(clinical_dir))
	        preprocess_clinical_data(clinical_file, can_dir)
        else:
            print('Tar download is missing. Skipping %s' %c)

        continue
        # RNA-Seq
        rnaseq_dir = input_dir + os.sep + 'stddata__2014_12_06' + os.sep + \
                    c + os.sep + '20141206' + os.sep + GDAC_PREFIX + c + \
                    '.mRNAseq_Preprocess.Level_3.2014120600.0.0'
        rnaseq_file = rnaseq_dir + os.sep + c + '.uncv2.mRNAseq_RSEM_all.txt'
        print("-"*40)
        print("\n\t mRNA-Seq")
        print("-"*40)

        if os.path.exists(rnaseq_dir + '.tar.gz'):
            tar_extract(rnaseq_dir +'.tar.gz', os.path.dirname(rnaseq_dir))
            preprocess_rnaseq_data(rnaseq_file, can_dir, c)
        else:
            print('Tar {} does not exist. Skipping..'.format(os.path.basename(rnaseq_dir + '.tar.gz')))


        # CNA
        cna_dir =   input_dir + os.sep + 'analyses__2014_10_17' + \
                    os.sep + c + os.sep + '20141017' + os.sep + GDAC_PREFIX + \
                    c + '-TP.CopyNumber_Gistic2.Level_4.2014101700.0.0'
        cna_file = cna_dir + os.sep + 'all_thresholded.by_genes.txt'

        print("-"*40)
        print("\n\t CNA")
        print("-"*40)

        if os.path.exists(cna_dir + '.tar.gz'):
            tar_extract(cna_dir + '.tar.gz', os.path.dirname(cna_dir))
            preprocess_cna_data(cna_file, can_dir, c)
        else:
            print('Tar {} does not exist. Skipping..'.format(os.path.basename(cna_dir + '.tar.gz')))


        # Mutation
        print("-"*40)
        print("\n\t MUTATION")
        print("-"*40)

        mut_dir = input_dir + os.sep + 'analyses__2014_10_17' + \
                 os.sep + c + os.sep + '20141017' + os.sep + \
                 'MutSigNozzleReportCV' + os.sep

        mut_file = mut_dir  + c + '-' + TSS_SYM[TSS_CODE[c]] +'.final_analysis_set.maf'
        if os.path.exists(mut_file):
            process_mutation_data_mutsig(mut_file, can_dir, c)
        else:
            print('{} does not exist. Skipping..'.format(os.path.basename(mut_file)))





def tar_extract(tf, dest_dir):
    """ Extract a tar file. """
    if os.path.exists(tf):
        tarfile.open(tf, 'r').extractall(dest_dir)
    else:
        error('Tar %s not found' %(tf))


# Program entry point
if __name__ == "__main__":
    df = preprocess_gdac_data()
    #print total_samples
