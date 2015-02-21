# -*- coding: utf-8 -*-
"""
General purpose utility functions
Created on Fri Jul  4 12:01:05 2014

@author: Nitin Singh
"""
from __future__ import print_function
import sys, os
import operator
import inspect
import pandas as pd

# Write warning message to stderr
def warn(*objs):
    print("WARNING in function " + inspect.stack()[1][3] + ": ", *objs, file=sys.stderr)
#    sys.stderr.flush()

# Write error message to stderr and QUIT
def error(*objs):
    print("ERROR in function "+ inspect.stack()[1][3] + ": ", *objs, file=sys.stderr)
    sys.exit(1)

# Write warning message to stderr
def info(*objs):
    print("INFO from function " + inspect.stack()[1][3] + ": ", *objs, file=sys.stdout)

def makedirs(dir_name):
    if not os.path.exists(dir_name):
        os.makedirs(os.path.abspath(dir_name))
    
# A helper function that returns a file path name appended with "append" string
def make_filename(orig_name, append='cleaned'):
    # Now save the cleaned data
    new_name = os.path.dirname(orig_name) + os.sep  + \
            '.'.join(os.path.basename(orig_name).split('.')[0:-1]) + \
            '_' + append + \
            '.txt'
            
    if os.path.isfile(new_name):
           print("New file ", new_name, " exists!")
           
    return new_name

# Sort dictionary by decreasing order of values    
def sort_dict(x, decreasing=True):
    sorted_x = sorted(x.iteritems(), key=operator.itemgetter(1), reverse=decreasing)
    
    return sorted_x

# Flush stdout and stderr
def flush():
    sys.stdout.flush()
    sys.stderr.flush()    

# Check if any value is missing in a pandas dataframe
def is_missing(df):
    s = df.isnull().any(axis=0)
    if any(s.values): 
        return True 
    else:
        return False 

def dict_append(d, k, v):
    """ Append to dict list values """
    if not isinstance(d, dict):
        error("Must pass a dict as first arg")

    if k in d:
        d[k].append(v)
    else:
        d[k] = [v]

def pd_print_full(x):
    pd.set_option('display.max_rows', len(x))
    print(x)
    pd.reset_option('display.max_rows')


def save_df(df, filename, dest_dir='', prefix=''):
    # Save the desired dataset in matrix format   
    input_fpath = os.path.abspath(filename)
    if not prefix:
        prefix = 'matrix_'

    if not dest_dir:
        dest_dir = os.path.dirname(input_fpath)

    of_matrix = dest_dir + os.sep + prefix + os.path.basename(input_fpath) 
    
    df.to_csv(of_matrix, index=False, sep='\t') 
    print("Wrote %s %d x %d matrix" %(os.path.basename(of_matrix), len(df.index)+1, len(df.columns)))   
