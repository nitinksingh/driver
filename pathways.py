
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 18 16:27:20 2015

@author: Nitin Singh

Clean pathways from MSigDB
"""
from __future__ import division
import sys
import os
from utils import *
import pandas as pd


def preprocess_pathway_data():
    """Load GSEA MSigDB Broad Pathway DB"""
    input_dir = '../data/pathways'
    filename = input_dir + os.sep + 'kegg_biocarta_pid.txt'

    pathways = {}
    with open(filename, 'r') as f:
        for line in f:
            p = line.strip().split('\t')
            pathways[p[0]] = p[2:]

    df = pd.DataFrame.from_dict(pathways, orient='index').transpose()

    return df
