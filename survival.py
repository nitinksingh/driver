# -*- coding: utf-8 -*-
"""
Created on Sun Feb  8 13:18:50 2015

@author: nitin
"""
from __future__ import division
import os
import math
import numpy as np
import matplotlib.pyplot as plt
import pylab
from matplotlib import rc
import pandas as pd
import scipy.stats as stats


from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test



def survival_analysis(df, kmshow=False):
    """ Perform survival analysis based on groups in the population
    Parmeters
    ---------------
    df: The input dataframe must have time, sensor, groups columns
    
    Return
    ---------------
    KM Plot
    """
    # Time to death or time to last follow up, right sensor alive patients
    T = df['time']
    E = ~ df['right_sensor']
    G = df['groups']
    
    groups = list(G.unique())
    ix = (G == groups[0])
    
    # p-value of this group v/s others
    p_values = pd.Series(index=groups)
    spt = logrank_test(T[ix], T[~ix], E[ix], E[~ix], suppress_print=True )
    p_val = '%.2E' % spt[1]
    p_values[groups[0]] = p_val
    
    # Fit and plot Kaplan-Meir Plot
    if kmshow:
        kmf = KaplanMeierFitter()
        kmf.fit(T[ix], E[ix], label=str(groups[0]) + " (" + str(p_val) + ")")    
        ax = kmf.plot()

    for g in groups[1:]:
        ix = (G == g)
        # p-value of this group v/s others
        spt = logrank_test(T[ix], T[~ix], E[ix], E[~ix], suppress_print=True )
        p_val = '%.2E' % spt[1]
        p_values[g] = p_val
        
        if kmshow:
            kmf.fit(T[ix], E[ix], label=str(g) + " (" + str(p_val) + ")")
            kmf.plot(ax=ax)


    return p_values


