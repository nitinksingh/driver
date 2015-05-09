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



from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
# import sklearn.metrics as metrics
def survival_analysis(df, kmshow=True, title='', xlabel='', ylabel='', figsize=(12, 10)):
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

    p_val = '%.2E' % spt.p_value
    p_values[groups[0]] = p_val

    pylab.rcParams['figure.figsize'] = figsize
    # Fit and plot Kaplan-Meir Plot
    if kmshow:
        kmf = KaplanMeierFitter()
        kmf.fit(T[ix], E[ix], label='Cluster ' + str(groups[0]) + " (size: " + str(len(T[ix])) + ", p-val:" + str(p_val) + ")")
        ax = kmf.plot(ci_show=False, show_censors=True, flat=True)

    for g in groups[1:]:
        ix = (G == g)
        # p-value of this group v/s others
        spt = logrank_test(T[ix], T[~ix], E[ix], E[~ix], suppress_print=True )
        p_val = '%.2E' % spt.p_value
        p_values[g] = p_val

        if kmshow:
            kmf.fit(T[ix], E[ix], label='Cluster ' + str(g) + " (size: " + str(len(T[ix])) + ", p-val:" + str(p_val) + ")")
            kmf.plot(ax=ax, ci_show=False, show_censors=True, flat=True)

    plt.title('Survival Analysis: ' + title, fontsize=18)
    plt.xlabel('Time ' + xlabel)

    return p_values