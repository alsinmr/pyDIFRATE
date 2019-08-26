#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions for efficicient sampling of an MD trajectory.
(separated from Ct_fast due to circular import issues)

Created on Fri Aug 23 10:17:11 2019

@author: albertsmith
"""

import numpy as np

def trunc_t_axis(nt,n=100,nr=10,**kwargs):
    """
    Calculates a log-spaced sampling schedule for an MD time axis. Parameters are
    nt, the number of time points, n, which is the number of time points to 
    load in before the first time point is skipped, and finally nr is how many
    times to repeat that schedule in the trajectory (so for nr=10, 1/10 of the
    way from the beginning of the trajectory, the schedule will start to repeat, 
    and this will be repeated 10 times)
    
    """
    
    logdt0=np.log10(1.50000001)/n
    
    index=list()
    index.append(0)
    dt=0
    while index[-1]<nt:
        index.append(index[-1]+np.round(10**dt))
        dt+=logdt0
        
    index=np.array(index)

    index=np.repeat(index,nr,axis=0)+np.repeat([np.arange(0,nt,nt/nr)],index.size,axis=0).reshape([index.size*nr])
    
    index=index[index<nt]
    index=np.unique(index).astype('int')
    
    return index


def get_count(index):
    """
    Returns the number of averages for each time point in the sparsely sampled 
    correlation function
    """
    N=np.zeros(index[-1]+1)
    n=np.size(index)
   
    for k in range(n):
        N[index[k:]-index[k]]+=1
        
    return N