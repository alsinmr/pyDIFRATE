#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 17:02:56 2019

@author: albertsmith

Functions for accelarated calculation of correlation functions

Some notes on these functions: Let's take an example. Say we have a trajectory,
sampled every 1 ps, out to 500 ns (500,000 pts). The first 1000 points each are
calculated from about 5e5 pairs of time points. However, we learn about 
approximately 3 orders of magnitude of dynamics from those first 1000 points. By
comparison, starting from 250 ns, the following 1000 pts of the correlation 
function are calculated from 2.5e5 time point pairs (similar accuracy). However,
there is virtually no new information between the time point at 250 ns and
251 ns. There is almost no decay, because all the correlation times comparable to
1 ns have already decayed. It stands to reason that we don't really need so 
many time points from later in the correlation function. In fact, it would make
sense to only calculate the correlation function on a log-spaced time axis. 

This is problematic, however, because we still need to load all time points to 
get the log-spaced time points. On the other hand, we could load log-spaced time
points from the trajectory, and calculate the correlation function for all 
possible time points available based on the spacing of the trajectory. Then,
long correlation times will still be common in the correlation function, but they
will not be accurately calculated. Hopefully, we can still successfully fit them
with detectors, and recover the information based on the number of time points 
instead of the accuracy
"""

import numpy as np
import multiprocessing as mp
import os
os.chdir('../data')
from data_class import data
os.chdir('../iRED')


def trunc_t_axis(nt,n=100,nr=10):
    """
    Calculates a log-spaced sampling schedule for an MD time axis. Parameters are
    nt, the number of time points, n, which is the number of time points to 
    load in before the first time point is skipped, and finally nr is how many
    times to repeat that schedule in the trajectory (so for nr=10, 1/10 of the
    way from the beginning of the trajectory, the scehdule will start to repeat, 
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
    index=np.sort(index).astype('int')
    
    return index
        

def get_trunc_vec(molecule,n=100,nr=10):
    """
    vec=get_trunc_vec(molecule,n=100,nr=10)
    

    
    """
    pass