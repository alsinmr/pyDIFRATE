#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 10 15:58:49 2019

@author: albertsmith
"""

"""
We do basic analysis of correlation functions, without iRED analysis, to 
compare to iRED results
"""

import numpy as np
import multiprocessing as mp
import os
os.chdir('../data')
from data_class import data
os.chdir('../iRED')
from iRED_ana import get_vec
from iRED_ana import align_vec

#%% Create a data object from the Correlation function results
def Ct2data(molecule,**kwargs):
    
    
    vec=get_vec(molecule.sel1,molecule.sel2,**kwargs)
    
    if 'align' in kwargs and kwargs.get('align').lower()[0]=='y':
        vec=align_vec(vec)
    
    ct=Ct(vec,**kwargs)
    
    S2=S2calc(vec)
    
    Ctdata=data(molecule=molecule,Ct=ct,S2=S2)

    return Ctdata

#%% Calculate correlation functions
def Ct(vec,**kwargs):    
    if 'dt' in kwargs:
        dt=kwargs.get('dt')
        nt=vec.get('t').size
        t=np.arange(0,dt*nt,dt)
    else:
        t=vec.get('t')      
    
    nb=vec.get('X').shape[0]
 
    "Prepare the data needed for each correlation function"    
    v1=list()
    for k in range(0,nb):
        v1.append(np.array([vec.get('X')[k,:],vec.get('Y')[k,:],vec.get('Z')[k,:]]))
    
    "Run in series or in parallel"
    if 'parallel' in kwargs and kwargs.get('parallel').lower()[0]=='n':
        ct0=list()
        for k in range(0,nb):
            ct0.append(Ct_parfun(v1[k]))         
    else:             
        nc=mp.cpu_count()
        if'n_cores' in kwargs:
            nc=np.min([kwargs.get('n_cores'),nc])
            
        with mp.Pool(processes=nc) as pool:
            ct0=pool.map(Ct_parfun,v1)
    

    ct={'t':t,'Ct':np.array(ct0)}
    
    return ct
        
                
           
def Ct_parfun(v):
    nt=np.shape(v)[1]
    for m in range(0,nt):
        v0=np.repeat(np.transpose([v[:,m]]),nt-m,axis=1)
        if m==0:
            ct=(3*np.sum(v0*v[:,m:],axis=0)**2-1)/2
        else:
            ct[0:-m]+=(3*np.sum(v0*v[:,m:],axis=0)**2-1)/2
            
    ct=ct/np.arange(nt,0,-1)
            
    return ct

#%% Calculate the order parameter
def S2calc(vec):
    v=[vec.get('X'),vec.get('Y'),vec.get('Z')]
    S2=np.zeros(np.shape(vec.get('X'))[0])
    for k in v:
        for m in v:
            S2+=np.mean(k*m,axis=1)**2
    
    S2=3/2*S2-1/2
    
    return S2        