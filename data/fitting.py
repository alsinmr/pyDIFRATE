#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  8 12:26:05 2019

@author: albertsmith
"""

import numpy as np
import data_class as dc
from scipy.optimize import lsq_linear as lsq
from scipy.stats import norm
import os
import multiprocessing as mp
os.chdir('../r_class')
from detectors import detect as dt
os.chdir('../data')

def fit_data(data,detect=None,**kwargs):
    if detect is None:
        if data.detect is None:
            print('A detect object must be provided in the input or as part of the data object')
            return
        else:
            detect=data.detect
    
    if detect.r(bond=0) is None:
        print('First optimize a set of detectors for analysis')
        return
    
    nb=data.R.shape[0]  #number of data points to fit (n_bonds)
    
    "Output object"
    out=dc.data()
    "The new sensitivities of the output data are the detectors used"
    out.sens=detect
    
    nd=detect.r(bond=0).shape[1]    #Number of detectors
    out.R=np.zeros([nb,nd])
    
    "Set some defaults for error analysis"
    if 'ErrorAna' in kwargs:
        ea=kwargs.get('ErrorAna')
        if ea.lower()[0:2]=='mc':
            if len(ea)>2:
                nmc=int(ea[2:])
            else:
                nmc=100
        else:
            nmc=0
    else:
        nmc=0
            
    if 'Conf' in kwargs:
        conf=kwargs.get('Conf')
    else:
        conf=0.68
    out.conf=conf
    
    
    if data.S2 is not None and not('subS2' in kwargs and kwargs.get('subS2').lower()[0]=='n'):
        print('Subtracting S2')
    
    "Here, we buildup up X with all the information required for each fit"
    "required: normalized data, normalized r, upper and lower bounds"
    X=list()
    for k in range(0,nb):
        rhoz=detect.rhoz(bond=k)
        UB=rhoz.max(axis=1)
        LB=rhoz.min(axis=1)
        r=detect.r(bond=k)
        
        if data.S2 is not None and not('subS2' in kwargs and kwargs.get('subS2').lower()[0]=='n'):
            R=(data.R[k,:]-data.S2[k]-detect.R0in(k))/data.R_std[k,:]
        else:
            R=(data.R[k,:]-detect.R0in(k))/data.R_std[k,:]
        r=r/np.repeat(np.transpose([data.R_std[k,:]]),r.shape[1],axis=1)
        
        X.append((r,R,LB,UB,conf,nmc))
    
    "Set up parallel processing"    
    if 'parallel' in kwargs and kwargs.get('parallel').lower()[0]=='y':
        nc=mp.cpu_count()
        if 'n_cores' in kwargs:
            nc=np.min([kwargs.get('n_cores'),nc])
        
        with mp.Pool(processes=nc) as pool:
            Y=pool.map(para_fit,X)
            
    else:
        "Else, run in a normal for loop"
        Y=list()
        for k in range(0,nb):
            Y.append(para_fit(X[k]))
            
    "Options to not save the input- possibly useful for MD data analysis"
    if 'save_input' in kwargs and kwargs.get('save_input').lower()[0]=='n':
        sv_in=False
    else:
        sv_in=True
       
    if sv_in:
        out.Rc=np.zeros(data.R.shape)

    nd=detect.r(bond=0).shape[1]
    out.R=np.zeros([nb,nd])
    out.R_std=np.zeros([nb,nd])
    out.R_l=np.zeros([nb,nd])
    out.R_u=np.zeros([nb,nd])        
    for k in range(0,nb):
        out.R[k,:]=Y[k][0]
        out.R_std[k,:]=Y[k][1]
        out.R_l[k,:]=Y[k][2]
        out.R_u[k,:]=Y[k][3]
        if sv_in:
            out.Rc[k,:]=np.dot(detect.r(bond=k),out.R[k,:])
        
        
    if sv_in:
        out.Rin=data.R
        out.Rin_std=data.R_std
        
    out.detect=dt(detect)
    
    out.ired=data.ired
    out.label=data.label
    
    return out

def para_fit(X):
    "Function to calculate results in parallel"
    Y=lsq(X[0],X[1],bounds=(X[2],X[3]))
    rho=Y['x']
    Rc=Y['fun']+X[1]
    
    if X[5]==0:
        std=np.sqrt(np.sum(np.linalg.pinv(X[0])**2,axis=1))
        nstd=norm.ppf(1/2+X[4]/2)
        u=nstd*std
        l=nstd*std
    else:
        Y1=list()
        nmc=max([X[5],np.ceil(2/X[4])])
        for k in range(0,X[5]):
            Y0=lsq(X[0],Rc+np.random.normal(size=X[1].shape))
            Y1.append(Y0['x'])
        std=np.std(Y1,axis=0)
        Y1sort=np.sort(Y1,axis=0)
        in_l=np.round(nmc*(1/2-X[4]/2))
        in_u=np.round(nmc*(1/2+X[4]/2))
        l=rho-Y1sort[int(in_l)]
        u=Y1sort[int(in_u)]-rho
       
    return rho,std,u,l