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

def fit_data(data,detect=None,bounds=True,ErrorAna=None,save_input=True,parallel=True,**kwargs):
    """
    Subsequent fitting is currently failing (I think), because we are later trying to 
    fit the detectors that result from the R2 exchange correction. Should have an 
    automatic mechanism to discard these in later fits.
    """
    
    
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
    out.sens=detect.copy()
    out.sens._disable()    #Clear the input sensitivities (restricts ability to further edit sens)
    
    "Delete the estimation of R2 due to exchange if included in the data here"
    if hasattr(data.sens,'detect_par') and data.sens.detect_par['R2_ex_corr'][0].lower()=='y':
        R=data.R.copy()[:,:-1]
        R_std=data.R_std.copy()[:,:-1]
        data=data.copy('shallow')    #We don't want to edit the original data object by deleting some of the R data
        "The shallow copy alone would still edit the original R data"
        "Replacing the matrices, however, should leave the orignal matrices untouched"
        data.R=R
        data.R_std=R_std
        
    
    nd=detect.r(bond=0).shape[1]    #Number of detectors
    out.R=np.zeros([nb,nd])
    
    "Set some defaults for error analysis"
    if ErrorAna is not None:
        ea=ErrorAna
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
    

#    "Set up parallel processing"
#    if 'parallel' in kwargs:
#        if kwargs.get('parallel')[0].lower()=='y':
#            para=True
#        else:
#            para=False
#    elif not(bounds):
#        para=False
#    else:
#        if nmc==0:
#            para=True
#        else:
#            para=True
#    
    if not(bounds):
        Y=list()
        for k in range(nb):
            if data.S2 is not None and not('subS2' in kwargs and kwargs.get('subS2').lower()[0]=='n'):
                R=(data.R[k,:]-data.S2[k]-detect.R0in(k))/data.R_std[k,:]
            else:
                R=(data.R[k,:]-detect.R0in(k))/data.R_std[k,:]
            r=detect.r(bond=k)
            r=r/np.repeat(np.transpose([data.R_std[k,:]]),r.shape[1],axis=1)
                        
            nstd=norm.ppf(1/2+conf/2)
            std=np.sqrt(np.sum(np.linalg.pinv(r)**2,axis=1))
            u=nstd*std
            l=nstd*std
            rho=np.dot(np.linalg.pinv(r),R)
            Y.append((rho,std,u,l))
        
    elif not(parallel):
        "Series processing (only on specific user request)"
        Y=list()
        for k in range(nb):
            rhoz=detect.rhoz(bond=k)
            UB=rhoz.max(axis=1)
            LB=rhoz.min(axis=1)
            r=detect.r(bond=k)
            
            if data.S2 is not None and not('subS2' in kwargs and kwargs.get('subS2').lower()[0]=='n'):
                R=(data.R[k,:]-data.S2[k]-detect.R0in(k))/data.R_std[k,:]
            else:
                R=(data.R[k,:]-detect.R0in(k))/data.R_std[k,:]
            r=r/np.repeat(np.transpose([data.R_std[k,:]]),r.shape[1],axis=1)
            X=(r,R,LB,UB,conf,nmc)
            Y.append(para_fit(X))
            print(k)
    else:
        "Here, we buildup up X with all the information required for each fit"
        "required: normalized data, normalized r, upper and lower bounds"
        X0=list()
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
            
            X0.append((r,R,LB,UB,conf,nmc))
        
        "Parallel processing (default)"
        nc=mp.cpu_count()
        if 'n_cores' in kwargs:
            nc=np.min([kwargs.get('n_cores'),nc])
            
        with mp.Pool(processes=nc) as pool:
            Y=pool.map(para_fit,X0)



    Rc=np.zeros(data.R.shape)

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
        Rc[k,:]=np.dot(detect.r(bond=k),out.R[k,:])+detect.R0in(k)

    if save_input:
        out.Rc=Rc
        
    out.sens.info.loc['stdev']=np.median(out.R_std,axis=0)
        
    if save_input:
        out.Rin=data.R
        out.Rin_std=data.R_std
    
        
    out.detect=dt(detect)
    
    out.ired=data.ired
    out.label=data.label
    

    out.chi=np.sum((data.R-Rc)**2/(data.R_std**2),axis=1)
    
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
       
    return rho,std,l,u

#%% Function to force a data object to be fully consistent with a positive dynamics distribution
def opt2dist(data,sens=None,parallel=True,return_dist=False,in_place=False,**kwargs):
    """
    Takes a distribution and sensitivity object (usually contained in the data
    object, but can be provided separately), and for each bond/residue, optimizes
    a distribution that approximately yields the set of detectors, while requiring
    that the distribution itself only contains positive values and has an integral
    of 1 (or 1-S2, if S2 is stored in data). Note that the distribution itself 
    is not a good reporter on dynamics; it is neither regularized or a stable
    description of dynamics. However, its use makes the detector responses more
    physically consistent
    
    opt_data=opt2dist(data,sens=None,para=True,return_dist=False,in_place=False)
    
    returns 0, 1, or 2 values, depending on the setting of return_dist and in_place
    
    """
    
    nb=data.R.shape[0]
    
    if data.S2 is None:
        S2=np.zeros(nb)
    else:
        S2=data.S2
        
    if sens is None:
        sens=data.sens

    "data required for optimization"
    X=[(R,R_std,sens._rho(bond=k),S2r) for k,(R,R_std,S2r) in enumerate(zip(data.R,data.R_std,S2))]        
    
    if parallel:
        nc=mp.cpu_count()
        if 'n_cores' in kwargs:
            nc=np.min([kwargs.get('n_cores'),nc])
            
        with mp.Pool(processes=nc) as pool:
            Y=pool.map(dist_opt,X)
    else:
        Y=[dist_opt(X0) for X0 in X]
    
    out=data if in_place else data.copy() #We'll edit out, which might be the same object as data
    
    dist=list()
    for k,y in enumerate(Y):
        out.R[k]=y[0]
        dist.append(y[1])
        
    "Output"
    if in_place and return_dist:
        return dist
    elif in_place:
        return
    elif return_dist:
        return (out,dist)
    else:
        return out
    
def dist_opt(X):
    """
    Optimizes a distribution that yields detector responses, R, where the 
    distribution is required to be positive, and have an integral of 1-S2
    
    Ropt,dist=dist_opt((R,R_std,rhoz,S2,dz))
    
    Note- intput is via tuple
    """
    
    R,R_std,rhoz,S2=X
    total=np.atleast_1d(1-S2)
    """Later, we may need to play with the weighting here- at the moment, we
    fit to having a sum of 1, but in fact it is not forced....it should be
    """
    
    ntc=rhoz.shape[1]
    rhoz=np.concatenate((rhoz/np.repeat(np.atleast_2d(R_std).T,ntc,axis=1),
        np.atleast_2d(np.ones(ntc))),axis=0)
    Rin=np.concatenate((R/R_std,total))
    
    dist=0
    while np.abs(np.sum(dist)-total)>1e-3:  #This is a check to see that the sum condition has been satisfied
        dist=lsq(rhoz,Rin,bounds=(0,1))['x']
        Rin[-1]=Rin[-1]*10
        rhoz[-1]=rhoz[-1]*10
    Ropt=np.dot(rhoz[:-1],dist)*R_std
    
    return Ropt,dist