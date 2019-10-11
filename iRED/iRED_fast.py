#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 10:06:35 2019

@author: albertsmith
"""

import numpy as np
import multiprocessing as mp
import os
os.chdir('../data')
from data_class import data
os.chdir('../iRED')
from MDAnalysis.analysis.align import rotation_matrix
from psutil import virtual_memory
from fast_funs import trunc_t_axis,S2calc,Ct,get_trunc_vec,get_count
from par_iRED import par_class as ipc
from time import time

#%% Run the full iRED analysis
def iRED_full(mol,rank=2,n=100,nr=10,**kwargs):
    """
    Runs the full iRED analysis for a given selection (or set of vec_special functions)
    Arguments are the rank (0 or 1), the sampling (n,nr), whether to align the 
    vectors (align_iRED='y'/'n', and refVecs, which may be a dict containing 
    a vector, created by DIFRATE, a tuple of strings selecting two sets of atoms
    defining bonds), or simply 'y', which will default to using the N-CA bonds in
    a protein.
    
    ired=iRED_full(mol,rank=2,n=100,nr=10,align_iRED='n',refVecs='n',**kwargs)
    
    """

    if 'nt' in kwargs:
        nt=np.min([mol.mda_object.trajectory.n_frames,kwargs.get('nt')])
    else:
        nt=mol.mda_object.trajectory.n_frames
    index=trunc_t_axis(nt,n,nr)
    vec=get_trunc_vec(mol,index,**kwargs)
    
    if 'align_iRED' in kwargs and kwargs.get('align_iRED').lower()[0]=='y':
        if 'refVecs' in kwargs:
            vec0=kwargs['refVecs']
            if isinstance(vec0,dict):
                pass
            elif len(vec0)==2 and isinstance(vec0[0],str) and isinstance(vec0[1],str):
                mol1=mol.copy()
                mol1.select_atoms(sel1=vec0[0],sel2=vec0[1])
                vec0=get_trunc_vec(mol1,index)
            elif isinstance(vec0,str) and vec0.lower()[0]=='y':
                s1='protein and name CA and around 1.6 N'
                s2='protein and name N and around 1.6 CA'
                mol1=mol.copy()
                mol1.select_atoms(sel1=s1,sel2=s2)
                vec0=get_trunc_vec(mol1,index)
            else:
                print('Warning: refVecs entry not valid, using input vectors as reference (without aligning)')
                vec0=vec
        else:
            vec0=vec
        vec=align_mean(vec)
        
        

        n_added_vecs=vec0.get('X').shape[1]
        for k in ['X','Y','Z']:
            vec[k]=np.concatenate((vec.get(k),vec0.get(k)),axis=1)        
        aligned=True
    else:
        aligned=False
        n_added_vecs=0
    
 
    M=Mmat(vec,rank)
    Yl=Ylm(vec,rank)
    aqt=Aqt(Yl,M)

    "Parallel calculation of correlation functions"
    ct=Cqt(aqt)
    ctinf=CtInf(aqt)
    dct=DelCt(ct,ctinf)
    ired={'rank':rank,'M':M['M'],'lambda':M['lambda'],'m':M['m'],'t':ct['t'],\
          'N':ct['N'],'index':ct['index'],'DelCt':dct['DelCt'].T,'CtInf':ctinf,\
          'Aligned':aligned,'n_added_vecs':n_added_vecs}
            
    return ired

#%% Generate a data object with iRED results
def iRED2data(molecule,rank=2,**kwargs):
    """Input a molecule object with selections already made, to get a full iRED 
    analysis, moved into a data object
    """
    
    
    ired=iRED_full(molecule,**kwargs)
    
    Ctdata=data(iRED=ired,molecule=molecule)
    Ctdata.sens.molecule=molecule
#    Ctdata.sens.molecule.set_selection()
    Ctdata.detect.molecule=Ctdata.sens.molecule
    
    return Ctdata
#%% Calculate the iRED M matrix
def Mmat(vec,rank=2):
    """Calculates the iRED M-matrix, yielding correlation of vectors at time t=0
    M = Mmat(vec,rank=2)
    M is returned as dictionary object, including the matrix itself, and also
    the 
    """
    
    X=vec['X'].T
    Y=vec['Y'].T
    Z=vec['Z'].T
    
    nb=X.shape[0]
    
    M=np.eye(nb)
    
    for k in range(0,nb-1):
        "These are the x,y,z positions for one bond"
        x0=np.repeat([X[k,:]],nb-k-1,axis=0)
        y0=np.repeat([Y[k,:]],nb-k-1,axis=0)
        z0=np.repeat([Z[k,:]],nb-k-1,axis=0)
        
        "We correlate those positions with all bonds having a larger index (symmetry of matrix allows this)"
        dot=x0*X[k+1:,:]+y0*Y[k+1:,:]+z0*Z[k+1:,:]
        
        if rank==1:
            val=np.mean(dot,axis=1)
        elif rank==2:
            val=np.mean((3*dot**2-1)/2,axis=1)
            
        M[k,k+1:]=val
        M[k+1:,k]=val
        
    Lambda,m=np.linalg.eigh(M)
    return {'M':M,'lambda':Lambda,'m':m,'rank':rank}

def Mlagged(vec,lag,rank=2):
    """Calculates the iRED M-matrix, with a lag time, which is provided by an 
    index or range of indices (corresponding to the separation in time points)
    M = Mlagged(vec,rank=2,lag)
    
    lag=10
    or 
    lag=[10,20]
    
    The first instance calculates M using time points separated by exactly the
    lag index. The second takes all time points separated by the first argument, 
    up to one less the last argument (here, separated by 10 up to 19)
    
    """
    
    X=vec['X'].T
    Y=vec['Y'].T
    Z=vec['Z'].T
    
    index0=vec['index']
    
    if np.size(lag)==1:
        lag=np.atleast_1d(lag)
    elif np.size(lag)==2:
        lag=np.arange(lag[0],lag[1])
    
    "Calculate indices for pairing time points separated within the range given in lag"
    index1=np.zeros(0,dtype=int)
    index2=np.zeros(0,dtype=int)
    for k in lag:
        i=np.isin(index0+k,index0)
        j=np.isin(index0,index0+k)
        index1=np.concatenate((index1,np.where(i)[0]))
        index2=np.concatenate((index2,np.where(j)[0]))
    
    nb=X.shape[0]
    M=np.eye(nb)
    
    for k in range(0,nb):
        "We correlate all times that have a second time within the lag range"
        x0=np.repeat([X[k,index1]],nb,axis=0)
        y0=np.repeat([Y[k,index1]],nb,axis=0)
        z0=np.repeat([Z[k,index1]],nb,axis=0)

        dot=x0*X[:,index2]+y0*Y[:,index2]+z0*Z[:,index2]
        
        if rank==1:
            val=np.mean(dot,axis=1)
        elif rank==2:
            val=np.mean((3*dot**2-1)/2,axis=1)
            
        M[k,:]=val
        
    return M
#%% Estimates cross-correlation of the eigenvectors of the M matrix
def Mrange(vec,rank,i0,i1):
    """Estimates the Mmatrix for frames offset by a minimum distance of i0 and
    a maximum distance of i1-1. All M-matrices are simply added together
    M=Mrange(vec,rank,i0,i1)
    """
    pass

#%% Calculates the spherical tensor components for the individual bonds
def Ylm(vec,rank=2):
    """
    Calculates the values of the rank-2 spherical components of a set of vectors
    Yl=Ylm(vec,rank)
    """
    X=vec.get('X')
    Y=vec.get('Y')
    Z=vec.get('Z')
    
    
    Yl=dict()
    if rank==1:
        c=np.sqrt(3/(2*np.pi))
        Yl['1,0']=c/np.sqrt(2)*Z
        a=(X+Y*1j)
#        b=np.sqrt(X**2+Y**2)
#        Yl['1,+1']=-c/2*b*a  #a was supposed to equal exp(i*phi), but wasn't normalized (should be normalized by b)
#        Yl['1,-1']=c/2*b*a.conjugate()  #Correction below
        Yl['1,+1']=-c/2*a
        Yl['1,-1']=c/2*a.conjugate()
    elif rank==2:
        c=np.sqrt(15/(32*np.pi))
        Yl['2,0']=c*np.sqrt(2/3)*(3*Z**2-1)
        a=(X+Y*1j)
        b=np.sqrt(X**2+Y**2)
        b2=b**2
        b[b==0]=1
#        Yl['2,+1']=2*c*Z*b*a
#        Yl['2,-1']=2*c*Z*b*a.conjugate()
        Yl['2,+1']=2*c*Z*a
        Yl['2,-1']=2*c*Z*a.conjugate()
#        a=np.exp(2*np.log(X+Y*1j))
#        b=b**2
#        Yl['2,+2']=c*b*a
#        Yl['2,-2']=c*b*a.conjugate()
        a=np.exp(2*np.log(a/b))
        Yl['2,+2']=c*b2*a
        Yl['2,-2']=c*b2*a.conjugate()
        
    Yl['t']=vec['t']
    Yl['index']=vec['index']
    return Yl

def Aqt(Yl,M):
    """
    Project the Ylm onto the eigenmodes
    aqt=Aqt(Yl,M)
    """
    aqt=dict()
    for k,y in Yl.items():
        if k!='t' and k!='index':
            aqt[k]=np.dot(M['m'].T,y.T).T
        else:
            aqt[k]=y
        
    return aqt

def Cqt(aqt,**kwargs):
    
    "Get number of cores"
    if 'parallel' in kwargs:
        p=kwargs.get('parallel')
        if isinstance(p,str) and p.lower()[0]=='n':
            nc=1
        elif isinstance(p,int):
            nc=p if p>0 else 1   #Check the # cores is bigger than 0
        else:                   #Default use parallel processing
            nc=mp.cpu_count()   #Use all available cores
    else:
        nc=mp.cpu_count()
        
    ref_num,v0=ipc.store_vecs(aqt,nc)
    try:
        t0=time()
        with mp.Pool(processes=nc) as pool:
            ct=pool.map(ipc.Ct,v0)
#            print('t={0}'.format(time()-t0))
        ct=ipc.returnCt(ref_num,ct)
    except:
        print('Error in calculating correlation functions')
    finally:
        ipc.clear_data(ref_num)
    
    index=aqt['index']
    N=get_count(index)
    dt=np.diff(aqt['t'][0:2])/np.diff(index[0:2])
    t=np.linspace(0,dt.squeeze()*np.max(index),index[-1]+1)
    i=N!=0
    N=N[i]
    t=t[i]
    ct=dict({'Ct':ct,'t':t,'index':index,'N':N})
    
    return ct

def Cij_t(aqt,i,j,**kwargs):
    """
    Calculates the cross correlation between modes in the iRED analysis, indexed
    by i and j 
    (this function should later be improved using parallel processing for multiple
    pairs of modes. Currently supports only one pair)
    c_ij=Cij_t(aqt,i,j,**kwargs)
    """
    
    index=aqt['index']
    n=np.size(index)
    
    
    for p,(name,a) in enumerate(aqt.items()):
        if p==0:
            ct=np.zeros(index[-1]+1)+0j
        if name!='index' and name!='t':
            for k in range(n):
                ct[index[k:]-index[k]]+=np.multiply(a[k:,i],a[k,j].conjugate())
    N0=get_count(index)
    nz=N0!=0
    N=N0[nz]
    dt=np.diff(aqt['t'][0:2])/np.diff(index[0:2])
    t=np.linspace(0,dt.squeeze()*np.max(index),index[-1]+1)
    t=t[nz]
    ct=np.divide(ct[nz].real,N)
    
    ct=dict({'Ct':ct,'t':t,'index':index,'N':N})
    
    return ct
    
#%% Estimate the correlation function at t=infinity
def CtInf(aqt):
    "Get final value of correlation function"
    ctinf=None
    for k in aqt.keys():
        if k!='t' and k!='index':
            a=aqt.get(k).mean(axis=0)
            if np.shape(ctinf)==():
                ctinf=np.real(a*a.conj())
            else:
                ctinf+=np.real(a*a.conj())
            
    return ctinf

#%% Estimate the correlation function at t=infinity
def Cij_Inf(aqt,i,j):
    "Get final value of correlation function"
    ctinf=None
    for k in aqt.keys():
        if k!='t' and k!='index':
            a=aqt.get(k)[:,i].mean()
            b=aqt.get(k)[:,j].mean()
            if np.shape(ctinf)==():
                ctinf=np.real(a*b.conj())
            else:
                ctinf+=np.real(a*b.conj())
            
    return ctinf

#%% Returns normalized correlation function
def DelCt(ct,ctinf):
    "Get a normalized version of the correlation function (starts at 1, decays to 0)"
    t=ct.get('t')
    ct=ct.get('Ct')
    nt=ct.shape[0]
    ctinf=np.repeat([ctinf],nt,axis=0)
    ct0=np.repeat([ct[0,:]],nt,axis=0)
    delCt={'t':t,'DelCt':(ct-ctinf)/(ct0-ctinf)}
    
    return delCt

def align_mean(vec0):
    """
    Aligns the mean direction of a set of vectors along the z-axis. This can be
    useful for iRED analysis, to mitigate the orientational dependence of the 
    iRED analysis procedure.
    
    vec = align_mean(vec0)
    """
    nt=vec0.get('X').shape[1]
    
    "Mean direction of the vectors"
    X0=vec0['X'].mean(axis=1)
    Y0=vec0['Y'].mean(axis=1)
    Z0=vec0['Z'].mean(axis=1)
    
    "Normalize the length"
    length=np.sqrt(X0**2+Y0**2+Z0**2)
    X0=np.divide(X0,length)
    Y0=np.divide(Y0,length)
    Z0=np.divide(Z0,length)
    
    "Angle away from the z-axis"
    beta=np.arccos(Z0)
    
    "Angle of rotation axis away from y-axis"
    "Rotation axis is at (-Y0,X0): cross product of X0,Y0,Z0 and (0,0,1)"
    theta=np.arctan2(-Y0,X0)
    
    
    xx=np.cos(-theta)*np.cos(-beta)*np.cos(theta)-np.sin(-theta)*np.sin(theta)
    yx=-np.cos(theta)*np.sin(-theta)-np.cos(-theta)*np.cos(-beta)*np.sin(theta)
    zx=np.cos(-theta)*np.sin(-beta)
    
    X=np.repeat(np.transpose([xx]),nt,axis=1)*vec0.get('X')+\
    np.repeat(np.transpose([yx]),nt,axis=1)*vec0.get('Y')+\
    np.repeat(np.transpose([zx]),nt,axis=1)*vec0.get('Z')
    
    xy=np.cos(-theta)*np.sin(theta)+np.cos(-beta)*np.cos(theta)*np.sin(-theta)
    yy=np.cos(-theta)*np.cos(theta)-np.cos(-beta)*np.sin(-theta)*np.sin(theta)
    zy=np.sin(-theta)*np.sin(-beta)
    
    Y=np.repeat(np.transpose([xy]),nt,axis=1)*vec0.get('X')+\
    np.repeat(np.transpose([yy]),nt,axis=1)*vec0.get('Y')+\
    np.repeat(np.transpose([zy]),nt,axis=1)*vec0.get('Z')

    xz=-np.cos(theta)*np.sin(-beta)
    yz=np.sin(-beta)*np.sin(theta)
    zz=np.cos(-beta)
    
    Z=np.repeat(np.transpose([xz]),nt,axis=1)*vec0.get('X')+\
    np.repeat(np.transpose([yz]),nt,axis=1)*vec0.get('Y')+\
    np.repeat(np.transpose([zz]),nt,axis=1)*vec0.get('Z')

    vec={'X':X,'Y':Y,'Z':Z,'t':vec0['t'],'index':vec0['index']}
    
    return vec