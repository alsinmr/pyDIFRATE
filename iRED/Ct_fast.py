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
from MDAnalysis.analysis.align import rotation_matrix
from psutil import virtual_memory

def Ct2data(molecule,n=100,nr=10,**kwargs):
    
    if molecule.sel1in is None:
        sel1in=np.arange(molecule.sel1.atoms.n_atoms)
    else:
        sel1in=molecule.sel1in
    if molecule.sel2in is None:
        sel2in=np.arange(molecule.sel2.atoms.n_atoms)
    else:
        sel2in=molecule.sel2in
        
        
    nt=molecule.mda_object.trajectory.n_frames
        
    index=trunc_t_axis(nt,n,nr)
        
    vec=get_trunc_vec(molecule.sel1,molecule.sel2,index,sel1in,sel2in,**kwargs)
    
    ct=Ct(vec,**kwargs)
    
    S2=S2calc(vec)
    
    Ctdata=data(molecule=molecule,Ct=ct,S2=S2,**kwargs)
    
    "Still need to make sure that the standard deviation is properly imported from here!!!"
    
    return Ctdata
    
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
    index=np.unique(index).astype('int')
    
    return index
        
def get_trunc_vec(sel1,sel2,index,sel1in=None,sel2in=None,**kwargs):
    """
    vec=get_trunc_vec(sel1,sel2,index,sel1in=None,sel2in=None)
    

    
    """
    
    "Indices to allow using the same atom more than once"
    if sel1in is None:
        sel1in=np.arange(sel1.n_atoms)
    if sel2in is None:
        sel2in=np.arange(sel2.n_atoms)
        
    if sel1.universe!=sel2.universe:
        print('sel1 and sel2 must be generated from the same MDAnalysis universe')
        return
        
    if np.size(sel1in)!=np.size(sel2in):
        print('sel1 and sel2 or sel1in and sel2in must have the same number of atoms')
        return
    
    nt=np.size(index) #Number of time steps
    na=np.size(sel1in) #Number of vectors
    
    X=np.zeros([nt,na])
    Y=np.zeros([nt,na])
    Z=np.zeros([nt,na])
    t=np.zeros([nt])

    traj=sel1.universe.trajectory
    if 'dt' in kwargs:
        dt=kwargs.get('dt')
    else:
        dt=traj.dt
        if traj.units['time']=='ps':    #Convert time units into ns
            dt=dt/1e3
        elif traj.units['time']=='ms':
            dt=dt*1e3

    ts=iter(traj)
    for k,t0 in enumerate(index):
        try:
            traj[t0]     #This jumps to time point t in the trajectory
        except:
            "Maybe traj[t] doesn't work, so we skip through the iterable manually"
            if k!=0:    
                for _ in range(index[k]-index[k-1]):
                    next(ts,None) 
        pos=sel1[sel1in].positions-sel2[sel2in].positions
#        pos=sel1.positions[sel1in]-sel2.positions[sel2in]
        X0=pos[:,0]
        Y0=pos[:,1]
        Z0=pos[:,2]
        
        length=np.sqrt(X0**2+Y0**2+Z0**2)
        
        X[k,:]=np.divide(X0,length)
        Y[k,:]=np.divide(Y0,length)
        Z[k,:]=np.divide(Z0,length)
        t[k]=dt*t0
        if k%int(nt/100)==0 or k+1==nt:
            printProgressBar(k+1, nt, prefix = 'Loading:', suffix = 'Complete', length = 50) 

    vec={'X':X,'Y':Y,'Z':Z,'t':t,'index':index}
    
    if not('alignCA' in kwargs and kwargs.get('alignCA').lower()[0]=='n'):
        "Default is always to align the molecule (usually with CA)"
        vec=align(vec,sel1.universe,**kwargs)
           
    return vec

def Ct(vec,**kwargs):    
    if'n_cores' in kwargs:
        nc=np.min([kwargs.get('n_cores'),nc])
    else:
        nc=mp.cpu_count()
        
    nb=vec['X'].shape[1]
    
    v0=list()   #Store date required for each core
    for k in range(nc):
        v0.append((vec['X'][:,range(k,nb,nc)],vec['Y'][:,range(k,nb,nc)],vec['Z'][:,range(k,nb,nc)],vec['index']))
    
    if 'parallel' in kwargs and kwargs.get('parallel').lower()[0]=='n':
        ct0=list()
        for v in v0:
            ct0.append(Ct_par(v))
    else:
        with mp.Pool(processes=nc) as pool:
            ct0=pool.map(Ct_par,v0)
    
    
    "Get the count of number of averages"
    index=vec['index']
    N=np.zeros(index[-1]+1)
    n=np.size(index)
    for k in range(n):
        N[index[k:]-index[k]]+=1
        
    i=N!=0
    N=N[i]
    
    ct=np.zeros([np.size(N),nb])
    N0=N
    N=np.repeat([N],np.shape(ct0)[2],axis=0).T
    for k in range(nc):
        ct[:,range(k,nb,nc)]=np.divide(ct0[k][i],N)
        
    
    dt=(vec['t'][1]-vec['t'][0])/(vec['index'][1]-vec['index'][0])
    t=np.linspace(0,dt*np.max(index),index[-1]+1)
    t=t[i]
    
    Ct={'t':t,'Ct':ct.T,'N':N0,'index':index}
    
    return Ct

def Ct_par(v):
    index=v[3]
    X=v[0]
    Y=v[1]
    Z=v[2]
    
    n=np.size(index)
    c=np.zeros([np.max(index)+1,np.shape(X)[1]])
    
    for k in range(n):
        c[index[k:]-index[k]]+=(3*(np.multiply(X[k:],X[k])+np.multiply(Y[k:],Y[k])\
             +np.multiply(Z[k:],Z[k]))**2-1)/2
#        if k%int(n/100)==0 or k+1==n:
#            printProgressBar(k+1, n, prefix = 'C(t) calc:', suffix = 'Complete', length = 50) 
    return c
    
def align(vec0,uni,**kwargs):
    """
    Removes overall rotation from a trajectory, by aligning to a set of reference
    atoms. Default is protein backbone CA. 
    """
    if 'align_ref' in kwargs:
        uni0=uni.select_atoms(kwargs.get('align_ref'))
    else:
        uni0=uni.select_atoms('name CA')
        if uni0.n_atoms==0:
            uni0=uni.select_atoms('name C11')   #Not sure about this. Alignment for lipids?
    
    ref0=uni0.positions-uni0.atoms.center_of_mass()
    
    SZ=np.shape(vec0.get('X'))
    index=vec0['index']
    "Pre-allocate the direction vector"
    vec={'X':np.zeros(SZ),'Y':np.zeros(SZ),'Z':np.zeros(SZ),'t':vec0.get('t'),'index':index} 

    nt=vec0['t'].size

    
    traj=uni.trajectory
    ts=iter(traj)
    for k,t0 in enumerate(index):
        try:
            traj[t0]     #This jumps to time point t in the trajectory
        except:
            "Maybe traj[t] doesn't work, so we skip through the iterable manually"
            if k!=0:    
                for _ in range(index[k]-index[k-1]):
                    next(ts,None) 
        "CA positions"
        pos=uni0.positions-uni0.atoms.center_of_mass()
        
        "Rotation matrix for this time point"
        R,_=rotation_matrix(pos,ref0)
        vec['X'][k,:]=vec0['X'][k,:]*R[0,0]+vec0['Y'][k,:]*R[0,1]+vec0['Z'][k,:]*R[0,2]
        vec['Y'][k,:]=vec0['X'][k,:]*R[1,0]+vec0['Y'][k,:]*R[1,1]+vec0['Z'][k,:]*R[1,2]
        vec['Z'][k,:]=vec0['X'][k,:]*R[2,0]+vec0['Y'][k,:]*R[2,1]+vec0['Z'][k,:]*R[2,2]
        if k%int(np.size(index)/100)==0 or k+1==nt:
            printProgressBar(k+1, np.size(index), prefix = 'Aligning:', suffix = 'Complete', length = 50) 
        
    return vec

def S2calc(vec):
    v=[vec.get('X'),vec.get('Y'),vec.get('Z')]
    S2=np.zeros(np.shape(vec.get('X'))[1])
    for k in v:
        for m in v:
            S2+=np.mean(k*m,axis=0)**2
    
    S2=3/2*S2-1/2
    
    return S2     

def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█'):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    # Print New Line on Complete
    if iteration == total: 
        print()