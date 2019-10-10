#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 10:24:19 2019

@author: albertsmith
"""

import numpy as np
import multiprocessing as mp
from MDAnalysis.analysis.align import rotation_matrix
from time import time
from parCt import par_class as pct

#%% Calculates a truncated time axis
def trunc_t_axis(nt,n=100,nr=10,**kwargs):
    """
    Calculates a log-spaced sampling schedule for an MD time axis. Parameters are
    nt, the number of time points, n, which is the number of time points to 
    load in before the first time point is skipped, and finally nr is how many
    times to repeat that schedule in the trajectory (so for nr=10, 1/10 of the
    way from the beginning of the trajectory, the schedule will start to repeat, 
    and this will be repeated 10 times)
    
    """
    
    """We set n=-1 to just get the full axis. We'll eventually eliminate the
    various non-fast functions and simply acquire them this way
    """
    if n==-1:
        index=np.arange(nt)
        return index
    
    "Step size: this log-spacing will lead to the first skip after n time points"
    logdt0=np.log10(1.50000001)/n
    
    index=list()
    index.append(0)
    dt=0
    while index[-1]<nt: #Step through until the last index is greater than nt
        index.append(index[-1]+np.round(10**dt))
        dt+=logdt0
        
    index=np.array(index)
    
    "Repeat this indexing nr times throughout the trajectory"
    index=np.repeat(index,nr,axis=0)+np.repeat([np.arange(0,nt,nt/nr)],index.size,axis=0).reshape([index.size*nr])
    
    "Eliminate indices >= nt, eliminate repeats, and sort the index"
    "(repeats in above line lead to unsorted axis, unique gets rid of repeats and sorts)"
    index=index[index<nt]
    index=np.unique(index.astype('int'))
    return index

#%% Determine how many frame pairs are averaged into each time point
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

#%% Estimate the order parameter
def S2calc(vec):
    """
    Calculates an estimate of the order parameter, according to
    3/2*(<x^2>+<y^2>+<z^2>+2<x*y>+2<x*z>+2<y*z>)-1/2 with averages performed
    over the complete vector
    
    S2=S2calc(vec)
    """
    v=[vec.get('X'),vec.get('Y'),vec.get('Z')]
    S2=np.zeros(np.shape(vec.get('X'))[1])
    for k in v:
        for m in v:
            S2+=np.mean(k*m,axis=0)**2
    
    S2=3/2*S2-1/2    
    return S2

#%% Returns the correlation function defined by vec
def Ct(vec,**kwargs):
    """
    Calculates the correlation functions for vectors with unequal spacing in the
    time axis. By default, uses parallel processing (using all available cores)
    Optional arguments are parallel, which determines whether or not to use 
    parallel processing ('y'/'n'), or optionally one may simply set parallel to
    the desired number of cores (parallel=4, for example)
    """    
    
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
           
    
    if 'n_cores' in kwargs:
        nc=np.min([kwargs.get('n_cores'),nc])
        print('Warning: n_cores argument will be removed in a later version. set parallel=n_cores')
        "Optional second argument. Not documented- possibly will be removed"
        
    nb=vec['X'].shape[1]
    
    
    if nc==1:
        "Might later replace this with the code in place"
        "But, should keep some variant- parallel version isn't necessarily stable"
        v0=list()   #Store date required for each core
        for k in range(nc):
            v0.append((vec['X'][:,range(k,nb,nc)],vec['Y'][:,range(k,nb,nc)],vec['Z'][:,range(k,nb,nc)],vec['index']))
    
    if nc==1:   #Series processing
        ct0=list()
        for v in v0:
            ct0.append(Ct_par(v))
    else:
        ref_num,v0=pct.store_vecs(vec,nc)
        try:
            with mp.Pool(processes=nc) as pool:
#                ct=pool.map(ctpar.Ct,v0)
                ct=pool.map(pct.Ct,v0)
            ct=pct.returnCt(ref_num,ct)
        finally:
            pct.clear_data(ref_num)
            
    "Get the count of number of averages"
    index=vec['index']
    N=get_count(index)
    
    "i finds all times points for which we have some averaging of the correlation function"    
    i=N!=0
    N=N[i]
    

    N0=N
    
    if nc==1:
        ct=np.zeros([np.size(N0),nb])
        for k in range(nc):
            N=np.repeat([N0],np.shape(ct0[k])[1],axis=0).T  #N with same shape as ct
            ct[:,range(k,nb,nc)]=np.divide(ct0[k][i],N) #Normalize correlation function based on how many averages
        
    
    dt=(vec['t'][1]-vec['t'][0])/(index[1]-index[0])
    t=np.linspace(0,dt*np.max(index),index[-1]+1)
    t=t[i]
    
    Ct={'t':t,'Ct':ct.T,'N':N0,'index':index}
    
    return Ct


#%% Parallel function to calculate correlation functions
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

#%% Load in the truncated vectors from the trajectory
def get_trunc_vec(molecule,index,**kwargs):
    """
    vec=get_trunc_vec(molecule,index,**kwargs)
    
    Returns time-dependent vectors defined in the molecule object. Usually this
    is vectors defined by atom selections in sel1 and sel2 (and possibly indexed
    by sel1in and sel2in). Alternatively, if function-defined vectors are stored
    in molecule._vf (molecule.vec_fun() returns vectors), then these will be 
    returned instead
    
    One must provide the molecule object, and an index determining which time 
    points to analyze.
    
    Optional arguments are dt, which re-defines the time step
    (vs. the time step returned by MDAnalysis), and align, which can be set to
    'y' and will remove overall motion by aligning all frames to a reference
    set of atoms. Default is CA in proteins. To change default, provide a second
    argument, align_ref, which is an MDAnalysis selection string. This string 
    will select from all atoms in the trajectory, and align them.
    
    
    """
    
    if molecule._vf is not None:
        vf=molecule.vec_fun
        special=True
    else:
        sel1=molecule.sel1
        sel2=molecule.sel2
        sel1in=molecule.sel1in
        sel2in=molecule.sel2in
        
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
        special=False
    
    nt=np.size(index) #Number of time steps
    if special:
        na=vf().shape[1]
    else:
        na=np.size(sel1in) #Number of vectors
    
    X=np.zeros([nt,na])
    Y=np.zeros([nt,na])
    Z=np.zeros([nt,na])
    t=np.zeros([nt])

    uni=molecule.mda_object
    traj=uni.trajectory
    if 'dt' in kwargs:
        dt=kwargs.get('dt')
    else:
        dt=traj.dt/1e3
#        if traj.units['time']=='ps':    #Convert time units into ns
#            dt=dt/1e3
#        elif traj.units['time']=='ms':
#            dt=dt*1e3
        

    ts=iter(traj)
    for k,t0 in enumerate(index):
        try:
            traj[t0]     #This jumps to time point t in the trajectory
        except:
            "Maybe traj[t] doesn't work, so we skip through the iterable manually"
            if k!=0:    
                for _ in range(index[k]-index[k-1]):
                    next(ts,None)
                    
        if special:
            "Run the function to return vector"
            X0,Y0,Z0=vf()
        else:
            "Else just get difference in atom positions"
            pos=sel1[sel1in].positions-sel2[sel2in].positions
    #        pos=sel1.positions[sel1in]-sel2.positions[sel2in]
            X0=pos[:,0]
            Y0=pos[:,1]
            Z0=pos[:,2]
        
        "Make sure length is one"
        length=np.sqrt(X0**2+Y0**2+Z0**2)
        X[k,:]=np.divide(X0,length)
        Y[k,:]=np.divide(Y0,length)
        Z[k,:]=np.divide(Z0,length)
        "Keep track of the time axis"
        t[k]=dt*t0
        if k%int(nt/100)==0 or k+1==nt:
            printProgressBar(k+1, nt, prefix = 'Loading:', suffix = 'Complete', length = 50) 

    vec={'X':X,'Y':Y,'Z':Z,'t':t,'index':index}
    
    "Re-align vectors to some set of reference atoms"
    if 'align' in kwargs and kwargs.get('align').lower()[0]=='y':
        "Default does not align molecule"
        vec=align(vec,uni,**kwargs)
        
#    "Re-align vectors so they all point along z"
#    if 'align_iRED' in kwargs and kwargs.get('align_iRED').lower()[0]=='y':
#        vec=align_mean(vec)
           
    return vec

#%% Removes 
def align(vec0,uni,**kwargs):
    """
    Removes overall rotation from a trajectory, by aligning to a set of reference
    atoms. Default is protein backbone CA. If no CA found, try C11 for lipids 
    (possibly this isn't standard- shouldn't create problems for the time being).
    Next try all carbons, and finally all atoms)
    """
    if 'align_ref' in kwargs:
        uni0=uni.select_atoms(kwargs.get('align_ref'))
    else:
        uni0=uni.select_atoms('name CA')    #Standard alignment for proteins
        if uni0.n_atoms==0:
            uni0=uni.select_atoms('name C11')   #Not sure about this. Alignment for lipids?
        if uni0.n_atoms==0:
            uni0=uni.select_atoms('type C') #Try for all carbons
        if uni0.n_atoms==0:
            uni0=uni.atoms #Take all atoms
    
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
        "Ref positions"
        pos=uni0.positions-uni0.atoms.center_of_mass()
        
        "Rotation matrix for this time point"
        R,_=rotation_matrix(pos,ref0)
        "Apply the rotation matrix to the input vector"
        vec['X'][k,:]=vec0['X'][k,:]*R[0,0]+vec0['Y'][k,:]*R[0,1]+vec0['Z'][k,:]*R[0,2]
        vec['Y'][k,:]=vec0['X'][k,:]*R[1,0]+vec0['Y'][k,:]*R[1,1]+vec0['Z'][k,:]*R[1,2]
        vec['Z'][k,:]=vec0['X'][k,:]*R[2,0]+vec0['Y'][k,:]*R[2,1]+vec0['Z'][k,:]*R[2,2]
        "Print out progress"
        if k%int(np.size(index)/100)==0 or k+1==nt:
            printProgressBar(k+1, np.size(index), prefix = 'Aligning:', suffix = 'Complete', length = 50) 
        
    return vec



#%% Progress bar for loading/aligning
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