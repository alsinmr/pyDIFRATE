#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 11:07:25 2020

@author: albertsmith

Module for performing Principal Component Analysis of an MD trajectory. 

PCA may be used to re-construct dynamics of only certain components

Eventually, we'll implement means of filtering the components to highlight
dynamics relevant to specific detectors
"""

import numpy as np
import os
curdir=os.getcwd()
os.chdir('../r_class')
from Ctsens import Ct
from detectors import detect
os.chdir('../data')
import data_class as dc
import fitting
os.chdir('../')
import MDAnalysis.analysis.pca as pca
import MDAnalysis as mda
os.chdir(curdir)


def runPCA(uni,select,start=None,stop=None,step=None):
    """
    Creates the PCA object, and runs the calculation. 
    """
    
    PSF_pca=pca.PCA(uni,select,start,stop,step)
    PSF_pca.run()
    return PSF_pca


def Ct_ft(x,y=None):
    """
    Calculates a time correlation for x with itself or with y, using Fourier
    transforms. For multi-dimensional data, the transform is taken along the 
    first dimension
    """
    if y is None:y=x
    X=np.fft.fft(np.concatenate((x,x[::-1]),axis=0),axis=0)
    Y=np.fft.fft(np.concatenate((y,y[::-1]),axis=0),axis=0)
    norm=np.sqrt(np.mean(x**2,axis=0)*np.mean(y**2,axis=0))
    ct=(np.fft.ifft(X*np.conjugate(Y),axis=0).real/x.shape[0])/norm/2
    return ct[:len(x)]

def pca2data(PSF_pca,n_components=None,pca_space=None,dt=None):
    """
    Extracts principle components from a trajectory and calculates their correlation
    functions, returning a data object that may be further analyzed
    """
    
    uni=PSF_pca.mean_atoms.universe
    atomgroup=uni.atoms[PSF_pca.mean_atoms.indices]        
    
    
    #Project positions onto PC vectors
    if pca_space is None:
        pca_space=PSF_pca.transform(atomgroup,n_components)
    #Calculate correlation functions for each PC
    ct0=Ct_ft(pca_space)
    #Time axis
    t=np.arange(0,PSF_pca.stop-PSF_pca.start,PSF_pca.step)
    if dt is None:
        dt=uni.trajectory.dt
    t=t*dt
    
    Ct={'Ct':ct0.T,'t':t}
    
    data=dc.data(Ct=Ct)
    data.sens.molecule.mda_object=uni
    data.detect.molecule=data.sens.molecule
    
    return data

def pca2traj(PSF_pca,filename,component_index,linear=False,start=None,stop=None,step=None,pca_space=None):
    """
    Creates an MD trajectory from selected components of a principle component 
    analysis. The result of principle component analysis should be provided
    (PSF_pca from runPCA or from MDAnalysis's pca module, after being run), 
    a filename is required (ending will determine file type, .xtc, .pdb, etc).
    A component index should  be provided, which determines which components will
    be included in the trajectory. Finally, if only a single component is requested,
    one may specify linear as True, which, instead of taking the path of the
    original trajectory for that component, will simply run the component from 
    -2x to 2x its standard deviation (101 steps)
    
    Start, stop, and step indicate which frames to use to back-calculate the 
    trajectory.
    """
    
    uni=PSF_pca.mean_atoms.universe
    atomgroup=uni.atoms[PSF_pca.mean_atoms.indices]
    
    nc0=PSF_pca.variance.size #Total number of components
    
    #Generate a logical index
    component_index=np.atleast_1d(component_index)
    if len(component_index)==nc0:
        ci=np.array(component_index,dtype=bool)
        n=np.argwhere(component_index).squeeze().max()+1
    else:
        ci=np.zeros(nc0,dtype=bool)
        ci[component_index]=True
        n=np.max(component_index)+1
    
    #Get the vector for the pca reconstruction
    if not(linear) or np.sum(ci)>1:
        if pca_space is None:
            pca_space=PSF_pca.transform(atomgroup,n)    
        pca_space=pca_space[:,ci[:n]]
        linear=False
    else:
        stdev=np.sqrt(PSF_pca.variance[ci])
        x=np.linspace(-2*stdev,2*stdev,101)
        pca_space=np.concatenate((x[:-1],x[:0:-1]))
        
    V=PSF_pca.p_components
    
    with mda.Writer(os.path.splitext(filename)[0]+'.pdb',atomgroup.n_atoms) as W:
        W.write(atomgroup)
    
    with mda.Writer(filename, atomgroup.n_atoms) as W:
        if linear:
            for A in pca_space:
                xyz=(np.dot(V[:,ci],A)+PSF_pca.mean).reshape(atomgroup.positions.shape)
                atomgroup.positions=xyz
                W.write(atomgroup)
        else:
            for ts,A in zip(uni.trajectory,pca_space):                
                xyz=(np.dot(V[:,ci],A)+PSF_pca.mean).reshape(atomgroup.positions.shape)
                atomgroup.positions=xyz
                W.write(atomgroup)

def get_dominant_modes(S2c,index,frac=0.9):
    """
    Finds the modes which constitute a certain fraction of the total motion of a
    given bond. Input is the S2c parameter (order parameter matrix for all
    components), the index of the bond of interest, and the fraction of the total
    motion to be reproduced:
        1-S2c[index,i].prod(axis=1)>=frac*(1-S2[index])
        
    i=get_dominant_modes(S2c,index,frac=0.9)
        
    """                

    i0=S2c[index].argsort()
    cp=(S2c[index][i0]).cumprod()
    
    ii=np.argwhere((1-cp)>=(1-cp[-1])*frac)[0][0]
    
    return i0[:ii]

def rho2traj(PSF_pca,rhoz,det_num,filename,pca_space=None):
    """
    Generates a trajectory from PCA results for a given detector sensitivity. 
    One must provide the PSF_pca object after being run, a sensitivity matrix 
    (for example, from detect.rhoz()), and a detector number. Correlation times
    are determined for each principle component. Components are then sorted into
    which detector they influence the most (usually, we recommend max-normalized
    detecotrs). Components sorted into the bin of the requested det_num
    """
    pass

def rho_sorter(PSF_pca,rhoz,pca_space=None):
    """
    Given the results of a PCA and a sensitivity matrix (rhoz), rho_sorter 
    determines which detector is most influenced by each mode.
    """
    
    if hasattr(PSF_pca,'R'):
        data=PSF_pca
    else:
        data=pca2data(PSF_pca,pca_space=pca_space)
    data.detect.r_target(rhoz,NT='M')
    fit=data.fit()
    
    return fit.R.argmax(axis=1)


def pca2zS2(PSF_pca,mol=None,nd=6,pca_space=None):
    """
    Calculates the correlation times for a set of principle components, and if
    mol is provided, also the order parameters for each bond and each principle
    component
    """
    
    uni=PSF_pca.mean_atoms.universe
    atomgroup=uni.atoms[PSF_pca.mean_atoms.indices]
    if pca_space is None:
        pca_space=PSF_pca.transform(atomgroup)

    data=pca2data(PSF_pca,pca_space=pca_space)
    
    data.detect.r_auto(nd)
    fit=data.fit()
    
    z,_,_,_=fitting.fit2tc(fit,df=1)
    
    
    
    if mol is None:
        return z #Just return z if no selection provided
    
    "Get indices to find the atoms of interest in the pca"
    i0=PSF_pca.mean_atoms.indices
    i1=np.array([np.argwhere(i==i0)[0,0]*3+np.arange(3) for i in mol.sel1.indices])
    i2=np.array([np.argwhere(i==i0)[0,0]*3+np.arange(3) for i in mol.sel2.indices])
    
    "Only keep components required for selection i1, i2"
    n=3*mol.sel1.n_atoms
    V0=PSF_pca.p_components
    V1=V0[i1,:]
    V2=V0[i2,:]
    
    del_mean=(PSF_pca.mean[i1]-PSF_pca.mean[i2]) #Mean difference in position
    d=np.sqrt((del_mean**2).sum(axis=1))
    
    eavg=np.repeat(np.atleast_3d(del_mean),pca_space.shape[1],axis=-1).swapaxes(0,1)
    
    e=(V1-V2).swapaxes(0,1)
    
    S20=[[np.zeros(e.shape[1:]) for _ in range(3)] for _ in range(3)]
    
    for k,p in enumerate(pca_space[::10]):
        vec=(e*p)+eavg
        d=np.sqrt((vec**2).sum(axis=0))
        vec=vec/d
        if k/100==round(k/100):
            print(k)
        for k in range(3):
            for j in range(3):
                S20[k][j]+=vec[k]*vec[j]

    n=pca_space[::10].shape[0]
    S2=-1/2*np.ones(e.shape[1:])
    for k in range(3):
        for j in range(3):
            S2+=3/2*(S20[k][j]/n)**2
            
    return z,S2

def pca2totalS2(PSF_pca,mol,pca_space=None):
    """
    Indirectly calculates the total S2 for each selection, based on the variances
    and principal components in the results of a PCA
    
    S2 = pca2totalS2(PSF_pca,mol,pca_space=None)
    """
    
    uni=PSF_pca.mean_atoms.universe
    atomgroup=uni.atoms[PSF_pca.mean_atoms.indices]
    if pca_space is None:
        pca_space=PSF_pca.transform(atomgroup)
    
    
    
    if mol is None:
        return z #Just return z if no selection provided
    
    "Get indices to find the atoms of interest in the pca"
    i0=PSF_pca.mean_atoms.indices
    i1=np.array([np.argwhere(i==i0)[0,0]*3+np.arange(3) for i in mol.sel1.indices])
    i2=np.array([np.argwhere(i==i0)[0,0]*3+np.arange(3) for i in mol.sel2.indices])
    
    "Only keep components required for selection i1, i2"
    V0=PSF_pca.p_components
    V1=V0[i1,:]
    V2=V0[i2,:]
    
    del_mean=(PSF_pca.mean[i1]-PSF_pca.mean[i2]) #Mean difference in position
    d=np.sqrt(((mol.sel1.positions-mol.sel2.positions)**2).sum(axis=1))
    
    eavg=del_mean.T/d
    
    e=((V1-V2).T/d).T
    
    S2=-1/2*np.ones(e.shape[0])
    
    for k in range(3):
        for j in range(3):
            S2+=3/2*(eavg[k]*eavg[j]+np.dot(e[:,k,:]*e[:,j,:],PSF_pca.variance))**2
        
    return S2  

def z_hist(z,S2c,z0=None):
    """
    Calculates histograms of the distributions of motion given the log-
    correlation times (z) and a matrix of S2 for each component. Note, this
    distribution is the real correlation times, not the effective correlation
    times (which is the observable). z0 defines the bins, which are otherwise 
    calculated automatically (z0 may also be set to just a number of bins, and
    the bounds are then automatically calculated).
    
    Notes: 
        1) The returned z-axis is the middle of each bin, not the edges. But, 
           the input bins, if user provided, should be the edges!
        2) In each bin, we return 1-prod(S2) for all S2 in that bin, not the sum
    
    out = z_hist(z,S2c,z0=None)
    """
    
    z=np.array(z)
    S2c=np.array(S2c)
    
    "Make sure z0 is set"
    if z0 is None:
        z0=np.linspace(z.min(),z.max()+1e-12,101)   #Capture the largest value
    elif np.size(z0)==1:
        z0=np.linspace(z.min(),z.max()+1e-12,z0+1)
    else:
        z0=np.array(z0)
        
    i=np.digitize(z,z0,right=False)
    
    I=np.array([1-(S2c[:,i==k]).prod(axis=1) for k in range(i.min(),i.max()+1)])
    
    z0=(z0[:-1]+z0[1:])/2
    
    return {'z':z0,'I':I}
    
def z_hist_eff(z,S2c,z0=None,norm=False):
    """
    Calculates histograms of the distributions of motion given the log-
    correlation times (z) and a matrix of S2 for each component. Note, this
    distribution is the effective correlation times (the observable z_eff). 
    z0 defines the bins, which are otherwise calculated automatically 
    (z0 may also be set to just a number of bins, and the bounds are then 
    automatically calculated)
    
    out = z_hist_eff(z,S2c,z0=None,norm=False)
    """
    
    "Make sure z0 is set"
    if z0 is None:
        z0=np.linspace(np.log10(1/(10**-z).sum())-.1,z.max()+.1,101)   #Capture the largest value
    elif np.size(z0)==1:
        z0=np.linspace(np.log10(1/(10**-z).sum())-.1,z.max()+.1,z0)
    else:
        z0=np.array(z0)
        
    
    I=np.zeros([np.size(z0),S2c.shape[0]])
    
    for z1,S20 in zip(z,S2c.T):
        
        I0=I*S20 #Scaling of everything that is already binned by S20

        i1=np.digitize(z1,z0,right=False) #Index of new contribution
        Is=(1-S20)*(1-I.sum(axis=0))    #Scaled amplitude of new contribution
        
        "Re-binning of everything already binned, for effective correlation times"
        ze=-np.log10(10**-z1+10**-z0)
        i2=np.digitize(ze,z0,right=False)
        I2=np.array([(I[i==i2,:]).sum(axis=0) for i in np.unique(i2)])*(1-S20)
    
    
        I=I0;
        I[i1-1]+=Is
        I[np.unique(i2)-1]+=I2
        
        I[I>1]=1    #Tiny cleanup to avoid roundoff errors
    
    z0=(z0[:-1]+z0[1:])/2
    
    return {'z':z0,'I':I[:-1].T}

def z_hist2det(hist,rhoz,z=None):
    """
    Calculates detector responses from a histogram of amplitudes (output of
    z_hist_eff). Input is the histogram and a matrix of sensitivities (or a 
    detector object). If a matrix is given, the corresponding log-correlation 
    times for those sensitivities must also be given
    
    rho = z_hist2det(hist,rhoz,z)
    """
    
    if z is None:
        z,rhoz=rhoz.z(),rhoz._rho()
    
    I=np.concatenate((np.atleast_2d(1-hist['I'].sum(axis=1)).T,hist['I']),axis=1)
    zh=np.concatenate(([100],hist['z']))  #Put the order parameter into motion at long tc
    
    rho=np.zeros([I.shape[0],rhoz.shape[0]])
    
    for z0,I0 in zip(zh,I.T):
        i=np.argmin(np.abs(z0-z))
        rho+=np.dot(np.atleast_2d(I0).T,[rhoz[:,i]])
        
    return rho