#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 14:23:32 2019

@author: albertsmith
"""

import numpy as np
import matplotlib.pyplot as plt

def plot_cc(Rcc,lbl=None,ax=None,norm='y',**kwargs):
    """"2D plot of the cross-correlation, given by a square matrix, and an axis label
    plot_cc(Rcc,lbl,ax=None,norm='y',**kwargs)
    """
    if ax==None:
        fig=plt.figure()
        ax=fig.add_subplot(111)
    else:
        fig=ax.figure
        
        
    if norm.lower()[0]=='y':
        dg=np.sqrt([np.diag(Rcc)])
        x=Rcc/np.dot(dg.T,dg)
    else:
        x=Rcc
        
    if lbl is not None:
        if isinstance(lbl[0],str):
            xaxis_lbl=lbl.copy()
            lbl=np.arange(np.size(lbl))
        else:
            xaxis_lbl=None
    else:
        lbl=np.arange(0,Rcc.shape[0])
        xaxis_lbl=None
    
    sz=(np.max(lbl)+1)*np.array([1,1])
    mat=np.zeros(sz)
    mat1=np.zeros([sz[0],sz[1],4])
    mat2=np.ones([sz[0],sz[1],4])*0.75
    mat2[:,:,3]=1
    
    for i,k in enumerate(lbl):
        mat[k][np.array(lbl)]=x[i,:]
        mat1[k,k,3]=1
        mat2[k,np.array(lbl),3]=0
        
#        mat1[:,:,3]=-(mat1[:,:,3]-1)
    
    if 'cmap' in kwargs:
        cmap=kwargs.get('cmap')
    else:
        cmap='Blues'

    cax=ax.imshow(np.abs(mat),interpolation=None,cmap=cmap)
    if norm.lower()[0]=='y':
        ax.imshow(mat1,interpolation=None)
    ax.imshow(mat2,interpolation=None)
    fig.colorbar(cax)

    if 'axis_label' in kwargs:
        axlbl=kwargs.get('axis_label')
    else:
        axlbl='Residue'
    
    ax.set_xlabel(axlbl)
    ax.set_ylabel(axlbl)
    
    if xaxis_lbl is not None:
        ax.set_xticks(lbl)
        ax.set_xticklabels(xaxis_lbl,rotation=90)
        ax.set_yticks(lbl)
        ax.set_yticklabels(xaxis_lbl,rotation=0)
    
    fig.show()
    
    return ax


def plot_rho(lbl,R,R_std=None,style='plot',color=None,ax=None,split=True,**kwargs):
    """
    Plots a set of rates or detector responses. 
    """
    
    if ax==None:
        ax=plt.figure().add_subplot()
    
    "We divide the x-axis up where there are gaps between the indices"
    lbl1=list()
    R1=list()
    R_u1=list()
    R_l1=list()
    
    if split:
        s0=np.where(np.concatenate(([True],np.diff(lbl)>1,[True])))[0]
    else:
        s0=np.array([0,np.size(R1)])
    
    for s1,s2 in zip(s0[:-1],s0[1:]):
        lbl1.append(lbl[s1:s2])
        R1.append(R[s1:s2])
        if R_std is not None:
            if np.ndim(R_std)==2:
                R_l1.append(R_std[0][s1:s2])
                R_u1.append(R_std[1][s1:s2])
            else:
                R_l1.append(R_std[s1:s2])
                R_u1.append(R_std[s1:s2])
        else:
            R_l1.append(None)
            R_u1.append(None)
    
    "Plotting style (plot or scatter, scatter turns the linestyle to '' and adds a marker)"
    if style.lower()[0]=='s':
        if 'marker' not in kwargs:
            kwargs['marker']='o'
        if 'linestyle' not in kwargs:
            kwargs['linestyle']=''
        
    
    for lbl,R,R_u,R_l in zip(lbl1,R1,R_u1,R_l1):
        if R_l is None:
            ax.plot(lbl,R,color=color,**kwargs)
        else:
            ax.errorbar(lbl,R,[R_l,R_u],color=color,**kwargs)
        if color is None:
            color=ax.get_children()[0].get_color()
                
    return ax    
    
    
    
    