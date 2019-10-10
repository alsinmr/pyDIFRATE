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