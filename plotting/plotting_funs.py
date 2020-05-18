#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 14:23:32 2019

@author: albertsmith
"""

import numpy as np
import matplotlib.pyplot as plt

def plot_cc(Rcc,lbl=None,ax=None,norm=True,index=None,**kwargs):
    """"2D plot of the cross-correlation, given by a square matrix, and an axis label
    plot_cc(Rcc,lbl,ax=None,norm='y',**kwargs)
    """
    if ax==None:
        fig=plt.figure()
        ax=fig.add_subplot(111)
    else:
        fig=ax.figure
        
        
    if norm:
        dg=np.sqrt([np.diag(Rcc)]) 
        "Should we use abs here or not?"
        x=np.abs(Rcc)/np.dot(dg.T,dg)
    else:
        x=Rcc
        
    if index is None:
        index=np.arange(x.shape[0])
    x=x[index][:,index]

        
    if lbl is not None and len(lbl)==x.shape[0]:
        lbl=np.array(lbl)[index]
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
    elif mat.min()<0:
        cmap='RdBu_r'
        if norm:mat[0,0],mat[-1,-1]=1,-1
    else:
        cmap='Blues'

    cax=ax.imshow(mat,interpolation=None,cmap=cmap)
    if norm:ax.imshow(mat1,interpolation=None)
    ax.imshow(mat2,interpolation=None)
    fig.colorbar(cax)

    if 'axis_label' in kwargs:
        axlbl=kwargs.get('axis_label')
    else:
        axlbl='Residue'
    
    ax.set_xlabel(axlbl)
    ax.set_ylabel(axlbl)
    
    "Limit to 50 axis labels"
    while xaxis_lbl is not None and len(lbl)>50:
        xaxis_lbl=np.array(xaxis_lbl)
        xaxis_lbl=xaxis_lbl[range(0,len(lbl),2)]
        lbl=lbl[range(0,len(lbl),2)] 
    
    if xaxis_lbl is not None:
        ax.set_xticks(lbl)
        ax.set_xticklabels(xaxis_lbl,rotation=90)
        ax.set_yticks(lbl)
        ax.set_yticklabels(xaxis_lbl,rotation=0)
    ax.invert_yaxis()
    fig.show()
    
    return ax


def plot_rho_series(data,fig=None,plot_sens=True,index=None,rho_index=None,errorbars=False,**kwargs):
    """
    Plots the full series of detector response (or a limited set, specified by rho_index)
    """


    if fig is None:
        fig=plt.figure()
        
    nd=data.R.shape[1]
    
    rho_index=np.atleast_1d(np.arange(nd) if rho_index is None else np.array(rho_index))
    
    if hasattr(data.sens,'detect_par') and data.sens.detect_par['R2_ex_corr'] and\
        nd-1 in rho_index:
        R2ex=True
    else:
        R2ex=False
    
    if plot_sens and data.sens is not None:
        nplts=np.size(rho_index)+2
        ax0=fig.add_subplot(int(nplts/2)+1,1,1)
        
        temp=data.sens._rho(rho_index,bond=None)
        if R2ex:
            temp[-1][:]=0
            
        hdl=ax0.plot(data.sens.z(),temp.T)
        ax0.set_xlabel(r'$\log_{10}(\tau$ / s)')
        ax0.set_ylabel(r'$\rho(z)$')
        ax0.set_xlim(data.sens.z()[[0,-1]])
        mini=np.min(temp)
        maxi=np.max(temp)
        ax0.set_ylim([mini-(maxi-mini)*.05,maxi+(maxi-mini)*.05])

            
        color=[h.get_color() for h in hdl]
    else:
        nplts=np.size(rho_index)
        color=plt.rcParams['axes.prop_cycle'].by_key()['color']
        
    ax=list()
    
    if index is not None:
        index=np.atleast_1d(index).astype(int)
    else:
        index=np.arange(data.R.shape[0]).astype(int)
    
    if np.size(data.label)==data.R.shape[0]:
        lbl=np.array(data.label)[index]
        if isinstance(lbl[0],str):
            xaxis_lbl=lbl.copy()
            lbl=np.arange(np.size(lbl))
        else:
            xaxis_lbl=None
    else:
        lbl=np.arange(np.size(index))
        xaxis_lbl=None
    
    for k,ri in enumerate(rho_index):
        if k==0:
            ax.append(fig.add_subplot(nplts,1,k+nplts-np.size(rho_index)+1))
        else:
            ax.append(fig.add_subplot(nplts,1,k+nplts-np.size(rho_index)+1,sharex=ax[0]))
        
                    
        if errorbars:
            if data.R_l is None:
                plot_rho(lbl,data.R[index,ri],data.R_std[:,ri],ax=ax[-1],\
                  color=color[k],**kwargs)
            else:
                plot_rho(lbl,data.R[index,ri],[data.R_l[index,ri],data.R_u[index,ri]],ax=ax[-1],\
                  color=color[k],**kwargs)
        else:
            plot_rho(lbl,data.R[index,ri],ax=ax[-1],color=color[k],**kwargs)
                             
        
        
        ax[-1].set_ylabel(r'$\rho_'+str(k)+'^{(\\theta,S)}$')
        
        yl=ax[-1].get_ylim()
        ax[-1].set_ylim([np.min([yl[0],0]),yl[1]])
        
         
        
        if k<np.size(rho_index)-1:
            if xaxis_lbl is not None:
                ax[-1].set_xticklabels(xaxis_lbl)
            plt.setp(ax[-1].get_xticklabels(),visible=False)
            "Limit to 50 axis labels"
        else:
            while xaxis_lbl is not None and len(lbl)>50:
                xaxis_lbl=xaxis_lbl[range(0,len(lbl),2)]
                lbl=lbl[range(0,len(lbl),2)]
            if xaxis_lbl is not None:
                ax[-1].set_xticks(lbl)
                ax[-1].set_xticklabels(xaxis_lbl,rotation=90)
            if R2ex:
                ax[-1].set_ylabel(r'$R_2^{ex} / s^{-1}$')
            
    fig.subplots_adjust(hspace=0.25)    
    fig.show()    
    return ax

def plot_rho(lbl,R,R_std=None,style='plot',color=None,ax=None,split=True,**kwargs):
    """
    Plots a set of rates or detector responses. 
    """
    
    if ax is None:
        ax=plt.figure().add_subplot(111)
    
    "We divide the x-axis up where there are gaps between the indices"
    lbl1=list()
    R1=list()
    R_u1=list()
    R_l1=list()
    
    lbl=np.array(lbl)   #Make sure this is a np array
    if not(np.issubdtype(lbl.dtype,np.number)):
        split=False
        lbl0=lbl.copy()
        lbl=np.arange(len(lbl0))
    else:
        lbl0=None
    
    if split:
        s0=np.where(np.concatenate(([True],np.diff(lbl)>1,[True])))[0]
    else:
        s0=np.array([0,np.size(R)])
    
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
    
    "Plotting style (plot,bar, or scatter, scatter turns the linestyle to '' and adds a marker)"
    if style.lower()[0]=='s':
        if 'marker' not in kwargs:
            kwargs['marker']='o'
        if 'linestyle' not in kwargs:
            kwargs['linestyle']=''
    elif style.lower()[0]=='b':
        if 'linestyle' not in kwargs:
            kwargs['linestyle']=''
        
    
    for lbl,R,R_u,R_l in zip(lbl1,R1,R_u1,R_l1):
        if R_l is None:
            ax.plot(lbl,R,color=color,**kwargs)
        else:
            ax.errorbar(lbl,R,[R_l,R_u],color=color,**kwargs)
        if style.lower()[0]=='b':
            kw=kwargs.copy()
            if 'linestyle' in kw: kw.pop('linestyle')
            ax.bar(lbl,R,color=color,**kw)
        if color is None:
            color=ax.get_children()[0].get_color()
    
    if lbl0 is not None:
        ax.set_xticks(lbl)
        ax.set_xticklabels(lbl0,rotation=90)
        
            
    return ax    
    
    
    
    