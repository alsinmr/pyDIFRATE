#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 14:23:32 2019

@author: albertsmith
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

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


def plot_rho_series(data,fig=None,plot_sens=True,index=None,rho_index=None,errorbars=False,style='plot',**kwargs):
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
            if data.R_l is None or data.R_u is None:
                plot_rho(lbl,data.R[index,ri],data.R_std[:,ri],ax=ax[-1],\
                  color=color[k],style=style,**kwargs)
            else:
                plot_rho(lbl,data.R[index,ri],[data.R_l[index,ri],data.R_u[index,ri]],ax=ax[-1],\
                  color=color[k],style=style,**kwargs)
        else:
            plot_rho(lbl,data.R[index,ri],ax=ax[-1],color=color[k],style=style,**kwargs)
                             
        
        
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
        ebar_clr=color
    elif style.lower()[0]=='b':
        if 'linestyle' not in kwargs:
            kwargs['linestyle']=''
        ebar_clr='black'
    else:
        ebar_clr=color
    
    for lbl,R,R_u,R_l in zip(lbl1,R1,R_u1,R_l1):
        if R_l is None:
            ax.plot(lbl,R,color=color,**kwargs)
        else:
            ax.errorbar(lbl,R,[R_l,R_u],color=ebar_clr,capsize=3,**kwargs)
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
 
#%% Plot the data fit
def plot_fit(lbl,Rin,Rc,Rin_std=None,info=None,index=None,exp_index=None,fig=None,ax=None):
    """
    Plots the fit of experimental data (small data sizes- not MD correlation functions)
    Required inputs are the data label, experimental rates, fitted rates. One may
    also input the standard deviation of the experimental data, and the info
    structure from the experimental data.
    
    Indices may be provided to specify which residues to plot, and which 
    experiments to plot
    
    A figure handle may be provided to specifiy the figure (subplots will be
    created), or a list of axis handles may be input, although this must match
    the number of experiments
    
    plot_fit(lbl,Rin,Rc,Rin_std=None,info=None,index=None,exp_index=None,fig=None,ax=None)
    
    one may replace Rin_std with R_l and R_u, to have different upper and lower bounds
    """
    
    "Apply index to all data"
    if index is not None:
        lbl=lbl[index]
        Rin=Rin[index]
        Rc=Rc[index]
        if Rin_std is not None: Rin_std=Rin_std[index]
        
    "Remove experiments if requested"
    if exp_index is not None:
        if info is not None: 
            info=info.loc[:,exp_index].copy
            info.columns=range(Rin.shape[0])
            
        Rin=Rin[:,exp_index]
        Rc=Rc[:,exp_index]
        if Rin_std is not None: Rin_std=Rin_std[:,exp_index]
    
    nexp=Rin.shape[1]    #Number of experiments
    
    
    #Make sure ax is a list
    if ax is not None and not(isinstance(ax,list)):ax=[ax]
    
    #Make sure there are the correct number of elements in ax
    if ax is not None and length(ax)!=nexp:
        print('Warning: Number of axes provided does not match number of experiments')
        ax=None
    
    #Generate a set of axes if they aren't provided
    if ax is None:  
        if fig is None: #Generate a figure if not provided
            fig=plt.figure()    
        SZ=np.sqrt(nexp)
        SZ=[np.ceil(SZ).astype(int),np.floor(SZ).astype(int)]
        if np.prod(SZ)<nexp: SZ[1]+=1
        
        ntop=np.mod(nexp,SZ[1]) #Number of plots to put in the top row
        if ntop==0:
            ntop=SZ[1]
            skip=0
        else:
            skip=np.floor((SZ[1]-ntop)/2).astype(int)   #Number of plots to skip at beginning
        ax=[fig.add_subplot(SZ[0],SZ[1],k+skip+1) for k in range(ntop)]  #Top row plots
        ax.extend([fig.add_subplot(SZ[0],SZ[1],k+1+SZ[1]) for k in range(nexp-ntop)]) #Other plots
        

    #Make sure the labels are set up
    """Make lbl a numpy array. If label is already numeric, then we use it as is.
    If it is text, then we replace lbl with a numeric array, and store the 
    original lbl as lbl0, which we'll label the x-axis with.
    """
    lbl=np.array(lbl)   #Make sure this is a np array
    if not(np.issubdtype(lbl.dtype,np.number)):
        split=False
        lbl0=lbl.copy()
        lbl=np.arange(len(lbl0))

                    
    else:
        lbl0=None
    
    #Sweep through each experiment
    clr=[k for k in colors.TABLEAU_COLORS.values()]     #Color table
    for k,a in enumerate(ax):
        a.bar(lbl,Rin[:,k],color=clr[np.mod(k,len(clr))])       #Bar plot of experimental data
        if Rin_std is not None:             
            a.errorbar(lbl,Rin[:,k],Rin_std[:,k],color='black',linestyle='',\
                       capsize=3) #Errorbar
        a.plot(lbl,Rc[:,k],linestyle='',marker='o',color='black',markersize=3)
        if k>=nexp-SZ[1]:   #Are we in the last row?
            if lbl0 is not None:
                "Shorten lbl0 if necessary"
                if len(lbl0)>50/SZ[0]:  #Let's say we can fit 50 labels in one figure
                    nlbl=np.floor(50/SZ[0])
                    space=np.floor(len(lbl0)/nlbl).astype(int)
                    i=range(0,len(lbl0),space)
                    a.set_xticks(i)
                    a.set_xticklabels(lbl0[i],rotation=90)
                else:
                    a.set_xticks(range(len(lbl0)))
                    a.set_xticklabels(lbl0,rotation=90)
        else:
            plt.setp(a.get_xticklabels(),visible=False)
        if k==0 or np.mod(k-ntop,SZ[1])==0:
            a.set_ylabel(r'R / s$^{-1}$')
        
        #Apply labels to each plot if we find experiment type in the info array
        if info is not None and 'Type' in info.index.to_numpy():
            if info[k]['Type'] in {'R1','NOE','R2'}:
                a.set_ylim(np.min(np.concatenate(([0],Rin[:,k],Rc[:,k]))),\
                   np.max(np.concatenate((Rin[:,k],Rc[:,k])))*1.25)
                i=info[k]
                string=r'{0} {1}@{2:.0f} MHz'.format(i['Nuc'],i['Type'],i['v0'])
                a.text(np.min(lbl),a.get_ylim()[1]*0.88,string,FontSize=8)
            else:
                a.set_ylim(np.min(np.concatenate(([0],Rin[:,k],Rc[:,k]))),\
                   np.max(np.concatenate((Rin[:,k],Rc[:,k])))*1.45)
                i=info[k]
                string=r'{0} {1}@{2:.0f} MHz'.format(i['Nuc'],i['Type'],i['v0'])
                a.text(np.min(lbl),a.get_ylim()[1]*0.88,string,FontSize=8)
                string=r'$\nu_r$={0} kHz, $\nu_1$={1} kHz'.format(i['vr'],i['v1'])
                a.text(np.min(lbl),a.get_ylim()[1]*0.73,string,FontSize=8)
                
    return ax        
            

def plot_Ct_fit(t,Ct_in,Ct_fit):
    """
    Placeholder for plots of the correlation function fits
    """
    return
   
#%% Plot the rate constant sensitivites
def plot_rhoz(sens,index=None,ax=None,bond=None,norm=False,mdl_num=None,**kwargs):
    """
    Plots the sensitivities found in any sensitivity (child of the model superclass)
    Input is the object, an index of the sensitivities to be plotted, whether
    or not to normalize the plots to 1, an axis to plot on, and any plot settings
    desired. For bond-specific sensitivities, the bond may also be provided.
    
    Except sens, all arguments are optional.
    
    plot_rhoz(sens,index=None,norm=False,ax=None,bond=None,**kwargs)
    """
    if index is None:
        index=sens.info.columns.values
        
    a,_=sens._rho_eff(exp_num=index,bond=bond,mdl_num=mdl_num)
    a=a.T
    
    if hasattr(sens,'detect_par') and sens.detect_par['R2_ex_corr']:
        a=a[:,:-1]

    if norm:
        norm_vec=np.max(np.abs(a),axis=0)
        a=a/np.tile(norm_vec,[np.size(sens.tc()),1])      
    
    if ax is None:
        fig=plt.figure()
        ax=fig.add_subplot(111)
        hdl=ax.plot(sens.z(),a)
#            ax=hdl[0].axes
    else:
        hdl=ax.plot(sens.z(),a)
    
    _set_plot_attr(hdl,**kwargs)
    
        
    ax.set_xlabel(r'$\log_{10}(\tau$ / s)')
    if norm:
        ax.set_ylabel(r'$R$ (normalized)')
    else:
        ax.set_ylabel(r'$R$ / s$^{-1}$')
    ax.set_xlim(sens.z()[[0,-1]])
    ax.set_title('Sensitivity (no model)')
    
    return hdl    

#%% Plot the quality of the rate constant reproduction
def plot_Rc(sens,exp_num=None,norm=True,bond=None,ax=None):
    """
    Plots the input sensitivities compared to their reproduction by fitting to
    detectors. Options are to specifiy experiments (exp_num), to normalize (norm),
    to specify a specific bond (bond), and a specific axis to plot onto (ax). 
    
    plot_Rc(sens,exp_num=None,norm=True,bond=None,ax=None)
    """
    
    nb=sens._nb()
    if nb==1:
        bond=0
                
    a=sens.Rin(bond).T
    b=sens.Rc(bond).T
    
    
    
    if np.size(exp_num)>1 or exp_num!=None:
        a=a[:,exp_num]
        b=b[:,exp_num]
    
    if norm:
        N=np.max(np.abs(a),axis=0)
        a=a/np.tile(N,[np.size(sens.tc()),1]) 
        b=b/np.tile(N,[np.size(sens.tc()),1])
    
    
    if ax is None:
        fig=plt.figure()
        ax=fig.add_subplot(111)
        hdl1=ax.plot(sens.z(),a,'k')
        hdl2=ax.plot(sens.z(),b,'r--')
#            hdl1=plt.plot(self.z(),a,'k')
#            hdl2=plt.plot(self.z(),b,'r--')
#            ax=hdl1[0].axes
    else:
        hdl1=ax.plot(sens.z(),a,'k')
        hdl2=ax.plot(sens.z(),b,'r--')
        
    ax.set_xlabel(r'$\log_{10}(\tau$ / s)')
    if norm:
        ax.set_ylabel(r'$R(z)$ (normalized)')
    else:
        ax.set_ylabel(r'$R(z) / $s$^{-1}$')
    ax.set_xlim(sens.z()[[0,-1]])
    ax.set_title('Rate Constant Reproduction')
    
    hdl=hdl1+hdl2
    
    return hdl

#%% Sets plot attributes from kwargs
def _set_plot_attr(hdl,**kwargs):
    """
    Get properties for a list of handles. If values in kwargs are found in props,
    then that attribute is set (ignores unmatched values)
    """
    if not(hasattr(hdl,'__len__')): #Make sure hdl is a list
        hdl=[hdl]
    
    props=hdl[0].properties().keys()
    for k in kwargs:
        if k in props:
            for m in hdl:
                getattr(m,'set_{}'.format(k))(kwargs.get(k))