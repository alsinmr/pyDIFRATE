#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  7 16:51:28 2019

@author: albertsmith
"""

import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.patches as patch
from matplotlib.collections import PatchCollection
os.chdir('../r_class')
from Ctsens import Ct
from detectors import detect
os.chdir('../chimera')
from chimera_funs import plot_cc as plt_cc3D
from chimera_funs import plot_rho
os.chdir('../data')
from fitting import fit_data

class data(object):
#%% Basic container for experimental and MD data
    def __init__(self,**kwargs):
        
        
        self.label=None
        
        self.R=None
        self.R_std=None
        self.R_u=None
        self.R_l=None
        self.conf=0.68
        
        self.S2=None
        self.tot_cc=None
        self.tot_cc_norm=None
        
        self.Rcc=list()
        self.Rcc_norm=list()
        
        self.Rin=None
        self.Rin_std=None
        self.Rc=None
        
        self.ired=None
        self.type=None
        
        self.sens=None
        self.detect=None
        
        
        
        self.load(**kwargs)
#%% Some options for loading in data        
    def load(self,**kwargs):
        "Load in correlation functions from an iRED calculation"
        if 'iRED' in kwargs:
            self.ired=kwargs.get('iRED')
            self.R=self.ired.get('DelCt')
            nt=self.ired.get('t').size
            self.sens=Ct(t=self.ired.get('t'),**kwargs)
            self.detect=detect(self.sens)
#            self.R_std=np.repeat([self.sens.info.loc['stdev']],self.R.shape[0],axis=0)
            
            norm=1/(self.ired.get('Ct')[:,0]-self.ired.get('CtInf'))
            norm=np.transpose([norm/norm[0:-(self.ired.get('rank')*2+1)].mean()])

#            self.std=np.dot(norm,[self.sens.info.loc['stdev']])
        

            self.R_std=np.dot(norm,[self.sens.info.loc['stdev']])
            
        if 'Ct' in kwargs:
            ct=kwargs.get('Ct')
            self.R=ct.get('Ct')
            self.sens=Ct(t=ct.get('t'),**kwargs)
            self.detect=detect(self.sens)
            self.R_std=np.repeat([self.sens.info.loc['stdev']],self.R.shape[0],axis=0)
            if 'S2' in kwargs:
                self.S2=kwargs.get('S2')
            molecule=kwargs.get('molecule')
            
        if 'molecule' in kwargs:
            mol=kwargs.get('molecule')
            if self is not None:
                self.sens.molecule=mol
                self.detect.molecule=mol
                self.sens.molecule.set_selection()
            self.label=mol.label
             
        if 'EstErr' in kwargs and kwargs.get('EstErr').lower()[0]=='y':
            self.detect.r_nopt(np.min([15,nt]))
            

    def new_detect(self,**kwargs):
        if 'sens' in kwargs:
            self.detect=detect(kwargs.get('sens'),**kwargs)
        else:
            self.detect=detect(self.sens,**kwargs)
            
#%% Run fit_data from the object     
    def fit(self,detect=None,**kwargs):
        if detect is None:
            detect=self.detect
            
        return fit_data(self,detect,**kwargs)
            
#%% Convert iRED data types into normal detector responses and cross-correlation matrices            
    def iRED2rho(self):
        if self.ired is None or not isinstance(self.sens,detect):
            print('Function only applicable to iRED-derived detector responses')
            return
        
        out=data()

        nd=self.sens.r(bond=0).shape[1]
        nb=self.R.shape[0]
        

        
        rank=self.ired.get('rank')
        ne=2*rank+1
        
        
        if self.sens.molecule.sel1in is not None:
            nb0=np.size(self.sens.molecule.sel1in)
        elif self.sens.molecule.sel2in is not None:
            nb0=np.size(self.sens.molecule.sel2in)
        else:
            nb0=self.sens.molecule.sel1.n_atoms
        
        out.R=np.zeros([nb0,nd])
        out.R_std=np.zeros([nb0,nd])
        out.R_l=np.zeros([nb0,nd])
        out.R_u=np.zeros([nb0,nd])
        
        for k in range(0,nd):
            lambda_rho=np.repeat([self.ired.get('lambda')[0:-ne]*self.R[0:-ne,k]],nb0,axis=0)
            out.R[:,k]=np.sum(lambda_rho*self.ired.get('m')[0:nb0,0:-ne]**2,axis=1)
             
            lambda_rho=np.repeat([self.ired.get('lambda')[0:-ne]*self.R_std[0:-ne,k]],nb0,axis=0)
            out.R_std[:,k]=np.sum(lambda_rho*self.ired.get('m')[0:nb0,0:-ne]**2,axis=1)
            
            lambda_rho=np.repeat([self.ired.get('lambda')[0:-ne]*self.R_l[0:-ne,k]],nb0,axis=0)
            out.R_l[:,k]=np.sum(lambda_rho*self.ired.get('m')[0:nb0,0:-ne]**2,axis=1)
            
            lambda_rho=np.repeat([self.ired.get('lambda')[0:-ne]*self.R_u[0:-ne,k]],nb0,axis=0)
            out.R_u[:,k]=np.sum(lambda_rho*self.ired.get('m')[0:nb0,0:-ne]**2,axis=1)
            
            
        
        out.sens=self.sens
        
        "Pre-allocate nd matrices for the cross-correlation calculations"
        out.tot_cc=np.zeros([nb0,nb0])
        for k in range(0,nd):
            out.Rcc.append(np.zeros([nb0,nb0]))
        "Loop over all eigenmodes"
        for k in range(0,nb-ne): #We exclude the 2*rank+1 global motions
            m=self.ired.get('m')[0:nb0,k]
            mat=self.ired.get('lambda')[k]*np.dot(np.transpose([m]),[m])
            "Calculate total correlation"
            out.tot_cc+=mat
            "Loop over all detectors"
            for m in range(0,nd):
                out.Rcc[m]+=mat*self.R[k,m]
                
        "Calculate the normalized correlation"
        dg=np.sqrt([np.diag(out.tot_cc)])
        out.tot_cc_norm=out.tot_cc/np.dot(np.transpose(dg),dg)
        for k in range(0,nd):
            dg=np.sqrt([np.diag(out.Rcc[k])])
            out.Rcc_norm.append(out.Rcc[k]/np.dot(np.transpose(dg),dg))
        
        if self.label is not None:
            out.label=self.label    
        elif self.sens is not None and np.size(self.sens.molecule.label)!=0:
            out.label=self.sens.molecule.label
            
        "Calculate the order parameters"

        lda=np.repeat([self.ired.get('lambda')[0:-ne]],nb0,axis=0)
        out.S2=1-np.sum(lda*self.ired.get('m')[0:nb0,0:-ne]**2,axis=1)

        return out
    
    def plot_rho(self,fig=None,plot_sens='y',index=None,**kwargs):
        if fig is None:
            fig=plt.figure()
            
        nd=self.R.shape[1]
        
        if plot_sens.lower()[0]=='y' and self.sens is not None:
            nplts=nd+2
            ax0=fig.add_subplot(int(nplts/2)+1,1,1)
            hdl=ax0.plot(self.sens.z(),np.transpose(self.sens._rho(np.arange(0,nd),bond=None)))
            ax0.set_xlabel(r'$\log_{10}(\tau$ / s)')
            ax0.set_ylabel(r'$\rho(z)$')
            ax0.set_xlim(self.sens.z()[[0,-1]])
        else:
            nplts=nd
            
        ax=list()
        
        if index is not None:
            index=np.atleast_1d(index).astype(int)
        else:
            index=np.arange(self.R.shape[0]).astype(int)
        
        if np.size(self.label)!=0:
            lbl=np.array(self.label)[index]
            if isinstance(lbl[0],str):
                xaxis_lbl=lbl.copy()
                lbl=np.arange(np.size(lbl))
            else:
                xaxis_lbl=None
        else:
            lbl=np.arange(np.size(index))
            xaxis_lbl=None
        
        
        for k in range(0,nd):
            ax.append(fig.add_subplot(nplts,1,k+nplts-nd+1))
            if self.R_l is None:
                ax[k].errorbar(lbl,self.R[index,k],self.R_std[:,k],color=hdl[k].get_color(),\
                  **kwargs)
            else:
                ax[k].errorbar(lbl,self.R[index,k],[self.R_l[index,k],self.R_u[index,k]],color=hdl[k].get_color(),\
                  **kwargs)
            ax[k].set_ylabel(r'$\rho_'+str(k)+'^{(\\theta,S)}$')
            
            if k<nd-1:
                ax[k].set_xticklabels([])
            elif xaxis_lbl is not None:
                ax[k].set_xticks(lbl)
                ax[k].set_xticklabels(xaxis_lbl,rotation=90)
        fig.subplots_adjust(hspace=0.25)
        
            
    def plot_cc(self,det_num,cutoff=None,ax=None,norm='y',**kwargs):
        if np.size(self.Rcc)==0:
            print('Data object does not contain cross-correlation data')
            print('First, create a data object from iRED analysis (data=iRED2data(...))')
            print('Then, analyze with detectors, data.r_auto(...);fit0=data.fit(...)')
            print('Finally, convert fit into normal detector responses, fit=fit0.iRED2rho()')
            return
        if ax==None:
            fig=plt.figure()
            ax=fig.add_subplot(111)
            
        if norm.lower()[0]=='y':
            if det_num is None:
                x=self.tot_cc_norm
            else:
                x=self.Rcc_norm[det_num]
        else:
            if det_num is None:
                x=self.tot_cc
            else:
                x=self.Rcc[det_num]
            
            
        if self.label is not None:
            lbl=self.label
            if isinstance(lbl[0],str):
                xaxis_lbl=lbl.copy()
                lbl=np.arange(np.size(lbl))
            else:
                xaxis_lbl=None
        else:
            lbl=np.arange(0,self.R.shape[0])
        
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
        ax.imshow(mat1,interpolation=None)
        ax.imshow(mat2,interpolation=None)
        fig.colorbar(cax)

        if 'axis_label' in kwargs:
            axlbl=kwargs.get('axis_label')
        else:
            axlbl='Residue'
        
        ax.set_xlabel(axlbl)
        ax.set_ylabel(axlbl)
        if det_num is None:
            ax.set_title('Total cross correlation')
        else:
            ax.set_title(r'Cross correlation for $\rho_{' + '{}'.format(det_num) + '}$')
        
        
        if xaxis_lbl is not None:
            ax.set_xticks(lbl)
            ax.set_xticklabels(xaxis_lbl,rotation=90)
            ax.set_yticks(lbl)
            ax.set_yticklabels(xaxis_lbl,rotation=0)
        
        
        
        return ax
       
        
                
    def draw_cc3D(self,bond,det_num=None,chain=None,fileout=None,scaling=None,norm='y',**kwargs):
        "bond is the user-defined label! Not the absolute index..."

        if self.label is None:
            print('User has not defined any bond labels, bond will now refer to the absolute index')
            index=bond
        elif any(np.atleast_1d(self.label)==bond):
            index=np.where(np.array(self.label)==bond)[0][0]
            print(index)

        if norm.lower()[0]=='y':
            if det_num is None:
                values=self.tot_cc_norm[index,:]
            else:
                values=self.Rcc_norm[det_num][index,:]
        else:
            if det_num is None:
                values=self.tot_cc[index,:]
            else:
                values=self.Rcc[det_num][index,:]
        
        "Take absolute value- I'm not convinced about this yet..."
        values=np.abs(values)

        if scaling is None:
#            "Default is to scale to the maximum of all correlations"
#            scale0=0
#            for k in range(0,np.shape(self.Rcc_norm)[0]):
#                a=self.Rcc_norm[k]-np.eye(np.shape(self.Rcc_norm)[1])
#                scale0=np.max([scale0,np.max(np.abs(a))])
            if norm.lower()[0]=='y':
#                if det_num is None:
#                    scale0=np.max(np.abs(self.tot_cc_norm)-np.eye(np.shape(self.tot_cc_norm)[0]))
#                else:
#                    scale0=np.max(np.abs(self.Rcc_norm[det_num]-np.eye(np.shape(self.Rcc_norm)[1])))
                scale0=1
            else:
                scale0=np.max(np.abs(values))
            scaling=1/scale0

        res1=self.sens.molecule.sel1.resids
        chain1=self.sens.molecule.sel1.segids
        res2=self.sens.molecule.sel2.resids
        chain2=self.sens.molecule.sel2.segids

        if np.all(res1==res2) and np.all(chain1==chain2):
            resi=res1
            chain=chain1
            chain[chain=='PROA']='p'
            plt_cc3D(self.sens.molecule,resi,values,resi0=bond,chain=chain,chain0=chain[index],\
                     fileout=fileout,scaling=scaling,color_scheme='blue',**kwargs)
        else:
            plt_cc3D(self.sens.molecule,None,values,resi0=bond,scaling=scaling,color_scheme='red',**kwargs)
#            print('Selections over multiple residues/chains- not currently implemented')
        
        
        
    def draw_rho3D(self,det_num=None,resi=None,fileout=None,scaling=None,**kwargs):
        
        if det_num is None:
            values=1-self.S2
        else:
            values=self.R[:,det_num]
        
      
        res1=self.sens.molecule.sel1.resids
        chain1=self.sens.molecule.sel1.segids
        res2=self.sens.molecule.sel2.resids
        chain2=self.sens.molecule.sel2.segids
        

        if np.size(res1)==np.size(res2) and (np.all(res1==res2) and np.all(chain1==chain2)):
            resi0=resi
            resi=res1
            chain=chain1
            chain[chain=='PROA']='p'
            
            
            
            if resi0 is not None:
                index=np.in1d(resi,resi0)
                resi=resi[index]
                chain=chain[index]
                values=values[index]
              
            if scaling is None:
                scale0=np.max(values)
                scaling=1/scale0
                
            plot_rho(self.sens.molecule,resi,values,chain=chain,\
                     fileout=fileout,scaling=scaling,**kwargs)
                
        else:
            if scaling is None:
                scale0=np.max(values)
                scaling=1/scale0
            
            plot_rho(self.sens.molecule,None,values,scaling=scaling,**kwargs)
#            print('Selections over multiple residues/chains- not currently implemented')
        
        