#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Copyright 2021 Albert Smith-Penzel

This file is part of pyDIFRATE

pyDIFRATE is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

pyDIFRATE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with pyDIFRATE.  If not, see <https://www.gnu.org/licenses/>.


Questions, contact me at:
albert.smith-penzel@medizin.uni-leipzig.de



Created on Wed Mar 31 09:30:03 2021

@author: albertsmith
"""

"""
Module for calculating various relaxation rate constants as a function of 
correlation time (tc). Amplitudes are corresponding to (1-S2)=1. 

Some rules for adding functions:
    1) The first argument is always tc for sensitivity functions
    2) Functions that you don't want to show up in the sens_obj as sensitivity
        functions should start with _ . (these may have any initial argument)
    3) Try not to make new variables for every new function. Use a constant set
        of variable names wherever possible
"""

from tools.DRtools import NucInfo
import numpy as np
import os
import pandas as pd
from plots.plotting_funs_new import plot_rhoz
        

#%% These functions handle setting up of the info object: sorting, repeats, defaults
def _defaults(info,**kwargs):
    """
    Populates the info array with default values, followed by input values    
    """
    ne=np.size(list(kwargs.values())[0])
    keys=info.T.columns.to_numpy()
    pd_dict={k:([None for _ in range(ne)] if (k=='Nuc' or k=='Nuc1') else np.zeros(ne)) for k in keys}
        
    info_new=info.__class__(pd_dict).T
    
    "End-of-file function"
    def eof(f):
        "Determines if we are at the end of the file"
        pos=f.tell()    #Current position in the file
        f.readline()    #Read out a line
        if pos==f.tell(): #If position unchanged, we're at end of file
            return True
        else:       #Otherwise, reset pointer, return False
            f.seek(pos)
            return False
    
    "Load the defaults into a dictionary"
    default_file=os.path.join(os.path.dirname(os.path.abspath(__file__)),'NMR_defaults.txt')
    defaults=dict()
    current=dict()
    with open(default_file,'r') as f:
        reading_default=False
        while not(eof(f)):
            a=f.readline().strip()
            if reading_default:   
                if a=='END DEFAULT':    #End of current default settings reached
                    reading_default=False
                else:
                    k,v=a.split(maxsplit=1) #Split into variable name and values
                    try:
                        v=np.array(v.split(),dtype=float)   #Split values if more than one
                    except:
                        v=np.array(v.split())   #If conversion to float fails, leave as string
                    if v.size==1:v=v[0] #Don't keep as an array if only one value
                    current[k]=v    #Store result
            else:
                if a=='BEGIN DEFAULT':
                    reading_default=True #Start reading a set of defaults
                    current=dict()
                    keys=f.readline().strip().split()
                    defaults.update({k:current for k in keys})
    
    "Now we populate our dataframe with the defaults"     
    if 'Nuc' in kwargs:
        Nuc=np.atleast_1d(kwargs.pop('Nuc'))
        for k,Nuc0 in enumerate(Nuc):
            if Nuc0 in defaults.keys():
                for key,value in defaults[Nuc0].items():
                    try:
                        info_new.at[key,k]=value
                    except:
                        print(key,k,value)
    
    "Finally, we override the defaults with our input values"
    for key,values in kwargs.items():
        for k,value in enumerate(np.atleast_1d(values)):
            info_new[k][key]=value

    return pd.concat([info,info_new],axis=1,ignore_index=True)


def _repeats(**kwargs):
    """
    Handles the treatment of arguments where multiple experiments may be specified,
    for example, where multiple experiments may be at v0=X, but each has a 
    different MAS speed.    
    """
    "Figure out how many experiments there are"
    ne=1
    for k,val in kwargs.items():
        if k!='dXY' and k!='Nuc1':
            ne=np.max([ne,np.array(val).size])
    
    out=dict()
    for k,val in kwargs.items():
        if k=='dXY' or k=='Nuc1':
            out[k]=np.array([val]).repeat(ne) if np.size(val)==1 else np.atleast_2d(val).repeat(ne,axis=0)
        else:
            rep=np.ceil(ne/np.size(val))
            val1=np.repeat(val,rep)[:ne]
            out[k]=val1
            
    return out

def _sorting():
    """
    Specifies the sorting of variables used in NMR_funs. Arguments found in this
    file that are not found
    """
    
    return ('Type','v0','v1','vr','offset','stdev','Nuc','Nuc1','dXY','CSA','eta','CSoff','QC','etaQ')

#%% This function used for adding some plotting options into sensitivity classes
"""
We can do this for importing arbitrary plotting functions into the different
object types. Might vary, for example, for a correlation function. However, first
argument should be the 'self' argument in all instances (here called obj)
"""

def _setup_functions():
    "Defines the functions to be added to the NMR sensitivity object"    
    def plot_R(self,exp_num=None,norm=False,ax=None,**kwargs):
        """
        Plots the sensitivites of the experiments. Default plots all experiments 
        without normalization. Set norm=True to normalize all experiments to 1. 
        Specify exp_num to only plot selected experiments. Set ax to specify the
        axis on which to plot
        
        plot_R(exp_num=None,norm=False,ax=None,**kwargs)
        """
        hdl=plot_rhoz(self,index=exp_num,norm=norm,ax=ax,**kwargs)
        ax=hdl[0].axes
        ax.set_ylabel(r'$R / s^{-1}$')
        ax.set_title('Experimental Sensitivities')
        return hdl
        "Defines functions for calling sensitivities of the NMR type"
    def R(self):
        return self._rho()     
    return {'plot_R':plot_R,'R':R}


#%% Functions for calculation of relaxation
def J(tc,v0):
    """
    Returns the spectral density at frequency v for the correlation time axis tc
    """
    return 2/5*tc/(1+(2*np.pi*v0*tc)**2)

def R1(tc,Nuc,v0,Nuc1=None,CSA=0,dXY=0,eta=0,vr=0,CSoff=0,QC=0,etaQ=0):
    """
    Returns the T1 (as R1) relaxation rate constant. Sources of relaxation are:
        
    Quadrupole coupling: Provide Nuc, v0, QC, and etaQ
    Dipole coupling (heteronuclear): Provide Nuc, Nuc1, v0, dXY
    Dipole coupling (homonuclear): Provide Nuc, Nuc1, v0, dXY, vr, CSoff
    CSA: Provide Nuc, v0, CSA, eta
    
    CSA and quadrupole relaxation use the same parameter, eta
    
    Defaults of most parameters are 0 (Nuc1 is None). 
    Multiple dipole couplings may be considered. Provide dXY and Nuc1 as a list
    
    All provided contributions will be included in the total rate constant.
    """

    v0=np.array(v0)*1e6     #1H resonance frequency (convert MHz to Hz)
    vr=np.array(vr)*1e3     #MAS frequency (convert kHz to Hz)
    dXY=np.atleast_1d(dXY)
    Nuc1=np.atleast_1d(Nuc1)
    assert Nuc1.size==dXY.size,"Nuc1 and dXY must have the same size"
    vX=NucInfo(Nuc)/NucInfo('1H')*v0
    CSA=CSA*vX/1e6
    Delv=np.array(CSoff)*vX/1e6
    R=np.zeros(tc.shape)

    "Dipole relaxation"
    for N1, dXY1 in zip(Nuc1,dXY):
        if N1 is not None:
            vY=NucInfo(N1)/NucInfo('1H')*v0
            S=NucInfo(N1,'spin')
            sc=S*(S+1)*4/3 # Scaling factor depending on the spin, =1 for spin 1/2            
            if vX==vY:  #Homonuclear
                R+=sc*(np.pi*dXY1/2)**2*(1/6*J(tc,Delv+2*vr)+1/6*J(tc,Delv-2*vr)\
                   +1/3*J(tc,Delv+vr)+1/3*J(tc,Delv-vr)+3*J(tc,vX)+6*J(tc,2*vX))
            else:       #Heteronuclear
                R+=sc*(np.pi*dXY1/2)**2*(J(tc,vX-vY)+3*J(tc,vX)+6*J(tc,vY+vX))
    
    print(R)
    "Quadrupole Relaxation"
    """
    Note that we calculate the orientationally averaged, initial relaxation
    rate constant. We do not account for multi-exponentiality as occuring from
    either different orientations or from relaxation through multiple spin
    states.
    
    Use R1Q for relaxation of quadrupolar order
    """
    S=NucInfo(Nuc,'spin')     
    if S>=1:      
        deltaQ=1/(2*S*(2*S-1))*QC*2*np.pi
        C=(deltaQ/2)**2*(1+etaQ**2/3) #Constant that scales the relaxation        
        if S==1:
            R+=C*(3*J(tc,vX)+12*J(tc,2*vX))
        elif S==1.5:
            R+=C*(36/5*J(tc,vX)+144/5*J(tc,2*vX))
        elif S==2.5:
            R+=C*(96/5*J(tc,vX)+384/5*J(tc,2*vX))
        else:
            print('Spin={0} not implemented for quadrupolar relaxation'.format(S))
    "CSA relaxation"
    R+=3/4*(2*np.pi*CSA)**2*(1+eta**2/3)*J(tc,vX)
    return R

def plot_R(self,exp_num=None,norm=False,ax=None,**kwargs):
    """
    Plots the sensitivites of the experiments. Default plots all experiments 
    without normalization. Set norm=True to normalize all experiments to 1. 
    Specify exp_num to only plot selected experiments. Set ax to specify the
    axis on which to plot
    
    plot_R(exp_num=None,norm=False,ax=None,**kwargs)
    """
    hdl=plot_rhoz(self,index=exp_num,norm=norm,ax=ax,**kwargs)
    ax=hdl[0].axes
    ax.set_ylabel(r'$R / s^{-1}$')
    ax.set_title('Experimental Sensitivities')
    return hdl
