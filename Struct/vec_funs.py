#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 13:21:49 2019

@author: albertsmith
"""

import MDAnalysis as mda
import numpy as np
from scipy.linalg import svd
import vf_tools as vft
import select_tools as selt


def new_fun(Type,molecule,**kwargs):
    """
    Creates a function to calculate a particular vector(s) from the MD trajectory.
    Mainly responsible for searching the vec_funs file for available functions and
    returning the appropriate function if found (Type determines which function to use)
    
    Required arguments are Type (string specifying the function to be used) and
    a molecule object (contains the MDAnalysis object)
    
    fun=new_fun(Type,molecule,**kwargs)
    
    """
    if Type in globals() and globals()[Type].__code__.co_varnames[0]=='molecule': #Determine if this is a valid vec_fun
        fun0=globals()[Type]
    else:
        raise Exception('Vector function "{0}" was not recognized'.format(Type))
    
    fun=fun0(molecule,**kwargs)
    
    return fun


#%% Frames
"""
Here, we define various functions that define the frames of different motions
in an MD trajectory. Each function should return another function that will
produce one or two vectors defining the frame (without arguments). Those vectors
should have X,Y,Z as the first dimension (for example, such that we can apply
X,Y,Z=v). Note this is the transpose of the outputs of MDanalysis positions
"""    

def peptide_plane(molecule,resids=None,segids=None,filter_str=None):
    "Peptide plane motion, defined by C,N,O positions"
    selN,selC,selO=selt.peptide_plane(molecule,resids,segids,filter_str)
    uni=molecule.mda_object
    def sub():
        box=uni.dimensions[0:3]
        v1=selO.positions-selC.positions
        v2=selN.positions-selC.positions
        v1=vft.pbc_corr(v1.T,box)
        v2=vft.pbc_corr(v2.T,box)
    
        return v1,v2
    return sub
    
def bond(molecule,sel1=None,sel2=None,Nuc=None,resids=None,segids=None,filter_str=None):
    "Bond defines the frame"
    if Nuc is not None:
        sel1,sel2=selt.protein_defaults(Nuc,molecule,resids,segids,filter_str)
    else:
        sel1=selt.sel_simple(sel1,molecule,resids,segids,filter_str)
        sel2=selt.sel_simple(sel2,molecule,resids,segids,filter_str)
    uni=molecule.mda_object
    def sub():
        box=uni.dimensions[0:3]
        v=sel1.positions-sel2.positions
        v=vft.pbc_corr(v.T,box)
        return v
    return sub

def superimpose(molecule,sel=None,resids=None,segids=None,filter_str=None):
    """
    Superimposes a selection of atoms to a reference frame (the first frame)
    
    Note that we may have multiple selections. In this case, then at least some
    of the arguments will be lists or higher dimensional. For this purpose, the
    sel_lists function is used (in select_tools.py)
    
    f=superimpose(molecule,sel=None,resids,None,segids=None,filter_str=None)
    
    f() returns vectors representing the rotation matrix
    """
    
    sel=selt.sel_lists(molecule,sel,resids,segids,filter_str)    
    uni=molecule.mda_object
    "Calculate the reference vectors"
    uni.trajectory.rewind()
    vref=list()
    i0=list()
    for s in sel:
        vr=s.positions
        i0.append(vft.sort_by_dist(vr))
        vref.append(np.diff(vr[i0[-1]],axis=0).T)
       
    def sub():
        R=list()
        box=uni.dimensions[:3]
        for s,vr,i in zip(sel,vref,i0):
            v=vft.pbc_corr(np.diff(s.positions[i],axis=0).T,box)   #Calculate vectors, periodic boundary correction
            R.append(vft.RMSalign(v,vr))    #Get alignment to reference vector
            
        return vft.R2vec(R)     #This converts R back into two vectors
    return sub        
            
            
            
    
    
    