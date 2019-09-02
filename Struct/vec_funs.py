#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 13:21:49 2019

@author: albertsmith
"""

import MDAnalysis as mda
import numpy as np

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
        print('Vector function "{0}" was not recognized'.format(Type))
        return
    
    fun=fun0(molecule,**kwargs)
    
    return fun
    
def moment_of_inertia(molecule,sel=None,select=None,**kwargs):
    """
    fun=moment_inertia(molecule,sel=None,sel0=None,**kwargs)
    Returns the x, y, and z components of the largest component of the moment 
    of inertia for the given selection.
    
    Note: runs selections function, with arguments sel and select to filter what
    atom groups have the moment of inertial calculated
    """

    sel0=selections(molecule,sel0=sel,select=select)

    def sub():
        vec=list()
        for sel in sel0:
            UT=sel.principal_axes()
            I=sel.moment_of_inertia()
            i=np.diag(UT.dot(I.dot(UT.T))).argmin()
            vec.append([UT[i]])
        return np.concatenate(vec,axis=0).T
    
    return sub


def removeMOI(molecule,sel=None,select=None,**kwargs):
    """
    Factors motion of the moment of inertia out of individual bond motions 
    (effectively removing overall reorienational motion, but not rotational motion)
    By default, factors out motion from all atoms that have the same residue as 
    the first atom in the selection (usually, this should also be the same 
    residue for the second atom)
    
    sel1 and sel2 define the bond direction. These may be previously defined in
    the molecule object, or given here (should be strings that select a set of
    atoms using MDAnalysis style). Nuc may also be used, to select default type of bonds
    
    select is a string that initially filters the selection (for example, 
    select='resname POPC' might pre-select only lipids in a simulation). Argument
    is optional
    
    sel defines a single selection, applied to all bonds, which will define the 
    moment of inertia. sel may be a string, which selects certain atoms in the
    MDanalysis universe, or can be an MDanalysis atom selection (the object itself
    as opposed to the selection string)
    """
    
    mol=molecule

    
    a,b=selections(molecule,**kwargs)
    sel1,sel1in=a
    sel2,sel2in=b
    
    if sel is None:
        resi0=np.zeros(sel1in.shape,dtype=int)
        for k,q in enumerate(sel1in):
            resi0[k]=sel1[q].resid
            
        sel=list()
        resi,index=np.unique(resi0,return_index=True) #all residues, and index for each selection
        for r in resi:
            sel.append(mol.mda_object.residues[r].atoms)
    else:
        a=selections(molecule,sel0=sel,select=select)
        sel,_=a
        index=np.zeros(np.shape(sel1in))
        resi0=np.zeros(sel1in.shape,dtype=int)
        resi=np.zeros(1)

    def sub():
        vec=list()
        for s in sel:
            UT=s.principal_axes()
            I=s.moment_of_inertia()
            i=np.diag(UT.dot(I.dot(UT.T))).argmin()
            vec.append(UT[i])
            
        vec0=sel1[sel1in].positions-sel2[sel2in].positions
        
        vec_out=np.zeros(vec0.shape[-1::-1])
        
        for resi1,vec1 in zip(resi,vec):
            i=resi1==resi0
            vec_out[:,i]=remove(vec0[i],vec)
        
        return vec_out                
    return sub

def rot_axis(molecule,axis=[0,0,1],**kwargs):
    """
    Calculates the rotation of a vector (defined by a pair of atoms) around a 
    fixed axis (for example, the z-axis of the MD axis system). Must provide the
    molecule object and the axis. The selection of atoms may be defined here 
    (using the same format as molecule.select_atoms), or if a selection is already
    made with molecule.select_atoms, one may simply give the molecule and axis here
    """    
    
#    if 'sel1' in kwargs or 'sel2' in kwargs or 'Nuc' in kwargs:
#        molecule.select_atoms(**kwargs) #Run selection if appropriate arguments are provided
#    
#    axis=np.atleast_1d(axis)
#    axis=axis/np.sum(axis**2)
#    
#    sel1=molecule.sel1
#    sel2=molecule.sel2
#    if molecule.sel1in is None:
#        sel1in=np.arange(sel1.n_atoms)
#    else:
#        sel1in=molecule.sel1in
#    if molecule.sel2in is None:
#        sel2in=np.arange(sel2.n_atoms)
#    else:
#        sel2in=np.molecule.sel2in
#   
    if 'sel1' not in kwargs:
        kwargs.update({'sel1':True})
    if 'sel2' not in kwargs:
        kwargs.update({'sel2':True})
    a,b=selections(molecule,**kwargs)
    sel1,sel1in=a
    sel2,sel2in=b
    
    def sub():
        vec0=sel1[sel1in].positions-sel2[sel2in].positions
        vec=project(vec0,axis,Type='plane')
        return vec
    
    return sub

def project(vec0,vec1,Type='plane'):
    """
    Projects a vector (vec0) onto either another vector or a plane (default). If
    projecting onto another vector, vec1 defines that vector. If projecting onto
    a plane, vec1 defines a vector normal to the plane. Set Type='plane' for a
    plane, and Type='vector' for projection onto the vector. vec0 may be an array
    of vectors (Nx3). vec1 may be a single vector, or an array of vectors the
    same size as vec0
    
    vec=project(vec0,vec1,Type='normal')
    """
    
    vec0=np.atleast_2d(vec0)
    if np.ndim(vec1)==1:
        a=np.dot(vec0,vec1)
        para_vec=np.dot(np.transpose([a]),[vec1])
    elif np.dim(vec1)==2 and np.shape(vec1)[0]==np.shape(vec0)[0]:
        a=np.multiply(vec0,vec1).sum(axis=1)
        para_vec=np.multiply(np.transpose([a]),vec1)
    else:
        print('vec1 must have one vector or the same number of vectors as vec0')
        return
    
    if Type.lower()[0]=='p':
        vec=vec0-para_vec
    elif Type.lower()[0]=='v':
        vec=para_vec
        
    return vec.T
    
def remove(vec0,vec1):
    """
    Removes motion of one vector from the motion of the other. That is, we find
    a transformation to rotate vec1 to the z-axis, and apply the same transformation
    to vec0. vec0 may be an array of vectors, and vec1 may be a single vector or
    an array of vectors the same size as vec0

    vec=remove(vec0,vec1)
    """    
    vec0=np.array(vec0)
    vec1=np.array(vec1)
    
    if vec0.shape[0]!=3:
        vec0=vec0.T
    if vec1.shape[0]!=3:
        vec1=vec1.T
    
    vec1=vec1/np.sqrt(np.sum(vec1**2,axis=0))
    
    theta=np.arctan2(vec1[0],vec1[1])
    phi=np.arccos(vec1[2])
    
    vec=np.zeros(vec0.shape)
    
    x=vec0[0]*np.cos(theta)-vec0[1]*np.sin(theta)
    y0=vec0[0]*np.sin(theta)+vec0[1]*np.cos(theta)
    z0=vec0[2]
    z=z0*np.cos(phi)+y0*np.sin(phi)
    y=-z0*np.sin(phi)+y0*np.cos(phi)

    return np.array([x,y,z])

def selections(molecule,sel0=None,sel1=None,sel2=None,select=None,Nuc=None,sel1in=None,sel2in=None):
    """
    General function for returning a set of selections (does not edit the molecule
    object)
    One may return up to three selections, specified by sel0, sel1, and sel2.
    Each argument may be given as a string (valid for atom selection in MDAnalysis),
    as an MDAnalysis selection, or as a list of MDAnalysis selections. sel1 and sel2
    may alternatively be set to True, which will then return the selection already
    stored in molecule.sel1 or molecule.sel2 (depending on if sel1 or sel2 is set
    to True). All selections will be filtered with the selection string (valid
    for MDAnalysis) given in select (if not set to None, which is default)
    
    (s0,s1,s2)=selections(molecule,sel0=None,sel1=None,sel2=None,select=None)
    
    Note, the size of the tuple returned by selections depends on whether sel0,
    sel1, and/or sel2 are set to None (for example, if only sel0 is not None, then
    only one selection is returned)
    
    sel0 is returned as a list of selections (in some cases, multiple selections are needed)
    sel1 and sel2 are returned as single selections, but are a tuple including the
    selection and an index. The index is usually just np.arange(sel1.n_atoms), a
    list numbering from 0 up to the number of atoms in the selection, but if an
    index is provided in mol.sel1in, this may be passed, or one may provide sel1in and
    sel2in directly to selection (simply passed to output)
    
    Finally, note that if sel0, sel1, and sel2 are all None, then a set of selections
    will be returned, where all atoms in each segment or each molecule (if only one
    segment in simulation, will then try residues)
    """
    
    mol=molecule
    
    "Get sel1 and/or sel2 out of molecule in case"
    if isinstance(sel1,bool) and sel1:
        sel1=mol.sel1
        if mol.sel1in is not None:
            sel1in=mol.sel1in
        else:
            sel1in=np.arange(sel1.n_atoms)
    if isinstance(sel2,bool) and sel2:
        sel2=mol.sel2
        if mol.sel2in is not None:
            sel2in=mol.sel2in
        else:
            sel2in=np.arange(sel2.n_atoms)
     
    if Nuc is not None:
        sel1,sel2=nuc_defaults(Nuc)
    
    if sel0 is not None:
        sel0=sel_helper(molecule,sel0,select)
    if sel1 is not None:
        sel1=sel_helper(molecule,sel1,select)
        if not(isinstance(sel1,list)):
            if sel1in is None:
                sel1=(sel1,np.arange(sel1.n_atoms))
            else:
                sel1=(sel1,sel1in)    
    if sel2 is not None:
        sel2=sel_helper(molecule,sel2,select)
        if not(isinstance(sel2,list)):
            if sel2in is None:
                sel2=(sel2,np.arange(sel2.n_atoms))
            else:
                sel2=(sel2,sel2in)   
            
    if sel0 is None and sel1 is None and sel2 is None:
        sel0=sel_helper(molecule,None,select)
        
    if sel0 is None and sel1 is None:
        return sel2
    elif sel0 is None and sel2 is None:
        return sel1
    elif sel1 is None and sel2 is None:
        return sel0
    elif sel0 is None:
        return (sel1,sel2)
    elif sel1 is None:
        return (sel0,sel2)
    elif sel2 is None:
        return (sel0,sel1)
            
def sel_helper(molecule,sel,select):
    uni=molecule.mda_object
    if isinstance(sel,list):
        "If list provided, just run over all values, and drop into a new list"
        sel_out=list()
        for sel1 in sel:
            sel_out.append(sel_helper(molecule,sel1,select))
    elif sel is None:
        "If no selection provided, break up by segments or residues"
        if select is None:
            sel0=uni.atoms
        else:
            sel0=uni.select_atoms(select)
        sel_out=list()
        if sel0.n_segments>1:
            for a in sel0.segments:
                sel_out.append(a.atoms)
        elif sel0.n_residues>1:
            for a in sel0.residues:
                sel_out.append(a.atoms)
        else:
            sel_out.append(sel0)
    elif isinstance(sel,str):
        "Make selection with a string"
        if select is None:
            sel0=uni.atoms
        else:
            sel0=uni.select_atoms(select)
            
        sel_out=sel0.select_atoms(sel)
    elif isinstance(sel.atoms,mda.AtomGroup):
        if select is None:
            sel_out=sel
        else:
            sel_out=sel.select_atoms(select)
            
    return sel_out

def nuc_defaults(Nuc):
    if Nuc.lower()=='15n' or Nuc.lower()=='n' or Nuc.lower()=='n15':
            sel1='(name H or name HN) and around 1.1 name N'
            sel2='name N and around 1.1 (name H or name HN)'
    elif Nuc.lower()=='co':
        sel1='name C and around 1.4 name O'
        sel2='name O and around 1.4 name C'
    elif Nuc.lower()=='ca':
        sel1='name CA and around 1.4 (name HA or name HA2)'
        sel2='(name HA or name HA2) and around 1.4 name CA'
        print('Warning: selecting HA2 for glycines. Use manual selection to get HA1 or both bonds')
        
    return (sel1,sel2)