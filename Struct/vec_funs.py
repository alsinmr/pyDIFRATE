#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 13:21:49 2019

@author: albertsmith
"""

import MDAnalysis as mda
import numpy as np

def new_fun(Type,molecule,**kwargs):
    try:
        fun0=globals()[Type]
    except:
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

    sel0=selections(molecule,sel,select)
    
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
    
    if 'sel1' in kwargs or 'sel2' in kwargs or 'Nuc' in kwargs:
        sel01=mol.sel1;sel02=mol.sel2;sel1in0=mol.sel1in;sel2in0=mol.sel2in
        mol.select_atoms(**kwargs) #Run selection if appropriate arguments are provided
        sel1=mol.sel1;sel2=mol.sel2;sel1in=mol.sel1in;sel2in=mol.sel2in;
        mol.sel1=sel01;mol.sel2=sel02;mol.sel1in=sel1in0;mol.sel2in=sel2in0
        """We're resetting molecule here. Maybe we should also implement this elsehwere
        (see rot_axis, moment of inertia)"""
    else:
        sel1=mol.sel1;sel2=mol.sel2;sel1in=mol.sel1in;sel2in=mol.sel2in
    
    if sel is not None:
        if isinstance(sel,str):
            sel=molecule.mda_object.select_atoms(sel)
        elif hasattr(sel,'atoms'):
            sel=sel.atoms

    if select is not None:
        "Apply the additional selection string if provided"
        if sel is not None:
            sel=sel.select_atoms(select)
        sel1=sel1.select_atoms(select)
        sel2=sel2.select_atoms(select)
            
    if sel1in is None:
        sel1in=np.arange(sel1.n_atoms)
    if sel2in is None:
        sel2in=np.arange(sel2.n_atoms)
        
    if np.size(sel1in)!=np.size(sel2in):
        print('Both selections must have the same size')
        return
    
    if sel is None:
        resi=np.zeros(sel1in.shape)
        for k,q in enumerate(sel1in):
            resi[k]=sel1[q]
            
        sel=list()
        resi,index=np.unique(resi,return_index=True) #all residues, and index for each selection
        for r in resi:
            sel.append(mol.mda_object.residues[r].atoms)
    else:
        sel=[sel]
        index=np.zeros(np.shape(sel1in))

    def sub():
        R=list()
        for s in sel:
            UT=sel.principal_axes()
            I=sel.moment_of_inertia()
            i=np.diag(UT.dot(I.dot(UT.T))).argmin()
            vec=UT[i]
            """Now we need to figure out how to rotate this vector to z, and then apply
            the same rotation to all vectors for that selection"""
            
        return

def rot_axis(molecule,axis=[0,0,1],**kwargs):
    """
    Calculates the rotation of a vector (defined by a pair of atoms) around a 
    fixed axis (for example, the z-axis of the MD axis system). Must provide the
    molecule object and the axis. The selection of atoms may be defined here 
    (using the same format as molecule.select_atoms), or if a selection is already
    made with molecule.select_atoms, one may simply give the molecule and axis here
    """    
    
    if 'sel1' in kwargs or 'sel2' in kwargs or 'Nuc' in kwargs:
        molecule.select_atoms(**kwargs) #Run selection if appropriate arguments are provided
    
    axis=np.atleast_1d(axis)
    axis=axis/np.sum(axis**2)
    
    sel1=molecule.sel1
    sel2=molecule.sel2
    if molecule.sel1in is None:
        sel1in=np.arange(sel1.n_atoms)
    else:
        sel1in=molecule.sel1in
    if molecule.sel2in is None:
        sel2in=np.arange(sel2.n_atoms)
    else:
        sel2in=np.molecule.sel2in
    
    def sub():
        vec0=sel1[sel1in].positions-sel2[sel2in].positions
        a=np.dot(vec0,axis)
        vec=vec0-np.dot(np.transpose([a]),[axis])
        return vec.T
    
    return sub
    
def selections(molecule,sel=None,select=None):
    """
    Function to select only some atoms on which to perform special vector calculations
    mda_selection=selections(molecule,sel=None,select=None)
    
    Note that sel ashould be a selection created from an MDanalysis universe, 
    or be a list of MDanalysis selections.
    
    select is a selection string (valid for MDAnalysis) that may be applied after 
    the initial selection in sel
    """    
    
    uni=molecule.mda_object
    sel_out=list()
    
    
    if sel is None:
        if molecule.sel1 is not None:
            sel_out.append(molecule.sel1)
            if select is not None:
                try:
                    sel_out[0]=sel_out[0].atoms.select_atoms(select)
                except:
                    print('select is not a valid selection string')
        else:
            if select is None:
                sel0=uni.atoms
            else:
                try:
                    sel0=uni.select_atoms(select)
                except:
                    print('select is not a valid selection string')
            if np.size(sel0.segments)>1:
                for sel in sel0.segments:
                    sel_out.append(sel.atoms)
            elif np.size(sel0.residues)>1:
                for sel in sel0.residues:
                    sel_out.append(sel.atoms)
            else:
                sel_out.append(sel0.atoms)
    elif isinstance(sel,mda.AtomGroup):
        sel_out.append(sel)
    else:
        if np.size(sel)>1:
            for sel1 in sel:
                sel_out.append(sel1.atoms)
        else:
            try:
                sel_out.append(sel.atoms)
            except:
                print('sel is not a valid MDanalysis atom group')
    
#    if sel is None:
#        if select is None:
#            sel0=uni.atoms
#        else:
#            try:
#                sel0=uni.select_atoms(select)
#            except:
#                print('select is not a valid selection string')
#                return
#        if molecule.sel1 is not None:
#            sel_out.append(molecule.sel1)
#            if select is not None:
#                sel_out[0].atoms.select_atoms(select)
#        elif np.size(sel0.segments)>1:
#            for sel in sel0.segments:
#                sel_out.append(sel.atoms)
#        elif np.size(sel0.residues)>1:
#            for sel in sel0.residues:
#                sel_out.append(sel.atoms)
#        else:
#            sel_out.append(sel0.atoms)
#    elif isinstance(sel,mda.AtomGroup):
#        sel_out.append(sel)
#    else:
#        for sel1 in sel:
#            sel_out.append(sel1)
#
#    """Here we run some cleanup on the selections, in case the user has not properly
#    selected atom groups (for example, if they enumerate segments, these aren't
#    quite atom groups yet, or if they take the whole universe)
#    """
#    for k,sel in enumerate(sel_out):
#        if not isinstance(sel,mda.AtomGroup):
#            if hasattr(sel,'atoms'):
#                sel=sel.atoms
#                sel_out[k]=sel
#
#                
#    if select is not None:
#        try:
#            keep=[True]*np.size(sel_out)
#            for k,sel in enumerate(sel_out):
#                sel_out[k]=sel.select_atoms(select)
#                if sel_out[k].n_atoms==0:
#                    keep[k]=False
#            sel_out=np.ndarray.tolist(np.array(sel_out)[keep])
##            return sel_out,keep
#        except:
#            print('select is not a valid selection string')
#            return


    return sel_out                
        