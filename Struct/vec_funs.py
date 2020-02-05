#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 13:21:49 2019

@author: albertsmith
"""

import numpy as np
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

def return_frame_info(Type=None):
    """
    Provides information as to what frames are available, and what arguments they
    take.
    
    frames=return_frame_info()  Returns list of the frames
    
    args=return_frame_info(Type)    Returns argument list and help info for Type
    """
    
    if Type is None:
        fun_names=list()
        for key in globals().keys():
            if hasattr(globals()[key],'__code__') and \
            globals()[key].__code__.co_varnames[0]=='molecule' and \
            key!='new_fun':
                fun_names.append(key)
        return fun_names
    else:
        if Type in globals():
            return globals()[Type].__code__.co_varnames[1:]
        else:
            print('Frame "{0}" was not recognized'.format(Type))
            return
            
def print_frame_info():
    """
    Prints out some information about the possible frames
    """
    fun_names=return_frame_info()
    print('Implemented frames are:')
    for f in fun_names:
        args=return_frame_info(f)
        print('{0} with arguments {1}'.format(f,args))
    
#%% Frames
"""
Here, we define various functions that define the frames of different motions
in an MD trajectory. Each function should return another function that will
produce one or two vectors defining the frame (without arguments). Those vectors
should have X,Y,Z as the first dimension (for example, such that we can apply
X,Y,Z=v). Note this is the transpose of the outputs of MDanalysis positions
"""    

def peptide_plane(molecule,resids=None,segids=None,filter_str=None,full=True):
    """
    Aligns the peptide plane motion. Two options exist, full=True performs an
    RMS alignment of the N,H,CA of the given residue and C',O,CA of the previous
    residue. 
    full=False uses only the positions of the N of the given residue and C',O
    of the previous.
    
    The former is notably slower, but also performs better when separating
    librational motion
    """
    "Peptide plane motion, defined by C,N,O positions"
    if full:
        "Get selections" 
        selCA,selH,selN,selCm1,selOm1,selCAm1=selt.peptide_plane(molecule,resids,segids,filter_str)
        
        "Get universe, reset time"
        uni=molecule.mda_object
        uni.trajectory.rewind()
        
        "Define function to calculate the vectors defining the plane"
        def vfun():
            v=list()
            for CA,H,N,Cm1,Om1,CAm1 in zip(selCA,selH,selN,selCm1,selOm1,selCAm1):
                v0=np.array([CA.position-N.position,
                            H.position-N.position,
                            N.position-Cm1.position,
                            Cm1.position-Om1.position,
                            Cm1.position-CAm1.position])
                box=uni.dimensions[:3]
                v.append(vft.pbc_corr(v0.T,box))
            return v
        
        "Get the reference vectors (at t=0)"
        vref=vfun()        
        
        def sub():
            R=list()
            vecs=vfun()
            R=[vft.RMSalign(vr,v) for v,vr in zip(vecs,vref)]
            return vft.R2vec(R)
        return sub
    else:
        "Peptide plane motion, defined by C,N,O positions"
        selN,selC,selO=selt.peptide_plane(molecule,resids,segids,filter_str,full)
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
        sel1=selt.sel_simple(molecule,sel1,resids,segids,filter_str)
        sel2=selt.sel_simple(molecule,sel2,resids,segids,filter_str)
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
            R.append(vft.RMSalign(vr,v))    #Get alignment to reference vector
            
        return vft.R2vec(R)     #This converts R back into two vectors
    return sub        
            

def chain_rotate(molecule,sel=None,Nuc=None,resids=None,segids=None,filter_str=None):
    """
    Creates a frame for which a chain of atoms (usually carbons) is aligned
    such that the vector formed by the previous and next heteroatom (not 1H)
    are aligned along z.
    
    Note that the frame is selected with a single initial selection, and the
    function automatically searches for the surrounding atoms. In case a methyl
    carbon is included, the rotation is defined by the carbon itself and its
    nearest neighbor, instead of the surrounding two atoms (which would then
    have to include a methyl proton)
    """

    uni=molecule.mda_object

    "Get the initial selection"
    if Nuc is not None:
        sel,_=selt.protein_defaults(Nuc,molecule,resids,segids,filter_str)
    else:
        sel=selt.sel_simple(molecule,sel,resids,segids,filter_str)
    
    "Get all atoms in the residues included in the initial selection"
    resids=np.unique(sel.resids)
    sel0=uni.residues[np.isin(uni.residues.resids,resids)].atoms
    
    "Get bonded"
    sel1,sel2=selt.find_bonded(sel,sel0=sel0,n=2,sort='cchain')
    
    "Replace 1H with the original selection"
    i=sel2.types=='H'
    
    sel20=sel2
    sel2=uni.atoms[:0]
    for s2,s,i0 in zip(sel20,sel,i):
        if i0:
            sel2+=s
        else:
            sel2+=s2
            
    
    def sub():
        box=uni.dimensions[0:3]
        v=sel2.positions-sel1.positions
        v=vft.pbc_corr(v.T,box)
        return v
    return sub
    
def methylCC(molecule,Nuc=None,resids=None,segids=None,filter_str=None):
    """
    Superimposes the C-X bond attached to a methyl carbon, and can separate
    methyl rotation from re-orientation of the overall methyl group
    
    Note- we only return one copy of the C–C bond, so a frame index is necessary
    """            
    if not(Nuc.lower()=='ivl' or Nuc.lower()=='ivla' or Nuc.lower()=='ch3'):
        print('Nuc must be ivl, ivla, or ch3 for the methylCC frame')
        return
    elif Nuc is None:
        Nuc='ch3'
    selC1,_=selt.protein_defaults(Nuc,molecule,resids,segids,filter_str)  
    selC1=selC1[::3]    #Above line returns 3 copies of each carbon      
    
    sel0=molecule.mda_object.atoms
    sel0=sel0.residues[np.isin(sel0.residues.resids,selC1.resids)].atoms
    
    selC2=sum([sel0.select_atoms('not name H* and around 1.6 atom {0} {1} {2}'\
                                 .format(s.segid,s.resid,s.name)) for s in selC1])
    
    def sub():
        box=molecule.mda_object.dimensions[:3]
        v=selC1.positions-selC2.positions
        v=vft.pbc_corr(v.T,box)
        return v
    return sub

def librations(molecule,sel1=None,sel2=None,Nuc=None,resids=None,segids=None,filter_str=None,full=True):
    """
    Defines a frame for which librations are visible. That is, for a given bond,
    defined by sel1 and sel2, we search for two other atoms bound to the 
    heteroatom (by distance). The reference frame is then defined by the 
    heteroatom and the additional two atoms, leaving primarily librational
    motion of the bond. We preferentially select the other two atoms for larger
    masses, but they may also be protons (for example, a methyl H–C bond will 
    be referenced to the next carbon but also another one of the protons of 
    the methyl group)
    
    In case the heteroatom only has two bound partners, the second atom in the
    bond will also be used for alignment, reducing the effect motion
    (not very common in biomolecules)
    
    librations(sel1,sel2,Nuc,resids,segids,filter_str)
    """
    if Nuc is not None:
        sel1,sel2=selt.protein_defaults(Nuc,molecule,resids,segids,filter_str)
    else:
        sel1=selt.sel_simple(molecule,sel1,resids,segids,filter_str)
        sel2=selt.sel_simple(molecule,sel2,resids,segids,filter_str)
        
    if sel1.masses.sum()<sel2.masses.sum():
        sel1,sel2=sel2,sel1 #Make sel1 the heteroatom
    
    resids=np.unique(sel1.resids)
    i=np.isin(sel1.universe.residues.resids,resids)    #Filter for atoms in the same residues
    sel0=sel1.universe.residues[i].atoms
    if full:
        "Slightly improved performance if we align all 4 bonds to the carbon"
        "Note, "
        sel2,sel3,sel4,sel5=selt.find_bonded(sel1,sel0,n=4,sort='mass')
        
        def vfun():
            v=list()
            for v1,v2,v3,v4,v5 in zip(sel1,sel2,sel3,sel4,sel5):
                v0=np.array([v2.position-v1.position,
                            v3.position-v1.position,
                            v4.position-v1.position,
                            v5.position-v1.position])
                box=uni.dimensions[:3]
                v.append(vft.pbc_corr(v0.T,box))
            return v
        
        uni=molecule.mda_object
        uni.trajectory.rewind()
        
        vref=vfun()
        
        def sub():
            R=list()
            vecs=vfun()
            R=[vft.RMSalign(vr,v) for v,vr in zip(vecs,vref)]
            return vft.R2vec(R)
    else:
        sel2,sel3=selt.find_bonded(sel1,sel0,n=2,sort='mass')
        
        uni=molecule.mda_object
        def sub():
            box=uni.dimensions[0:3]
            v1=sel2.positions-sel1.positions
            v2=sel1.positions-sel3.positions
            v1=vft.pbc_corr(v1.T,box)
            v2=vft.pbc_corr(v2.T,box)
            return v1,v2
    return sub
    
def librations0(molecule,sel1=None,sel2=None,Nuc=None,resids=None,segids=None,filter_str=None):
    """
    Defines a frame for which librations are visible. That is, for a given bond,
    defined by sel1 and sel2, we search for two other atoms bound to the 
    heteroatom (by distance). The reference frame is then defined by the 
    heteroatom and the additional two atoms, leaving primarily librational
    motion of the bond. We preferentially select the other two atoms for larger
    masses, but they may also be protons (for example, a methyl H–C bond will 
    be referenced to the next carbon but also another one of the protons of 
    the methyl group)
    
    In case the heteroatom only has two bound partners, the second atom in the
    bond will also be used for alignment, reducing the effect motion
    (not very common in biomolecules)
    
    librations(sel1,sel2,Nuc,resids,segids,filter_str)
    """
    if Nuc is not None:
        sel1,sel2=selt.protein_defaults(Nuc,molecule,resids,segids,filter_str)
    else:
        sel1=selt.sel_simple(molecule,sel1,resids,segids,filter_str)
        sel2=selt.sel_simple(molecule,sel2,resids,segids,filter_str)
        
    if sel1.masses.sum()<sel2.masses.sum():
        sel1,sel2=sel2,sel1 #Make sel1 the heteroatom
    
    resids=np.unique(sel1.resids)
    i=np.isin(sel1.universe.residues.resids,resids)    #Filter for atoms in the same residues
    sel0=sel1.universe.residues[i].atoms
    sel2,sel3,sel4,sel5=selt.find_bonded(sel1,sel0,n=4,sort='mass')
    
    def vfun():
        v=list()
        for v1,v2,v3,v4,v5 in zip(sel1,sel2,sel3,sel4,sel5):
            v0=np.array([v2.position-v1.position,
                        v3.position-v1.position,
                        v4.position-v1.position,
                        v5.position-v1.position])
            box=uni.dimensions[:3]
            v.append(vft.pbc_corr(v0.T,box))
        return v
    
    uni=molecule.mda_object
    uni.trajectory.rewind()
    
    vref=vfun()
    
    def sub():
        R=list()
        vecs=vfun()
        R=[vft.RMSalign(vr,v) for v,vr in zip(vecs,vref)]
        return vft.R2vec(R)
        
    return sub        