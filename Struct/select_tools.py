#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 13:06:44 2019

@author: albertsmith
"""

"""
Library of selection tools, to help define the selections for correlation
function calculation, frame definition, etc.
"""
import MDAnalysis as mda
import numpy as np
import numbers

  
def sel0_filter(mol,resids=None,segids=None,filter_str=None):
    """
    Performs initial filtering of all atoms in an MDA Universe. Filtering may
    be by resid, segid, and/or a selection string. Each selector will default
    to None, in which case it is not applied
    
    sel0=sel0_filter(mol,resids,segids,filter_str)
    """
    if hasattr(mol,'mda_object'):
        sel0=mol.mda_object.atoms
    elif hasattr(mol,'atoms'):
        sel0=mol.atoms
    else:
        print('mol needs to be a molecule object or an atom group')
        return
    
    if segids is not None:
        segids=np.atleast_1d(segids)
        i=np.isin(sel0.segments.segids,segids)
        sel_si=sel0.segments[np.argwhere(i).squeeze()].atoms
        sel0=sel0.intersection(sel_si)
    if resids is not None:
        resids=np.atleast_1d(resids)
        i=np.isin(sel0.residues.resids,resids)
        sel_ri=sel0.residues[np.argwhere(i).squeeze()].atoms
        sel0=sel0.intersection(sel_ri)

    if filter_str is not None:
        sel_fs=sel0.select_atoms(filter_str)
        sel0=sel0.intersection(sel_fs)

    return sel0

#%% Simple selection 

def sel_simple(mol,sel=None,resids=None,segids=None,filter_str=None):
    """
    Takes a selection out of the molecule object, where that selection may
    be an atom group produced by the user, may be a selection string, or
    may be extracted from mol.sel1/mol.sel2. In each case, the output is a MDAnalysis
    atom group.
    
    One may further filter the given group with a valid MDAnalysis selection
    string, with a specification of specific residues or specific segments
    
    sel = sel_simple(sel,mol=None,resids,segids,filter_str=None)
    
    set 'sel' to be a string, an atom group, or simply 1 or 2 to use the selections
    stored in mol
    """
    
    """If sel has atoms as an attribute, we make sure it's an atom group
    (could be a residue group or universe, for example)
    """
    if hasattr(mol,'atoms'): sel=mol.atoms 
    
    
    if sel is None:
        if not(isinstance(sel,mda.AtomGroup)):
            print('If the molecule object is not provided, then sel must be an atom group')
            return
        sel=sel0_filter(mol,resids,segids,filter_str)
        return sel
    
    if isinstance(sel,str):
        sel0=sel0_filter(mol,resids,segids,filter_str)
        sel=sel0.select_atoms(sel)
    elif isinstance(sel,numbers.Real) and sel==1:
        sel=sel0_filter(mol.sel1,resids,segids,filter_str)
    elif isinstance(sel,numbers.Real) and sel==2:
        sel=sel0_filter(mol.sel2,resids,segids,filter_str)
    elif isinstance(sel,mda.AtomGroup):
        sel=sel0_filter(sel,resids,segids,filter_str)
    else:
        print('sel is not an accepted data type')
        return
        
    return sel


def sel_lists(mol,sel=None,resids=None,segids=None,filter_str=None):
    """
    Creates multiple selections from single items or lists of sel, resids,
    segids, and filter_str.
    
    Each argument (sel,resids,segids,filter_str) may be None, may be a single
    argument (as for sel_simple), or may be a list of arguments. If more than
    one of these is a list, then the lists must have the same length. Applies
    sel_simple for each item in the list. The number of selections returns is 
    either one (no lists used), or the length of the lists (return will always
    be a list)
    
    sel_list=sel_lists(mol,sel=None,resids=None,segids=None,filter_str=None)
    """

    "First apply sel, as a single selection or list of selections"
    if hasattr(sel,'atoms') or isinstance(sel,str) or sel==1 or sel==2:
        sel=sel_simple(mol,sel)
        n=1
    elif isinstance(sel,list):
        sel=[sel_simple(mol,s) for s in sel]
        n=len(sel)
    elif sel is None:
        sel=mol.mda_object.atoms
        n=1
    else:
        print('sel data type was not recognized')
        return
    
    "Apply the resids filter"
    if resids is not None:
        if hasattr(resids,'__iter__') and hasattr(resids[0],'__iter__'):
            if n==1:
                n=len(resids)
                sel=[sel_simple(sel,resids=r) for r in resids]
            elif len(resids)==n:
                sel=[sel_simple(s,resids=r) for s,r in zip(sel,resids)]
            else:
                print('Inconsistent sizes for selections (resids)')
        else:
            if n==1: 
                sel=sel_simple(sel,resids=resids) 
            else:
                sel=[sel_simple(s,resids=resids) for s in sel]
            
    "Apply the segids filter"
    if segids is not None:
        if not(isinstance(segids,str)) and hasattr(segids,'__iter__') and hasattr(segids[0],'__iter__'):
            if n==1:
                n=len(segids)
                sel=[sel_simple(sel,segids=si) for si in segids]
            elif len(segids)==n:
                sel=[sel_simple(s,segids=si) for s,si in zip(sel,segids)]
            else:
                print('Inconsistent sizes for selections (segids)')
        else:
            if n==1: 
                sel=sel_simple(sel,segids=segids) 
            else:
                sel=[sel_simple(s,segids=segids) for s in sel]
                
    "Apply the filter_str"
    if filter_str is not None:
        if np.ndim(filter_str)>0:
            if n==1:
                n=len(filter_str)
                sel=[sel_simple(sel,filter_str=f) for f in filter_str]
            elif len(filter_str)==n:
                sel=[sel_simple(s,filter_str=f) for s,f in zip(sel,filter_str)]
            else:
                print('Inconsistent sizes for selections (filter_str)')
        else:
            if n==1:
                sel=sel_simple(sel,filter_str=filter_str)
            else:
                sel=[sel_simple(s,filter_str=filter_str) for s in sel]
                
    if n==1:
        sel=[sel]
        
    return sel

#%% Specific selections for proteins
def protein_defaults(Nuc,mol,resids=None,segids=None,filter_str=None):
    """
    Selects pre-defined pairs of atoms in a protein, usually based on nuclei that
    are observed for relaxation. One may also select specific residues, specific
    segments, and apply a filter string
    
    sel1,sel2=protein_defaults(Nuc,mol,resids,segids,filter_str)
    
    Nuc is a string and can be:
    N (15N,n,n15,N15), CA (13CA,ca,ca13,CA13), C (CO, 13CO, etc.)
    """
    
    sel0=sel0_filter(mol,resids,segids,filter_str)
        
    if Nuc.lower()=='15n' or Nuc.lower()=='n' or Nuc.lower()=='n15':       
        sel1=sel0.select_atoms('name N and around 1.1 (name H or name HN)')                 
        sel2=sel0.select_atoms('(name H or name HN) and around 1.1 name N')        
    elif Nuc.lower()=='co' or Nuc.lower()=='13co' or Nuc.lower()=='co13' or Nuc.lower()=='c':
        sel1=sel0.select_atoms('name C and around 1.4 name O')
        sel2=sel0.select_atoms('name O and around 1.4 name C')
    elif Nuc.lower()=='ca' or Nuc.lower()=='13ca' or Nuc.lower()=='ca13':
        sel1=sel0.select_atoms('name CA and around 1.5 (name HA or name HA2)')
        sel2=sel0.select_atoms('(name HA or name HA2) and around 1.5 name CA')
        print('Warning: selecting HA2 for glycines. Use manual selection to get HA1 or both bonds')
    elif Nuc[:3].lower()=='ivl' or Nuc[:3].lower()=='ch3':
        if Nuc[:4].lower()=='ivla':
            fs0='resname ILE Ile ile VAL val Val LEU Leu leu ALA Ala ala'
            Nuc0=Nuc[4:]
        elif Nuc[:3].lower()=='ivl':
            fs0='resname ILE Ile ile VAL val Val LEU Leu leu'
            Nuc0=Nuc[3:]
        else:
            fs0=None
            Nuc0=Nuc[3:]
        filter_str=filter_str if fs0 is None else (fs0 if filter_str is None else \
                            '('+filter_str+') and ('+fs0+')')
        select=None
        if 't' in Nuc0.lower() or 'l' in Nuc0.lower():select='l'
        if 'r' in Nuc0.lower():select='r'
        
        sel1,sel2=find_methyl(mol,resids,segids,filter_str,select=select)
        
#        
#        if Nuc[-1].lower()=='t' or Nuc[-2].lower()=='t':    #Truncated list- only one C per residue
#            sel0=sel0-sel0.select_atoms('(resname VAL val Val and name CG2) or \
#                                         (resname ILE ile Ile and name CG2) or \
#                                         (resname LEU leu Leu and name CD1)')
#        if Nuc[-1].lower()=='r' or Nuc[-1].lower()=='l':
#            sel0=sel0-sel0.select_atoms('(resname VAL val Val and name CG{0}) or \
#                                         (resname ILE ile Ile and name CG2) or \
#                                         (resname LEU leu Leu and name CD{0})'\
#                                            .format('2' if Nuc[-1].lower()=='l' else '1'))
#            
#        if Nuc[:4].lower()=='ivla':
#            sel0C=sel0.select_atoms('resname ILE Ile ile VAL val Val LEU Leu leu ALA Ala ala and name C*')
#            sel0H=sel0.select_atoms('resname ILE Ile ile VAL val Val LEU Leu leu ALA Ala ala and name H*')
#        elif Nuc[:3].lower()=='ivl':
#            sel0C=sel0.select_atoms('resname ILE Ile ile VAL val Val LEU Leu leu and name C*')
#            sel0H=sel0.select_atoms('resname ILE Ile ile VAL val Val LEU Leu leu and name H*')
#        else:
#            sel0C=sel0.select_atoms('resname ILE Ile ile VAL val Val LEU Leu leu ALA Ala ala MET Met met THR Thr thr and name C*')
#            sel0H=sel0.select_atoms('resname ILE Ile ile VAL val Val LEU Leu leu ALA Ala ala MET Met met THR Thr thr and name H*')
#        ids=list()
#        for s in sel0C:
#            if (sel0H+sel0C).select_atoms('name H* and around 1.15 atom {0} {1} {2}'.format(s.segid,s.resid,s.name)).n_atoms==3:
#                ids.append(s.id)
#        sel1=sel0[np.isin(sel0.ids,ids)]
#        sel1=sel1[np.repeat([np.arange(sel1.n_atoms)],3,axis=1).reshape(sel1.n_atoms*3)]
#        sel2=(sel1+sel0H).select_atoms('name H* and around 1.15 name C*')
        
        if '1' in Nuc:
            sel1=sel1[::3]
            sel2=sel2[::3]
          
    return sel1,sel2

def find_methyl(mol,resids=None,segids=None,filter_str=None,select=None):
    """
    Finds methyl groups in a protein for a list of residues. Standard selection
    options are used. 
    
    select may be set to 'l' or 'r' (left or right), which will select one of the
    two methyl groups on valine or leucine, depending on their stereochemistry. In
    this mode, only the terminal isoleucine methyl group will be returned.
    
    To just get rid of the gamma methyl on isoleucine, set select to 'ile_d'
    """
    mol.mda_object.trajectory[0]
    sel0=sel0_filter(mol,resids,segids,filter_str)
    selC0,selH0=sel0.select_atoms('name C*'),sel0.select_atoms('name H*')
    index=np.array([all(b0 in selH0 for b0 in b)\
           for b in np.array(find_bonded(selC0,selH0,n=3,d=1.5,sort='massi')).T])
#    index=np.all([b.names[0]=='H' for b in find_bonded(selC0,selH0,n=3,d=1.5)],axis=0)
    selH=find_bonded(selC0[index],sel0=selH0,n=3,d=1.5)

    selH=np.array(selH).T
    selC=selC0[index]
    
    "First, we delete the gamma of isoleucine if present"
    if select is not None:
        ile=[s.resname.lower()=='ile' for s in selC] #Find the isoleucines
        if any(ile):
            exclude=[s.sum() for s in selH[ile]]
            nxt=find_bonded(selC[ile],sel0=sel0,exclude=exclude,n=1,sort='cchain')[0]
            keep=np.array([np.sum([b0.name[0]=='H' for b0 in b])==2 \
                    for b in np.array(find_bonded(nxt,sel0=sel0,exclude=selC[ile],n=2,sort='massi')).T])
    #        keep=np.sum([b.types=='H' for b in find_bonded(nxt,sel0=sel0,exclude=selC[ile],n=2)],axis=0)==2
            index=np.ones(len(selC),dtype=bool)
            index[ile]=keep
            selC,selH=selC[index],selH[index]
    
    if select is not None and (select[0].lower() in ['l','r']):
        val_leu=[s.resname.lower()=='val' or s.resname.lower()=='leu' for s in selC]
        if any(val_leu):
            exclude=[s.sum() for s in selH[val_leu][::2]]
            nxt0=find_bonded(selC[val_leu][::2],sel0=sel0,exclude=exclude,n=1,sort='cchain')[0]
            exclude=np.array([selC[val_leu][::2],selC[val_leu][1::2]]).T
            exclude=[e.sum() for e in exclude]
            nxt1=find_bonded(nxt0,sel0=sel0,exclude=exclude,n=1,sort='cchain')[0]
            nxtH=find_bonded(nxt0,sel0=sel0,exclude=exclude,n=1,sort='massi')[0]
            
            cross=np.cross(nxtH.positions-nxt0.positions,nxt1.positions-nxt0.positions)
            dot0=(cross*selC[val_leu][::2].positions).sum(1)
            dot1=(cross*selC[val_leu][1::2].positions).sum(1)
            keep=np.zeros(np.sum(val_leu),dtype=bool)
            keep[::2]=dot0>=dot1
            keep[1::2]=dot0<dot1
            if select[0].lower()=='l':keep=np.logical_not(keep)
            
            index=np.ones(len(selC),dtype=bool)
            index[val_leu]=keep
            selC,selH=selC[index],selH[index]
        
    selH=np.concatenate(selH).sum()
    selC=np.sum(np.repeat(selC,3))    
    
    return selC,selH
    
    

def find_bonded(sel,sel0=None,exclude=None,n=3,sort='dist',d=1.65):
    """
    Finds bonded atoms for each input atom in a given selection. Search is based
    on distance. Default is to define every atom under 1.65 A as bonded. It is 
    recommended to also provide a second selection (sel0) out of which to search
    for the bound atoms. If not included, the full MD analysis universe is searched.
    
    Note- a list of selections is returned. Sorting may be determined in one
    of several ways (set sort)
        'dist':     Sort according to the nearest atoms
        'mass':     Sort according to the largest atoms first
        'massi':    Sort according to smallest atoms first (H first)
        'cchain':   Sort, returing C atoms preferentially (followed by sorting by mass)
    
    One may also exclude a set of atoms (exclude), which then will not be returned
    in the list of bonded atoms. Note that exclude should be a list the same
    size as sel (either a selection the same size as sel, or a list of selections
    with a list length equal to the number of atoms in sel)
    """
    
    if not(hasattr(sel,'__len__')):sel=[sel]
    
    out=[sel[0].universe.atoms[0:0] for _ in range(n)]  #Empty output
    
    if sel0 is None:
        sel0=sel[0].universe
    
    for m,s in enumerate(sel):
        sel01=sel0.select_atoms('point {0} {1} {2} {3}'.format(*s.position,d))
        sel01=sel01-s #Exclude self
        if exclude is not None:
            sel01=sel01-exclude[m]
        if sort[0].lower()=='d':
            i=np.argsort(((sel01.positions-s.position)**2).sum(axis=1))
        elif sort[0].lower()=='c':
            C=np.array([s.name[0]=='C' for s in sel01])
#            C=sel01.type=='C'
            nC=np.logical_not(C)
            i1=np.argsort(sel01[nC].masses)[::-1]
            C=np.argwhere(C)[:,0]
            nC=np.argwhere(nC)[:,0]
            i=np.concatenate((C,nC[i1]))
        elif sort.lower()=='massi':
            i=np.argsort(sel01.masses)
        else:
            i=np.argsort(sel01.masses)[::-1]
        sel01=sel01[i]
        for k in range(n):
            if len(sel01)>k:
                out[k]+=sel01[k]
            else:
                out[k]+=s
                
    return out        
        
    
#%% This allows us to use a specific keyword to make an automatic selection
"""
Mainly for convenience, cleanliness in code construction
To add a new keyword, simply define a function of the same name, that returns
the desired selections.
Note that mol must always be an argument (the molecule object)
resids,segids,and filter_str must also be arguments, or **kwargs must be included
"""
def keyword_selections(keyword,mol,resids=None,segids=None,filter_str=None,**kwargs):
    if keyword in globals() and globals()[keyword].__code__.co_varnames[0]=='mol': #Determine if this is a valid vec_fun
        fun0=globals()[keyword]
    else:
        raise Exception('Keyword selection "{0}" was not recognized'.format(keyword))
    
    fun=fun0(mol=mol,resids=resids,segids=segids,filter_str=filter_str,**kwargs)
    
    return fun

def peptide_plane(mol,resids=None,segids=None,filter_str=None,full=True):
    """
    Selects the peptide plane. One may also provide resids, segids,
    and a filter string. Note that we define the residue as the residue containing
    the N atom (whereas the C, O, and one Ca of the same peptide plane are actually in
    the previous residue).
    
    returns 6 selections:
    selCA,selH,selN,selCm1,selOm1,selCAm1   
    (selCA, selH, and selN are from residues in resids, and 
    selCm1, selOm1, selCAm1 are from residues in resids-1)
    
    or if full = False, returns 3 selections
    selN,selCm1,selOm1
    
    Note that peptide planes for which one of the defining atoms is missing will
    be excluded
    """
    sel0=sel0_filter(mol,resids,segids,filter_str)
    if resids is None:
        resids=sel0.resids
    selm1=sel0_filter(mol,np.array(resids)-1,segids,filter_str)
    
    if full:
#        selN=(sel0.union(selm1)).select_atoms('protein and (name N and around 1.5 name HN H CD) and (around 1.4 (name C and around 1.4 name O))')
        selN=(sel0.union(selm1)).select_atoms('protein and (name N and around 1.7 name HN H CD) and (around 1.7 (name C and around 1.7 name O))')
    else:  #We don't need the HN to be present in this case  
#        selN=(sel0.union(selm1)).select_atoms('protein and (name N and around 1.4 (name C and around 1.4 name O))')
        selN=(sel0.union(selm1)).select_atoms('protein and (name N and around 1.7 (name C and around 1.7 name O))')

    i=np.isin(selN.resids,resids)
    selN=selN[i]    #Maybe we accidently pick up the N in the previous plane? Exclude it here
    resids=selN.resids
    "Re-filter the original selection for reduced resid list"
    sel0=sel0_filter(sel0,resids)
    selm1=sel0_filter(selm1,np.array(resids)-1)
    if full:
#        selH=sel0.residues.atoms.select_atoms('protein and (name H HN CD and around 1.5 name N)')
        selH=sel0.residues.atoms.select_atoms('protein and (name H HN CD and around 1.7 name N)')
        selCA=sel0.residues.atoms.select_atoms('protein and (name CA and around 1.7 name N)')
    
#    i=np.argwhere(np.isin(sel0.residues.resids,sel1.residues.resids-1)).squeeze()
#    selCm1=selm1.residues.atoms.select_atoms('protein and (name C and around 1.4 name O)')
    selCm1=selm1.residues.atoms.select_atoms('protein and (name C and around 1.7 name O)')
#    selOm1=selm1.residues.atoms.select_atoms('protein and (name O and around 1.4 name C)')
    selOm1=selm1.residues.atoms.select_atoms('protein and (name O and around 1.7 name C)')
    if full:
#        selCAm1=selm1.residues.atoms.select_atoms('protein and (name CA and around 1.6 name C)')
        selCAm1=selm1.residues.atoms.select_atoms('protein and (name CA and around 1.7 name C)')
    
    if full:
        return selCA,selH,selN,selCm1,selOm1,selCAm1
    else:
        return selN,selCm1,selOm1
    

def get_chain(atom,sel0,exclude=None):
    if exclude is None:exclude=[]
    '''searching a path from a methyl group of a residue down to the C-alpha of the residue
    returns a list of atoms (MDA.Atom) beginning with the Hydrogens of the methyl group and continuing
    with the carbons of the side chain
    returns empty list if atom is not a methyl carbon'''
    final=False
    def get_bonded():
        '''it happens, that pdb files do not contain bond information, in that case, we switch to selection
        by string parsing'''
        return np.sum(find_bonded([atom],sel0,n=4,d=1.7))
    
    a_name=atom.name.lower()
    a_type=atom.name[0].lower()
    if 'c'==a_name and len(exclude):
      return [atom]
    elif a_name == "n":
      return []
    connected_atoms = []
    bonded = get_bonded()
    if len(exclude)==0:
      if np.sum(np.fromiter(["h"==a.type.lower() for a in bonded],dtype=bool)) == 3:
        final=True  
        for a in bonded:
          if "h"==a.name[0].lower():
            connected_atoms.append(a)
        if not "c"==a_type:
          return []
      else:
        return []
    connected_atoms.append(atom)
    exclude.append(atom)
    for a in bonded:
      if not a in exclude:
        nxt = get_chain(a,sel0,exclude)
        for b in nxt:
           connected_atoms.append(b)
    if len(connected_atoms)>1:
      if final:
          return np.sum(connected_atoms)
      else:
          return connected_atoms
    else:
      return []

def search_methyl_groups(residue):
    methyl_groups = []
    for atom in residue.atoms:
        chain = get_chain(atom,residue,[])
        if len(chain):
            methyl_groups.append(chain)
    return methyl_groups