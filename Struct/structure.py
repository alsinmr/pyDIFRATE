#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 15:05:19 2019

@author: albertsmith
"""

import MDAnalysis as mda
import numpy as np
import os

class molecule(object):
    def __init__(self,*args):
        self.mda_object=None
        self.sel1=None
        self.sel2=None
        self.sel1in=None
        self.sel2in=None
        self.label_in=list()
        self.label=list()
        self.vXY=np.array([])
        self.vCSA=np.array([])
        self.Ralign=list()

        self.pdb=None #Container for a pdb extracted from the mda_object
        self.pdb_id=None
        "We might want to delete this pdb upon object deletion"
        
        if np.size(args)>0:
            self.load_struct(*args)

    def load_struct(self,*args):   
        self.mda_object=mda.Universe(*args)
        
    def select_atoms(self,sel1=None,sel2=None,index1=None,index2=None,Nuc=None,sel1in=None,sel2in=None,**kwargs):
        if Nuc!=None:
            if Nuc.lower()=='15n' or Nuc.lower()=='n':
                if 'select' in kwargs:
                    sel=self.mda_object.select_atoms(kwargs.get('select'))
                else:
                    sel=self.mda_object.select_atoms('name *')
                
                self.sel1=sel.select_atoms('(name H or name HN) and around 1.1 name N')
                self.sel2=sel.select_atoms('name N and around 1.1 (name H or name HN)')
                
                self.label_in=self.sel1.resids                                
        else:
            if sel1!=None:
                if index1!=None:
                    self.sel1=self.mda_object.select_atoms(sel1)[index1]
                else:
                    self.sel1=self.mda_object.select_atoms(sel1)
            if sel2!=None:
                if index2!=None:
                    self.sel2=self.mda_object.select_atoms(sel2)[index2]
                else:
                    self.sel2=self.mda_object.select_atoms(sel2)
              
            if sel1in is not None:
                self.sel1in=sel1in
            if sel2in is not None:
                self.sel2in=sel2in
        
    
    def add_label(self,label=None):
        self.label_in=label
        
    def set_selection(self,**kwargs):
        
        if self.sel1in is None and self.sel2in is None:
            nr=np.size(self.sel1.resids)
        elif self.sel1in is None:
            nr=np.size(self.sel2in)
        else:
            nr=np.size(self.sel1in)

        
        vec=np.zeros([nr,3])
        
    #    for k in range(0, nt-2):
    
        if 'tstep' in kwargs:
            tstep=kwargs.get('tstep')
        else:
            tstep=int(self.mda_object.trajectory.n_frames/100)
            if tstep==0:
                tstep=1
    
        nt=self.mda_object.trajectory.n_frames
    
        for k in range(0,nt,tstep):
            try:
                self.mda_object.trajectory[k]
            except:
                if k!=0:
                    for _ in range(0,tstep):
                        self.mda_object.trajectory.next()
    
            if self.sel1in is None and self.sel2in is None:
                vec+=self.sel1.positions-self.sel2.positions
            elif self.sel1in is None:
                for m,q in enumerate(self.sel2in):
                    vec[m,:]+=self.sel1.positions[m]-self.sel2.positions[q]
            elif self.sel2in is None:
                for m,q in enumerate(self.sel1in):
                    vec[m,:]+=self.sel1.positions[q]-self.sel2.positions[m]
            else:
                for m,q in enumerate(self.sel1in):
                    vec[m,:]+=self.sel1.positions[q]-self.sel2.positions[self.sel2in[m]]
                    
        len=np.sqrt(np.sum(np.power(vec,2),axis=1))
        vec=np.divide(vec,np.reshape(np.repeat(len,3),vec.shape)) #Normalize the vector
        
        self.vXY=vec
        if np.shape(self.label_in)[0]==nr:
            self.label=self.label_in
        
    def MDA2pdb(self,tstep=None,selection=None,**kwargs):
        "Provide a molecule, print a certain frame to pdb for later use in chimera"
        

        uni=self.mda_object

        
        if tstep is None:
            tstep=int(uni.trajectory.n_frames/2)
            
        dir_path = os.path.dirname(os.path.realpath(__file__))
    
        full_path=os.path.join(dir_path,os.path.basename(self.mda_object.filename)+'_{0}'.format(tstep)+'.pdb')
        
        try:
            uni.trajectory[tstep]
        except:
            uni.trajectory.rewind()
            for k in range(0,tstep):
                uni.trajectory.next()
        
        if selection is not None:
            a=uni.select_atoms(selection)
        else:
            a=uni.select_atoms('protein')
            if a.n_atoms==0:
                a=uni.select_atoms('name *')
                
        a.write(full_path)
        
        self.pdb=full_path
        self.pdb_id=np.array(a.ids)
        
        return full_path