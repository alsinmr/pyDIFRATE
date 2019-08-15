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
        
        self.__MDA_info=None
        
        if np.size(args)>0:
            self.load_struct(*args)

    def load_struct(self,*args):   
        self.mda_object=mda.Universe(*args)
        
    def select_atoms(self,sel1=None,sel2=None,sel1in=None,sel2in=None,index1=None,index2=None,Nuc=None,resi=None,select=None,**kwargs):
        
        if select is not None:
            sel=self.mda_object.select_atoms(kwargs.get('select'))
        else:
            sel=self.mda_object.select_atoms('name *')
            
        if resi is not None:
            string=''
            for res in resi:
                string=string+'resid {0:.0f} or '.format(res)
            string=string[0:-4]
            sel=sel.select_atoms(string)
        
        if Nuc!=None:
            if Nuc.lower()=='15n' or Nuc.lower()=='n':                    
                self.sel1=sel.select_atoms('(name H or name HN) and around 1.1 name N')
                self.sel2=sel.select_atoms('name N and around 1.1 (name H or name HN)')
            elif Nuc.lower()=='CO':
                self.sel1=sel.select_atoms('name C and around 1.4 name O')
                self.sel2=sel.select_atoms('name O and around 1.4 name C')
            elif Nuc.lower()=='CA':
                self.sel1=sel.select_atoms('name CA and around 1.4 (name HA or name HA2)')
                self.sel2=sel.select_atoms('(name HA or name HA2) and around 1.4 name CA')
                print('Warning: selecting HA2 for glycines. Use manual selection to get HA1 or both bonds')
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
            """
            I think that this averaging over the trajectory should probably include
            a re-alignment of each molecule with itself at each time point. Consider
            the problems otherwise- for a lipid rotating in a membrane, we'd end
            up with all average bonds pointing approximately along the lipid normal.
            """
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
    
    def del_MDA_object(self):
        """
        In some cases, it is necessary to delete the MD analysis objects 
        (for example, when saving, we can't pickle the MD object). This function
        deletes the object after first saving information required to reload
        it and the atom selections
        """
        if self.mda_object is None:
            "Do nothing if no universe is stored"
            return
        else:
            uni=self.mda_object
            info=dict()
            self.__MDA_info=info
        "Save the filenames used for the universe"
        info.update({'filename':uni.filename})
        if hasattr(uni.trajectory,'filenames'):
            info.update({'filenames':uni.trajectory.filenames})
            
        "Save the id numbers of the selections"
        if self.sel1 is not None:
            info.update({'sel1':self.sel1.ids})
        if self.sel2 is not None:
            info.update({'sel2':self.sel2.ids})
        
        "Set the MD analysis objects to None"
        self.mda_object=None
        self.sel1=None
        self.sel2=None
        
    def reload_MDA(self):
        if self.__MDA_info is None:
            "Do nothing if MD analysis object hasn't been deleted"
            return
        info=self.__MDA_info
        if 'filenames' in info:
            uni=mda.Universe(info['filename'],info['filenames'].tolist())
        else:
            uni=mda.Universe(info['filename'])
        self.mda_object=uni
        
        sel0=uni.select_atoms('name *')
        if 'sel1' in info:
            self.sel1=sel0[info['sel1']]
        if 'sel2' in info:
            self.sel2=sel0[info['sel2']]
        
        self.__MDA_info=None
        