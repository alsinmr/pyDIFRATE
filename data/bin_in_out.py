#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 16:37:08 2019

@author: albertsmith
"""

import pickle


def save_bin(filename,obj):
    """
    |save_bin saves a python object. 
    |
    |save_bin(filename,obj)
    |
    |Fails if that object contains an MDanalysis object
    """
    
    with open(filename,'wb') as f:
        pickle.dump(obj,f)
        
def load_bin(filename):
    """
    |Loads a python object
    |
    |obj = load_bin(filename)
    |
    |If object saved with save_DIFRATE, reload with load_DIFRATE
    """
    with open(filename,'rb') as f:
        obj=pickle.load(f)
        
    return obj


def save_DIFRATE(filename,obj):
    """
    |save_DIFRATE saves a DIFRATE object. 
    |
    |save_DIFRATE(filename,obj)
    |
    |Deletes the MDanalysis object before saving- this object otherwise creates
    |a pickling error 
    """
    
    """Note- I don't understand why this function is necessary. The MDAnalysis
    universe exists in the atom selections, and can recovered from these. 
    Nonetheless, pickling fails if we don't first remove the saved universe.
    """
    
    if hasattr(obj,'sens') and hasattr(obj,'detect'):
        if obj.sens is not None and obj.sens.molecule is not None:
            obj.sens.molecule.mda_object=None
        if obj.detect is not None and obj.detect.molecule is not None:
            obj.detect.molecule.mda_object=None
    elif hasattr(obj,'molecule'):
        obj.molecule.mda_object=None
    elif hasattr(obj,'mda_object'):
        obj.mda_object=None
        
    
    save_bin(filename,obj)
    
    
def load_DIFRATE(filename):
    """
    |load_DIFRATE loads a DIFRATE object from a file
    |
    |obj = load_DIFRATE(filename)
    |
    |Replaces the mda_object in the various DIFRATE objects
    """
    
    obj=load_bin(filename)
    
    if hasattr(obj,'sens') and hasattr(obj,'detect'):
        if obj.sens is not None and obj.sens.molecule is not None and obj.sens.molecule.sel1 is not None:
            obj.sens.molecule.mda_object=obj.sens.molecule.sel1.universe
        if obj.detect is not None and obj.detect.molecule is not None and obj.detect.molecule.sel1 is not None:
            obj.detect.molecule.mda_object=obj.detect.molecule.sel1.universe
    elif hasattr(obj,'molecule') and obj.molecule.sel1 is not None:
        obj.molecule.mda_object=obj.molecule.sel1.universe
    elif hasattr(obj,'mda_object') and obj.sel1 is not None:
        obj.mda_object=obj.sel1.universe
        
    return obj