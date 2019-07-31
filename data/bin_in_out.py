#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 16:37:08 2019

@author: albertsmith
"""

import pickle

def save_bin(filename,obj):
    """
    |Save data produced in pyDIFRATE. If fails for given object, try the specific
    |input and output methods for that object
    """
    
    a=pickle.dumps(obj)
    with open(filename,'wb') as f:
        pickle.dump(obj,f)
        
def load_bin(filename):
    """
    |Loads
    """
    with open(filename,'rb') as f:
        obj=pickle.load(f)
        
    return obj