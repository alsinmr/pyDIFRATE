#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions for loading experimental relaxation rate constants from NMR

Created on Mon Jul 29 14:14:04 2019

@author: albertsmith
"""

import os
import numpy as np

os.chdir('../r_class')
from sens import rates
os.chdir('../data')

def load_from_file(filename):
    keys0=np.array(['info','data','model'])
    
    rate_args=None
    R=None
    Rstd=None
    mdl_args=None
    
    with open(filename,'r') as f:
        while not eof(f):
            a=f.readline()
            if a.strip().lower()=='info':
                rate_args=read_info(f,keys0)
            elif a.strip().lower()=='data':
                R,Rstd=read_data(f,keys0)
            elif a.strip().lower()=='model':
                mdl_args=read_model(f,keys0)
                
    print(rate_args)
    return rate_args

def read_data(f,keys0):
    pass

def read_model(f,keys0):
    pass

def read_info(f,keys0):

    rate_args=list()
    rat=rates()
    keywords=np.concatenate((rat.retExper(),rat.retSpinSys())) #These are the possible variables to load
    
    cont=True
    
    while not eof(f) and cont:
        pos=f.tell()
        a=f.readline()
        print(a)
        args=dict()
        used=list()
        if np.isin(a.strip().lower(),keywords):
            name=a.strip().lower()
            if name in used:    #We reset to a new set of experiments if a parameter is repeated (usually 'Type')
                rate_args.append(args)
                used=list()
                args=dict()
            else:
                used.append(name)
                
            val=np.array(f.readline().strip().split())
            print(val)
            try:
                val=val.astype('float')
            except:
                pass
            args.update({name:val})
        
        elif np.isin(a.strip().lower(),keys0):
            cont=False
        
        rate_args.append(args)
    return rate_args

def eof(f):
    "Determines if we are at the end of the file"
    pos=f.tell()    #Current position in the file
    f.readline()    #Read out a line
    if pos==f.tell(): #If position unchanged, we're at end of file
        return True
    else:       #Otherwise, reset pointer, return False
        f.seek(pos)
        return False
        