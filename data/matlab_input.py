#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Copyright 2021 Albert Smith-Penzel

This file is part of pyDIFRATE

pyDIFRATE is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

pyDIFRATE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with pyDIFRATE.  If not, see <https://www.gnu.org/licenses/>.


Questions, contact me at:
albert.smith-penzel@medizin.uni-leipzig.de


Created on Mon Sep  2 13:03:49 2019

@author: albertsmith
"""

from scipy.io import loadmat as LM
import numpy as np
from pyDIFRATE.r_class.sens import rates as Rates
#import os
#os.chdir('../r_class')
#from sens import rates as Rates
#os.chdir('../data')


def mat2py(filename):
    out=dict()  #Output will be a dictionary containing each saved variable
    x=LM(filename,squeeze_me=True)  #The initial data load
    var_names=get_names(x)    #Get name of each variable
    
    for name in var_names:
        out[name]=load_from_dict(x,name)
    
    return out
    
def get_names(x):
    """
    Returns the names of the variables stored in the .mat file from the initial\
    result of LM(filename)
    """
    
    names=list()
    skipped=['__header__','__version__','__globals__']
    for key in x.keys():
        if key not in skipped:
            names.append(key)
    return names
        
def load_from_dict(x,name):
    """
    First step in converting from scipy.io.loadmat: load from a dictionary
    """
    var=x[name]
    out=load_var(var)
    return out
    
def load_var(var):
    """
    Try to determine variable type and load. Works recursively to descend into
    structures from matlab
    """
    
    if not isinstance(var,np.ndarray):  #If not a numpy array, must be string, number, etc
        out=var
#    elif var.dtype.names is not None:
    elif var.dtype.kind=='V':       #Not sure if correct yet. Seems void is used for named structures
        out=dict()
        for n in var.dtype.names:
            try:
                out[n]=load_var(var[n].item())  #This works unless data coming from a cell in matlab
            except:
                out[n]=load_cell(var[n])    #Has a size, but no item() function. Should be a cell
    elif var.shape!=(): #Has a shape
        if var.dtype.kind not in ['O','V']:  #Just numbers- pass directly
            out=var
        else:   #Otherwise probably a list
            out=list()
            for v in var:
                out.append(load_var(v))
    else:
        print('Was not able to convert variable')
        return None
    
    return out


def load_cell(var):
    """
    Load cells into a list or a list or lists
    """ 
    out=list()       
    if var.ndim>1:  #Just use load_cell recursively
        for v in var:
            out.append(load_cell(v))
    else:   #Otherwise just load each element in the list
        for v in var:
            out.append(load_var(var))
            
    return out

def mat2rates(rates_in):
    """
    Load a rates object from a dictionary generated by mat2py
    First:
        matlab_vars=mat2py(filename)
        rates_in=matlab_vars['rates_loc']
        rates=mat2rates(rates_in)
    """
    ri=rates_in
    rates=Rates(tc=ri['tc'])
    info=ri['info'].copy()

    
    Rin=list()
    for Type,in0 in info.items():
        names_py=['stdev','dXY','offset']
        names_mat=['std','dHX','off']
        for m,n in zip(names_py,names_mat):
            if n in in0:
                in0[m]=in0.pop(n)   #Some of the variable names differ in matlab/python. Update here
                
        for k in range(in0['v0'].size):
            in1=dict()
            for m,val in in0.items():
                if val[k]==val[k]:   #This is horrifying. I can't figure out what nan is being used here, but it's not equal to itself
                    in1[m]=val[k]
            if 'Nuc1' not in in1:
                in1['Nuc1']='1H'
            try:
                rates.new_exp(Type=Type,**in1)
                Rin.append([ri[Type][k]])
            except:
                print('Failed to load {0} type experiments with index {1} from matlab'.format(Type,k))
            
    Rin=np.concatenate(Rin,axis=0)
    Rc=rates.R()
    norm=np.max(np.abs(Rin),axis=1)
    test=np.abs(Rin-Rc).max(axis=1)/norm>5e-4
    for k,t in enumerate(test):
        if t:
            print('Calculated sensitivities for exp.#{0} do not agree with those from MATLAB'.format(k))
    
    return rates

def mat2r(r_in):
    """
    Load a detector object from matlab. If a cell is loaded, this will be loaded
    as a single detector object 
    """
    pass

def mat2data(data_in):
    """
    Loads a data set from matlab (either NMR data or results of a detectors fit)
    """
    pass