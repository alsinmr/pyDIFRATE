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
<<<<<<< HEAD
os.chdir('../data')
=======
from detectors import detect
os.chdir('../data')
import data_class as dc
>>>>>>> data_input_updates

def load_from_file(filename):
    keys0=np.array(['info','data','model'])
    
<<<<<<< HEAD
    rate_args=None
    R=None
    Rstd=None
    mdl_args=None
=======
    data=dc.data()
    data.sens=rates()
    
    rate_args=list()
    mdl_args=list()
>>>>>>> data_input_updates
    
    with open(filename,'r') as f:
        while not eof(f):
            a=f.readline()
            if a.strip().lower()=='info':
                rate_args=read_info(f,keys0)
<<<<<<< HEAD
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
=======
                for k in rate_args:
                    data.sens.new_exp(**k)
            elif a.strip().lower()=='data':
                R,Rstd=read_data(f,keys0)
                data.R=R
                data.R_std=Rstd
            elif a.strip().lower()=='model':
                mdl_args.append(read_model(f,keys0))
    
    mdl=False
    for mdls in mdl_args:
        data.sens.new_mdl(**mdls)
        mdl=True
    
    if mdl:
        data.detect=detect(data.sens,mdl_num=0)
    else:
        data.detect=detect(data.sens)
            
    return data

def read_data(f,keys0):    
    cont=True
    R=list()
    Rstd=list()
    ne=0
    while not(eof(f)) and cont:
        pos=f.tell()
        a=f.readline()
        
        
        if np.isin(a.strip(),['R','Rstd']):
            if a.strip()=='R':
                R.append(read_lines(f,np.concatenate((keys0,['R','Rstd','label']))))
            elif a.strip()=='Rstd':
                Rstd.append(read_lines(f,np.concatenate((keys0,['R','Rstd','label']))))
            elif a.strip()=='label':
                label=read_label(f,np.concatenate((keys0,['R','Rstd','label'])))
        elif np.isin(a.strip(),keys0):
            cont=False
            f.seek(pos)
    
    if np.size(R)!=0:
        R=np.concatenate(R)
    else:
        R=None
        print('Warning: no data found in data entry')
        return None,None

    if np.size(Rstd)!=0:
        Rstd=np.concatenate(Rstd)
    else:
        Rstd=None
    

    if Rstd is None:
        print('Warning: Standard deviations are not provided')
        print('Standard deviations set equal to 1/10 of the median of the rate constants')
        ne=R.shape[0]
        Rstd=np.repeat([np.median(R,axis=0)],ne,axis=0)
    elif np.any(R.shape!=Rstd.shape):
        print('Warning: Shape of standard deviation does not match shape of rate constants')
        print('Standard deviations set equal to 1/10 of the median of the rate constants')
        ne=R.shape[0]
        
        Rstd=np.repeat([np.median(R,axis=0)],ne,axis=0)
        
    
    
    return R,Rstd

def read_lines(f,keys0):
    
    R=list()
    ne=0
    cont=True
    
    while not(eof(f)) and cont:
        pos=f.tell()
        a=f.readline()
        if np.isin(a.strip(),keys0):
            cont=False
            f.seek(pos)
        else:
            try:
                R0=np.atleast_1d(a.strip().split()).astype('float')
            except:
                print('Warning: Could not convert data in file into float')
                return None
            if R0.size>0:
                if ne==0:
                    ne=R0.size
                elif ne!=R0.size:
                    print('Inconsistent column sizes, data input aborted')
                    return None
                
                R.append(R0)
                
    return np.atleast_2d(R)

def read_label(f,keys0):
    
    label=list()
    
                
def read_model(f,keys0):
    mdl_pars=dict()
    cont=True
    
    while not(eof(f)) and cont:
        pos=f.tell()
        a=f.readline()
        if np.isin(a.strip(),keys0):
            cont=False
            f.seek(pos)
        else:
            name=a.strip()
            a=f.readline()
            val=np.atleast_1d(a.strip().split())
            try:
                val=val.astype('float')
            except:
                pass
            if val.size==1:
                val=val[0]
            mdl_pars.update({name:val})
    
    return mdl_pars
def read_info(f,keys0):
    temp=rates()
    keywords=np.concatenate((temp.retExper(),temp.retSpinSys())) #These are the possible variables to load

    rate_args=list()
    args=dict()
    used=list()
    cont=True
    while not eof(f) and cont:
        pos=f.tell()
        a=f.readline()

        if np.isin(a.strip(),keywords):
            name=a.strip()
            if name in used:    #We reset to a new set of experiments if a parameter is repeated (usually 'Type')
                rate_args.append(args)
                used=list()
#                print(args)
>>>>>>> data_input_updates
                args=dict()
            else:
                used.append(name)
                
<<<<<<< HEAD
            val=np.array(f.readline().strip().split())
            print(val)
            try:
                val=val.astype('float')
            except:
                pass
=======
            val=f.readline().strip().split()
            try:
                val=np.array(val).astype('float')
            except:
                pass
            
>>>>>>> data_input_updates
            args.update({name:val})
        
        elif np.isin(a.strip().lower(),keys0):
            cont=False
<<<<<<< HEAD
        
        rate_args.append(args)
=======
            f.seek(pos)
        
    if args:
        rate_args.append(args)
        
>>>>>>> data_input_updates
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
        