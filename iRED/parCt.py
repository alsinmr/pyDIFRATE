#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 12:36:42 2019

@author: albertsmith
"""
import numpy as np
#import traceback
from iRED.fast_index import get_count

#%% Parallel class for fast parallel calculation
class par_class():
    """
    We create this class to perform parallel calculation of correlation functions.
    The variables required for calculation are stored as *class* variables, which
    are dictionaries. Because they are class variables, these are available to
    any instance of par_class (that is, if multiple par_class *objects* are 
    created by different processes, they will all have access to these dicts. To
    make problems unlikely, we assign the dictionary keys using a random number.
    That number is passed out to the parent process, and then used to find
    the correct data later. In principle, although different processes are 
    looking into the same dictionaries, they should avoid using/editing the same
    data)
    """
    X=dict()
    Y=dict()
    Z=dict()
    
    ct=dict()
    
    keys=dict()
    nb=dict()
    nk=dict()
    index=dict()
        
    @classmethod    
    def Ct(cls,v):
        i,ref_num=v
        index=cls.index[ref_num]
        X=cls.X[i]
        Y=cls.Y[i]
        Z=cls.Z[i]
        n=np.size(index)

        ct=np.zeros([index[-1]+1,np.shape(X)[1]])
        
        for k in range(n):
            ct[index[k:]-index[k]]+=(3*(np.multiply(X[k:],X[k])+np.multiply(Y[k:],Y[k])\
                 +np.multiply(Z[k:],Z[k]))**2-1)/2
        
        "Store results of correlation function calculation"
#        cls.storeCt(i,ct)
#        cls.ct[i]=ct
        return ct
    
    @classmethod
    def CtFT(cls,v):
        i,ref_num=v
        index=cls.index[ref_num]
        X=cls.X[i]
        Y=cls.Y[i]
        Z=cls.Z[i]
        
        
#    
    @classmethod
    def store_vecs(cls,vec,nc):
        """Responsible for sorting out the vectors for each process.
        Uses class variables, which are effectively global, but indexes them randomly
        so that we shouldn't end up accessing the same variables in multiple processes
        """
        
        nk=nc   #Maybe we should change this to reduce memory usage. Right now just nc steps
        
        """nc is the number of cores to be used, and nk the number of chunks to
        do the calculation in. Currently equal.
        """
        
        ref_num=np.random.randint(0,1e9) 
        
        cls.keys[ref_num]=ref_num+np.arange(nk) #Keys where the data is stored
        cls.nb[ref_num]=vec['X'].shape[1]   #Number of correlation functions (n bonds)
        cls.nk[ref_num]=nk  #Number of chunks
        cls.index[ref_num]=vec['index'] #Index of frames taken
        nb=cls.nb[ref_num]  
        for k,i in enumerate(cls.keys[ref_num]):     #Separate and store parts of the vector
            cls.X[i]=vec['X'][:,range(k,nb,nk)]
            cls.Y[i]=vec['Y'][:,range(k,nb,nk)]
            cls.Z[i]=vec['Z'][:,range(k,nb,nk)]
            
        v0=list()    
        for i in cls.keys[ref_num]:
            v0.append((i,ref_num))
        
        return ref_num,v0
    
    @classmethod
    def returnCt(cls,ref_num,ct):
        nk=cls.nk[ref_num]
        index=cls.index[ref_num]
        N0=get_count(index)
        nz=N0!=0
        N0=N0[nz]
        nb=cls.nb[ref_num]
        
#        ct=list()
#        for i in cls.keys[ref_num]:
#            ct.append(cls.ct[i])
        
        ct0=np.zeros([np.size(N0),nb])
        for k,c in enumerate(ct):
            N=np.repeat([N0],np.shape(c)[1],axis=0).T
            ct0[:,range(k,nb,nk)]=np.divide(c[nz],N)
    
        return ct0
    
    @classmethod
    def clear_data(cls,ref_num):
        locs=['X','Y','Z','ct']
        if ref_num in cls.keys:
            for ref0 in cls.keys[ref_num]:
                for loc in locs:
                    if ref0 in getattr(cls,loc):
                        del getattr(cls,loc)[ref0]
        else:
            print('Data already deleted')
            
        locs=['keys','nb','nk','index']
        for loc in locs:
            if ref_num in getattr(cls,loc):
                del getattr(cls,loc)[ref_num]

    @classmethod
    def _clear_all(cls):
        locs=['X','Y','Z','ct']
        for loc in locs:
            while len(getattr(cls,loc).keys())!=0:
                try:
                    k=list(getattr(cls,loc).keys())
                    cls.clear_data(k[0])
                except:
                    pass
                
