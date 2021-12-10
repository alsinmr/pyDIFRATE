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



Created on Wed Mar 31 09:37:47 2021

@author: albertsmith
"""


"""
Object for generating sensitivities of arbitrary experiment type. 
"""

import numpy as np
import Sensitivities.mdl_sens as mdl
import importlib
import pandas as pd
from types import MethodType

class Info(pd.DataFrame):
    """
    This is essentially just a pandas DataFrame, with a small modification: it
    keeps track of if a change has been made to the values. To determine if a
    change has been made, we run info.Updated, which returns True if an update
    has been made. In this case, we need to refresh the sensitivities stored in
    self.__R and self.__RCSA. We do this lazily– only check when they've been 
    called.
    """
    def __init__(self,data=None,index=None,columns=None,dtype=None,copy=False):
        "Initialize, and add the updated variable"
        super().__init__(data,index,columns,dtype,copy)
        self.__oldvalue=None #If initialized, assume that calcs need updated
    def calcs_updated(self):
        "Indicate that __R and __RCSA objects have been recalculated"
        self.__oldvalue=self.copy() #Copy of self. Values not updated
    @property
    def T(self):
        "Make sure this remains an instance of Info if transpose is used"
        out=super().T
        return Info(out)
    def copy(self,deep=True):
        return Info(super().copy(deep))
    def drop(self,*args,**kwargs):
        "Make sure the result of drop yields an Info object"        
        return Info(super().drop(*args,**kwargs))
    @property
    def Updated(self):
        "Returns True if the dataframe has been updated"
        return not(self.equals(self.__oldvalue))

    

class Sens(mdl.Model):
    """
    Sens is a general class for generating and storing sensitivity object. 
    """
    def __init__(self,method,tc=None,z=None,Type=None,**kwargs):
        """
        Initiate an instance of the Sens object, where Sens is the general 
        sensitivity object. Can be used as a parent class for creating other
        sensitivity objects. The argument data_type determines where Sens
        searches for functions for calculating sensitivities. For example,
        setting to "NMR" will then search in the NMR_funs.py library.
        
        One may initialize Sens with arguments for generating sensitivities of
        specific experiments. Give as kwargs (all kwargs passed to new_exp 
        function)
        
        In principle, Sens may be used directly as the sensitivity object. 
        Classes inheriting from Sens may have additional data-specific functions.
        
        Sens is a child of mdl_sens.model class.
        """
        
        super().__init__(tc,z)  #Initialize the model class
        
        self.method=method
        
        "This is a module, containing all the functions for calculating sensitivites"
        self._funs=importlib.import_module('Sensitivities.functions.'+method+'_funs')
        
        "Initialize info object"
        args=np.insert(self.fun_args(display=False),0,'Type')
        if '_sorting' in dir(self._funs):
            "If sorting present, it allows us to enforce a preferred variable order"
            args0=self._funs._sorting()
            for a in args:
                if a not in args0:
                    args0.append(a)
            args=args0
            
        
        self.info=Info(columns=args,data=np.zeros([0,len(args)]),index=[]).T #Initialize info object 
        
        "Initialize storage for rate constant calculations"
        self.__R=list()
        self.__RCSA=list()  #This is for separating sensitivities from a second contribution
                            #In NMR, this second contribution is a non-colinear CSA
        
        "Add plotting/calling functions"
        if '_setup_functions' in dir(self._funs):
            plots=self._funs._setup_functions()
            for key,fun in plots.items():
                setattr(self,key,MethodType(fun,self))
       
        self.new_exp(Type,**kwargs)

#%% Get information about available sensitivities
    def fun_names(self,display=True):
        """
        Returns a list of function names. Set display=False to return the list
        as output
        """
        fun_name=list()
        for f in dir(self._funs):
            fun=getattr(self._funs,f)
            if hasattr(fun,'__code__'):
                code=fun.__code__
                if code.co_argcount!=0 and code.co_varnames[0]=='tc' and f[0]!='_':
                    fun_name.append(f)
            
        if not(display):return fun_name
        print(fun_name)
        
    def fun_args(self,fun_name=None,display=True):
        """
        Prints or returns a list of arguments for a given function, or if no 
        function is specified, returns all possible arguments
        
        """
        
        if fun_name is None:
            fun_name=self.fun_names(display=False)

            if display:
                arg_list=list()
                for f in fun_name:
                    fun=getattr(self._funs,f)
                    arg_list.append(fun.__code__.co_varnames[1:fun.__code__.co_argcount])
                for a,f in zip(arg_list,fun_name):
                    print(f+' has arguments: ',a)
            else:
                arg_list=list()
                for f in fun_name:
                    fun=getattr(self._funs,f)
                    arg_list.extend(fun.__code__.co_varnames[1:fun.__code__.co_argcount])
                return np.unique(arg_list)
        else:
            assert fun_name in dir(self._funs),'{0} does not exist for {1} sensitivities'.format(fun_name,self.method)
            fun=getattr(self._funs,fun_name)
            arg_list=fun.__code__.co_varnames[1:fun.__code__.co_argcount]
            if not(display):return arg_list
            print(fun_name+' has arguments: ',arg_list)
                
    def fun_help(self,fun_name):
        """
        Returns help for input function, along with the list of arguments that
        are accepted
        """
        assert fun_name in dir(self._funs),'{0} does not exist for {1} sensitivities'.format(fun_name,self.method)
        help(getattr(self._funs,fun_name))
            
#%% Load a new experiment            
    def new_exp(self,Type=None,info=None,**kwargs):
        """
        Adds an experiment to the sensitivity object. Provide list of parameters
        or Dataframe containing parameters for new experiment. Note the dataframe
        must contain ALL required parameters for this method (dataframe in info)
        """
        
        if info:    #Concatenate a new info frame
            self.info=pd.concat([self.info,info],axis=1,ignore_index=True)
        
        if Type:kwargs['Type']=Type     #Put the required experiment Type into kwargs
        
        if kwargs:                      
            if '_repeats' in dir(self._funs):   #Repeat kwargs the appropriate number of times
                kwargs=self._funs._repeats(**kwargs)
            
            out=self._funs._defaults(self.info,**kwargs)
            self.info=out
    #%% Delete experiment
    def del_exp(self,exp_num):
        """
        Deletes an experiment(s) by giving the experiment number (or a list of
        numbers)
        """
        
        if np.size(exp_num)>1:
            exp_num=np.array(np.atleast_1d(exp_num))
            exp_num[::-1].sort()    #Crazy, but this sorts exp_num in descending order
            "delete largest index first, because otherwise the indices will be wrong for later deletions"
            for m in exp_num:
                self.del_exp(m)
        else:
            if np.ndim(exp_num)>0:
                exp_num=np.array(exp_num[0])
            self.info=self.info.drop(exp_num,axis=1)
#            del self.__R[exp_num]
#            del self.__RCSA[exp_num]
                
            self.info.set_axis(np.arange(np.size(self.info.axes[1])),axis=1,inplace=True)
            self._clear_stored(exp_num)
   
    def _rho(self,exp_num=None,bond=None):
        """
        Calculates the sensitivities of a given experiment or experiments, and
        returns as a numpy array
        """

        
        if self.info.Updated:
            self.__R=list()
            for k,exper in self.info.iteritems():
                pd_dict=exper.to_dict()
                Type=pd_dict.pop('Type')
                args=self.fun_args(Type,display=False)
                kwargs=dict()
                for key,value in pd_dict.items():
                    if key in args:
                        kwargs[key]=value
                    
                self.__R.append(getattr(self._funs,Type)(self.tc,**kwargs))
                
            self.info.calcs_updated()
        
        if exp_num is None:exp_num=self.info.axes[1]
        exp_num=np.atleast_1d(exp_num)
        
        return np.atleast_2d(np.array(self.__R)[exp_num])
        
        