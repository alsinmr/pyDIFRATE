#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 21:41:57 2019

@author: albertsmith
"""

import numpy as np
import pandas as pd
import DIFRATE_funs as dff
import matplotlib.pyplot as plt
import mdl_sens as mdl
import sens

class rates(mdl.model):  
    def __init__(self,tc=None,**kwargs):
        
        
        """Probably a better way to do this, but I need to identify which
        child of mdl_sens is which later. Using isinstance requires me to 
        import the children into mdl_sens, but also import mdl_sens into its
        children. This seems to create some strange dependence so that I can't
        actually load any of the classes any more"""
        
        self._class='rates'
        self._origin='rates'
        """The detectors class may have bond-specific sensitivities in _rho. We
        need to know if this is the case for the mdl_sens class to work 
        properly
        """
        self._BondSpfc='no'
        
        """Get user defined tc if provided. Options are to provide the tc 
        vector directly, to provide the log directly (define z instead of tc), 
        or to specify it as a start, end, and number of steps, which are then 
        log-spaced (3 entries for tc or z)
        """
        if tc is None:
            if 'z' in kwargs:
                z=kwargs.get('z')
                if np.size(z)==3:
                    self.__tc=np.logspace(z[0],z[1],z[2])
                else:
                    self.__tc=np.power(10,kwargs.get('z'))
                "Allow users to input z instead of tc"
            else:
                self.__tc=np.logspace(-14,-3,200)
        elif np.size(tc)==3:
            self.__tc=np.logspace(np.log10(tc[0]),np.log10(tc[1]),tc[2])
        else:
            self.__tc=tc
        """We don't allow editing of the tc vector; you must initialize a new 
        instance of rates if you want to change it"""
        
        
        """If you want to edit the code to include new experiments, and these 
        require new variables, they MUST be added to one of these lists
        """
        
        "Names of the experimental variables that are available"
        self.__exper=['Type','v0','v1','vr','offset','stdev']
        "Names of the spin system variables that are available"
        self.__spinsys=['Nuc','Nuc1','dXY','CSA','QC','eta','theta']
    
        "We need to initialize self.info"
        self.info=None  
        "Here we feed in all the information on the experiments"
        self.new_exp(**kwargs)

        "Initialize some storage for rate constant calculation"
        self.__R=np.zeros([0,np.size(self.__tc)])
        self.__info=pd.DataFrame(index=self.__exper+self.__spinsys)
        
        super().__init__()
    
    def new_exp(self,**kwargs):
        "Count how many experiments are given"
        ne=0
        for k in self.__exper:
            if k in kwargs:
                ne=np.max([ne,np.size(kwargs.get(k))])
                
        "Move all input variables to __sys and __exp"
        "Defaults that don't depend on the observed nucleus can be set here"
        self.__exp=dict()
        for k in self.__exper:
            if k in kwargs:
                self.__exp.update({k : kwargs.get(k)})
            else:
                self.__exp.update({k : None})
        
        self.__sys=dict()
        for k in self.__spinsys:
            if k in kwargs:
                self.__sys.update({k : kwargs.get(k)})
            elif k=='Nuc':
                self.__sys.update({k : '15N'})
            else:
                self.__sys.update({k : None})

        
        self.__cleanup(ne) 
        self.__set_defaults(ne)
                
        "Copy all experiments from a previous rates variable, input as info"
        if 'info' in kwargs and isinstance(kwargs.get('info'),sens.rates):
            info0=kwargs.get('info').info
        else:
            info0=kwargs.get('info')
            

        

        "Create the new pandas array"
        info=pd.concat([pd.DataFrame.from_dict(self.__exp),pd.DataFrame.from_dict(self.__sys)],axis=1)
        
        try: #I don't love using try....but not sure what to do here
            info=pd.concat([info0,info.T],axis=1,ignore_index=True)
        except:
            info=info.T
        
        if not isinstance(self.info,pd.DataFrame):
            self.info=info
        else:
            self.info=pd.concat([self.info,info],axis=1,ignore_index=True)
         
#%% Make sure inputs all are the correct type (numpy arrays)         
    "Function to make sure all inputs are arrays, and have the correct sizes"    
    def __cleanup(self,ne):
        
        "Check that all experimental variables can be treated as arrays"
        
        for k in self.__exper:
            a=self.__exp.get(k)
            if not isinstance(a,(list,np.ndarray,pd.DataFrame)):
                a=[a]*ne
            elif np.size(a)!=ne:
                "We tile the output if the number of experiments doesn't match up"
                a=a*int(np.ceil(ne/np.size(a)))
                a=a[0:ne]
            if not isinstance(a,np.ndarray):
                a=np.array(a)
            self.__exp.update({k:a})
                
        for k in self.__spinsys:
            a=self.__sys.get(k)
            if not isinstance(a,(list,np.ndarray)):
                a=[a]*ne
            elif k=='dXY' or k=='Nuc1':
                b=np.array([None]*ne)
                for m in range(0,ne):
                    b[m]=a
                a=b
            if not isinstance(a,np.ndarray):
                a=np.array(a)
            
            if a.dtype.str[0:2]=='<U':
                a=a.astype('<U6')
                
            self.__sys.update({k:a})
            
        
#%% Setting defaults             
    "Function for setting defaults"
    def __set_defaults(self,ne):
        Nuc=self.__sys.get('Nuc')
        for k in range(0,ne):
            if Nuc[k].upper()=='15N' or Nuc[k].upper()=='N15' or Nuc[k].upper()=='N':
                self.__N15_def(k)
            elif Nuc[k].upper()=='13CA' or Nuc[k].upper()=='CA':
                self.__CA_def(k)
            elif Nuc[k].upper()=='13CO' or Nuc[k].upper()=='CO':
                self.__CO_def(k)
            elif Nuc[k].upper()=='CD2H' or Nuc[k].upper()=='13CD2H' or Nuc[k].upper()=='CHD2' or Nuc[k].upper()=='13CHD2':
                self.__CD2H_def(k)
            elif Nuc[k].upper()=='2H' or Nuc[k].upper()=='D':
                self.__2H_def(k)
        
            
            
            for m in self.__spinsys:
                a=self.__sys.get(m)[k]
                if a is None and m!='Nuc1':
                    a=0
                    self.__sys.get(m)[k]=a
                    
            for m in self.__exper:
                a=self.__exp.get(m)[k]
                if np.size(a)==1 and a is None and m!='stdev':
                    a=0
                    self.__exp.get(m)[k]=a
                
    "Function called to set N15 defaults"
    def __N15_def(self,k):
        self.__sys.get('Nuc')[k]='15N'
        if np.size(self.__sys.get('Nuc1')[k])==1 and self.__sys.get('Nuc1')[k] is None:
            self.__sys.get('Nuc1')[k]='1H'
        if np.size(self.__sys.get('dXY')[k])==1 and any([self.__sys.get('dXY')[k] is None]):
            self.__sys.get('dXY')[k]=dff.dipole_coupling(.102,'15N',self.__sys.get('Nuc1')[k])
        if self.__sys.get('CSA')[k] is None:
            self.__sys.get('CSA')[k]=113
        if self.__sys.get('theta')[k] is None:
            self.__sys.get('theta')[k]=23
        
                
    "Function called to set 13CA defaults"
    def __CA_def(self,k):
        self.__sys.get('Nuc')[k]='13C'
        if self.__sys.get('Nuc1')[k] is None:
            self.__sys.get('Nuc1')[k]='1H'
        if self.__sys.get('dXY')[k] is None:
            self.__sys.get('dXY')[k]=dff.dipole_coupling(.1105,'13C',self.Nuc1[k])
        if self.__sys.get('CSA')[k] is None:
            self.__sys.get('CSA')[k]=20
        if self.__sys.get('theta')[k] is None:
            self.__sys.get('theta')[k]=0
            
    "Function called to set 13CO defaults"
    def __CO_def(self,k):
        self.__sys.get('Nuc')[k]='13C'

        if self.__sys.get('CSA')[k] is None:
            self.__sys.get('CSA')[k]=155
        if self.__sys.get('theta')[k] is None:
            self.__sys.get('theta')[k]=0
        self.__sys.get('Nuc1')[k]=None
        self.__sys.get('dXY')[k]=0.0   
         
    "Function called to set Methyl CD2H defaults"
    def __CD2H_def(self,k):
        self.__sys.get('Nuc')[k]='13C'
        if self.__sys.get('Nuc1')[k] is None:
            self.__sys.get('Nuc1')[k]=['1H','2H','2H']
        if self.__sys.get('dXY')[k] is None:
            self.__sys.get('dXY')[k]=[dff.dipole_coupling(.1115,'1H','13C'),
                      dff.dipole_coupling(.1110,'2H','13C'),
                      dff.dipole_coupling(.1110,'2H','13C')]
        if self.__sys.get('CSA')[k] is None:
            self.__sys.get('CSA')[k]=16.6667  
        if self.__sys.get('theta')[k] is None:
            self.__sys.get('theta')[k]=0
            
    "Function called to set 2H defaults (Quadrupole only)"     
    def __2H_def(self,k):
        self.__sys.get('Nuc')[k]='2H'
        if self.__sys.get('dXY')[k] is None:
            self.__sys.get('dXY')[k]=0.0
        if self.__sys.get('CSA')[k] is None:
            self.__sys.get('CSA')[k]=0
        if self.__sys.get('theta')[k] is None:
            self.__sys.get('theta')[k]=0
        if self.__sys.get('QC')[k] is None:
            self.__sys.get('QC')[k]=170e3
            
        "Function called to set O17 in water defaults (Quadrupole only)"     
    def __H217O_def(self,k):
        self.__sys.get('Nuc')[k]='17O'
        if self.__sys.get('dXY')[k] is None:
            self.__sys.get('dXY')[k]=0.0
        if self.__sys.get('CSA')[k] is None:
            self.__sys.get('CSA')[k]=0
        if self.__sys.get('theta')[k] is None:
            self.__sys.get('theta')[k]=0
        if self.__sys.get('QC')[k] is None:
            self.__sys.get('QC')[k]=8.2e6

#%% Delete experiment
    def del_exp(self,exp_num):
        if np.size(exp_num)>1:
            for k in exp_num:
                self.info=self.info.drop(k,axis=1)
        else:
            self.info=self.info.drop(exp_num,axis=1)
            
#%% Adjust a parameter
    "We can adjust all parameters of a given type, or just one with the experiment index"
    def set_par(self,type,index=None,value=None):
        if index is None:
            self.info.at[type,:]=value
        else:
            self.info.at[type,index]=value

#%% Correlation time axes   
    "Return correlation times or log of correlation times"        
    def tc(self):
        return self.__tc
    
    def z(self):
        return np.log10(self.__tc)
    
#%% Rate constant calculations
    "Calculate rate constants for given experiment"
    def R(self,exp_num):
        "Make sure we're working with numpy array"
        exp_num=np.atleast_1d(exp_num)
            
        ntc=self.__tc.size
        ne=exp_num.size
        R=np.zeros([ne,ntc])
        for k in range(0,ne):
            "Look to see if we've already calculated this sensitivity, return it if so"
            count=0
            test=False
            n=self.__R.shape[0]
            while count<n and not test:
                if self.__info.iloc[:,count].eq(self.info.loc[:,exp_num[k]]).all():
                    test=True
                else:
                    count=count+1
                    
            if test:
                R[k,:]=self.__R[count,:]
            else:
                "Otherwise, calculate the new sensitivity, and store it"
                R[k,:]=dff.rate(self.__tc,self.info.loc[:,exp_num[k]])
                self.__R=np.vstack([self.__R,R[k,:]])
                self.__info=pd.concat([self.__info,self.info.loc[:,exp_num[k]]],axis=1,ignore_index=True)
        return R
    
    def _rhoCSA(self,exp_num,bond=None):
        """Calculates relaxation due to CSA only. We need this function to 
        allow application of anisotropic models, which then have different 
        influence depending on the direction of the interaction tensor. CSA
        points in a different direction (slighlty) than the dipole coupling
        """
        "Make sure we're working with numpy array"
        exp_num=np.atleast_1d(exp_num)            
            

        ntc=self.__tc.size
        ne=exp_num.size
        R=np.zeros([ne,ntc])
        for k in range(0,ne):
            "Get the information for this experiment"
            exper=self.info.loc[:,exp_num[k]].copy()
            "Turn off other interactions"
            exper.at['Nuc1']=None
            exper.at['dXY']=0
            exper.at['QC']=0
            
            "Look to see if we've already calculated this sensitivity, return it if so"
            count=0
            test=False
            n=self.__R.shape[0]
            while count<n and not test:
                if self.__info.iloc[:,count].eq(exper).all():
                    test=True
                else:
                    count=count+1
                    
            if test:
                R[k,:]=self.__R[count,:]
            else:
                "Otherwise, calculate the new sensitivity, and store it"
                R[k,:]=dff.rate(self.__tc,exper)
                self.__R=np.vstack([self.__R,R[k,:]])
                self.__info=pd.concat([self.__info,exper],axis=1,ignore_index=True)
        return R
    
    
#%% Plot the rate constant sensitivites
    def plot_R(self,exp_num=None,ax=None,**kwargs):
        
        if exp_num is None:
            exp_num=self.info.columns.values
            
        a=self.R(exp_num).T
        if 'norm' in kwargs and kwargs.get('norm')[0].lower()=='y':
            norm=np.max(a,axis=0)
            a=a/np.tile(norm,[np.size(self.tc()),1])      
        
        if ax is None:
            fig=plt.figure()
            ax=fig.add_subplot(111)
            hdl=ax.plot(self.z(),a)
#            ax=hdl[0].axes
        else:
            hdl=ax.plot(self.z(),a)
        
        self._set_plot_attr(hdl,**kwargs)
        
            
        ax.set_xlabel(r'$\log_{10}(\tau$ / s)')
        ax.set_ylabel(r'$R$ / s$^{-1}$')
        ax.set_xlim(self.z()[[0,-1]])
        ax.set_title('Rate Constant Sensitivity (no model)')
        
        return hdl

#%% Return the names of the experiment and sys variables
    def retSpinSys(self):
        return self.__spinsys
    def retExper(self):
        return self.__exper
        
#%% Hidden output of rates (semi-hidden, can be found if the user knows about it ;-) )
    def _rho(self,exp_num,bond=None):
        """The different children of mdl_sens will have different names for 
        their sensitivities. For example, this class returns R, which are the 
        rate constant sensitivities, but the correlation function class returns
        Ct, and the detector class returns rho. Then, we have a function, 
        _rho(self), that exists and functions the same way in all children
        """
        return self.R(exp_num)
    
    def Reff(self,exp_num=None,mdl_num=0,bond=None):
        R,R0=self._rho_eff(exp_num,mdl_num,bond)
        
        return R,R0