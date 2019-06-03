#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 22:07:08 2019

@author: albertsmith
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import DynamicModels as dm
import os
os.chdir('../Struct')
import structure
os.chdir('../r_class')
from scipy.interpolate import interp1d as interp


class model(object):
    def __init__(self):
        
        if self._class=='Ct':
            self.__Reff=np.zeros([0,np.size(self.t()),np.size(self.tc())])
            self.__R0=np.zeros([0,np.size(self.t())])
        else:
            self.__Reff=np.zeros([0,np.size(self.tc())])
            self.__R0=np.zeros([0,1])
        self.__mdlinfo=pd.DataFrame(index=self.retExper()+self.retSpinSys())
        self.__tMdl=list()
        self.__AMdl=list()
        self.MdlPar=list()
        self.tMdl=list()
        self.AMdl=list()
        self.molecule=structure.molecule()
        
            
    def new_mdl(self,tMdl=None,AMdl=None,Model=None,**kwargs):
        
        if tMdl!=None and AMdl!=None:
            if not isinstance(tMdl,np.ndarray):
                tMdl=np.array(tMdl)
            if tMdl.shape==():
                tMdl=np.array([tMdl])
            if not isinstance(AMdl,np.ndarray):
                AMdl=np.array(AMdl)
            if AMdl.shape==():
                AMdl=np.array([AMdl])
                
            
            MdlPar=dict(Model='Direct',BondSpfc='no')
            
            self.MdlPar.append(MdlPar)
            self.tMdl.append(tMdl)
            self.AMdl.append(AMdl)
        elif Model=='Combined' and 'mdl_nums' in kwargs:
            mdl_nums=kwargs.get('mdl_nums')
            if not isinstance(mdl_nums,np.ndarray):
                mdl_nums=np.array(mdl_nums)
            if mdl_nums.shape==():
                mdl_nums=np.array([mdl_nums])
            
            BndSpfc='no'
            Models=list()
            for k in mdl_nums:
                Models.append(self.MdlPar[k])
                if self.MdlPar[k]['BondSpfc']=='yes':
                    BndSpfc='yes'
            
            MdlPar=dict(Model='Combined',BondSpfc=BndSpfc,SubModels=Models)
            
            
            tMdl=np.array([]);     #Empty model
            AMdl=np.array([]);   #
            for k in mdl_nums:
                tMdl,AMdl,_=dm.ModelSel('Combined',tMdl1=tMdl,AMdl1=AMdl,tMdl2=self.tMdl[k],AMdl2=self.AMdl[k])
            
            self.MdlPar.append(MdlPar)
            self.tMdl.append(tMdl)
            self.AMdl.append(AMdl)
            
        else:
            if dm.ModelBondSpfc(Model) and self.molecule.vXY.size==0:
                print('Before defining an model with anisotropic motion, import a structure and select the desired bonds')
            else:
                tMdl,AMdl,BndSp=dm.ModelSel(Model,'dXY',self.molecule,**kwargs)
                if BndSp=='yes' and self._class!='Ct':
                    _,A,_=dm.ModelSel(Model,'dCSA',self.molecule,**kwargs)
                    AMdl=[AMdl,A]
                    AMdl=np.swapaxes(AMdl,0,1)
                    
                MdlPar=dict(Model=Model,BondSpfc=BndSp,**kwargs)
                
                self.MdlPar.append(MdlPar)
                self.tMdl.append(tMdl)
                self.AMdl.append(AMdl)
                
    def del_mdl(self,mdl_num):
        del self.AMdl[mdl_num]
        del self.tMdl[mdl_num]
        del self.MdlPar[mdl_num]
        
    
    def _rho_eff(self,exp_num=None,mdl_num=0,bond=None,**kwargs):
        """This function is mostly responsible for searching for a pre-existing
        calculation of the model and experiment
        """
        
        if np.size(exp_num)==1 and exp_num==None:
            exp_num=self.info.columns.values
        
        "Make sure we're working with numpy array for exp_num"
        if not isinstance(exp_num,np.ndarray):
            exp_num=np.array(exp_num)
        if exp_num.shape==():
            exp_num=np.array([exp_num])
        
        
        "Get tMdl and AMdl from the model array"
        if self.MdlPar[mdl_num].get('BondSpfc')=='yes':
            A=self.AMdl[mdl_num][bond]
        else:
            A=self.AMdl[mdl_num]
        """A bond-specific model still has all the same correlation times for
        every motion. Only the influence of each correlation time changes for 
        different bonds (that is, the amplitude changes, but not the 
        correlation time). Note, that for some orientations, an amplitude can
        become 0."""
        tMdl=self.tMdl[mdl_num]
        

        if self._class=='Ct':
            
            count=0
            test=False
            n=self.__Reff.shape[0]
            while count<n and not test:
                if np.ndarray.tolist(tMdl)==np.ndarray.tolist(self.__tMdl[count]) \
                and np.ndarray.tolist(A)==np.ndarray.tolist(self.__AMdl[count]):
                    test=True
                else:
                    count=count+1
                    
            if test:
                R=self.__Reff[count,:,:]
                R0=self.__R0[count,:]
            else:
                print('It\'s slow to apply models directly to Ct objects- instead, calculate detectors first')
                "Otherwise, calculate the new sensitivity, and store it"
                R,R0=self.__apply_mdl(exp_num,tMdl,A)
                self.__Reff=np.vstack([self.__Reff,[R]])
                self.__R0=np.vstack([self.__R0,R0])
                self.__tMdl.append(tMdl)
                self.__AMdl.append(A)
                        
            
        else:
            ntc=np.size(self.tc())
            ne=exp_num.size
            R=np.zeros([ne,ntc])
            R0=np.zeros([ne,1])
            for k in range(0,ne):
                "Look to see if we've already calculated this sensitivity, return it if so"
                count=0
                test=False
                n=self.__Reff.shape[0]
    
                while count<n and not test:
                    if np.ndarray.tolist(tMdl)==np.ndarray.tolist(self.__tMdl[count]) \
                    and np.ndarray.tolist(A)==np.ndarray.tolist(self.__AMdl[count]):
                        if self.__mdlinfo.iloc[:,count].eq(self.info.loc[:,exp_num[k]]).all():
                            test=True
                        else:
                            count=count+1
                    else:
                        count=count+1
                        
                if test:
                    R[k,:]=self.__Reff[count,:]
                    R0[k]=self.__R0[count]
                else:
                    "Otherwise, calculate the new sensitivity, and store it"
                    a,b=self.__apply_mdl(exp_num[k],tMdl,A)
                    R[k,:]=a
                    R0[k]=b
                    self.__Reff=np.vstack([self.__Reff,R[k,:]])
                    self.__R0=np.vstack([self.__R0,R0[k]])
                    self.__mdlinfo=pd.concat([self.__mdlinfo,self.info.loc[:,exp_num[k]]],axis=1,ignore_index=True,sort=True)
                    self.__tMdl.append(tMdl)
                    self.__AMdl.append(A)
        return R,R0
    
    def _rho_effCSA(self,exp_num=None,mdl_num=0,bond=None):
        """Same as above, but only for the CSA interaction
        """

        
        if np.size(exp_num)==1 and exp_num==None:
            exp_num=self.info.columns.values
            
        "Make sure we're working with numpy array for exp_num"
        if not isinstance(exp_num,np.ndarray):
            exp_num=np.array(exp_num)
        if exp_num.shape==():
            exp_num=np.array([exp_num])
            
        ntc=np.size(self.tc())
        ne=exp_num.size
        
            
        if self._origin=='Ct':
            """Check if this is the correlation function class, where we do
            not need to account for CSA separately"""
            R=np.zeros([ne,ntc])
            R0=np.zeros([ne,1])
        else:
            "Get tMdl and AMdl from the model array"
            if self.MdlPar[mdl_num].get('BondSpfc')=='yes':
                A=self.AMdl[mdl_num][bond]
            else:
                A=self.AMdl[mdl_num]
            """A bond-specific model still has all the same correlation times for
            every motion. Only the influence of each correlation time changes for 
            different bonds (that is, the amplitude changes, but not the 
            correlation time). Note, that for some orientations, an amplitude can
            become 0."""
            tMdl=self.tMdl[mdl_num]
            
    
    

            R=np.zeros([ne,ntc])
            R0=np.zeros([ne,1])
            for k in range(0,ne):
                "Look to see if we've already calculated this sensitivity, return it if so"
                count=0
                test=False
                n=self.__Reff.shape[0]
    
                exper=self.info.loc[:,exp_num[k]].copy()
                exper.at['dXY']=0
                exper.at['QC']=0
    
                while count<n and not test:
                    if np.ndarray.tolist(tMdl)==np.ndarray.tolist(self.__tMdl[count]) \
                    and np.ndarray.tolist(A)==np.ndarray.tolist(self.__AMdl[count]):
                        if self.__mdlinfo.iloc[:,count].eq(exper).all():
                            test=True
                        else:
                            count=count+1
                    else:
                        count=count+1
                        
                if test:
                    R[k,:]=self.__Reff[count,:]
                    R0[k]=self.__R0[count]
                else:
                    "Otherwise, calculate the new sensitivity, and store it"
                    """We'll do this by simply adding the new experiment to the 
                    list, calculating the effective sensitivity, and then deleting
                    the new experiment again. 
                    """
                    self.new_exp(info=exper)  #Add the new experiment
                    n=self.info.columns.values[-1] #Index of the new experiment
    
                    
    
                    a,b=self.__apply_mdl(n,tMdl,A)
                    R[k,:]=a
                    R0[k]=b
                    self.__Reff=np.vstack([self.__Reff,R[k,:]])
                    self.__R0=np.vstack([self.__R0,R0[k]])
                    self.__mdlinfo=pd.concat([self.__mdlinfo,exper],axis=1,ignore_index=True)
                    self.__tMdl.append(tMdl)
                    self.__AMdl.append(A)
                    
                    self.del_exp(n) #Delete the experiment to hide this operation from the user
                
        return R,R0
    
    def __apply_mdl(self,exp_num,tMdl,A):
        "tMdl is a list of correlation times in the model, and A the amplitudes"
        "Note that if A does not add to 1, we assume that S2 is non-zero"
        
        
        if np.size(np.shape(A))==2:
            for m in range(0,2):
                if m==0:
#                    Ri=self.__temp_exper(exp_num,'vCSA')
                    Ri=self._rhoCSA(exp_num)
                    nA=np.size(A[1])
                    S2=1-sum(A[1])
                    R0=0
                    R=S2*Ri #The part of the sensitivity unaffected by the applied model (0's for isotropic motion)
                else:
                    """"self._rho() contains all interactions, so we just 
                    subtract away the CSA from the first step"""
                    Ri=self._rho(exp_num)-Ri
                    nA=np.size(A[0])
                    S2=1-sum(A[0])

                    R+=S2*Ri #The part of the sensitivity unaffected by the applied model (0's for isotropic motion)    
                    
                for k in range(0,nA):
                    R0+=A[m][k]*np.interp(np.log10(tMdl[k]),self.z(),np.ndarray.flatten(Ri))
                    
                
                for k in range(0,nA):
                    R+=A[m][k]*np.interp(self.zeff(tMdl[k]),self.z(),np.ndarray.flatten(Ri))
                
                    
                R+=-R0
            
        else:
            Ri=self._rho(exp_num)
            Rinterp=interp(self.z(),Ri,fill_value='extrapolate')
            nA=np.size(A)
            S2=1-sum(A)
            R0=np.zeros(exp_num.shape)
            R=S2*Ri #The part of the sensitivity unaffected by the applied model (0's for isotropic motion)
            for k in range(0,nA):
                R0+=np.squeeze(A[k]*Rinterp(np.log10(tMdl[k])))
            for k in range(0,nA):
                R+=A[k]*Rinterp(self.zeff(tMdl[k]))
            if np.size(np.shape(R0))==0:
                R+=-R0
            else:
                R+=-np.repeat(np.transpose([R0]),self.tc().size,axis=1)
        
        return R,R0
    
    def __temp_exper(self,exp_num,inter):
        """When we calculate dipole/CSA relaxation under a bond-specific model 
        (ex. Anisotropic diffusion), we actually need to apply a different 
        model to the motion of the CSA and dipole. To do this, we create a new 
        experiment without the dipole, or without the CSA, calculate its 
        sensitivity, and then delete the experiment from the users's scope 
        after we're done with it. This gets passed back to _rho_eff, where the 
        new model is applied to experiments with CSA and dipole separately"""
        
        exper=self.info.loc[:,exp_num].copy()
        if inter=='dXY':
            exper.at['CSA']=0
        else:
            exper.at['dXY']=0
            exper.at['QC']=0    
            """We should make sure the quadrupole doesn't count twice. I guess
            this shouldn't matter, because we usually neglect dipole and CSA
            relaxation when a quadrupole is present, but if the user puts them
            in for some reason, it would result in a double-counting of the
            quadrupole relaxation"""
        self.new_exp(info=exper)  #Add the new experiment
        n=self.info.columns.values[-1] #Index of the new experiment
        R=self._rho(n)
        
        self.del_exp(n) #Delete the experiment to hide this operation from the user
        
        return R
    
    def clear_stored(self):
        """mdl_sens objects store results of model calculations. However, in 
        some situations, this may take up too much memory. We can clear them 
        with this function"""
        self.__Reff=np.zeros([0,np.size(self.tc())])
        self.__R0=np.zeros([0,1])
        self.__mdlinfo=pd.DataFrame(index=self.retExper()+self.retSpinSys())
        self.__tMdl=list()
        self.__AMdl=list()
        
    def zeff(self,t,tau=None):
        if tau==None:
            return self.z()+np.log10(t)-np.log10(10**self.z()+t)
        else:
            return np.log10(t)+np.log10(tau)-np.log10(t+tau)
        
    def plot_eff(self,exp_num=None,mdl_num=0,bond=None,ax=None,**kwargs):
        
            
        a,b=self._rho_eff(exp_num,mdl_num,bond)
        a=a.T
        if 'norm' in kwargs and kwargs.get('norm')[0].lower()=='y':
            norm=np.max(np.abs(a),axis=0)
            a=a/np.tile(norm,[np.size(self.tc()),1])      
        
        if ax==None:
            fig=plt.figure()
            ax=fig.add_subplot(111)
            hdl=ax.plot(self.z(),a)
#            hdl=plt.plot(self.z(),a)
#            ax=hdl[0].axes
        else:
            hdl=ax.plot(self.z(),a)
            
        ax.set_xlabel(r'$\log_{10}(\tau$ / s)')
        ax.set_ylabel(r'$R$ / s$^{-1}$')
        ax.set_xlim(self.z()[[0,-1]])
        ax.set_title('Sensitivity for Model #{0}'.format(mdl_num))
        
        return hdl
    
    def _set_plot_attr(self,hdl,**kwargs):
        props=hdl[0].properties().keys()
        for k in kwargs:
            if k in props:
                for m in hdl:
                    getattr(m,'set_{}'.format(k))(kwargs.get(k))