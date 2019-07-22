#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 22:07:08 2019

@author: albertsmith
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import DynamicModels as dm
import os
os.chdir('../Struct')
import structure
os.chdir('../r_class')
from scipy.interpolate import interp1d as interp


class model(object):
    def __init__(self):
        
#        if self._class=='Ct':
#            self.__Reff=np.zeros([0,np.size(self.t()),np.size(self.tc())])
#            self.__R0=np.zeros([0,np.size(self.t())])
#        else:
#            self.__Reff=np.zeros([0,np.size(self.tc())])
#            self.__R0=np.zeros([0,1])
    
#        self.__mdlinfo=pd.DataFrame(index=self.retExper()+self.retSpinSys())
#        self.__tMdl=list()
#        self.__AMdl=list()
        
        self.__Reff=list()
        self.__R0=list()
        self.__ReffCSA=list()
        self.__R0CSA=list()
        
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
    
        self.__Reff.append(None)
        self.__R0.append(None)
        self.__ReffCSA.append(None)
        self.__R0CSA.append(None)
        
            
    def del_mdl(self,mdl_num):
        del self.AMdl[mdl_num]
        del self.tMdl[mdl_num]
        del self.MdlPar[mdl_num]
        del self.__Reff[mdl_num]
        del self.__R0[mdl_num]
        del self.__ReffCSA[mdl_num]
        del self.__R0CSA[mdl_num]
        
    def del_mdl_calcs(self):
        self.__Reff=list(np.repeat(None,np.size(self.MdlPar)))
        self.__R0=list(np.repeat(None,np.size(self.MdlPar)))
        self.__ReffCSA=list(np.repeat(None,np.size(self.MdlPar)))
        self.__R0CSA=list(np.repeat(None,np.size(self.MdlPar)))
    
    def _rho_eff(self,exp_num=None,mdl_num=0,bond=None,**kwargs):
        """This function is mostly responsible for searching for a pre-existing
        calculation of the model and experiment
        """
        
        if self.__Reff[mdl_num] is None:
            Reff,R0,ReffCSA,R0CSA=self.__apply_mdl(self.tMdl[mdl_num],self.AMdl[mdl_num])
            self.__Reff[mdl_num]=Reff
            self.__R0[mdl_num]=R0
            self.__ReffCSA[mdl_num]=ReffCSA
            self.__R0CSA[mdl_num]=R0CSA
            
        
        if exp_num is None and bond is None:
            R=self.__Reff[mdl_num]
            R0=self.__R0[mdl_num]
        elif exp_num is None:
            R=self.__Reff[mdl_num][bond,:,:]
            R0=self.__R0[mdl_num][bond,:]
        else:
            R=self.__Reff[mdl_num][bond,exp_num,:]
            R0=self.__R0[mdl_num][bond,exp_num]
            
        return R,R0
        
#        if np.size(exp_num)==1 and exp_num==None:
#            exp_num=self.info.columns.values
#        
#        "Make sure we're working with numpy array for exp_num"
#        if not isinstance(exp_num,np.ndarray):
#            exp_num=np.array(exp_num)
#        if exp_num.shape==():
#            exp_num=np.array([exp_num])
#        
#        
#        "Get tMdl and AMdl from the model array"
#        if self.MdlPar[mdl_num].get('BondSpfc')=='yes':
#            A=self.AMdl[mdl_num][bond]
#        else:
#            A=self.AMdl[mdl_num]
#        """A bond-specific model still has all the same correlation times for
#        every motion. Only the influence of each correlation time changes for 
#        different bonds (that is, the amplitude changes, but not the 
#        correlation time). Note, that for some orientations, an amplitude can
#        become 0."""
#        tMdl=self.tMdl[mdl_num]
#        
#
#        if self._class=='Ct':
#            
#            count=0
#            test=False
#            n=self.__Reff.shape[0]
#            while count<n and not test:
#                if np.ndarray.tolist(tMdl)==np.ndarray.tolist(self.__tMdl[count]) \
#                and np.ndarray.tolist(A)==np.ndarray.tolist(self.__AMdl[count]):
#                    test=True
#                else:
#                    count=count+1
#                    
#            if test:
#                R=self.__Reff[count,:,:]
#                R0=self.__R0[count,:]
#            else:
#                print('It\'s slow to apply models directly to Ct objects- instead, calculate detectors first')
#                "Otherwise, calculate the new sensitivity, and store it"
#                R,R0=self.__apply_mdl(exp_num,tMdl,A)
#                self.__Reff=np.vstack([self.__Reff,[R]])
#                self.__R0=np.vstack([self.__R0,R0])
#                self.__tMdl.append(tMdl)
#                self.__AMdl.append(A)
#                        
#            
#        else:
#            ntc=np.size(self.tc())
#            ne=exp_num.size
#            R=np.zeros([ne,ntc])
#            R0=np.zeros([ne,1])
#            for k in range(0,ne):
#                "Look to see if we've already calculated this sensitivity, return it if so"
#                count=0
#                test=False
#                n=self.__Reff.shape[0]
#    
#                while count<n and not test:
#                    if np.ndarray.tolist(tMdl)==np.ndarray.tolist(self.__tMdl[count]) \
#                    and np.ndarray.tolist(A)==np.ndarray.tolist(self.__AMdl[count]):
#                        if self.__mdlinfo.iloc[:,count].eq(self.info.loc[:,exp_num[k]]).all():
#                            test=True
#                        else:
#                            count=count+1
#                    else:
#                        count=count+1
#                        
#                if test:
#                    R[k,:]=self.__Reff[count,:]
#                    R0[k]=self.__R0[count]
#                else:
#                    "Otherwise, calculate the new sensitivity, and store it"
#                    a,b=self.__apply_mdl(exp_num[k],tMdl,A)
#                    R[k,:]=a
#                    R0[k]=b
#                    self.__Reff=np.vstack([self.__Reff,R[k,:]])
#                    self.__R0=np.vstack([self.__R0,R0[k]])
#                    self.__mdlinfo=pd.concat([self.__mdlinfo,self.info.loc[:,exp_num[k]]],axis=1,ignore_index=True,sort=True)
#                    self.__tMdl.append(tMdl)
#                    self.__AMdl.append(A)
#        return R,R0
    
    def _rho_effCSA(self,exp_num=None,mdl_num=0,bond=None):
        """Same as above, but only for the CSA interaction
        """
        
        if self.__ReffCSA[mdl_num] is None:
            Reff,R0,ReffCSA,R0CSA=self.__apply_mdl(self.tMdl[mdl_num],self.AMdl[mdl_num])
            self.__Reff[mdl_num]=Reff
            self.__R0[mdl_num]=R0
            self.__ReffCSA[mdl_num]=ReffCSA
            self.__R0CSA[mdl_num]=R0CSA
            
        
        if exp_num is None and bond is None:
            R=self.__ReffCSA[mdl_num]
            R0=self.__R0CSA[mdl_num]
        elif exp_num is None:
            R=self.__ReffCSA[mdl_num][bond,:,:]
            R0=self.__R0CSA[mdl_num][bond,:]
        else:
            R=self.__ReffCSA[mdl_num][bond,exp_num,:]
            R0=self.__R0CSA[mdl_num][bond,exp_num]
            
        return R,R0
#
#        
#        if np.size(exp_num)==1 and exp_num==None:
#            exp_num=self.info.columns.values
#            
#        "Make sure we're working with numpy array for exp_num"
#        if not isinstance(exp_num,np.ndarray):
#            exp_num=np.array(exp_num)
#        if exp_num.shape==():
#            exp_num=np.array([exp_num])
#            
#        ntc=np.size(self.tc())
#        ne=exp_num.size
#        
#            
#        if self._origin=='Ct':
#            """Check if this is the correlation function class, where we do
#            not need to account for CSA separately"""
#            R=np.zeros([ne,ntc])
#            R0=np.zeros([ne,1])
#        else:
#            "Get tMdl and AMdl from the model array"
#            if self.MdlPar[mdl_num].get('BondSpfc')=='yes':
#                A=self.AMdl[mdl_num][bond]
#            else:
#                A=self.AMdl[mdl_num]
#            """A bond-specific model still has all the same correlation times for
#            every motion. Only the influence of each correlation time changes for 
#            different bonds (that is, the amplitude changes, but not the 
#            correlation time). Note, that for some orientations, an amplitude can
#            become 0."""
#            tMdl=self.tMdl[mdl_num]
#            
#    
#    
#
#            R=np.zeros([ne,ntc])
#            R0=np.zeros([ne,1])
#            for k in range(0,ne):
#                "Look to see if we've already calculated this sensitivity, return it if so"
#                count=0
#                test=False
#                n=self.__Reff.shape[0]
#    
#                exper=self.info.loc[:,exp_num[k]].copy()
#                exper.at['dXY']=0
#                exper.at['QC']=0
#    
#                while count<n and not test:
#                    if np.ndarray.tolist(tMdl)==np.ndarray.tolist(self.__tMdl[count]) \
#                    and np.ndarray.tolist(A)==np.ndarray.tolist(self.__AMdl[count]):
#                        if self.__mdlinfo.iloc[:,count].eq(exper).all():
#                            test=True
#                        else:
#                            count=count+1
#                    else:
#                        count=count+1
#                        
#                if test:
#                    R[k,:]=self.__Reff[count,:]
#                    R0[k]=self.__R0[count]
#                else:
#                    "Otherwise, calculate the new sensitivity, and store it"
#                    """We'll do this by simply adding the new experiment to the 
#                    list, calculating the effective sensitivity, and then deleting
#                    the new experiment again. 
#                    """
#                    self.new_exp(info=exper)  #Add the new experiment
#                    n=self.info.columns.values[-1] #Index of the new experiment
#    
#                    
#    
#                    a,b=self.__apply_mdl(n,tMdl,A)
#                    R[k,:]=a
#                    R0[k]=b
#                    self.__Reff=np.vstack([self.__Reff,R[k,:]])
#                    self.__R0=np.vstack([self.__R0,R0[k]])
#                    self.__mdlinfo=pd.concat([self.__mdlinfo,exper],axis=1,ignore_index=True)
#                    self.__tMdl.append(tMdl)
#                    self.__AMdl.append(A)
#                    
#                    self.del_exp(n) #Delete the experiment to hide this operation from the user
#                
#        return R,R0
    
    def __apply_mdl(self,tMdl,A):
        "tMdl is a list of correlation times in the model, and A the amplitudes"
        "Note that if A does not add to 1, we assume that S2 is non-zero"
        
        
        "Get the experimental sensitivities"
        R=self._rho(self.info.axes[1])
        RCSA=self._rhoCSA(self.info.axes[1])
        R+=-RCSA #We operate on relaxation from dipole and CSA separately
        
        "Shapes of matrices, preallocation"
        SZA=np.shape(A)
        if np.size(SZA)>1:
            SZA=SZA[0]
            iso=False
        else:
            iso=True
            
        SZR=R.shape
        SZeff=np.concatenate([np.atleast_1d(SZA),np.atleast_1d(SZR)])
        SZ0=np.concatenate([np.atleast_1d(SZA),[SZR[0]]])
        
        "Contributions to relaxation coming from model with non-zero S2"
        if np.size(A)>1:
            S2=1-np.sum(A[:,0,:],axis=1)
            S2CSA=1-np.sum(A[:,1,:],axis=1)
        else:
            S2=1-np.sum(A)
            S2CSA=S2
        Reff=np.reshape(np.dot(np.transpose([S2]),np.reshape(R,[1,np.prod(SZR)])),SZeff)
        ReffCSA=np.reshape(np.dot(np.transpose([S2CSA]),np.reshape(RCSA,[1,np.prod(SZR)])),SZeff)
        
        R0=np.zeros(SZ0)
        R0CSA=np.zeros(SZ0)
        
        "Loop over all correlation times in model"
        for k,tc in enumerate(tMdl):
            "Matrix to transform from z to zeff (or simply to evaluate at z=log10(tc) with M0)"
            M,M0=self.z2zeff(tc)
            
            
            R00=np.dot(M0,R.T)
            R0CSA0=np.dot(M0,RCSA.T)
            Reff0=np.reshape(np.dot(M,R.T).T-np.transpose([R00]),[1,np.prod(SZR)])
            ReffCSA0=np.reshape(np.dot(M,RCSA.T).T-np.transpose([R0CSA0]),[1,np.prod(SZR)])

            if iso:
                Reff+=A[k]*np.reshape(Reff0,SZeff)
                R0+=A[k]*np.reshape(R00,SZ0)
                ReffCSA+=A[k]*np.reshape(ReffCSA0,SZeff)
                R0CSA+=A[k]*np.reshape(R0CSA0,SZ0)
            else:
                A0=A[:,0,k]
                Reff+=np.reshape(np.dot(np.transpose([A0]),Reff0),SZeff)
                R0+=np.reshape(np.dot(np.transpose([A0]),[R00]),SZ0)
                A0=A[:,1,k]
                ReffCSA+=np.reshape(np.dot(np.transpose([A0]),ReffCSA0),SZeff)
                R0CSA+=np.reshape(np.dot(np.transpose([A0]),[R0CSA0]),SZ0)
            
        if iso:
            Reff=Reff[0]
            R0=R0[0]
            ReffCSA=ReffCSA[0]
            R0CSA=R0CSA[0]
        
        return Reff,R0,ReffCSA,R0CSA
        
    def z2zeff(self,tc):
        
        z=self.z()
        zeff=z+np.log10(tc)-np.log10(tc+10**z)  #Calculate the effective log-correlation time
        zeff[zeff<z[0]]=z[0]                    #Cleanup: no z shorter than z[0]
        zeff[zeff>=z[-1]]=z[-1]-1e-12           #Cleanup: no z longer than z[-1]
        i=np.digitize(zeff,z,right=False)-1     #Index to find longest z such that z<zeff
        sz=np.size(z)
        M=np.zeros([sz,sz])                     #Pre-allocate Matrix for rho->rho_eff transform
        
        dz=z[1:]-z[0:-1]
        wt=(z[i+1]-zeff)/dz[i]
        M[np.arange(0,sz),i]=wt
        M[np.arange(0,sz),i+1]=1-wt
        
        zi=np.log10(tc)                        #Calculate the log of input tc
        i=np.digitize(zi,z,right=False)-1     #Index to find longest z such that z<zeff
        M0=np.zeros([sz])                     #Pre-allocate Matrix for rho->rho_eff transform
        
        wt=(z[i+1]-zi)/dz[i]
        M0[i]=wt
        M0[i+1]=1-wt
        
        return M,M0
        
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
        
        if bond is None and np.size(a.shape)==3:
            maxi=np.max(a,axis=0)
            mini=np.min(a,axis=0)
            a=np.mean(a,axis=0)
            pltrange=True
        else:
            pltrange=False
        
        a=a.T
        maxi=maxi.T
        mini=mini.T
        if 'norm' in kwargs and kwargs.get('norm')[0].lower()=='y':
            norm=np.max(np.abs(a),axis=0)
            a=a/np.tile(norm,[np.size(self.tc()),1])      
            
            if pltrange:
                maxi=maxi/np.tile(norm,[np.size(self.tc()),1])
                mini=mini/np.tile(norm,[np.size(self.tc()),1])
        
        if ax==None:
            fig=plt.figure()
            ax=fig.add_subplot(111)
            hdl=ax.plot(self.z(),a)
#            hdl=plt.plot(self.z(),a)
#            ax=hdl[0].axes
        else:
            hdl=ax.plot(self.z(),a)
            
        if pltrange:
            x=np.concatenate([self.z(),self.z()[-1::-1]],axis=0)
            print(a.shape[0])
            for k in range(0,a.shape[1]):
                y=np.concatenate([mini[:,k],maxi[-1::-1,k]],axis=0)
                xy=np.concatenate(([x],[y]),axis=0).T
                patch=Polygon(xy,facecolor=hdl[k].get_color(),edgecolor=None,alpha=0.5)
                ax.add_patch(patch)
            
            
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