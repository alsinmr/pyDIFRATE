#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 20:32:25 2019

@author: albertsmith
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import mdl_sens as mdl
from numpy.linalg import svd
#from scipy.sparse.linalg import svds
from scipy.sparse.linalg import eigs
from scipy.optimize import linprog
from scipy.optimize import lsq_linear as lsqlin
import multiprocessing as mp

class detect(mdl.model):
    def __init__(self,sens,exp_num=None,mdl_num=None):
        """ We initiate the detectors class by giving it a sens/Ctsens class, from
        which it extracts the specified experiments and models for each 
        experiment.
        """
        
        self.n=None;
        self._class='detector'
        self._origin=sens._origin
        
        self.__tc=sens.tc()
        ntc=np.size(self.__tc)
        """This is the maximum number of singular values to return if not
        specified. Probably its higher than necessary.
        """
        self.__maxN=20;
        
        
        if np.size(exp_num)==1 and exp_num is None:
            exp_num=sens.info.columns.values
                    
        "Make sure we're working with numpy array for exp_num"
        exp_num=np.atleast_1d(exp_num)
        
        ne=np.size(exp_num)
        
        if mdl_num is None:
            mdl_num=-1

        "Make sure we're working with numpy array for mdl_num"
        mdl_num=np.atleast_1d(mdl_num)

#        "If all mdl_num are the same, replace with a single entry"
#        if np.size(mdl_num)>1 and np.all(mdl_num[0]==mdl_num):
#            mdl_num=mdl_num[0]
#        
        
        "Store all the experiment and model information"
        self.info_in=sens.info.loc[:,exp_num].copy()
        self.MdlPar_in=sens.MdlPar.copy()
        self.mdl_num=mdl_num.copy()
        
        
        k=0
        nM=np.size(self.MdlPar_in)
        while k<nM:
            if not(np.any(self.mdl_num==k)):
                del self.MdlPar_in[k]
                in1=np.where(np.array(self.mdl_num)!=None)[0]
                in2=np.where(self.mdl_num[in1]>k)[0]
                self.mdl_num[in1[in2]]+=-1
                nM=nM-1
            else:
                k=k+1  
        
        if np.all(self.mdl_num==-1):
            self.mdl_num=[]
        
        "Determine if any models are bond specific"
        self.BondSpfc='no'
        if sens._rho(bond=-1).ndim==3:
            self.BondSpfc='yes'     #If the previously applied models are bond-specific, we need to maintain bond speciicity
        else:
            for k in self.MdlPar_in:
                if k.get('BondSpfc')=='yes':
                    self.BondSpfc='yes'
        

        "Pass the molecule object"
        """Note that this is not a copy, but rather a pointer to the same object.
        If you edit the object, it will be changed in the original sens object
        Best to set the selections first with the sens object, and leave alone
        here"""
        self.molecule=sens.molecule                       

        
        "How many bonds are there?"
        if self.BondSpfc=='yes':
            nb=self.molecule.vXY.shape[0]
        else:
            nb=1
        
        "Storage for the input rate constants"
        self.__R=list()         #Store experimental sensitivities
        self.__R0=list()
        self.__RCSA=list()      #Store experimental sensitivities for CSA only
        self.__R0CSA=list()

        "Load in the sensitivity of the selected experiments"
        if np.size(mdl_num)==1:
            if self.BondSpfc=='yes':
#                for k in range(0,nb):
#                    a,b=sens._rho_eff(exp_num=exp_num,mdl_num=mdl_num[0],bond=k)
#                    self.__R.append(a)
#                    self.__R0.append(b)
#                    a,b=sens._rho_effCSA(exp_num=exp_num,mdl_num=mdl_num[0],bond=k)
#                    self.__RCSA.append(a)
#                    self.__R0CSA.append(b)
                a,b=sens._rho_eff(exp_num=exp_num,mdl_num=mdl_num[0],bond=-1)
                c,d=sens._rho_effCSA(exp_num=exp_num,mdl_num=mdl_num[0],bond=-1)
                for k in range(0,nb):
                    self.__R.append(a[k])
                    self.__R0.append(b[k])
                    self.__RCSA.append(c[k])
                    self.__R0CSA.append(d[k])
            elif mdl_num!=-1:
                a,b=sens._rho_eff(exp_num=exp_num,mdl_num=mdl_num[0])
                self.__R.append(a)
                self.__R0.append(b)
                a,b=sens._rho_effCSA(exp_num=exp_num,mdl_num=mdl_num[0])
                self.__RCSA.append(a)
                self.__R0CSA.append(b)
            else:
                self.__R.append(sens._rho(exp_num))
                self.__R0.append(np.zeros(self.__R[0].shape[0]))
                self.__RCSA.append(sens._rhoCSA(exp_num))
                self.__R0CSA.append(np.zeros(self.__RCSA[0].shape[0]))
        else:
            "In this case, we have to get the experiments one at a time"
            if self.BondSpfc=='yes':
                for k in range(0,nb):
                    self.__R.append(np.zeros([ne,ntc]))
                    self.__R0.append(np.zeros(ne))
                    self.__RCSA.append(np.zeros([ne,ntc]))
                    self.__R0CSA.append(np.zeros(ne))
                    for m in range(0,ne):
                        a,b=sens._rho_eff(exp_num=exp_num[m],mdl_num=mdl_num[m],bond=k)
                        self.__R[k][m,:]=a
                        self.__R0[k][m]=b
                        a,b=sens._rho_effCSA(exp_num=exp_num[m],mdl_num=mdl_num[m],bond=k)
                        self.__RCSA[k][m,:]=a
                        self.__R0CSA[k][m]=b
            else:
                self.__R.append(np.zeros([ne,ntc]))
                self.__R0.append(np.zeros(ne))
                self.__RCSA.append(np.zeros([ne,ntc]))
                self.__R0CSA.append(np.zeros(ne))
                for m in range(0,ne):
                    a,b=sens._rho_eff(exp_num=exp_num[m],mdl_num=mdl_num[m])
                    self.__R[0][m,:]=a
                    self.__R0[0][m]=b
                    a,b=sens._rho_effCSA(exp_num=exp_num[m],mdl_num=mdl_num[m])
                    self.__RCSA[0][m,:]=a
                    self.__R0CSA[0][m]=b
       
        "Names of the experimental variables that are available"
        self.__exper=['rho','z0','z0_std','Del_z','Del_z_std','stdev']
        "Names of the spin system variables that are available"
        self.__spinsys=[]
        "Initialize self.info"
        self.info=None
        
        self.detect_par=dict()
        self.detect_par['Normalization']='M'
        self.detect_par['inclS2']='no'
        self.detect_par['NegAllow']=0.5
        
        "Pass the normalization"
        a=self.info_in.loc['stdev'].to_numpy()   
        "Replace None with max of abs of sensitivity"
        index=a==None
        b=np.max(np.abs(np.mean(self.__R,axis=0)),axis=-1)
        a[index]=b[index]
        self.norm=np.divide(1,a).astype('float64')

        "Storage for the detection vectors"
        self.__r=list(np.zeros(nb))   #Store the detection vectors
        self.__rho=list(np.zeros(nb))   #Store the detector sensitivities
        self.__rhoAvg=None
        self.__rAvg=None
        self.__Rc=list(np.zeros(nb))        #Store the back-calculated sensitivities
        self.__RcAvg=None
        self.__rhoCSA=list(np.zeros(nb)) #CSA only sensitivities
        
        "Store SVD matrix for parallel function"
        self.__Vt=None

        self.z0=list(np.zeros(nb))
        self.Del_z=list(np.zeros(nb))
        
        self.SVD=list(np.zeros(nb))
        self.SVDavg=dict()
        
        "Store error for r_auto routine"
        self.__r_auto=dict()
        
        super().__init__()

        "The previous line clears self.molecule, so we have to load it again :-/"        
        self.molecule=sens.molecule        
        

#%% Performs and stores results of singular value decomposition of experimental sensitivities    
    def getSVD(self,bond=None,n=None):
        "Function to perform (and store) all singular value decomposition calculations"
        ne=np.shape(self.__R)[1]
        if n is None:
            n=np.min([np.shape(self.__R)[1],self.__maxN])
        
        if bond is None:
            if 'S' in self.SVDavg.keys() and self.SVDavg['S'].size>=n:
                U=self.SVDavg['U'][:,0:n]
                S=self.SVDavg['S'][0:n]
                Vt=self.SVDavg['Vt'][0:n,:]
                VtCSA=0
            else:
                norm=np.repeat(np.transpose([self.norm]),np.size(self.__tc),axis=1)
                U,S,Vt=svd0(np.multiply(np.mean(self.__R,axis=0),norm),n)
                
                self.SVDavg['U']=U
                self.SVDavg['Vt']=Vt
                self.SVDavg['S']=S
                VtCSA=0
                
                U=U[:,0:n]
                S=S[0:n]
                Vt=Vt[0:n,:]
        else:
            if self.SVD[bond]!=0 and self.SVD[bond]['S'].size>=n:
                U=self.SVD[bond]['U'][:,0:n]
                S=self.SVD[bond]['S'][0:n]
                Vt=self.SVD[bond]['Vt'][0:n,:]
                VtCSA=self.SVD[bond]['VtCSA'][0:n,:]
                
            else:                    
                norm=np.repeat(np.transpose([self.norm]),np.size(self.__tc),axis=1)

                U,S,Vt=svd0(np.multiply(self.Rin(bond),norm),n)
                U=U[:,0:np.size(S)]
                    
                
                VtCSA=np.dot(np.diag(1/S),np.dot(U.T,np.multiply(self._RCSAin(bond),norm)))
                
                if self.SVD[bond]==0:
                    self.SVD[bond]=dict()
                
                self.SVD[bond]['U']=U
                self.SVD[bond]['S']=S
                self.SVD[bond]['Vt']=Vt
                self.SVD[bond]['VtCSA']=VtCSA
                
                U=U[:,0:n]
                S=S[0:n]
                Vt=Vt[0:n,:]
                VtCSA=VtCSA[0:n,:]
    
        return U,S,Vt,VtCSA
#%% Generate r matrix for fitting tests (detector sensitivies are not optimized- and not well-separated)
    def r_no_opt(self,n,bond=None,**kwargs):
        
        self.n=n
        nb=np.shape(self.__R)[0]

        if nb==1:
            bond=0
            
        if bond is not None and np.size(bond)==1 and np.atleast_1d(bond)[0]==-1:
            bond=np.arange(0,nb)
            
        if bond is None:
            U,S,Vt,VCSA=self.getSVD(None,n)
            self.__rAvg=np.multiply(np.repeat(np.transpose([1/self.norm]),n,axis=1),np.dot(U,np.diag(S)))
            self.__rhoAvg=Vt
            norm=np.repeat(np.transpose([self.norm]),np.size(self.__tc),axis=1)
            self.__RcAvg=np.divide(np.dot(U,np.dot(np.diag(S),Vt)),norm)
            self.SVDavg['T']=np.eye(n)
            self.SVDavg['stdev']=1/S
            if 'sort_rho' not in kwargs:
                kwargs['sort_rho']='n'

            self.__r_info(None,**kwargs)
        else:
            
            bond=np.atleast_1d(bond)

            for k in bond:
                U,S,Vt,VCSA=self.getSVD(k,n)
                self.__r[k]=np.multiply(np.repeat(np.transpose([1/self.norm]),n,axis=1),np.dot(U,np.diag(S)))
                self.__rho[k]=Vt
                self.__rhoCSA[k]=VCSA
                norm=np.repeat(np.transpose([self.norm]),np.size(self.__tc),axis=1)
                self.__Rc[k]=np.divide(np.dot(U,np.dot(np.diag(S),Vt)),norm)
                self.SVD[k]['T']=np.eye(n)
                self.__r_info(k,**kwargs)
                
            if 'sort_rho' not in kwargs:
                kwargs['sort_rho']='n'
            self.__r_info(None,**kwargs)

#%% Automatic generation of detectors from a set of sensitivities                 
    def r_auto(self,n,bond=None,**kwargs):
        self.n=n
        "Get input or defaults"
        if 'inclS2' in kwargs:
            inclS2=kwargs.get('inclS2')
            self.detect_par['inclS2']=inclS2
        else:
            inclS2=self.detect_par.get('inclS2')
        if 'Neg' in kwargs:
            Neg=kwargs.get('Neg')
            self.detect_par['NegAllow']=Neg
        else:
            Neg=self.detect_par.get('NegAllow')
          
            
        nb=np.shape(self.__R)[0]
        "If bond set to -1, run through all orientations."
        if bond is None:
            bonds=np.zeros(0)
        elif np.size(bond)==1 and np.atleast_1d(bond)[0]==-1:
            bond=None
            bonds=np.arange(0,nb)
        else:
            bond=np.atleast_1d(bond)
            bonds=bond[1:]
            bond=bond[0]
            
            
        if nb==1:
            "If we only have one set of sensitivities (that is, no orientation dependence, then don't use averages"
            bond=0
            
        if bond is None:
            "Here we operate on the average sensitivities"
            U,S,Vt,VCSA=self.getSVD(None,n)
            norm=np.repeat(np.transpose([self.norm]),np.size(self.__tc),axis=1)
            self.__RcAvg=np.divide(np.dot(U,np.dot(np.diag(S),Vt)),norm)
        else:                
            "We work on the first bond given, and use r_target for the remaining bonds"
            U,S,Vt,VCSA=self.getSVD(bond,n)
            norm=np.repeat(np.transpose([self.norm]),np.size(self.__tc),axis=1)
            self.__Rc[bond]=np.divide(np.dot(U,np.dot(np.diag(S),Vt)),norm)
     
        ntc=np.size(self.__tc) #Number of correlation times
        err=np.zeros(ntc)       #Error of fit
        

        "Prepare data for parallel processing"
        Y=list()
        for k in range(0,ntc):
            Y.append((Vt,k))

        "Default is parallel processing"
        if 'parallel' in kwargs and kwargs.get('parallel').lower()[0]=='n':
            X=list()
            for k in range(0,ntc):
                X.append(linprog_par(Y[k]))
        else:
            with mp.Pool() as pool:
                X=pool.map(linprog_par,Y)
            
        """We optimized detectors at every correlation time (see __linprog_par),
        which have values at 1 for the given correlation time. We want to keep 
        those detectors where the maximum is closest to the correlation time set
        to 1. We search for those here:
        """
        for k in range(0,ntc):
            err[k]=np.abs(np.argmax(np.dot(Vt.T,X[k]))-k)
        
        
        """Ideally, the number of detectors equals the number of minima in err,
        however, due to calculation error, this may not always be the case. We
        start searching for values where err=0. If we don't have enough, we
        raise this value in steps (looking for err=1, err=2), until we have
        enough. If, in one step, we go from too few to too many, we eliminate, 
        one at a time, the peak that is closest to another peak.
        """
        test=True
        thresh=0
        while test:
            pks=np.where(err<=thresh)[0]
            if pks.size==n:
                test=False
            elif pks.size<n:
                thresh=thresh+1
            elif pks.size>n:

                while pks.size>n:
                    a=np.argsort(np.diff(pks))
#                    pks=np.concatenate([pks[a[np.size(pks)-n:]],[pks[-1]]])
                    pks=np.concatenate([pks[a[1:]],[pks[-1]]])
                    pks.sort()
                test=False
        
        "Save the linear combinations for the best detectors"
        T=np.zeros([n,n])
        for k in range(0,n):
            T[k,:]=X[pks[k]]
        
        rhoz=np.dot(T,Vt)
        
        plt.plot(self.z(),rhoz[-1])
        """Detectors that are not approaching zero at the end of the range of
        correlation times tend to oscillate where they do approach zero. We want
        to push that oscillation slightly below zero
        """
        
        for k in range(0,n):
            if (rhoz[k,0]>0.95*np.max(rhoz[k,:]) or rhoz[k,-1]>0.95*np.max(rhoz[k,:])) and Neg!=0:

                reopt=True #Option to cancel the re-optimization in special cases
                
                if rhoz[k,0]>0.95*np.max(rhoz[k,:]):
                    pm=1;
                else:
                    pm=-1;                        
                
                temp=rhoz[k,:]
                "Locate maxima and minima in the detector"
                mini=np.where((temp[2:]-temp[1:-1]>=0) & (temp[1:-1]-temp[0:-2]<=0))[0]+1
                maxi=np.where((temp[2:]-temp[1:-1]<=0) & (temp[1:-1]-temp[0:-2]>=0))[0]+1
                
                """Filter out minima that occur at more than 90% of the sensitivity max,
                since these are probably just glitches in the optimization.
                """
                if np.size(mini)>=2 and np.size(maxi)>=2:
                    mini=mini[(temp[mini]<.9) & (temp[mini]<.05*np.max(-pm*np.diff(temp[maxi])))]
                elif np.size(mini)>=2:
                    mini=mini[temp[mini]<.9]
                    
                if np.size(maxi)>=2:
                    maxi=maxi[(temp[maxi]<.9) & (temp[maxi]>0.0)]
#                    maxi=maxi[(temp[maxi]<.9) & (temp[maxi]>0.0*np.max(-pm*np.diff(temp[maxi])))]
                
                
                if rhoz[k,0]>0.95*np.max(rhoz[k,:]):
                    "Calculation for the first detection vector"

                    if np.size(maxi)>=2 & np.size(mini)>=2:
                        step=int(np.round(np.diff(mini[0:2])/2))
                        slope2=-(temp[maxi[-1]]-temp[maxi[0]])*Neg/(maxi[-1]-maxi[0])
                    elif np.size(maxi)==1 and np.size(mini)>=1:
                        step=maxi[0]-mini[0]
                        slope2=temp[maxi[0]]*Neg/step
                    else:
                        reopt=False
                        
                    if reopt:
                        a=np.max([1,mini[0]-step])
                        slope1=-temp[maxi[0]]/step*Neg
                        line1=np.arange(0,-temp[maxi[0]]*Neg-1e-12,slope1)
                        line2=np.arange(-temp[maxi[0]]*Neg,1e-12,slope2)
                        try:
                            target=np.concatenate((np.zeros(a),line1,line2,np.zeros(ntc-a-np.size(line1)-np.size(line2))))
                        except:
                            reopt=False
                                
                else:
                    "Calculation otherwise (last detection vector)"
                    if np.size(maxi)>=2 & np.size(mini)>=2:
                        step=int(np.round(np.diff(mini[-2:])/2))
                        slope2=-(temp[maxi[0]]-temp[maxi[-1]])*Neg/(maxi[0]-maxi[-1])
                    elif np.size(maxi)==1 and np.size(mini)>=1:
                        step=mini[-1]-maxi[0]
                        slope2=-temp[maxi[0]]*Neg/step
                    else:
                        reopt=False
                        
                    if reopt:
                        a=np.min([ntc,mini[-1]+step])
                        slope1=temp[maxi[-1]]/step*Neg
    
                        line1=np.arange(-temp[maxi[-1]]*Neg,1e-12,slope1)
                        line2=np.arange(0,-temp[maxi[-1]]*Neg-1e-12,slope2)                    
                        target=np.concatenate((np.zeros(a-np.size(line1)-np.size(line2)),line2,line1,np.zeros(ntc-a)))
                    

                if reopt:
                    Y=(Vt,pks[k],target)
                    
                    X=linprog_par(Y)
                    T[k,:]=X
                    rhoz[k,:]=np.dot(T[k,:],Vt)
                 
        "Save the results into the detect object"
#        self.r0=self.__r
        if bond is None:
            self.__rAvg=np.multiply(np.repeat(np.transpose([1/self.norm]),n,axis=1),\
                        np.dot(U,np.linalg.solve(T.T,np.diag(S)).T))
            self.__rhoAvg=rhoz
            self.SVDavg['T']=T
            self.__r_info(bond,**kwargs)
            self.__r_auto={'Error':err,'Peaks':pks,'rho_z':self.__rhoAvg}
            self.__r_norm(bond,**kwargs)
            if np.size(bonds)>0:
                if 'NT' in kwargs: #We don't re-normalize the results of detectors obtained with r_target
                    kwargs.pop('NT')
                if 'Normalization' in kwargs:
                    kwargs.pop('Normalization')
                self.r_target(n,self.rhoz(None),bonds,**kwargs)
        else:

            """This isn't correct yet- if more than one bond, we want to 
            use the result for the average calculation as a target for the 
            individual bonds, not loop over all bonds with the result here
            """
            self.__r[bond]=np.multiply(np.repeat(np.transpose([1/self.norm]),n,axis=1),\
                        np.dot(U,np.linalg.solve(T.T,np.diag(S)).T))
            self.__rho[bond]=rhoz
            self.__rhoCSA[bond]=np.dot(T,VCSA)
            self.SVD[bond]['T']=T
            self.__r_info(bond,**kwargs)
            self.__r_auto={'Error':err,'Peaks':pks,'rho_z':self.__rho[bond]}
            self.__r_norm(bond,**kwargs)
            if np.size(bonds)>0:
                if 'NT' in kwargs: #We don't re-normalize the results of detectors obtained with r_target
                    kwargs.pop('NT')
                if 'Normalization' in kwargs:
                    kwargs.pop('Normalization')
                self.r_target(n,self.__rho[bond],bonds,**kwargs)

    def r_target(self,n,target,bond=None,**kwargs):
        "Set sensitivities as close to some target function as possible"
        
        target=np.atleast_2d(target)
        
        nb=np.shape(self.__R)[0]
        
        "If bond set to -1, run through all orientations."
        if bond is not None and np.size(bond)==1 and np.atleast_1d(bond)[0]==-1:
            bond=np.arange(0,nb)
                
        if nb==1:
            bond=0
            

        if bond is None:
            "Here we operate on the average sensitivities"
            U,S,Vt,VCSA=self.getSVD(None,n)
            norm=np.repeat(np.transpose([self.norm]),np.size(self.__tc),axis=1)
            self.__RcAvg=np.divide(np.dot(U,np.dot(np.diag(S),Vt)),norm)
            
            T=lsqlin_par((Vt,target))
            
            rhoz=np.dot(T,Vt)
            self.__rAvg=np.multiply(np.repeat(np.transpose([1/self.norm]),n,axis=1),\
                        np.dot(U,np.linalg.solve(T.T,np.diag(S)).T))
            self.__rhoAvg=rhoz
            self.SVDavg['T']=T
            if 'NT' in kwargs:
                self.__r_norm(None,**kwargs)            
            self.__r_info(bond,**kwargs)
        else:
            Y=list()
            bond=np.atleast_1d(bond)
            
            for k in bond:
                U,S,Vt,VCSA=self.getSVD(k,n)
                norm=np.repeat(np.transpose([self.norm]),np.size(self.__tc),axis=1)
                self.__Rc[k]=np.divide(np.dot(U,np.dot(np.diag(S),Vt)),norm)
                Y.append((Vt,target))
                
            "Default is parallel processing"
            if 'parallel' in kwargs and kwargs.get('parallel').lower()[0]=='n':
                T=list()
                for k in Y:
                    T.append(lsqlin_par(k))
            else:
                with mp.Pool() as pool:
                    T=pool.map(lsqlin_par,Y)
                    
            for index,k in enumerate(bond):
                U,S,Vt,VCSA=self.getSVD(k,n)
                self.SVD[k]['T']=T[index]
                self.__r[k]=np.multiply(np.repeat(np.transpose([1/self.norm]),n,axis=1),\
                    np.dot(U,np.linalg.solve(T[index].T,np.diag(S)).T))
                self.__rho[k]=np.dot(T[index],Vt)
                self.__rhoCSA[k]=np.dot(T[index],VCSA)
                self.SVD[k]['T']=T[index]
                if 'NT' in kwargs:
                    self.__r_norm(None,**kwargs)   
                
        if 'sort_rho' not in kwargs:
            kwargs['sort_rho']='n'
        self.__r_info(bond,**kwargs)
    
    
    def __addS2(self,bond=None,**kwargs):
        if 'NT' in kwargs:
            NT=kwargs.get('NT')
        elif 'Normalization' in kwargs:
            NT=kwargs.get('Normalization')
        else:
            NT=self.detect_par.get('Normalization')
            
        pass
    
    def __R2_ex_corr(self,bond=None):
        index=self.info_in.loc['Type']=='R2'
        
        pass
    
    def __r_norm(self,bond=None,**kwargs):
        "Applies equal-max or equal-integral normalization"
        if 'NT' in kwargs:
            NT=kwargs.get('NT')
            self.detect_par['Normalization']=NT
        elif 'Normalization' in kwargs:
            NT=kwargs.get('Normalization')
            self.detect_par['Normalization']=NT
        else:
            NT=self.detect_par.get('Normalization')
            
        nb=np.shape(self.__R)[0]
        if nb==1:
            bond=0
        
        rhoz=self.rhoz(bond)
        
        NT=self.detect_par.get('Normalization')
        nd=self.n
        
        dz=np.diff(self.z()[0:2])
        
        for k in range(0,nd):
            if NT.upper()[0]=='I':
                sc=np.sum(rhoz[k,:])*dz
            elif NT.upper()[0]=='M':
                sc=np.max(rhoz[k,:])
            else:
                print('Normalization type not recognized (use "N" or "I"). Defaulting to equal-max')
                sc=np.max(rhoz[k,:])                

            if bond is None:
                self.__rhoAvg[k,:]=rhoz[k,:]/sc
                self.SVDavg['T'][k,:]=self.SVDavg['T'][k,:]/sc
                self.__rAvg[:,k]=self.__rAvg[:,k]*sc
            else:
                self.__rho[bond][k,:]=rhoz[k,:]/sc
                self.__rhoCSA[bond][k,:]=self.__rhoCSA[bond][k,:]/sc
                self.SVD[bond]['T'][k,:]=self.SVD[bond]['T'][k,:]
                self.__r[bond][:,k]=self.__r[bond][:,k]*sc
                
                
            
    
    def __r_info(self,bond=None,**kwargs):
        """Calculates some parameters related to the detectors generates, z0,
        Del_z, and standard deviation of resulting detectors. Also resorts the
        detectors according to z0
        """
        nb=np.shape(self.__R)[0]
        

        match=True
        if self.__r[0].ndim==2:
            nd0=np.shape(self.__r[0])[1]
        else:
            match=False
            nd0=0
            
        stdev=np.zeros(nd0)
        
        if bond is None:
#            a=np.arange(0,nb)
            a=np.arange(0,0)
            match=False
        else:
            a=np.atleast_1d(bond)
            

        for k in a:
            if self.__r[0].ndim==2:
                nd=np.shape(self.__r[k])[1]
            else:
                nd=0
                            
            
            if nd0!=nd:
                match=False
            
            if nd>0:
                z0=np.divide(np.sum(np.multiply(self.__rho[k],\
                        np.repeat([self.z()],nd,axis=0)),axis=1),\
                        np.sum(self.__rho[k],axis=1))
    
                if 'sort_rho' in kwargs and kwargs.get('sort_rho').lower()[0]=='n':
                    i=np.arange(0,np.size(z0))
                else:
                    i=np.argsort(z0)
                
                self.z0[k]=z0[i]
                self.SVD[k]['T']=self.SVD[k]['T'][i,:]
                self.__r[k]=self.__r[k][:,i]
                self.__rho[k]=self.__rho[k][i,:]
                self.Del_z[k]=np.diff(self.z()[0:2])*np.divide(np.sum(self.__rho[k],axis=1),
                          np.max(self.__rho[k],axis=1))
                self.SVD[k]['stdev']=np.sqrt(np.dot(self.SVD[k]['T']**2,1/self.SVD[k]['S'][0:nd]**2))
                
                if match:
                    stdev+=self.SVD[k]['stdev']
                
        if match:
            a=dict()
            a.update({'z0' : np.mean(self.z0,axis=0)})
            a.update({'Del_z' : np.mean(self.Del_z,axis=0)})
            a.update({'stdev' : stdev/nb})
            if nb>1:
                a.update({'z0_std' : np.std(self.z0,axis=0)})
                a.update({'Del_z_std': np.std(self.Del_z,axis=0)})
            else:
                a.update({'z0_std' : np.zeros(nd)})
                a.update({'Del_z_std' : np.zeros(nd)})
                
            self.info=pd.DataFrame.from_dict(a)
            self.info=self.info.transpose()
            
            
            
    
    def r(self,bond=None):
        nb=np.shape(self.__R)[0]
        if nb==1:
            bond=0
            
        if bond is None:
            if self.__rAvg is None:
                print('First generate the detectors for the average sensitivities')
            else:
                return self.__rAvg
        else:
            if np.size(self.__r[bond])==1:
                print('First generate the detectors for the selected bond')
            else:
                return self.__r[bond]
    
    def rhoz(self,bond=None):
        nb=np.shape(self.__R)[0]
        if nb==1:
            bond=0
            
        if bond is None:
            if self.__rAvg is None:
                print('First generate the detectors for the average sensitivities')
            else:
                return self.__rhoAvg
        else:
            if np.size(self.__rho[bond])==1:
                print('First generate the detectors for the selected bond')
            else:
                if bond==-1:
                    return np.array(self.__rho)
                else:
                    return self.__rho[bond]
            
    def Rc(self,bond=None):
        nb=np.shape(self.__R)[0]
        if nb==1:
            bond=0
            
        if bond is None:
            if self.__RcAvg is None:
                print('First generate the detectors to back-calculate rate constant sensitivities')
            else:
                return self.__RcAvg
        else:
            if np.size(self.__Rc[bond])==1:
                print('First generate the detectors for the selected bond')
            else:
                return self.__Rc[bond]
            
    
    
    def Rin(self,bond=0):
        nb=np.shape(self.__R)[0]
        if nb==1:
            bond=0
        return self.__R[bond]
    
    def R0in(self,bond=0):
        nb=np.shape(self.__R)[0]
        if nb==1:
            bond=0
        return self.__R0[bond]
    
    def _RCSAin(self,bond=0):
        nb=np.shape(self.__R)[0]
        if nb==1:
            bond=0
        return self.__RCSA[bond]
    
    def _R0CSAin(self,bond=0):
        nb=np.shape(self.__R)[0]
        if nb==1:
            bond=0
        return self.__R0CSA[bond]
    
    def retExper(self):
        return self.__exper
    
    def retSpinSys(self):
        return self.__spinsys
        
    def tc(self):
        return self.__tc
    
    def z(self):
        return np.log10(self.__tc)
    
    def _rho(self,exp_num=None,bond=None):
        """The different children of mdl_sens will have different names for 
        their sensitivities. For example, this class returns rho_z, which are the 
        rate constant sensitivities, but the correlation function class returns
        Ct, and the detector class returns rho. Then, we have a function, 
        _rho(self), that exists and functions the same way in all children
        """
        

             
        if bond is None:
            bond=0
        
        if np.size(self.__rho[bond])==1:
            print('First generate the detectors for the selected bond')
            return
                
        if exp_num is None:
            exp_num=self.info.columns
        
        exp_num=np.atleast_1d(exp_num)
        
        
        if bond==-1:
            rhoz=self.rhoz(bond)
            if rhoz.ndim==3:
                rhoz=rhoz[:,exp_num,:]
            elif rhoz.ndim==2:
                rhoz=rhoz[exp_num,:]
        else:
            rhoz=self.rhoz(bond)[exp_num,:]
                
        return rhoz
    
    def _rhoCSA(self,exp_num=None,bond=None):
        """The different children of mdl_sens will have different names for 
        their sensitivities. For example, this class returns R, which are the 
        rate constant sensitivities, but the correlation function class returns
        Ct, and the detector class returns rho. Then, we have a function, 
        _rho(self), that exists and functions the same way in all children
        """
        
        if bond is None:
            bond=0
        
        if np.size(self.__rhoCSA[bond])==1:
            print('First generate the detectors for the selected bond')
            return
        
        
        if exp_num is None:
            exp_num=self.info.columns
        
        exp_num=np.atleast_1d(exp_num)
        
        
        if bond==-1:
            rhoz=np.array(self.__rhoCSA)
            if rhoz.ndim==3:
                rhoz=rhoz[:,exp_num,:]
            elif rhoz.ndim==2:
                rhoz=rhoz[exp_num,:]
            
            if rhoz.shape[0]==1:
                rhoz=rhoz[0]
        else:
            rhoz=self.__rhoCSA[bond]
            rhoz=rhoz[exp_num,:]
            
        
        return rhoz
    
    def plot_rhoz(self,bond=None,rho_index=None,ax=None,**kwargs):
        "Create an axis if not given"
        if ax is None:
            fig=plt.figure()
            ax=fig.add_subplot(111)
            
        nb=np.shape(self.__R)[0]
        if nb==1:
            bond=0
            
        if rho_index is None:
            "Get all rho"
            if np.size(bond)==1:
                nd=self.rhoz(bond).shape[-2]
            else:
                nd=self.rhoz(bond[0]).shape[-2]
            rho_index=np.arange(0,nd)
            

            
        if bond==-1:
            a=self.rhoz(bond=None).T
        else:                    
            a=self.rhoz(bond).T
        a=a[:,rho_index]

        hdl=ax.plot(self.z(),a)
        
        if bond is not None and np.size(bond)==1 and np.atleast_1d(bond)[0]==-1:
            bond=None
            ntc=self.tc().size
            maxi=np.zeros([nd,ntc])
            mini=1e12*np.ones([nd,ntc])
            for k in range(0,nb):
                maxi=np.max(np.concatenate(([maxi],[self.rhoz(k)]),axis=0),axis=0)
                mini=np.min(np.concatenate(([mini],[self.rhoz(k)]),axis=0),axis=0)

            x=np.concatenate((self.z(),self.z()[-1::-1]),axis=0)              
            for k in rho_index:
                y=np.concatenate((mini[k,:],maxi[k,-1::-1]),axis=0)
                xy=np.concatenate(([x],[y]),axis=0).T
                patch=Polygon(xy,facecolor=hdl[k].get_color(),alpha=0.5)
                
                ax.add_patch(patch)
        
        self._set_plot_attr(hdl,**kwargs)
        
        ax.set_xlabel(r'$\log_{10}(\tau$ / s)')
        ax.set_ylabel(r'$\rho_n(z)$')
        ax.set_xlim(self.z()[[0,-1]])
        
        
        return hdl
    
    def plot_r_opt(self,fig=None):
        if fig is None:
            fig=plt.figure()
        
        ax1=fig.add_subplot(211)
        ax1.plot(self.z(),self.__r_auto.get('Error'))
        max_err=self.__r_auto.get('Error').max()
        
        ax2=fig.add_subplot(212)
        hdls=ax2.plot(self.z(),self.__r_auto.get('rho_z').T)
        
        for index, k in enumerate(np.sort(self.__r_auto.get('Peaks'))):
            ax1.plot(np.repeat(self.z()[k],2),[-max_err/20,0],linewidth=.5,color=[0,0,0])
            ax1.text(self.z()[k],-max_err*1/10,r'$\rho_{'+str(index+1)+'}$',horizontalalignment='center',\
            verticalalignment='center',color=hdls[index].get_color())
            ax2.text(self.z()[k],self.__r_auto.get('rho_z')[index,:].max()+.05,\
            r'$\rho_{'+str(index+1)+'}$',horizontalalignment='center',\
            verticalalignment='center',color=hdls[index].get_color())
        
        ax1.set_xlabel(r'$\log_{10}(\tau$ / s)')
        ax1.set_ylabel(r'Opt. Error, $\Delta$(max)')
        ax1.set_xlim(self.z()[[0,-1]])
        ax1.set_ylim([-max_err*3/20,max_err*21/20])
        
        ax2.set_xlabel(r'$\log_{10}(\tau$ / s)')
        ax2.set_ylabel(r'$\rho_n(z)$')
        min_rho=self.__r_auto.get('rho_z').min()
        max_rho=self.__r_auto.get('rho_z').max()
        
        ax2.set_xlim(self.z()[[0,-1]])
        ax2.set_ylim([min_rho-.05,max_rho+.1])
        
        return hdls 
        
        
        
    def plot_Rc(self,bond=None,exp_num=None,ax=None,**kwargs):
        nb=np.shape(self.__R)[0]
        if nb==1:
            bond=0
                    
        a=self.Rin(bond).T
        b=self.Rc(bond).T
        
        
        
        if np.size(exp_num)>1 or exp_num!=None:
            a=a[:,exp_num]
            b=b[:,exp_num]
        
        if 'norm' in kwargs and kwargs.get('norm')[0].lower()=='y':
            norm=np.max(np.abs(a),axis=0)
            a=a/np.tile(norm,[np.size(self.tc()),1]) 
            b=b/np.tile(norm,[np.size(self.tc()),1])
        
        
        if ax is None:
            fig=plt.figure()
            ax=fig.add_subplot(111)
            hdl1=ax.plot(self.z(),a,'k')
            hdl2=ax.plot(self.z(),b,'r--')
#            hdl1=plt.plot(self.z(),a,'k')
#            hdl2=plt.plot(self.z(),b,'r--')
#            ax=hdl1[0].axes
        else:
            hdl1=ax.plot(self.z(),a,'k')
            hdl2=ax.plot(self.z(),b,'r--')
            
        ax.set_xlabel(r'$\log_{10}(\tau$ / s)')
        ax.set_ylabel(r'$R(z) / $s$^{-1}$')
        ax.set_xlim(self.z()[[0,-1]])
        ax.set_title('Rate Constant Reproduction')
        
        hdl=hdl1+hdl2
        
        return hdl

def svd0(X,n):
    if np.shape(X)[0]>np.shape(X)[1]:
#        U,S,Vt=svds(X,k=n,tol=0,which='LM')    #Large data sets use sparse svd to avoid memory overload
#        U=U[:,-1::-1]      #svds puts out eigenvalues in opposite order of svd
#        S=S[-1::-1]
#        Vt=Vt[-1::-1,:]
        S2,V=eigs(np.dot(np.transpose(X),X),k=n)
        S=np.sqrt(S2.real)
        U=np.dot(np.dot(X,V.real),np.diag(1/S))
        Vt=V.real.T
    else:
        U,S,Vt=svd(X)       #But, typically better results from full calculation
        U=U[:,0:np.size(S)] #Drop all the empty vectors
    
    return U,S,Vt
    
def linprog_par(Y):
    """This function optimizes a detector sensitivity that has a value of 1
    at correlation time k, and cannot go below some specific value at all
    other correlation times (usually 0). While satisfying these 
    requirements, the sensitivity is minimized.
    """
    Vt=Y[0]
    k=Y[1]
    ntc=np.shape(Vt)[1]
    
    if np.size(Y)==3:
        target=Y[2]
    else:
        target=np.zeros(ntc)
        
    try:
        x=linprog(np.sum(Vt,axis=1),-Vt.T,-target,[Vt[:,k]],1,bounds=(-500,500),method='interior-point',options={'disp' :False,})
        x=x['x']
        if np.any(np.dot(Vt.T,x)<(np.min(target)-.0001)):
            "This is sketchy. Linprog is return np.dot(Vt.T,x)<-1, but yields success. Shouldn't happend"
            x=np.ones(Vt.shape[0])
    except:
        x=np.ones(Vt.shape[0])
        
    return x

def lsqlin_par(Y):
    Vt=Y[0]
    target=Y[1]
    
    nSVD=np.shape(Vt)[0]
    n=np.shape(target)[0]
    
    T=np.zeros([nSVD,nSVD])
    

    for k in range(0,n):        
        x0=lsqlin(np.transpose(Vt),target[k,:],lsq_solver='exact')
        T[k,:]=x0['x']
        
    for k in range(n,nSVD):
#        a=np.argmin(np.sum(T**2,axis=0))
#        T[k,a]=1
        T[k,k]=1
    
    return T