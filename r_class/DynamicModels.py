#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 14:51:02 2019

@author: albertsmith
"""

"""Here we store all models for motion. We keep the basic input in the mdl_sens
class very general, so you can add whatever parameters desired here (pass with 
**kwargs)
Note that we furthermore can have access to the structure, imported via 
the mdanalysis module 
"""

import numpy as np
from numpy import inf

def ModelBondSpfc(Model):
    "Add models with bond-specific dynamics to this list"
    mdls=np.array(['AnisoDif']) 
    return np.any(Model==mdls)
    
    
def ModelSel(Model,direct='dXY',structure=None,**kwargs):
    if Model=='IsoDif':
        tMdl,AMdl=IsoDif(**kwargs)
        BndSpfc='no'
        """We must always say if the model is bond specific, so we know if an
        array of models is being returned"""
    elif Model=='AnisoDif':
        tMdl,AMdl=AnisoDif(structure,direct,**kwargs)
        BndSpfc='yes'
    elif Model=='Combined':
        tMdl,AMdl=Combined(tMdl1=kwargs.get('tMdl1'),AMdl1=kwargs.get('AMdl1'),\
                           tMdl2=kwargs.get('tMdl2'),AMdl2=kwargs.get('AMdl2'))
        BndSpfc='no'

    "Make sure we return np arrays with a dimension, and sort the correlation times"
    if not isinstance(tMdl,np.ndarray):
        tMdl=np.array(tMdl)
    if tMdl.shape==():
        tMdl=np.array([tMdl])
    if not isinstance(AMdl,np.ndarray):
        AMdl=np.array(AMdl)
    if AMdl.shape==():
        AMdl=np.array([AMdl])
        
    tMdl[tMdl==inf]=1000
    return tMdl,AMdl,BndSpfc

"Simple isotropic diffusion"
def IsoDif(**kwargs):

    if 'tM' in kwargs:        
        tMdl=kwargs.get('tM')
    elif 'tm' in kwargs:
        tMdl=kwargs.get('tm')
    elif 'tr' in kwargs:
        tMdl=kwargs.get('tr')
    elif 'tR' in kwargs:
        tMdl=kwargs.get('tR')
        
    AMdl=1
    
    return tMdl,AMdl
    
def AnisoDif(struct,direct='vXY',**kwargs):
   
    """First we get the diffusion tensor, and also Diso and D2. This can be 
    input either as the principle values, Dxx, Dyy, and Dzz, or as the trace of
    the tensor (the isotropic value, tM), plus optionally the anisotropy, xi, 
    and the asymmetry, eta
    """
    if 'Dxx' in kwargs and 'Dyy' in kwargs and 'Dzz' in kwargs:
        Dzz=kwargs.get('Dzz')
        Dxx=kwargs.get('Dxx')
        Dyy=kwargs.get('Dyy')
        Diso=1/3*(Dxx+Dyy+Dzz)
        Dsq=(Dxx*Dyy+Dyy*Dzz+Dzz*Dxx)/3;
    else:
        if 'tM' in kwargs:        
            tM=kwargs.get('tM')
        elif 'tm' in kwargs:
            tM=kwargs.get('tm')
        elif 'tr' in kwargs:
            tM=kwargs.get('tr')
        elif 'tR' in kwargs:
            tM=kwargs.get('tR')
            
        if 'xi' in kwargs:
            xi=kwargs.get('xi')
        else:
            xi=1
        if 'eta' in kwargs:
            eta=kwargs.get('eta')
        else:
            eta=0
            
        Diso=1/(6*tM);
        Dzz=3*Diso*xi/(2+xi);
        Dxx=(3*Diso-(2/3*eta*(xi-1)/xi+1)*Dzz)/2;
        Dyy=2/3*eta*Dzz*(xi-1)/xi+Dxx;
        Dsq=(Dxx*Dyy+Dyy*Dzz+Dzz*Dxx)/3;
        
    "We the relaxation rates"    
    D1=4*Dxx+Dyy+Dzz;
    D2=Dxx+4*Dyy+Dzz;
    D3=Dxx+Dyy+4*Dzz;
    D4=6*Diso+6*np.sqrt(Diso**2-Dsq);
    D5=6*Diso-6*np.sqrt(Diso**2-Dsq);
    


    dx=(Dxx-Diso)/np.sqrt(Diso**2-Dsq);
    dy=(Dyy-Diso)/np.sqrt(Diso**2-Dsq);
    dz=(Dzz-Diso)/np.sqrt(Diso**2-Dsq);
    
    
    "We rotate the vectors in structure"
    if 'euler' in kwargs and direct=='vXY':
        vec=RotVec(kwargs.get('euler'),struct.vXY)
    elif 'euler' in kwargs:
#        vec=RotVec(kwargs.get('euler'),struct.vCSA) 
        "Use the ABOVE LINE! We need to add support for calculating the CSA direction first...."
        vec=RotVec(kwargs.get('euler'),struct.vXY)
    else:
        print('Euler angles are required')
        return
        
    
    n=vec.shape[0]
    tM=np.zeros([5])
    A=np.zeros([n,5])
    
    for k in range(0,n):
        m=vec[k,:]
        res1=(1/4)*(3*(m[0]**4+m[1]**4+m[2]**4)-1)
        res2=(1/12)*(dx*(3*m[0]**4+6*m[1]**2*m[2]**2-1)\
        +dy*(3*m[1]**4+6*m[2]**2*m[0]**2-1)\
        +dz*(3*m[2]**4+6*m[0]**2*m[1]**2-1))
        
        A[k,0]=3*(m[1]**2)*(m[2]**2);
        A[k,1]=3*(m[0]**2)*(m[2]**2); 
        A[k,2]=3*(m[0]**2)*(m[1]**2); 
        A[k,3]=res1-res2;
        A[k,4]=res1+res2;
        
    tM[0]=1/D1
    tM[1]=1/D2
    tM[2]=1/D3
    tM[3]=1/D4
    tM[4]=1/D5
    
    return tM,A

def Combined(tMdl1,AMdl1,tMdl2,AMdl2):
    nt1=tMdl1.size
    nt2=tMdl2.size
    if np.size(tMdl1)==0:
        tMdl=tMdl2
        AMdl=AMdl2
    elif np.size(tMdl2)==0:
        tMdl=tMdl1
        AMdl=AMdl1
    else:
        tMdl=np.zeros((nt1+1)*(nt2+1)-1)
    
        tMdl[0:nt1]=tMdl1
        tMdl[nt1:nt1+nt2]=tMdl2
        
        
        for k in range(0,nt1):
            for m in range(0,nt2):
                tMdl[nt1+nt2+m+k*nt2]=tMdl1[k]*tMdl2[m]/(tMdl1[k]+tMdl2[m])
                
        AMdl1=np.swapaxes(AMdl1,0,-1)
        AMdl2=np.swapaxes(AMdl2,0,-1)

        
        if AMdl1.shape[1:]!=AMdl2.shape[1:]:
            if AMdl1.ndim>AMdl2.ndim:
                for k in range(1,AMdl1.ndim):
                    AMdl2=np.repeat(np.array([AMdl2.T]),AMdl1.shape[k],axis=k)
            else:
                for k in range(1,AMdl2.ndim):
                    AMdl1=np.repeat(np.array([AMdl1.T]),AMdl2.shape[k],axis=k) 
        
        S21=1-np.sum(AMdl1,axis=0)  
        S22=1-np.sum(AMdl2,axis=0)          

        AMdl=np.zeros(np.concatenate(([(nt1+1)*(nt2+1)-1],AMdl2.shape[1:])))
        AMdl[0:nt1]=np.multiply(np.repeat([S22],AMdl1.shape[0],axis=0),AMdl1)
        AMdl[nt1:nt1+nt2]=np.multiply(np.repeat([S21],AMdl2.shape[0],axis=0),AMdl2)
        
        for k in range(0,nt1):
            for m in range(0,nt2):
                AMdl[nt1+nt2+m+k*nt2]=np.multiply(AMdl1[k],AMdl2[m])
        
        AMdl=np.swapaxes(AMdl,0,-1)

    return tMdl,AMdl
    
def RotVec(euler,vec):
    def Rz(theta):
        return np.array([[np.cos(theta),np.sin(theta),0],[-np.sin(theta),np.cos(theta),0],[0,0,1]])
    def Ry(theta):
        return np.array([[np.cos(theta),0,-np.sin(theta)],[0,1,0],[np.sin(theta),0,np.cos(theta)]])
    
    
    
    return Rz(euler[2]).dot(Ry(euler[1]).dot(Rz(euler[0]).dot(vec.T))).T 
    
    