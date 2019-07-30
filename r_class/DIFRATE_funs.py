#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 09:32:49 2019

@author: albertsmith

Collection of useful functions for DIFRATE
"""

import pandas as pd
import re
import numpy as np
import os 

def NucInfo(Nuc=None,info='gyro'):
    """ Returns the gyromagnetic ratio for a given nucleus. Usually, should be 
    called with the nucleus and mass number, although will default first to 
    spin 1/2 nuclei if mass not specified, and second to the most abundant 
    nucleus. A second argument, info, can be specified to request the 
    gyromagnetic ratio ('gyro'), the spin ('spin'), the abundance ('abund'), or 
    if the function has been called without the mass number, one can return the 
    default mass number ('mass'). If called without any arguments, a pandas 
    object is returned containing all nuclear info ('nuc','mass','spin','gyro',
    'abund')
    """
    
    Nucs=[]
    MassNum=[]
    spin=[]
    g=[]
    Abund=[]

    dir_path = os.path.dirname(os.path.realpath(__file__))
    
    with open(dir_path+"/GyroRatio") as f:
        data=f.readlines()
        for line in data:
            line=line.strip().split()
            MassNum.append(int(line[1]))
            Nucs.append(line[3])
            spin.append(float(line[5]))
            g.append(float(line[6]))
            Abund.append(float(line[7]))
    
    NucData=pd.DataFrame({'nuc':Nucs,'mass':MassNum,'spin':spin,'gyro':g,'abund':Abund})
    
    
    if Nuc==None:
        return NucData
    else:
        
        if Nuc=='D':
            Nuc='2H'
        
        mass=re.findall(r'\d+',Nuc)
        if not mass==[]:
            mass=int(mass[0])
            
        
        Nuc=re.findall(r'[A-Z]',Nuc.upper())
        
        if np.size(Nuc)>1:
            Nuc=Nuc[0].upper()+Nuc[1].lower()
        else:
            Nuc=Nuc[0]
            
            
            
        NucData=NucData[NucData['nuc']==Nuc]
       
        if not mass==[]:    #Use the given mass number
            NucData=NucData[NucData['mass']==mass]
        elif any(NucData['spin']==0.5): #Ambiguous input, take spin 1/2 nucleus if exists
            NucData=NucData[NucData['spin']==0.5] #Ambiguous input, take most abundant nucleus
        elif any(NucData['spin']>0):
            NucData=NucData[NucData['spin']>0]
        
        NucData=NucData[NucData['abund']==max(NucData['abund'])]
            
        
        h=6.6260693e-34
        muen=5.05078369931e-27
        
        NucData['gyro']=float(NucData['gyro'])*muen/h
#        spin=float(NucData['spin'])
#        abund=float(NucData['abund'])
#        mass=float(NucData['spin'])
        if info=='allinfo':
            return NucData
        else:
            return float(NucData[info])

def dipole_coupling(r,Nuc1,Nuc2):
    """ Returns the dipole coupling between two nuclei ('Nuc1','Nuc2') 
    separated by a distance 'r' (in nm). Result in Hz (gives full anisotropy,
    not b12, that is 2x larger than b12)
    """
    
    gamma1=NucInfo(Nuc1)
    gamma2=NucInfo(Nuc2)
    
    h=6.6260693e-34 #Plancks constant in J s
    mue0 = 12.56637e-7  #Permeability of vacuum [T^2m^3/J]
    
    return h*2*mue0/(4*np.pi*(r/1e9)**3)*gamma1*gamma2
    
    
def J(tc,v):
    "Returns the spectral density"
    return 2/5*tc/(1+(2*np.pi*v*tc)**2)

def rate(tc,exper):
    """Returns the sensitivity of an experiment, specified by exper (one 
    column of a pandas array), for a given set of correlation times, tc.
    """
    if exper.loc['Type']=='R1':
        return R1(tc,exper)
    elif exper.loc['Type']=='R1p':
        return R1p(tc,exper)
    elif exper.loc['Type']=='NOE':
        return NOE(tc,exper)
    elif exper.loc['Type']=='R2':
        return R2(tc,exper)
    elif exper.loc['Type']=='cct':
        return cct(tc,exper)
    elif exper.loc['Type']=='ccl':
        return ccl(tc,exper)
    elif exper.loc['Type']=='Ct':
        return Ct(tc,exper)
    elif exper.loc['Type']=='R1Q':
        return R1Q(tc,exper)
    else:
        print('Experiment type {0} was not recognized'.format(exper.loc['Type']))
        return
    
def R1(tc,exper):
    v0=exper['v0']*1e6
    dXY=exper['dXY']
    Nuc=exper['Nuc']
    Nuc1=exper['Nuc1']
    QC=exper['QC']
    eta=exper['eta']
    vX=NucInfo(Nuc)/NucInfo('1H')*v0
    CSA=exper['CSA']*vX/1e6
    R=np.zeros(tc.shape)

    if Nuc1 is not None and dXY is not None:
        "Dipolar relaxation"
        if np.size(dXY)==1:
            vY=NucInfo(Nuc1)/NucInfo('1H')*v0
            S=NucInfo(Nuc1,'spin')
            sc=S*(S+1)*4/3 # Scaling factor depending on the spin, =1 for spin 1/2
            R+=sc*(np.pi*dXY/2)**2*(J(tc,vX-vY)+3*J(tc,vX)+6*J(tc,vY+vX))
        else:
            for k in range(0,np.size(dXY)):
                S=NucInfo(Nuc1[k],'spin')
                sc=S*(S+1)*4/3 # Scaling factor depending on the spin, =1 for spin 1/2
                vY=NucInfo(Nuc1[k])
                R+=sc*(np.pi*dXY[k]/2)**2*(J(tc,vX-vY)+3*J(tc,vX)+6*J(tc,vY+vX))
                
    "CSA relaxation"
    R+=3/4*(2*np.pi*CSA)**2*J(tc,vX)

    if QC!=0:
        "Quadrupolar relaxation"
        """Note that these formulas give the initial rate of relaxation, that 
        is, the average rate of relaxation for all orientations, and furthermore
        does not include deviations due to multi-exponential relaxation
        """
        S=NucInfo(Nuc,'spin')
        deltaQ=1/(2*S*(2*S-1))*QC*2*np.pi
        C=(deltaQ/2)**2*(1+eta**2/3) #Constant that scales the relaxation
        if S==0.5:
            print('No quadrupole coupling for S=1/2')
        elif S==1:
            R+=C*(3*J(tc,vX)+12*J(tc,2*vX))
        elif S==1.5:
            R+=C*(36/5*J(tc,vX)+144/5*J(tc,2*vX))
        elif S==2.5:
            R+=C*(96/5*J(tc,vX)+384/5*J(tc,2*vX))
        else:
            print('Spin not implemented')
            
    return R

def R1Q(tc,exper):
    """This function calculates the relaxation rate constant for relaxation of
    quadrupolar order
    """
    v0=exper['v0']*1e6
    Nuc=exper['Nuc']
    QC=exper['QC']
    eta=exper['eta']
    vX=NucInfo(Nuc)/NucInfo('1H')*v0
    
    S=NucInfo(Nuc,'spin')
    deltaQ=1/(2*S*(2*S-1))*QC*2*np.pi
    C=(deltaQ/2)**2*(1+eta**2/3)    #Constant scaling the relaxation
    if S==0.5:
        print('Quadrupolar relaxation experiment not possible for spin 1')
    elif S==1:
        R=C*9*J(tc,vX)
    elif S==1.5:
        R=C*(36*J(tc,vX)+36*J(tc,2*vX))
    elif S==2.5:
        R=C*(792/7*J(tc,vX)+972/7*J(tc,2*vX))
    else:
        print('Spin not implemented')
        
    return R

def R1p(tc,exper):
    v0=exper['v0']*1e6
    dXY=exper['dXY']
    Nuc=exper['Nuc']
    Nuc1=exper['Nuc1']
    QC=exper['QC']
    eta=exper['eta']
    vr=exper['vr']*1e3
    v1=exper['v1']*1e3
    off=exper['offset']*1e3
    vX=NucInfo(Nuc)/NucInfo('1H')*v0
    CSA=exper['CSA']*vX/1e6
    R=np.zeros(tc.shape)
    
    "Treat off-resonance spin-lock"
    ve=np.sqrt(v1**2+off**2)
    if ve==0:
        theta=np.pi/2
    else:
        theta=np.arccos(off/ve)
    
    R10=R1(tc,exper)    #We do this first, because it includes all R1 contributions
    "Start here with the dipole contributions"
    if Nuc1 is not None:
        if np.size(dXY)==1:
            vY=NucInfo(Nuc1)/NucInfo('1H')*v0
            S=NucInfo(Nuc1,'spin')
            sc=S*(S+1)*4/3 #Scaling depending on spin of second nucleus
            R1del=sc*(np.pi*dXY/2)**2*(3*J(tc,vY)+
                      1/3*J(tc,2*vr-ve)+2/3*J(tc,vr-ve)+2/3*J(tc,vr+ve)+1/3*J(tc,2*vr+ve))
        else:
            
            R1del=np.zeros(tc.shape)
            for k in range(0,np.size(dXY)):
                vY=NucInfo(Nuc1[k])/NucInfo('1H')*v0
                S=NucInfo(Nuc1[k],'spin')
                sc=S*(S+1)*4/3 #Scaling depending on spin of second nucleus
                R1del+=sc*(np.pi*dXY[k]/2)**2*(3*J(tc,vY)+
                          1/3*J(tc,2*vr-ve)+2/3*J(tc,vr-ve)+2/3*J(tc,vr+ve)+1/3*J(tc,2*vr+ve))
    else:
        R1del=np.zeros(tc.shape)
    "CSA contributions"
    R1del+=1/6*(2*np.pi*CSA)**2*(1/2*J(tc,2*vr-ve)+J(tc,vr-ve)+J(tc,vr+ve)+1/2*J(tc,vr+ve))
    
    "Here should follow the quadrupole treatment!!!"    
    
    "Add together R1 and R1p contributions, depending on the offset"
    R+=R10+np.sin(theta)**2*(R1del-R10/2) #Add together the transverse and longitudinal contributions    
    return R

def R2(tc,exper):
    exper['off']=0
    exper['v1']=0
    
    return R1p(tc,exper)   

def NOE(tc,exper):
    v0=exper['v0']*1e6
    dXY=exper['dXY']
    Nuc=exper['Nuc']
    Nuc1=exper['Nuc1']
    QC=exper['QC']
    eta=exper['eta']
    vr=exper['vr']*1e3
    v1=exper['v1']*1e3
    off=exper['offset']*1e3
    vX=NucInfo(Nuc)/NucInfo('1H')*v0
    CSA=exper['CSA']*vX/1e6
    R=np.zeros(tc.shape)
    
    if Nuc1!=None:
        vY=NucInfo(Nuc1)/NucInfo('1H')*v0
        S=NucInfo(Nuc1,'spin')
        sc=S*(S+1)*4/3 # Scaling factor depending on the spin, =1 for spin 1/2
        R+=sc*(np.pi*dXY/2)**2*(-J(tc,vX-vY)+6*J(tc,vY+vX))
        
    return R
        