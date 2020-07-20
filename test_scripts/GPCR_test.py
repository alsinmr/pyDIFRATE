#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 16:48:29 2019

@author: albertsmith
"""

import os
curdir=os.getcwd()
os.chdir('/Users/albertsmith/Documents/GitHub/pyDIFRATE')
import DIFRATE as DR
os.chdir('iRED')
import iRED_fast as irf
os.chdir(curdir)

# Data locations
folder='/Users/albertsmith/Documents/Dynamics/GPCR_MD_data/GPCR_run1/'
dcd=folder+'reduced_1ns.dcd'
psf=folder+'step5_assembly.xplor.psf'

#Load trajectory, select bond (N defaults to backbone H–N bond)
mol=DR.molecule()
mol.load_struct(psf,dcd)
mol.select_atoms(Nuc='N')
mol.set_selection()
mol.align('name CA')

#Calculates correlation functions of eigenmodes (correlation calculated at t=0)
data=irf.iRED2data(mol,n=5,dt=1,align_iRED='y')   #For faster calculation, we only sample some frames. Higher n = more frames. dt is timestep in ns

#Set up data analysis. 5 generates 5 detectors
data.detect.r_auto(5)

#ANalyze data (detector analysis of the correlation functions of the eigenmodes)
fit0=data.fit()
#Map results back to individual bonds
"Note that for dynamic regions, the eignemode analysis is not returning accurate detector responses- can be explained by cross-correlation of modes for t not equal 0"
fit=fit0.iRED2rho()

#Plot results
ax=fit.plot_rho()  #Plot detector responses
fit.plot_cc(2)  #Plot cross-correlation, for some detector (currently set to 1)


#Re-analyze detector responses without eigenmode analysis
mol1=mol.copy()
mol1.align(select='name CA')
data1=DR.Ct_fast.Ct2data(mol,n=5,dt=1) #Load correlation functions
data1.detect=data.detect    #Copy detectors from the above eigenmode analysis (note, n must be the same for both analyses)
fit1=data1.fit()        #Fit the data
fit1.plot_rho() #Plot the results

for k,a in enumerate(ax):
    a.plot(fit1.label,fit1.R[:,k],color='grey')
    a.autoscale(enable=True,tight=True)

"""
The results from these two analyses should match (or at least be pretty similar).
Differences probably results from incompleteness of eigenmode analysis, so that results in
fit1 should be accurate, but in fit, have some error
"""


#%% Check performance when truncating the trajectory
data2=DR.Ct_fast.Ct2data(mol,n=1,nt=329,dt=1)
data2.detect.r_auto(4)

data1=DR.Ct_fast.Ct2data(mol,n=1,dt=1) #Load correlation functions

data1.detect.r_target(data2.detect.rhoz(),n=5,NT='M')

fit1=data1.fit()
fit2=data2.fit()

ax=fit1.plot_rho()
for k,a in enumerate(ax):
    a.plot(fit2.label,fit2.R[:,k],color='grey')
    
#%% Check performance with pre-analysis (unoptimized detectors) and direct analysis
data1.detect.r_no_opt(10)

fit00=data1.fit()

fit00.detect.r_auto(5)
fit01=fit00.fit()

data1.detect.r_auto(5)
fit02=data1.fit()

ax=fit01.plot_rho()
for k,a in enumerate(ax):
    a.plot(fit02.label,fit02.R[:,k],color='grey')