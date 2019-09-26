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

#Calculates correlation functions of eigenmodes (correlation calculated at t=0)
data=irf.iRED2data(mol,n=15,dt=1)   #For faster calculation, we only sample some frames. Higher n = more frames. dt is timestep in ns

#Set up data analysis. 5 generates 5 detectors
data.detect.r_auto(5)

#ANalyze data (detector analysis of the correlation functions of the eigenmodes)
fit0=data.fit()
#Map results back to individual bonds
"Note that for dynamic regions, the eignemode analysis is not returning accurate detector responses- can be explained by cross-correlation of modes for t not equal 0"
fit=fit0.iRED2rho()

#Plot results
fit.plot_rho()  #Plot detector responses
fit.plot_cc(1)  #Plot cross-correlation, for some detector (currently set to 1)


#Re-analyze detector responses without eigenmode analysis
data1=DR.Ct_fast.Ct2data(mol,n=15,dt=1) #Load correlation functions
data1.detect=data.detect    #Copy detectors from the above eigenmode analysis (note, n must be the same for both analyses)
fit1=data1.fit()        #Fit the data
fit1.R_u[:]=0   #Set error bars to zero (I'm not sure how to estimate the error yet, but I guess right now it's being over estimated)
fit1.R_l[:]=0
fit1.plot_rho() #Plot the results

"""
The results from these two analyses should match (or at least be pretty similar).
Differences probably results from incompleteness of eigenmode analysis, so that results in
fit1 should be accurate, but in fit, have some error
"""
