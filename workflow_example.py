#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 29 16:51:35 2019

@author: albertsmith
"""

"""
The main idea is to have a program for comparison of MD and NMR relaxation data,
and use MD-derived correlation to improve the interpretation of the NMR parameters
"""

import DIFRATE as DR

#%% MD workflow

"File locations"
GPCR_dcd='/Users/albertsmith/Documents/Dynamics/GPCR_MD_data/GPCR_run1/reduced_1ns.dcd'
GPCR_psf='/Users/albertsmith/Documents/Dynamics/GPCR_MD_data/GPCR_run1/step5_assembly.xplor.psf'


"We open the trajectory"
mol=DR.molecule(GPCR_psf,GPCR_dcd)
"We are interested in HN motion- abbreviated here just as N (NMR relaxation is done on backbone 15N)"
mol.select_atoms(Nuc='N')

"This just produces a pdb onto which we can plot results"
mol.MDA2pdb(selection='protein')

"This calculates the correlation functions"
Ctdata=DR.iRED.iRED2data(mol,2,parallel='y',dt=1,align='y')
""" mol contains the trajectory object
2 is the rank of the tensors used for correlation (1 or 2, use 2 for comparison to NMR relaxation)
parallel uses parallel processing('y' is default)
dt is the timestep in nanoseconds (not sure why, but MDanalysis gets the wrong value, so we can override it here)
align removes orientation dependence from the cross-correlation calculation.
"""

"Optimize a set of detectors. Here we use 6 detectors"
Ctdata.detect.r_auto(6)
"Visualize the detector optimization"
Ctdata.detect.plot_r_opt()


"""Fit correlation functions, using the set of detectors just optimized
Note that this is a fit of the correlation function of the eigenmodes, so we 
need to transform it back into the basis of individual bonds
"""
fit=Ctdata.fit()

"""Converts the detector responses of the individual modes into auto- and cross-
correlated detector responses
"""
final=fit.iRED2rho()

"Plot some of the results"
final.plot_R()  #Detector responses
final.draw_rho3D(2) #Detector responses for detector 2 on molecule (requires Chimera)
final.plot_cc(3) #Plot correlation matrix
final.draw_cc3D(48.3)   #Plot row of correlation matrix on molecule (requires Chimera)

"""
There's quite a bit still missing- output to text files, and saving the critical
information. It would be nice if I could save the current state of the classes, 
to pick up work where left off (for large trajectories).
"""

#%% NMR workflow
""" I'll update this later- similar to MD, but currently it's not set up to have
a nice way to load data in from a text file or similar.
"""

#%% NMR and MD comparison
""" Include here methods for constructing detectors for direct NMR and MD comparison
"""