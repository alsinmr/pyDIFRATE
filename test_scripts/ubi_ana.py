#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 13:24:43 2019

@author: albertsmith
"""

import DIFRATE as DR

data=DR.io.load_NMR('/Users/albertsmith/Documents/GitHub/pyDIFRATE/test_scripts/ubi_data.txt')     #Experimental data and experiment info
data.sens.molecule.load_struct('/Users/albertsmith/Documents/GitHub/pyDIFRATE/test_scripts/1d3z.pdb')  #Structure of ubiquitin
data.sens.molecule.select_atoms(Nuc='15N',resi=data.label)  #Set atoms
data.sens.molecule.set_selection()
data.sens.new_mdl(Model='IsoDif',tM=4.84e-9)    #Isotropic tumbling model
data.sens.new_mdl(Model='AnisoDif',tM=4.84e-9,xi=1.18,euler=np.array([120.0,150.0,0.0])*np.pi/180.0)    #Anisotropic tumbling model


data.new_detect(mdl_num=0)  #Use model of isotropic tumbling

data.detect.r_auto(4,R2_ex_corr='yes')  #4 detectors, 5th detector controls for exchange contribution to R2
fit0=data.fit() #Fit data to detectors

data.new_detect(mdl_num=1)  #Use model of anisotropic tumbling
#data.detect.r_auto(4,-1,R2_ex_corr='yes')  #This line is supposed to perform the function of following two lines, but seems to be broken
data.detect.r_auto(4,R2_ex_corr='yes')  #Optimize detectors for average sensitivites
data.detect.r_target(bond=-1)   #Set individual sensitivities to have approximately same sensitivy as above optimized detectors
fit1=data.fit() #Fit data to detectors

fit0.plot_rho()
fit1.plot_rho()