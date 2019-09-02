#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 16:48:29 2019

@author: albertsmith
"""


import DIFRATE as DR
import iRED_fast as irf

folder='/Users/albertsmith/Documents/Dynamics/GPCR_MD_data/GPCR_run1/'
dcd=folder+'reduced_1ns.dcd'
psf=folder+'step5_assembly.xplor.psf'

mol=DR.molecule()
mol.load_struct(psf,dcd)
mol.select_atoms(Nuc='15N')
mol.set_selection()

data=irf.iRED2data(mol,n=15,dt=1)

data.detect.r_auto(5)
fit0=data.fit()
fit=fit0.iRED2rho()

fit.plot_rho()


data1=DR.Ct_fast.Ct2data(mol,n=15,dt=1)
data1.detect=data.detect
fit1=data1.fit()
fit1.R_u[:]=0
fit1.R_l[:]=0
fit1.plot_rho()