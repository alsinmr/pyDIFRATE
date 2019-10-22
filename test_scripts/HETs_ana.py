#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 11:19:48 2019

@author: albertsmith
"""

import os
curdir=os.getcwd()
os.chdir('/Users/albertsmith/Documents/GitHub/pyDIFRATE')
import DIFRATE as DR
os.chdir('iRED')
import iRED_fast as irf
os.chdir(curdir)

folder='/Users/albertsmith/Documents/Dynamics/HETs_reana/'


#%% Analyze experimental data
#Experimental data and experiment info
dataNMR=DR.io.load_NMR(folder+'HETs_data.txt')   
mol=dataNMR.sens.molecule
mol.load_struct(topo)
mol.select_atoms(Nuc='15N',select='segid B',resi=dataNMR.label-218+72)
mol.set_selection()

dataR1=dataNMR.copy()
dataR1.del_exp(range(3,8))
dataR1p=dataNMR.copy()
dataR1p.del_exp(range(3))

dataR1.detect.r_auto(2,NT='I')
dataR1p.detect.r_auto(2,NT='I')

fitR1=dataR1.fit()
fitR1p=dataR1p.fit()

fitR1.plot_rho()

#%% Load trajectory data
traj='/Volumes/My Book/HETs_AMBER/298K_protein/HETs_whole.xtc'
topo='/Volumes/My Book/HETs_AMBER/298K_protein/topo.pdb'

mol=DR.molecule()
mol.load_struct(topo,traj)
mol.select_atoms(Nuc='N')
mol.set_selection()
select='name CA and segid B and (resid 117-130 or resid 81-94)'
mol.align(select=select)

data=irf.iRED2data(mol,n=5,align_iRED='n')
data.detect.r_auto(6,NT='I')

fit0=data.fit()
fit=fit0.iRED2rho()

data_a=irf.iRED2data(mol,n=5,align_iRED='y')
data_a.detect=data.detect.copy()

fit0a=data_a.fit()
fit_a=fit0a.iRED2rho()

#%% Analyze without iRED
mol=DR.molecule()
mol.load_struct(topo,traj)
select='name CA and segid B and (resid 117-130 or resid 81-94)'
mol.align(select=select)
mol.select_atoms(Nuc='N')
mol.set_selection()


data1=DR.Ct_fast.Ct2data(mol,n=5)
data1.detect=data.detect.copy()

data1.detect.r_auto(5)

fit1=data1.fit()


#%% Check performance when truncating the trajectory
data2=DR.Ct_fast.Ct2data(mol,n=5,nt=int(77141/5))
data2.detect.r_target(data1.detect.rhoz()[0:4],NT='M')

data1.detect.r_target(data2.detect.rhoz(),n=5,NT='M')

fit1=data1.fit()
fit2=data2.fit()

ax=fit1.plot_rho()
for k,a in enumerate(ax):
    a.plot(fit2.label,fit2.R[:,k],color='grey')