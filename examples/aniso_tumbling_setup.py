#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 16:27:22 2021

@author: albertsmith
"""

import pyDIFRATE as DR
import matplotlib.pyplot as plt

nmr=DR.sens.NMR(Type='R1',v0=[100,200,400,500,600,800],Nuc='15N')
nmr.new_exp(Type=['R2','NOE'],v0=800,Nuc='15N')
nmr.new_exp(Type=['R2','NOE'],v0=600,Nuc='15N')

pdb='1d3z.pdb'
nmr.molecule.load_struct(pdb)
nmr.molecule.select_atoms(Nuc='15N')
nmr.new_mdl(Model='AnisoDif',tM=4.84e-9,xi=.2,eta=0,euler=[0,0,0])  

fig=plt.figure('Orientation Dependence for Anisotropic Tumbling')
fig.clear()
ax=[fig.add_subplot(2,5,k+1) for k in range(10)]
for k,a in enumerate(ax):
    for m in range(len(nmr.molecule.sel1)):
        a.plot(nmr.z(),nmr.Reff(exp_num=k,mdl_num=0,bond=m)[0])
        
r=DR.sens.detect(nmr,mdl_num=0)
r.r_auto(5)         #This optimizes the average detectors
r.r_target(bond=-1) #This optimizes the individual bonds

fig=plt.figure('Orientation Dependence for Detector Sensitivities')
fig.clear()
ax=[fig.add_subplot(2,3,k+1) for k in range(5)]
for k,a in enumerate(ax):
    for m in range(len(nmr.molecule.sel1)):
        a.plot(r.z(),r.rhoz(bond=m)[k])

