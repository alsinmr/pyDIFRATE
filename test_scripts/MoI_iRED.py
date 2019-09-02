#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 17:18:18 2019

@author: albertsmith
"""

#%% Required modules
import numpy as np
import MDAnalysis as mda
import os
curdir=os.getcwd()
os.chdir('/Users/albertsmith/Documents/GitHub/pyDIFRATE')
import DIFRATE as DR
os.chdir('iRED')
import iRED_fast as irf
os.chdir(curdir)
import matplotlib.pyplot as plt

#%% File locations, etc.
dir_str='/Users/albertsmith/Documents/Dynamics/Lipids/POPC/MD_ana/overall/'
file_str='oam{0}_nd{1}_n{2}_range{3}_{4}'

dcd='/Volumes/My Book/MD_{0}/step6.6_equilibration.gro'
psf0='/Volumes/My Book/MD_{0}/run1.part{1:04d}.xtc'

nr=256
nd=11

psf=list()
for k in range(2,140):
    psf.append(psf0.format(nr,k))

mol=DR.molecule(dcd.format(nr),psf[:])  #Create the molecule object

uni=mol.mda_object
sel0=uni.select_atoms('resname POPC')
res=10
nl=50

loc=sel0.select_atoms('resid {0}'.format(res)).center_of_mass()
dist=list()
for res in sel0.residues:
    dist.append(np.sqrt(np.sum((loc-res.atoms.center_of_mass())**2)))

x=np.argsort(dist)    
sel=list()
for k in x[:nl]:
    sel.append(sel0.residues[k].atoms)



res=10
select='resname POPC and (resid {0} or around 20 resid {0})'.format(res)
mol.vec_special(Type='moment_of_inertia',sel=sel)
mol.label=range(mol.vec_fun().shape[1])

data=irf.iRED2data(mol,n=1)
data.detect.r_auto(nd,NT='I')

fit0=data.fit()
fit=fit0.iRED2rho()


z0=fit.sens.info.loc['z0']

ax=plt.figure().add_subplot(1,1,1)
ax.plot(z0,fit.R.mean(axis=0))


data1=DR.Ct_fast.Ct2data(mol,n=1)
data1.detect=data.detect
fit1=data1.fit()

ax.plot(z0,fit1.R.mean(axis=0))