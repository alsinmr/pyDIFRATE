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
os.chdir(curdir)

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


start=1
stop=50
select='resname POPC and resid {0}-{1}'.format(start,stop)
mol.vec_special(Type='moment_of_inertia',select=select)


