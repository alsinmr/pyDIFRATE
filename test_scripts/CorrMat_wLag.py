#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 10:12:01 2019

@author: albertsmith
"""

import os
curdir=os.getcwd()
os.chdir('/Users/albertsmith/Documents/GitHub/pyDIFRATE')
import DIFRATE as DR
os.chdir('iRED')
import iRED_fast as irf
os.chdir('../plotting')
import plotting_funs as pf 
os.chdir(curdir)

import numpy as np


folder='/Users/albertsmith/Documents/Dynamics/GPCR_MD_data/GPCR_run1/'
dcd=folder+'reduced_1ns.dcd'
psf=folder+'step5_assembly.xplor.psf'

#Load trajectory, select bond (N defaults to backbone H–N bond)
mol=DR.molecule()
mol.load_struct(psf,dcd)
mol.select_atoms(Nuc='N')
mol.set_selection()


#%% Calculate a set of vectors and the M matrix
index=irf.trunc_t_axis(mol.mda_object.trajectory.n_frames,n=-1,dt=1)    #Truncated set of frames
vec=irf.get_trunc_vec(mol,index,dt=1)   #Load set of vectors
M=irf.Mmat(vec) #Correlation Matrix
Ylm=irf.Ylm(vec)    #Spherical Tensors
aqt=irf.Aqt(Ylm,M)  #Project spherical tensors into eigenbasis of M

#%% Calculate correlation matrix with a lag time
Mlag=irf.Mlagged(vec,600)

#%% Diagonalize Matrices (with normalization)
dg=np.array([M['lambda']])
M_dn=np.divide(np.dot(np.dot(M['m'].T,M['M']),M['m']),np.sqrt(np.dot(dg.T,dg)))

Mlad_d=np.dot(np.dot(M['m'].T,Mlag),M['m'])
Mlag_dn=np.divide(np.dot(np.dot(M['m'].T,Mlag),M['m']),np.sqrt(np.dot(dg.T,dg)))

#%% Plot the matrices
pf.plot_cc(M['M'],norm='n')
pf.plot_cc(Mlag,norm='n')
pf.plot_cc(M_dn[:-5,:-5],norm='n')
pf.plot_cc(Mlag_dn[:-5,:-5],norm='n')
pf.plot_cc(Mlag_d[:-5,:-5],norm='n')
pf.plot_cc(Mlag_dn,norm='n')

#%% Plot some cross correlation functions
Cqt=irf.Cqt(aqt)    #Correlation functions of modes
Ctinf=irf.CtInf(aqt)    #Estimate of "order parameter" for modes
DelCt=irf.DelCt(Cqt,Ctinf)  #Normalized correlation function

t=DelCt['t']

fig=plt.figure()
ax0=fig.add_subplot(2,1,1)
ax0.plot(t,DelCt['DelCt'][:,360])

ax1=fig.add_subplot(2,1,2)

C_360_365=irf.Cij_t(aqt,360,365)
Ci_360_365=irf.Cij_Inf(aqt,360,365)

ax1.plot(t,C_360_365['Ct'])