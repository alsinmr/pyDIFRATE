#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 13:42:37 2021

@author: albertsmith
"""

import numpy as np
import os
import sys
sys.path.append('/Users/albertsmith/Documents/GitHub/')
import pyDIFRATE as DR
import matplotlib.pyplot as plt
from pyDIFRATE.Struct import select_tools as selt

directory='/Users/albertsmith/Documents/Dynamics/HETs_Methyl_Loquet/procMD/full_noopt'
#%% Load MD simulation
MDdir='/Volumes/My Book/HETs/'
files0=['_3pw','_4pw','_MET_3pw','_MET_4pw','5_MET_4pw_cb10_anothertry']
files=[os.path.join(MDdir,'MDSimulation','HETs'+f+'.xtc') for f in files0]
sim=lambda k:DR.molecule(os.path.join(MDdir,'HETs_5chain_B.pdb' if k==4 else 'HETs_3chain.pdb'),files[k])


mol=sim(4)
nr=int(len(mol.mda_object.residues)/71)
mol.mda_object.residues.resids=np.atleast_2d(np.arange(219,290)).repeat(nr,axis=0).reshape(nr*71)
mol.select_atoms(Nuc='ivla',segids='B')
mol.clear_frames()
mol.tensor_frame(sel1=1,sel2=2)
mol.new_frame(Type='methylCC',Nuc='ivla',segids='B',frame_index=np.arange(40).repeat(3))

frames=DR.frames.frames2data(mol,n=-1,tf=2000,mode='sym')

r=frames[0].detect
r.r_auto(7,NegAllow=0)
for f in frames:f.detect=r
fr_fit=[f.fit() for f in frames]

ax=fr_fit[0].plot_rho()
for k,a in enumerate(ax):
    a.plot(fr_fit[1].R[:,k],color='grey')




#%% Try with 2 frames
mol.clear_frames()
mol.tensor_frame(sel1=1,sel2=2)
mol.new_frame(Type='hops_3site',Nuc='ivla',segids='B')
mol.new_frame(Type='methylCC',Nuc='ivla',segids='B')
#mol.new_frame(Type='side_chain_chi',Nuc='ivla',segids='B',n_bonds=0)
mol.new_frame(Type='side_chain_chi',Nuc='ivla',segids='B',n_bonds=1)
mol.new_frame(Type='side_chain_chi',Nuc='ivla',segids='B',n_bonds=2)


frames=DR.frames.frames2data(mol,n=-1,tf=2000,mode='auto')

r=frames[0].detect
r.r_auto(7,NegAllow=0)
for f in frames:f.detect=r
fr_fit=[f.fit() for f in frames]

ax=fr_fit[0].plot_rho()
for k,a in enumerate(ax):
    a.plot(fr_fit[1].R[:,k],color='grey')
  
#%% Try with isoleucine, chi1/chi2 rotation    
mol.select_atoms(Nuc='ivlt',segids='B',filter_str='resname ILE')
mol.clear_frames()
mol.tensor_frame(sel1=1,sel2=2)
resids=mol.sel1.residues.resids
frame_index=np.arange(resids.size).repeat(3)

cg1=selt.sel_simple(mol,sel='name CG1',resids=resids,segids='B')
cb=selt.sel_simple(mol,sel='name CB',resids=resids,segids='B')
ca=selt.sel_simple(mol,sel='name CA',resids=resids,segids='B')
co=selt.sel_simple(mol,sel='name C',resids=resids,segids='B')

mol.new_frame(Type='hops_3site',Nuc='ivlat',resids=resids,segids='B')
mol.new_frame(Type='methylCC',Nuc='ivlat',resids=resids,segids='B',frame_index=frame_index)
mol.new_frame(Type='bond',sel1=cg1,sel2=cb,sel3=ca,frame_index=frame_index)
mol.new_frame(Type='bond',sel1=cb,sel2=ca,sel3=co,frame_index=frame_index)

frames=DR.frames.frames2data(mol,n=-1,tf=2000,mode='auto')
t=frames[0].sens.info.loc['t'].to_numpy()
#stdev=np.diff(np.log10(t))
#stdev[0]=1e-6
#stdev=np.append(stdev,stdev[-1])
#
#
stdev=(t/t.max())**.5
stdev[0]=stdev[1]/100

frames[0].sens.info.loc['stdev']=stdev
frames[0].new_detect()
for f in frames:f.sens=frames[0].sens
for f in frames:f.R_std[:]=f.sens.info.loc['stdev'].to_numpy()

r=frames[0].detect
r.r_auto(7,NegAllow=0)
for f in frames:f.detect=r
fr_fit=[f.fit() for f in frames]

ax=fr_fit[0].plot_rho()
for k,a in enumerate(ax):
    a.plot(fr_fit[1].R[:,k],color='grey')
    a.set_ylim([0,1 if k==0 else np.max([.2,a.get_ylim()[1]])])
    
for f in fr_fit[2:]:
    ax=f.plot_rho()
    for k,a in enumerate(ax):
        *_,Rc=DR.fitting.fit2tc(f,df=3)
#        a.plot(range(f.R.shape[0]),Rc[:,k],color='black',linestyle=':')
        a.set_ylim([0,1 if k==0 else np.max([.2,a.get_ylim()[1]])])


#%% Plot selected correlation functions
i=11
fig=plt.figure()
ax=[fig.add_subplot(1,2,k) for k in range(1,3)]

ax[0].semilogx(t[1:],frames[0].R[i,1:])
ax[0].semilogx(t[1:],frames[1].R[i,1:],color='black',linestyle=':')
for f in frames[2:]:
    ax[1].semilogx(t[1:],f.R[i,1:])
ax[1].legend(['Libration','Methyl rot.',r'$\chi_1$',r'$\chi_2$','overall'])
for a in ax:a.set_ylim([0,1])


#%% Some 3D plotting
chimera_cmds=[]
directory='/Users/albertsmith/Documents/Jobs/Berlin2021/Figures/'
fileout=os.path.join(directory,'HETs_{0}_rho{1}')
scene=os.path.join(directory,'HETs.cxs')
for m,f in enumerate(fr_fit[2:4]):
    for k in range(5):
        f.draw_rho3D(k,scene=scene)