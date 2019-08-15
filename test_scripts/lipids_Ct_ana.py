#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 14:38:11 2019

@author: albertsmith
"""

import numpy as np
import MDAnalysis as mda


dcd='/Volumes/My Book/MD/step6.0_minimization.gro'

psf=list()
for k in range(2,80):
    psf.append('/Volumes/My Book/MD/run1.part{0:04d}.xtc'.format(k))    


uni=mda.Universe(dcd,psf[0:80])

#uni=mda.Universe(dcd,'/Volumes/My Book/MD/reduced_1ns_whole.xtc')

mol=DR.molecule()

mol.mda_object=uni

sel01=['C13','C13','C13','C14','C14','C14','C15','C15','C15',
       'C12','C12',
       'C11','C11',
       'C1','C1',
       'C2',
       'C3','C3',
       'C31',
       'C32','C32',
       'C33','C33',
       'C34','C34',
       'C35','C35',
       'C36','C36',
       'C37','C37',
       'C38','C38',
       'C39','C39',
       'C310','C310',
       'C311','C311',
       'C312','C312',
       'C313','C313',
       'C314','C314',
       'C315','C315',
       'C316','C316','C316',
       'C21',
       'C22','C22',
       'C23','C23',
       'C24','C24',
       'C25','C25',
       'C26','C26',
       'C27','C27',
       'C28','C28',
       'C29',
       'C210',
       'C211','C211',
       'C212','C212',
       'C213','C213',
       'C214','C214',
       'C215','C215',
       'C216','C216',
       'C217','C217',
       'C218','C218','C218'
       ]


sel02=['H13A','H13B','H13C','H14A','H14B','H14C','H15A','H15B','H15C',
       'H12A','H12B',
       'H11A','H11B',
       'HA','HB',
       'HS',
       'HX','HY',
       'O32',
       'H2X','H2Y',
       'H3X','H3Y',
       'H4X','H4Y',
       'H5X','H5Y',
       'H6X','H6Y',
       'H7X','H7Y',
       'H8X','H8Y',
       'H9X','H9Y',
       'H10X','H10Y',
       'H11X','H11Y',
       'H12X','H12Y',
       'H13X','H13Y',
       'H14X','H14Y',
       'H15X','H15Y',
       'H16X','H16Y','H16Z',
       'O22',
       'H2S','H2R',
       'H3S','H3R',
       'H4S','H4R',
       'H5S','H5R',
       'H6S','H6R',
       'H7S','H7R',
       'H8S','H8R',
       'H91',
       'H101',
       'H11S','H11R',
       'H12S','H12R',
       'H13S','H13R',
       'H14S','H14R',
       'H15S','H15R',
       'H16S','H16R',
       'H17S','H17R',
       'H18S','H18R','H18T'
       ]

label0=['gamma1A','gamma1B','gamma1C','gamma2A','gamma2B','gamma2C','gamma3A','gamma3B','gamma3C',
       'betaA','betaB',
       'alphaA','alphaB',
       'g3A','g3B',
       'g2',
       'g1X','g1Y',
       'C1_1',
       'C1_2X','C1_2Y',
       'C1_3X','C1_3Y',
       'C1_4X','C1_4Y',
       'C1_5X','C1_5Y',
       'C1_6X','C1_6Y',
       'C1_7X','C1_7Y',
       'C1_8X','C1_8Y',
       'C1_9X','C1_9Y',
       'C1_10X','C1_10Y',
       'C1_11X','C1_11Y',
       'C1_12X','C1_12Y',
       'C1_13X','C1_13Y',
       'C1_14X','C1_14Y',
       'C1_15X','C1_15Y',
       'C1_16X','C1_16Y','C1_16Z',
       'C2_1',
       'C2_2X','C2_2Y',
       'C2_3X','C2_3Y',
       'C2_4X','C2_4Y',
       'C2_5X','C2_5Y',
       'C2_6X','C2_6Y',
       'C2_7X','C2_7Y',
       'C2_8X','C2_8Y',
       'C2_9',
       'C2_10',
       'C2_11X','C2_11Y',
       'C2_12X','C2_12Y',
       'C2_13X','C2_13Y',
       'C2_14X','C2_14Y',
       'C2_15X','C2_15Y',
       'C2_16X','C2_16Y',
       'C2_17X','C2_17Y',
       'C2_18X','C2_18Y','C2_18Z',
       ]



res_in=[114,15,56,45,117]
res_in=range(1,21)

i1=np.zeros(np.unique(sel01).size*np.size(res_in),'int64')

a,b=np.unique(sel01,return_inverse=True)

for m in range(0,np.size(res_in)):
    sel0=uni.select_atoms('resid {0}'.format(res_in[m]))
    for k,s in enumerate(a):
        i1[k+m*a.size]=sel0.select_atoms('name {0}'.format(s)).indices[0]



resi_str=list()
for k in res_in:
    resi_str.append('resid {0}'.format(k))
    
sel0=uni.select_atoms(*resi_str)

mol.sel1=uni.select_atoms('name *')[i1]

mol.sel1in=list()
for k in range(0,np.size(res_in)):
    mol.sel1in.append(b+a.size*k)
mol.sel1in=np.concatenate(mol.sel1in)


i1=np.zeros(np.unique(sel02).size*np.size(res_in),'int64')

a,b=np.unique(sel02,return_inverse=True)

for m in range(0,np.size(res_in)):
    sel0=uni.select_atoms('resid {0}'.format(res_in[m]))
    for k,s in enumerate(a):
        i1[k+m*a.size]=sel0.select_atoms('name {0}'.format(s)).indices[0]


sel0=uni.select_atoms(*resi_str)
mol.sel2=uni.select_atoms('name *')[i1]

mol.sel2in=list()
for k in range(0,np.size(res_in)):
    mol.sel2in.append(b+a.size*k)
mol.sel2in=np.concatenate(mol.sel2in)



"Load in a pdb here"
mol.MDA2pdb(tstep=0,selection='resid 114 or resid 15 or resid 56 or resid 45 or resid 117')


label=list()
for k in res_in:
    for m in label0:    
        label.append('{0}_{1}'.format(m,k))

mol.label=np.array(label)


"Indices for first molecule"
in_head=np.arange(18)
in_SN1=np.arange(18,50)
in_SN2=np.arange(50,84)



"Calculate the experimental detectors"
rates=DR.rates(Type='R1',v0=[600,700])
rates.new_exp(Type='R1p',v0=600,vr=10,v1=[17,35])

r=DR.detect(rates)
r.r_auto(3,NT='M')



"iRED Analysis"
ired=DR.iRED.iRED2data(mol,2,tstep=10,dt=1,align='n',alignCA='y',align_ref='name C21 or name C31')

 


ired.detect.r_auto(7)
fit_ired=ired.fit()
fit_final=fit_ired.iRED2rho()

fit_final.plot_rho(index=in_head)
fit_final.plot_rho(index=in_SN1)
fit_final.plot_rho(index=in_SN2)

ired.detect.r_target(4,r.rhoz())
fit_ired1=ired.fit()
fit_final1=fit_ired.iRED2rho()

fit_final1.plot_rho(index=in_head)
fit_final1.plot_rho(index=in_SN1)
fit_final1.plot_rho(index=in_SN2)


"Normal Analysis"
#data=DR.Ct_ana.Ct2data(mol,tstep=1,dt=.005,align_ref='name C21 or name C31')

data=DR.Ct_ana.Ct2data(mol,tstep=1,dt=1,align_ref='name C21 or name C31')
data.detect.r_auto(7)


fit=data.fit()

fit.plot_rho(index=in_head)
fit.plot_rho(index=in_SN1)
fit.plot_rho(index=in_SN2)


data.detect.r_target(7,r.rhoz())
fit1=data.fit()

fit1.plot_rho(index=in_head)
fit1.plot_rho(index=in_SN1)
fit1.plot_rho(index=in_SN2)





