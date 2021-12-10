#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Copyright 2021 Albert Smith-Penzel

This file is part of pyDIFRATE

pyDIFRATE is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

pyDIFRATE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with pyDIFRATE.  If not, see <https://www.gnu.org/licenses/>.


Questions, contact me at:
albert.smith-penzel@medizin.uni-leipzig.de

Created on Mon May  6 12:41:02 2019

@author: albertsmith
"""

#%% iRED analysis without alignment
GPCR_dcd='/Users/albertsmith/Documents/Dynamics/GPCR_MD_data/GPCR_run1/reduced_1ns.dcd'
GPCR_psf='/Users/albertsmith/Documents/Dynamics/GPCR_MD_data/GPCR_run1/step5_assembly.xplor.psf'

GPCR_dcd2='/Users/albertsmith/Documents/Dynamics/GPCR_MD_data/GPCR_run2/reduced_1ns.dcd'
GPCR_psf2='/Users/albertsmith/Documents/Dynamics/GPCR_MD_data/GPCR_run2/step5_assembly.xplor.psf'

scene_file='/Users/albertsmith/Documents/Dynamics/Talks/170519_InstitutSeminar/GPCR_scene.py'
scene_file2='/Users/albertsmith/Documents/Dynamics/Talks/170519_InstitutSeminar/GPCR_scene2.py'

mol=DR.molecule(GPCR_psf,GPCR_dcd)
mol.select_atoms(Nuc='15N')
mol.MDA2pdb(selection='protein')

mol2=DR.molecule(GPCR_psf2,GPCR_dcd2)
mol2.select_atoms(Nuc='15N')
mol2.MDA2pdb(selection='protein')


Ctdata=DR.iRED.iRED2data(mol,2,parallel='y',dt=1,align='n')
Ctdata.detect.r_auto(5)

fit=Ctdata.fit(parallel='n',ErrorAna='')

final=fit.iRED2rho()


#%% iRED analysis with alignment

Ctalign=DR.iRED.iRED2data(mol,2,parallel='y',dt=1,align='y')
Ctalign.detect.r_auto(5)

fit_align=Ctalign.fit(parallel='n',ErrorAna='')

final_align=fit_align.iRED2rho()


Ctalign2=DR.iRED.iRED2data(mol2,2,parallel='y',dt=1,align='y')
Ctalign2.detect.r_auto(5)

fit2=Ctalign2.fit()
final2=fit2.iRED2rho()

#%% Normal analysis (no iRED)

Ctdata0=DR.Ct_ana.Ct2data(mol,parallel='y',dt=1,align='n')
Ctdata0.detect.r_auto(5)
fit0=Ctdata0.fit()



fig=plt.figure()

for k in range(0,4):
    fig.add_subplot(2,2,k+1)
    plt.spy(final.Rcc_norm[k]>.25)

final.plot_R()

fig1=plt.figure()

ax=fig1.add_subplot(111)
ax.plot(final.label,final.S2)

