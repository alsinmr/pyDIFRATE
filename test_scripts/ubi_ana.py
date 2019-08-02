#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 13:24:43 2019

@author: albertsmith
"""

data=DR.io.load_NMR('../test_scripts/ubi_data.txt')
data.sens.molecule.load_struct('../test_scripts/1d3z.pdb')
data.sens.molecule.select_atoms(Nuc='15N',resi=data.label)
data.sens.molecule.set_selection()
data.sens.new_mdl(Model='AnisoDif',tM=4.84e-9,xi=1.18,euler=np.array([120.0,150.0,0.0])*np.pi/180.0)


data.del_exp(range(13,16))

data.new_detect(mdl_num=0)

data.detect.r_auto(4)
fit0=data.fit()

data.new_detect(mdl_num=1)
data.detect.r_auto(4,-1)
fit1=data.fit()

fit0.plot_rho()
fit1.plot_rho()