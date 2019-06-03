#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 14:45:31 2019

@author: albertsmith
"""

file='../../../ssNMR/linear_approach/solution_state/aniso_proc/1d3z.pdb'

#Ct=Ctsens.Ct(t=[0,500,.005])
#Ct.molecule.load_struct(file)
#Ct.molecule.select_atoms(Nuc='15N')
#Ct.molecule.set_selection()
#Ct.new_mdl(Model='IsoDif',tM=4.84e-9)
#Ct.new_mdl(Model='AnisoDif',tM=100e-9,xi=5,euler=[0,0,0])
#Ct.new_mdl(Model='Combined',mdl_nums=[0,1])


rates=DR.rates(Type='R1',v0=[400,500,850])
r=DR.detect(rates)
r.r_auto(2)