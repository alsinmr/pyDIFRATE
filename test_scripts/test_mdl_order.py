#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 16:29:42 2019

@author: albertsmith
"""

file='/Users/albertsmith/Documents/Dynamics/pyDIFRATE/Struct/2kj3.pdb'
rates=DR.rates(Type='R1',v0=[400,500,850],Nuc='15N',stdev=2)
rates.new_exp(Type='R2',v0=[400,500,850],Nuc='15N',stdev=6)
rates.new_exp(Type='NOE',v0=[400,500,850],Nuc='15N',stdev=0.4)
rates.molecule.load_struct(file)
rates.molecule.select_atoms(Nuc='15N',select='segid B')
rates.molecule.set_selection()
#rates.new_mdl(Model='AnisoDif',tM=4.84e-9,xi=4,euler=[0,0,0])
rates.new_mdl(Model='IsoDif',tM=4.84e-9)


bond=0
nd=4


r0=DR.detect(rates)

r0.r_no_opt(7,bond)
#r0.new_mdl(Model='AnisoDif',tM=4.84e-9,xi=4,euler=[0,0,0])
r0.new_mdl(Model='IsoDif',tM=4.84e-9)


r01=DR.detect(r0,mdl_num=0)
#r01.r_auto(nd,bond)


r1=DR.detect(rates,mdl_num=0)

r1.r_auto(nd,bond)

hdl=r1.plot_rhoz(bond)

r01.r_target(nd,r1.rhoz(bond),bond)

r01.plot_rhoz(bond,ax=hdl[0].axes)


