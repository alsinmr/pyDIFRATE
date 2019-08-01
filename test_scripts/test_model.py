#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 23:12:00 2019

@author: albertsmith
"""

#importlib.reload(DynM)
#importlib.reload(sens)


file='/Users/albertsmith/Documents/Dynamics/pyDIFRATE/Struct/2kj3.pdb'
rates=DR.rates(Type='R1',v0=[400,500,850],Nuc='15N')
rates.new_exp(Type='R2',v0=850,Nuc='15N')
rates.new_exp(Type='NOE',v0=[400,500,850],Nuc='15N')
rates.molecule.load_struct(file)
rates.molecule.select_atoms(Nuc='15N',select='segid B')
rates.molecule.set_selection()

rates.new_mdl(Model='AnisoDif',tM=4.84e-9,xi=4,euler=[0,0,0])
