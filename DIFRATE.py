#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 13:55:59 2019

@author: albertsmith
"""

import os

curdir=os.getcwd()

DRloc=os.path.dirname(os.path.abspath(__file__))
print(DRloc)
os.chdir(DRloc)

os.chdir('r_class')
from sens import rates
from Ctsens import Ct
from detectors import detect
import DIFRATE_funs as funs

os.chdir('../data')
from data_class import data
from fitting import fit_data as fit
import in_out as io

os.chdir('../iRED')
import iRED_ana as iRED
import Ct_ana

os.chdir('../Struct')
from structure import molecule

os.chdir('../chimera')
import chimera_funs as chimera

os.chdir(curdir)



del(curdir)
del(DRloc)
del(os)