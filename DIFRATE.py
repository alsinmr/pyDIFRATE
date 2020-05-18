#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 13:55:59 2019

@author: albertsmith
"""

import os

curdir=os.getcwd()

DRloc=os.path.dirname(os.path.abspath(__file__))
os.chdir(DRloc)

os.chdir('r_class')
from sens import rates
from Ctsens import Ct
from detectors import detect
import DIFRATE_funs as funs

os.chdir('../data')
from data_class import data
import fitting
import in_out as io

os.chdir('../iRED')
import iRED_ana as iRED
import Ct_ana
import Ct_fast

os.chdir('../Struct')
from structure import molecule
import frame2traj

os.chdir('../chimera')
import chimera_funs as chimera

os.chdir('../tools')
import DRtools as tools

os.chdir('../plotting')

import plotting_funs as plot

os.chdir(curdir)



del(curdir)
del(DRloc)
del(os)