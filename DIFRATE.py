#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 13:55:59 2019

@author: albertsmith
"""
from r_class import all_sens as sens
from data.data_class import data
from Struct.structure import molecule
from data import fitting
from data import in_out as io
from plots import plotting_funs as plotting
from tools import DRtools as tools
from iRED import Ct_fast
from Struct import eval_fr as frames
from chimera import chimeraX_funs as chimeraX
from chimera import cmx_3D_plots as cmx_plots

"Temporary hack to recover old data objects that have been saved"
import DR_old
#from data import data_class







#from r_class import DIFRATE_funs



#import os

#curdir=os.getcwd()




#DRloc=os.path.dirname(os.path.abspath(__file__))
#os.chdir(DRloc)
#
#os.chdir('r_class')
#from sens import rates
#from Ctsens import Ct
#from detectors import detect
#import DIFRATE_funs as funs
#
#os.chdir('../data')
#from data_class import data
#import fitting
#import in_out as io
#
#os.chdir('../iRED')
#import iRED_ana as iRED
#import Ct_ana
#import Ct_fast
#
#os.chdir('../Struct')
#from structure import molecule
##import frame2traj

#os.chdir('../chimera')
#import chimera_funs as chimera
#
#os.chdir('../tools')
#import DRtools as tools
#
#os.chdir('../plotting')
#
#import plotting_funs as plot
#
#os.chdir(curdir)

#
#
#
#del(curdir)
#del(DRloc)
#del(os)