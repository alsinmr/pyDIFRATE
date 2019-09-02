#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 12:17:42 2019

@author: albertsmith
"""
import os
curdir=os.getcwd()
os.chdir('/Users/albertsmith/Documents/GitHub/pyDIFRATE/')
import DIFRATE as DR
os.chdir('iRED')
import iRED_fast as irf
os.chdir(curdir)
aqt=DR.io.load_DIFRATE('aqt_test')

ct=irf.Cqt(aqt)