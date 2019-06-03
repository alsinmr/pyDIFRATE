#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  6 12:34:29 2019

@author: albertsmith
"""

import importlib.util

spec=importlib.util.spec_from_file_location('structure','../r_class/structure.py')
structure=importlib.util.module_from_spec(spec)
spec.loader.exec_module(structure)

import MDanalysis as mda
import numpy as np

class traj(structure.molecule):
    def __init__(self):
        super().__init__()
        
    