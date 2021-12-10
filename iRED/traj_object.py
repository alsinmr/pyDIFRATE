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
        
    