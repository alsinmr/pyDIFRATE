#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  8 13:53:26 2019

@author: albertsmith
"""

import os
import importlib

def import_path(path2file,asname=None,package=None):
    curdir=os.getcwd()
    try:
        DR_dir=os.path.dirname(os.path.abspath('DIFRATE.py'))
        print(DR_dir)
        os.chdir(DR_dir+os.path.dirname(path2file))
        file=os.path.basename(path2file)
        print(os.getcwd())
        if asname is None:
            asname=file

        if package is None:
            asname=importlib.import_module(file)
        else:
            asname=importlib.import_module(path2file,package=package)
            
    except:
        print('Import from '+path2file+' failed')    
    finally:
        os.chdir(curdir)
        
