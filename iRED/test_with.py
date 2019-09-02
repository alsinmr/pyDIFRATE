#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 12:22:45 2019

@author: albertsmith
"""

class test_with():
    def __init__(self,arg):
        print('init run')
        self.arg=arg
    def __enter__(self):
        print('enter run')
        return self
    def __exit__(self, exc_type, exc_value, tb):
        print('exit run')
        if exc_type is not None:
            traceback.print_exception(exc_type, exc_value, tb)