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


Created on Tue Aug 10 17:11:14 2021

@author: albertsmith
"""

"""
Reference:
https://www.cgl.ucsf.edu/chimerax/docs/user/commands/remotecontrol.html
"""
import numpy as np
import os
from pyDIFRATE.chimera.chimeraX_funs import get_path,py_line,WrCC,chimera_path,run_command

class CMXRemote():
    ports=list()
    PIDs=list()
    closed=list()
    port0=60958
    
    @classmethod
    def launch(cls):
        if True in cls.closed:
            ID=np.argwhere(cls.closed)[0,0]
        else:
            ID=len(cls.ports)
            cls.PIDs.append(None)
            cls.closed.append(True)
            cls.ports.append(cls.port0+ID) 
            
        full_path=get_path('chimera_script{0:06d}.py'.format(ID))     #Location to write out chimera script
        with open(full_path,'w') as f:
            py_line(f,run_command())
            WrCC(f,'remotecontrol rest start port {0}'.format(cls.ports[ID]))
        
        cls.PIDs[ID]=os.spawnl(os.P_NOWAIT,chimera_path(),chimera_path(),full_path)
        cls.closed[ID]=False
        
        return ID
        
    @classmethod
    def run_command(cls,ID,string):
        string=string.replace(' ','+')
        return os.system('curl http://127.0.0.1:{0}/run?command={1}'.format(cls.ports[ID],string))
        
    @classmethod
    def kill(cls,ID):
        if ID=='all':
            for k,P in enumerate(cls.PIDs):
                os.system('kill {0}'.format(P))
        else:
            if not(cls.closed[ID]):
                os.system('kill {0}'.format(cls.PIDs[ID]))
    
    @classmethod
    def close(cls,ID):
        if ID=='all':
            for k in range(len(cls.ports)):
                if not(cls.closed[k]):
                    cls.run_command(k,'exit')
                    cls.closed[k]=True
        else:
            if not(cls.closed[ID]):
                cls.run_command(ID,'exit')
                cls.closed[ID]=True
    
