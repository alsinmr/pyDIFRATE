#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 12 17:26:01 2019

@author: albertsmith
"""

import os
import numpy as np
import MDAnalysis as md
os.chdir('../Struct')
from structure import molecule
os.chdir('../chimera') 

def chimera_path():
     return '/Applications/Chimera.app/Contents/MacOS/chimera'

def get_path(filename):
    """This function opens a file for writing, putting it in the same folder as 
    the chimera_funs.py script"""
    dir_path = os.path.dirname(os.path.realpath(__file__))
    
    full_path=os.path.join(dir_path,filename)
    
    return full_path


def plot_cc(mol,resi,values,resi0,chain=None,chain0=None,scaling=1,\
            filename=None,scene=None,fileout=None,**kwargs):
    if filename is None:
        if mol.pdb is not None:
            filename=mol.pdb
        else:
            filename=mol.MDA2pdb(selection='protein',**kwargs)
        
    
    "Scale the values for better viewing"
    if scaling is not None:
        values=values*scaling
    
    "Make sure the auto-correlated value is 1 no matter what"
    values[resi==resi0]=1
        
    chimera_setup(resi,values,fileout,chain=chain,filename=filename,resi0=resi0,chain0=chain0,scaling=scaling,scene=scene,**kwargs)
    
def plot_rho(mol,resi,values,chain=None,chain0=None,scaling=1,\
            filename=None,scene=None,fileout=None,**kwargs):
    if filename is None:
        if mol.pdb is not None:
            filename=mol.pdb
        else:
            filename=mol.MDA2pdb(selection='protein',**kwargs)
        
    
    "Scale the values for better viewing"
    if scaling is not None:
        values=values*scaling
        
    chimera_setup(resi,values,fileout,chain=chain,filename=filename,scaling=scaling,scene=scene,**kwargs)


def peptide_plane(attr,residue,value,chain=None):
    "Generate attribute files for a full peptide plane"
    full_path=get_path('attr_'+attr+'.txt')
    with open(full_path,'w+') as f:
        
        "File header"
        f.write('attribute: {0}\n'.format(attr))
        f.write('match mode: any\n')
        f.write('recipient: atoms\n')
        
        for k,res in enumerate(residue):
            if np.size(chain)==np.size(residue):
                f.write('\t:{0}.{1}'.format(res,chain[k]))
            elif chain is not None:
                f.write('\t:{0}.{1}'.format(res,chain))
            else:
                f.write('\:t{0}'.format(res))
            
            f.write('@N,H,HN,CA|')
            
            if np.size(chain)==np.size(residue):
                f.write(':{0}.{1}'.format(res-1,chain[k]))
            elif chain is not None:
                f.write(':{0}.{1}'.format(res-1,chain))
                
            f.write('@C,O\t{0}\n'.format(value[k]))
        
def pp_sel_string(residue,chain=None):
    residue=np.array(residue)
    
    if residue.shape==():
        residue=np.array([residue])
    
    res0=''
    res1=''
    for k in residue:
        res0+=',{0}'.format(k)
        res1+=',{0}'.format(k-1)
    res0=res0[1:]
    res1=res1[1:]
    
    if chain is not None:
#        if np.size(residue)==1:
#            string=':.{0}:{1}&@N,H,CA|:.{0}:.{2}&@C,0'.format(chain,residue,residue-1)
#        else:
        string=':{0}.{1}&@C,O'.format(res1,chain)
        string=string+'|:{0}.{1}@N,H,HN,CA'.format(res0,chain)
            
    else:
        string=':{0}&@C,O'.format(res1)
        string=string+'|:{0}@N,H,HN,CA'.format(res0)
        
    return string
#%% Writes out a script to open chimera with detector information plotted on the molecule
def chimera_setup(resi,value,fileout=None,style='pp',color_scheme=None,chain=None,filename=None,scene=None,\
                   sc_rad=3.2,resi0=None,chain0=None,chimera_commands=None,png_opts=None,edit=True,scaling=None,**kwargs):
    full_path=get_path('chimera_script.py')
    
    "Variable with all chains in it"
    if chain is not None:
        chain=np.array(chain)
        ch0=np.unique(chain)
    else:
        ch0=None
    
    "Default color schemes"
    if color_scheme is None:
        if resi0 is None:
            color_scheme='red'
        else:
            color_scheme='rb'
            
    
    "Create attribute files for radius size and coloring"
    if style.lower()=='pp':
        if sc_rad is not None:
            peptide_plane('radius',resi,np.max([np.power(np.abs(value),1/3)*sc_rad,np.zeros(np.size(value))*.001],axis=0),chain)
        peptide_plane('rho',resi,value,chain)
    else:
        print('Style not defined')
        return
    
    with open(full_path,'w+') as f:
        
        f.write('from chimera import runCommand as rc\n')
        "Load scene if given"
        
        if scene is not None:
            WrCC(f,'open {0}'.format(scene))
        elif filename is not None: 
            if 'model_id' in kwargs:
                WrCC(f,'open {0},{1}'.format(kwargs.get('model_id'),filename))
            else:
                WrCC(f,'open {0}'.format(filename))
            WrCC(f,'set bg_color white')
        else:
            print('You must either provide a filename (pdb), or a scene with molecule already loaded')
            return
        
        "Set everything to grey wires at the beginning"
        WrCC(f,'~ribbon')
        WrCC(f,'~display')
        WrCC(f,'display @N,C,CA,O,H,HN')
        WrCC(f,'represent wire')
        WrCC(f,'linewidth 5')
        WrCC(f,'color grey')
        
        
        "Define the attributes (color and radius)"
        WrCC(f,'defattr {0} raiseTool false'.format(get_path('attr_rho.txt')))
        if sc_rad is not None:
            WrCC(f,'defattr {0} raiseTool false'.format(get_path('attr_radius.txt')))
        
        
        "Define the color scheme (currently red/blue or red-shade)"
        if color_scheme.lower()=='rb':
            WrCC(f,'rangecolor rho -1.0 blue 0 tan 1.0 red')
        elif color_scheme.lower()=='red':
            WrCC(f,'rangecolor rho 0 yellow 0.25 orange 1.0 red')
        elif color_scheme.lower()=='blue':
            WrCC(f,'rangecolor rho 0 tan 0.33 sky blue 0.667 medium blue 1.0 navy blue')
        else:
            print('Color scheme not recognized')
            return
        
        "Color/scale whole peptide plane according to attributes"
        "Make whole peptide plane ball-and-stick style"
        if style.lower()=='pp':
            if chain is None:
                WrCC(f,'sel {0}'.format(pp_sel_string(resi)))
                WrCC(f,'display sel')
                WrCC(f,'represent bs sel')
            else:
                for k in np.unique(ch0):
                    index=chain==k
                    WrCC(f,'sel {0}'.format(pp_sel_string(resi[index],k)))
                    WrCC(f,'display sel')
                    WrCC(f,'represent bs sel')
            if resi0 is not None:          
                if chain0 is not None:
                    WrCC(f,'sel {0}'.format(pp_sel_string(resi0,chain0)))
                else:
                    WrCC(f,'sel {0}'.format(pp_sel_string(resi0)))
                "If resi0 given, color it black"
                WrCC(f,'color black sel')
        elif style.lower()=='bond':
            print('Style not defined')
            return
        
        
        
        
        "Background transparency, and de-select all"
        WrCC(f,'set bgTransparency')
        WrCC(f,'~sel')
        
        if not('label' in kwargs and kwargs.get('label').lower()[0]=='n'):
            print('')
#            WrCC(f,'2dlabels create maximum text scaled by {0:.3f} color black size 20 xpos .1 ypos .1'.format(scaling))
            
        "Pass user-specified commands to chimera"
        if chimera_commands is not None:
            if np.size(chimera_commands)==1:
                WrCC(f,chimera_commands)
            else:
                for k in chimera_commands:
                    WrCC(f,k)
        

        "Save a png if fileout is given"    
        if fileout is not None:
            if fileout[-4:].lower()!='.png':
                fileout=fileout+'.png'
            fileout=os.path.abspath(fileout)
                
            if png_opts is not None:
                "Take user-specified png settings"
                WrCC(f,'copy file {0} png {1}'.format(fileout,png_opts))
            else:    
                "Default settings"
                WrCC(f,'copy file {0} png width 2000 height 1000 dpi 300 supersample 3'.format(fileout))
                
        if not edit and fileout is not None:
            WrCC(f,'stop')
    """
    There's a problem here: we spawn the chimera process, which calls a few 
    different files that we've just created. They always have the same name, so
    if we perform this operation in a loop, the file creation and overwriting 
    is too fast, and all opened chimera instances end up with the same data.
    We should later try to have different file names for the generated files, 
    and figure out a way to delete those files on exit from chimera.
    """
    os.spawnl(os.P_NOWAIT,chimera_path(),chimera_path(),full_path)
    
def WrCC(f,command):
    "Function to print chimera commands correctly"
    f.write('rc("{0}")\n'.format(command))
                