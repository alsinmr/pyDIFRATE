#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 15:05:03 2020

@author: albertsmith

Collection of functions for plotting residual tensors onto molecules in Chimera
"""


import os
curdir=os.getcwd()
import numpy as np
from scipy.spatial import Delaunay,SphericalVoronoi
os.chdir('../Struct')
from vf_tools import d2
os.chdir('../chimera')
from chimera_funs import chimera_path,run_command,get_path,WrCC
from shutil import copyfile
os.chdir(curdir)

def tensor2xyz(delta,eta,alpha,beta,gamma,sc=2.18,q=8):
    """
    Calculates x,y, and z coordinates for a rank-2 tensor having delta and eta, 
    and euler angles alpha, beta, and gamma. The scaling is set so that an 
    unaveraged tensor spans the full length of an H–C bond. Can be adjusted for
    other bond lengths.
    """
    
    sc=np.sqrt(2/3)*sc
    
    a,b=sphere(q)
    
    A=[-1/2*delta*eta,0,np.sqrt(3/2)*delta,0,-1/2*delta*eta]   #Components in PAS
    
    #0 component after rotation by a and b
    A0=np.array([A[mp+2]*d2(b,m=0,mp=mp)*np.exp(1j*mp*a) for mp in range(-2,3)]).sum(axis=0).real
    
    #Coordinates before rotation by alpha, beta, gamma
    x0=np.cos(a)*np.sin(b)*np.abs(A0)*sc/2
    y0=np.sin(a)*np.sin(b)*np.abs(A0)*sc/2
    z0=np.cos(b)*np.abs(A0)*sc/2

    
    #Rotate by alpha
    x1,y1,z1=x0*np.cos(alpha)+y0*np.sin(alpha),-x0*np.sin(alpha)+y0*np.cos(alpha),z0
    #Rotate by beta
    x2,y2,z2=x1*np.cos(beta)-z1*np.sin(beta),y1,np.sin(beta)*x1+np.cos(beta)*z1
    #Rotate by gamma
    x,y,z=x2*np.cos(gamma)+y2*np.sin(gamma),-x2*np.sin(gamma)+y2*np.cos(gamma),z2
    
    """
    Here, we create a group index. The plots come out nicer if we separate positive
    and negative parts (A0) and separate different lobes of the tensor
    """
    
#    index=list()
#    i=np.argwhere(A0>=0).squeeze()
#    index.append(i)
#    index.append(i+n)
#    
#    if eta>0.15:
#        i=np.argwhere(np.logical_and(A0<=0,np.logical_or(a<=np.pi/2,a>3*np.pi/2))).squeeze()
#        index.append(i)
#        i=np.argwhere(np.logical_and(A0<=0,np.logical_or(a>np.pi/2,a<=3*np.pi/2))).squeeze()
#        index.append(i)
#    else:
#        i=np.argwhere(A0<=0).squeeze()
#        index.append(np.concatenate((i,i+n)))
    
#    index=[np.argwhere(A0>=0).squeeze(),np.argwhere(A0<=0).squeeze()]

    return x,y,z,A0

def sphere(q=8):
    """
    Returns a set of Euler angles (alpha,beta) for powder averaging. beta is
    limited to be less than pi/2 (results should be reflected over z)
    """
    
    NP=np.ceil(q*np.sin(np.arange(0.5,q+0.5)*np.pi/(2*q))).astype(int)
    NP41=4*NP+1
    
    nF0=NP.sum()
    nF = NP41.sum()
    
    alpha=np.zeros(nF)
    beta=np.zeros(nF)
    
    theta_j = 0;
    count = 0;

    for j in range(q):
        dtheta = np.arccos(np.cos(theta_j) - NP[j]/nF0 ) - theta_j
        beta[np.arange(count,count+NP41[j])] = theta_j + dtheta/2
        dphi = np.pi/(2*NP[j]);
        alpha[count:count+4*NP[j]] = np.linspace(0.5*dphi,(4*NP[j] - 0.5)*dphi,4*NP[j]);
        alpha[count+4*NP[j]]=alpha[count+4*NP[j]-1]-2*np.pi
        count = count + NP41[j];
        theta_j += dtheta;
    
    alpha=np.concatenate(([0],alpha,alpha[::-1],[0]))
    beta=np.concatenate(([0],beta,np.pi-beta[::-1],[np.pi]))
    
    return alpha,beta

def write_bild(filename,tensors,pos,sc=2.18,q=8):
    """
    Creates a text file (.bild) with instructions to draw the tensors. Should provide
    a filename, the tensors (as output from the frametools: a dictionary with
    fields delta,eta,euler, each of which is a numpy array (N,N,3xN)). Should also
    provide a list of positions (numpy array Nx3- output of MDAnalysis positions)
    
    """
    
    delta=tensors['delta']
    eta=tensors['eta']
    alpha,beta,gamma=tensors['euler']
    
    pstring='.polygon {0:.6f} {1:.6f} {2:.6f} {3:.6f} {4:.6f} {5:.6f} {6:.6f} {7:.6f} {8:.6f}\n'
    
    a0,b0=sphere(q)
    vert=Delaunay(np.vstack((a0,b0)).T).vertices
    
    
    with open(filename,'w') as f:
        for d,e,a,b,g,p in zip(delta,eta,alpha,beta,gamma,pos):
            x,y,z,A0=tensor2xyz(d,e,a,b,g,sc=sc,q=q)
            x,y,z=x+p[0],y+p[1],z+p[2]
            
            f.write('.color 1 0 0\n')
            for t in vert:
                if np.all(A0[t]>=0):
                    p=np.array([[x[t0],y[t0],z[t0]] for t0 in t]).reshape(9)
                    f.write(pstring.format(*p))
            f.write('.color 0 0 1\n')
            for t in vert:
                if np.all(A0[t]<=0):
                    p=np.array([[x[t0],y[t0],z[t0]] for t0 in t]).reshape(9)
                    f.write(pstring.format(*p))
            
def draw_tensors(data,sc=None,q=8,tensors='D2inf',index=None,fileout=None,scene=None,pos=None,png_opts=None):
    """
    Given a data object that contains averaged tensors, plot_tensors will plot
    those tensors onto the molecule in chimera. 
    
    Optional arguments:
        sc:     Deterines the length of tensors. If provided, this value specifies
                the length of a tensor with delta=1,eta=0. If not provide or None,
                the largest tensor is scaled to extend 2x the mean bond length,
                which is determined from the atom positions in data.sens.mol.sel1
                and data.sens.mol.sel2 (if pos is provided, then bond lengths 
                are instead assumed to be 1.09)
        q:      Quality of the tensor plot. Default is 10 (smoother tensors for
                higher values)
        index:  Index of the tensors to be plotted. Default plots all tensors
        pos:    Position of the tensors to plotted. Default position taken from
                mean position between mol.sel1 and mol.sel2
        tensors:Specify which tensor to plot or provide the tensors directly. 
                Options are "D2inf" (D2 evaluated at infinite time), "avg_tensor"
                which is the averaged rank-2 tensor, or one may provide a dict
                specifying the tensor directly (with delta (N,),eta(N,),euler(3,N)).
                Note, if tensors is provided directly, one may replace the data
                object with a molecule object
        fileout:Save the resulting figure to a png
        scene:  Load a scene saved previously in chimera
        png_opts:A string passed to chimera while saving the image
                (Command: copy file {fileout} png {png_opts})
    """  
    
    "Setup"
    if not isinstance(tensors,str) and hasattr(data,'pdb'):
        mol=data
    elif isinstance(tensors,str):
        if tensors[0].lower()=='d' and 'D2inf' in data.vars:
            tensors=data.vars['D2inf']
        elif tensors[0].lower()=='a' and 'avg_tensors' in data.vars:
            tensors=data.vars['avg_tensors']
        else:
            print('Tensors not found or not provided')
            return
        mol=data.sens.molecule
    
    
    if index is None:
        index=np.arange(tensors['delta'].size)
    else:
        index=np.array(index)
        tensors=tensors.copy()  #If we have to edit tensors, be careful to leave original untouched
        tensors['delta']=tensors['delta'][index]
        tensors['eta']=tensors['eta'][index]
        tensors['euler']=tensors['euler'][:,index]
        
    if pos is None:
        pos=(mol.sel1.positions[index]+mol.sel2.positions[index])/2
    
    
    
    "Check if tensors and positions have the same size"
    if pos.shape[0]!=tensors['delta'].size:
        print('Number of positions and number of tensors do not match')
        return
    
    if sc is None:
        if pos is None:
            sc=2*np.sqrt(((mol.sel1.positions[index]-mol.sel2.positions[index])**2).sum(axis=1)).mean/tensors['delta'].max()
        else:
            sc=2.18/tensors['delta'].max()
     
    rand_index=np.random.randint(1e6)       
    bild=get_path('tensors{0:06d}.bild'.format(rand_index))
    
    mol.mda_object.trajectory.rewind() #We'll always plot into the first frame!
    
    write_bild(bild,tensors,pos,sc=sc,q=q)

    if mol.pdb is not None:
        filename=mol.pdb
    else:
        mol.MDA2pdb()
        filename=mol.pdb

       
    full_path=get_path('chimera_script{0:06d}.py'.format(rand_index))
    
    uni=mol.mda_object
    if 'CA' in uni.atoms.names and 'C' in uni.atoms.names and 'N' in uni.atoms.names and 'O' in uni.atoms.names:
        "We guess this is probably a protein"
        protein=True
    else:
        protein=False    
    
    with open(full_path,'w+') as f:
        f.write('import os\n')
        f.write(run_command()) #imports runCommand into python
        f.write('try:\n')
        

        "Load scene if given"
        
        if scene is not None:
            WrCC(f,'open {0}'.format(scene))
        elif filename is not None: 
            WrCC(f,'open {0}'.format(filename))
            WrCC(f,'set bg_color white')
        else:
            print('You must either provide a filename (pdb), or a scene with molecule already loaded')
            return
        
        WrCC(f,'~ribbon')
        WrCC(f,'display')
        
        WrCC(f,'set bgTransparency')
        WrCC(f,'~sel')
        WrCC(f,'open bild:'+bild)
        WrCC(f,'setattr m stickScale .75')
        
        
        
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
                

            
        f.write('except:\n\tprint("Error occured in chimera script")\n')
        f.write('finally:\n')
        f.write('\tos.remove("'+bild+'")\n')
        f.write('\tos.remove("'+full_path+'")\n')
        f.write('\tos.remove("'+full_path+'c")\n') 
        
        
    "Copy the created chimera files to names in the chimera folder (ex. for debugging)"
    
    copyfile(full_path,full_path[:-9]+'.py')
    copyfile(bild,bild[:-11]+'.bild')
    os.spawnl(os.P_NOWAIT,chimera_path(),chimera_path(),full_path)
        