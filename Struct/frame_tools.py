#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 11:43:17 2019

@author: albertsmith
"""
import os
curdir=os.getcwd()
import numpy as np
os.chdir('../Struct')
import vf_tools as vft
os.chdir('../iRED')
from fast_index import trunc_t_axis
from Ct_fast import vec2data
from iRED_fast import vec2iRED
from fast_funs import get_trunc_vec
os.chdir(curdir)

#%% Load data from a set of frames into multiple data objects
def frames2data(mol,frame_funs=None,frame_index=None,ffavg='direction',n=100,nr=10,label=None,**kwargs):
    """
    Loads a set of frame functions into several data objects (the number of data
    objects produced is equal to the number of frames)
    
    Frames may be manually produced and included in a list (frame_funs), or 
    alternatively can be loaded into the molecule object itself (mol.new_frame)
    """
    
    vecs=get_vecs(mol,frame_funs,frame_index,ffavg,n,nr,**kwargs)
    if label is None and frame_funs is None:
        label=mol._frame_info['label']
    
    data=list()
    
    for v in vecs:
        data.append(vec2data(v,molecule=mol,**kwargs))
        if label is not None:
            data[-1].label=label
        
    return data

#%% Load data from a set of frames, and process with an iRED analysis
def frames2iRED(mol,frame_funs=None,frame_index=None,ffavg='direction',n=100,nr=10,\
                label=None,split_data=False,refVecs=False,**kwargs):
    """
    Loads a set of frame functions into a data object (or objects)
    
    Frames may be manually produced and included in a list (frame_funs), or 
    alternatively can be loaded into the molecule object itself (mol.new_frame)
    
    Currently, refVecs is only set to True or False, if True, mol.sel1 and 
    mol.sel2 will be used as reference vectors.
    """
    
    
        
    
    vecs=get_vecs(mol,frame_funs,frame_index,ffavg,n,nr,**kwargs)
    
    if refVecs:
        vf=mol._vf
        mol._vf=None
        kwargs['refVecs']=get_trunc_vec(mol,vecs[0]['index'])
        mol._vf=vf
        
    if label is None and frame_funs is None:
        label=mol._frame_info['label']
    
    if not(split_data):
        vec=vecs[0]
        vec['X']=np.concatenate([v['X'] for v in vecs],axis=1)
        vec['Y']=np.concatenate([v['Y'] for v in vecs],axis=1)
        vec['Z']=np.concatenate([v['Z'] for v in vecs],axis=1)
        data=vec2iRED(vec,molecule=mol,**kwargs)
        
        if label is not None:
            lbl=list()
            for k in range(len(vecs)):
                lbl.append(['f{0}_'.format(k)+'{0}'.format(l) for l in label])
                
            data.label=np.concatenate(lbl,axis=0)
    else:
        data=list()
        for v in vecs:
            data.append(vec2iRED(v,molecule=mol,**kwargs))
            if label is not None: data[-1].label=label
        
    return data
 
def get_vecs(mol,frame_funs=None,frame_index=None,ffavg='direction',n=100,nr=10,**kwargs):
    "Get the set of frame functions"
    if frame_funs is None:
        if mol._vf is not None:
            frame_funs=mol._vf
            if frame_index is None:
                frame_index=mol._frame_info['frame_index']
        else:
            print('No frame functions provided')
            return
    
    traj=mol.mda_object.trajectory
    index=trunc_t_axis(traj.n_frames,n,nr)
    
    dt=[kwargs['dt'] if 'dt' in kwargs else traj.dt/1e3]
    
    vec0=ini_vec_load(traj,frame_funs,frame_index,index,dt)
    vecs=applyFrame(vec0,ffavg)
    
    return vecs
    
#%% Load in vectors for each frame
def ini_vec_load(traj,frame_funs,frame_index=None,index=None,dt=None):
    """
    Loads vectors corresponding to each frame, defined in a list of frame functions.
    Each element of frame_funs should be a function, which returns one or two
    vectors.
    
    traj should be the trajectory iterable

    index (optional) is used for sparse sampling of the trajectory
    
    dt gives the step size of the MD trajectory (for correcting incorrect step sizes in trajectory)
    """
    
    if hasattr(frame_funs,'__call__'):frame_funs=[frame_funs]  #In case only one frame defined (unusual usage)

    nf=len(frame_funs)
    nt=traj.n_frames
    
    if index is None: index=np.arange(nt)
    if dt is None: dt=traj.dt
    
    t=index*dt
    v=[list() for _ in range(nf)]
      
    for c,i in enumerate(index):
        traj[i] #Go to current frame
        for k,f in enumerate(frame_funs):
            v[k].append(f())
        "Print the progress"
        if c%int(len(index)/100)==0 or c+1==len(index):
            printProgressBar(c+1, len(index), prefix = 'Loading Ref. Frames:', suffix = 'Complete', length = 50) 
    
    SZ=list()        
    
    for k,v0 in enumerate(v):
        v[k]=np.array(v0)
        """Put the vectors in order such that if two vectors given, each vector
        is an element of the first dimension, and the second dimension is X,Y,Z
        (If only one vector, X,Y,Z is the first dimension)
        """
        if v[k].ndim==4:
            v[k]=((v[k].swapaxes(0,1)).swapaxes(1,2)).swapaxes(2,3)
        else:
            v[k]=(v[k].swapaxes(0,1)).swapaxes(1,2)
            
        SZ.append(v[k].shape[-2])
        
    SZ=np.array(SZ)
    SZ=SZ[SZ!=1]

    "Below line causing problems- comment for the moment...means frame_index must be given elsewhere..."    
#    if np.all(SZ==SZ[0]):
#        frame_index=np.repeat([np.arange(SZ[0])],nf,axis=0)    
    
    "Somehow, the above line is required for proper iRED functioning...very strange"
    "!!!! Fix and understand above line's glitch!!!!"
    
    return {'n_frames':nf,'v':v,'t':t,'index':index,'frame_index':frame_index} 

#%% Apply the frames
def applyFrame(vecs,ffavg='direction',return_avgs=False):
    """
    Calculates vectors, which may then be used to determine the influence of
    each frame on the overall correlation function (via detector analysis)
    
    For each frame, we will calculate its trajectory with respect to the previous
    frame (for example, usually the first frame should be some overall motion
    and the last frame should be the bonds themselves. The overall motion is 
    calculated w.r.t. the lab frame)
    
    ffavg is "Faster-frame averaging", which determines how we take into account
    the influence of outer frames on the inner frames. 
     ffavg='off' :      Neglect influence of outer frames, only calculate behavior
                        of frames with respect to the previous frame
     ffavg='direction': Apply motion of a given frame onto a vector pointing
                        in the same direction as the tensor obtained by averaging
                        the tensor in the outer frame
     ffavg='full' :     Apply motion using direction, and also scale
                        (currently unavailable...)
    
    """

    v0=vecs['v']
    nf=len(v0)
    
    sincos=list()
    avgs=list()
    
    fi=vecs['frame_index']
    
    for k,v in enumerate(v0):
        "Extract vector or vectors of the current frame"
        if len(v)==2:
            v1,v2=v
        else:
            v1,v2=v,None
            
        "Apply previous frame"
        if k!=0: #Lab frame for first set of vectors
            """Unpack results from previous frame, reverse results for application
            to the next frame, and unpack those  results into the R function
            """
            sc=vft.pass2act(*sincos[-1])
            v1=vft.R(apply_index(v1,fi[k]),*sc)
            if v2 is not None:
                v2=vft.R(apply_index(v2,fi[k]),*sc)
        else:
            v1=apply_index(v1,fi[k])
            if v2 is not None:
                v2=apply_index(v2,fi[k])

        sincos.append(vft.getFrame(v1,v2)) #Store the result
        
        D=vft.D2(*sincos[-1])   #Calculate Tensor components
        Davg=D.mean(axis=-1)     #Take average of tensor components (over the time axis)
        out=vft.Spher2pars(Davg) #Convert these back into delta,eta,euler angles
        avgs.append({'delta':out[0],'eta':out[1],'euler':out[2:]})
    
    vec_out=[None for _ in range(nf)]
    nt=vecs['t'].size

    for k in range(nf):
        if k==nf-1:
            v=vft.R([0,0,1],*sincos[k])
        else:
            if ffavg.lower()[0]=='d':
                vm1=vft.R([0,0,1],*avgs[k+1]['euler'])
                vm1=np.atleast_3d(vm1).repeat(nt,axis=2)
            elif ffavg.lower()[0]=='o':
                vm1=[0,0,1]
            elif ffavg.lower()[0]=='f':
                print('Full averaging not implemented yet')
                return
            
            
            v=vft.R(vm1,*sincos[k])
        vec_out[k]={'X':v[0].T,'Y':v[1].T,'Z':v[2].T,'t':vecs['t'],'index':vecs['index'],'frame_index':fi[k]}
            
        
    if return_avgs:
        return vec_out,avgs
    else:
        return vec_out
            
    
    
def apply_index(v0,index):
    """
    Applies the frame index to a set of variables (the frame index is applied
    to the second dimension)
    """    
    SZ=[np.shape(v0)[0],np.size(index),np.shape(v0)[2]]
    vout=np.zeros(SZ)
    for k,v in enumerate(v0):
        vout[k]=v[np.array(index).astype(int)]
    return vout
    

#%% Progress bar for loading/aligning
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█'):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    # Print New Line on Complete
    if iteration == total: 
        print()