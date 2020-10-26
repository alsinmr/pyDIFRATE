#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 20:58:00 2020

@author: albertsmith
"""

import os
import cv2
import matplotlib.pyplot as plt
curdir=os.getcwd()
import numpy as np
os.chdir('../Struct')
from vf_tools import Spher2pars,norm,getFrame,Rspher,pbc_corr
os.chdir('../chimera')
import vis3D as v3D
from vis3D import py_line,WrCC
from chimera_funs import chimera_path,run_command,get_path
os.chdir(curdir)

def time_axis(nt0=1e5,nt=300,step='log',dt=0.005,fr=15,mode='time'):
    """
    Constructs a linear- or log-spaced time axis from an initial time axis. One
    may return the time axis itself, or an index for accessing the appropriate
    frames out of the original time axis.
    
        nt0     :   Number of (used) time points in original trajectory.
        nt      :   Number of time points in desired axis
        step    :   Type of step (linear or log)
        dt      :   Step size (in ns, 'time' and 'avg' only)
        fr      :   Frame rate (frames/s, '
        mode    :   Type of axis. 'time' will simply return the time axis itself
                    'index' will return an index for accessing the desired time
                    points from the trajectory. 'avg' will return the time elapsed
                    in the trajectory per second. 'avg_index' will return the
                    index corresponding specifically for 'avg'
                    
    """
    
    if mode=='index':dt=1 #If the mode is just to return the index, then we don't consider the real time step
    
    
    if step.lower()=='log':
        t=(np.logspace(0,np.log10(nt0),nt,endpoint=True)-1)*dt
    else:
        t=(np.linspace(0,nt0,nt,endpoint=False))*dt
    
    if mode.lower()=='time':return t
    
    if mode.lower()=='index':
        t=np.round(t).astype(int)
        if step=='log':t+=-1
        return t

    if mode.lower()=='avg' or 'avg_index':
        Dt=list()
        for k in range(len(t)):
            i1=np.max([0,k-int(fr/2)])
            i2=np.min([k+int(fr/2),len(t)-1])
            Dt.append(fr*(t[i2]-t[i1])/(i2-i1))
            
        Dt=np.array(Dt)
        
        if mode.lower()=='avg':return Dt
        t0=np.arange(nt0)*dt
        t=np.array([np.argmin(np.abs(Dt0-t0)).squeeze() for Dt0 in Dt])
        return t
    else:
        print('Unrecognized mode: used "time", "index", "avg", or "avg_index"')
        return    
        

def tensor_evol(ct_m0,mol,index=None,nt0=1e5,nt=300,step='log',sc=2.09,\
                file_template='images/tensor_{0:06d}.jpg',select='all',scene=None,\
                save_opts='width 1000 height 600 supersample 2',chimera_cmds=None,\
                marker=None):
    
    
    
    
    "Setup the residue index"
    ct_m0=np.array(ct_m0)
    nr=ct_m0.shape[1]
    if index is None:
        index=np.arange(nr,dtype=int)
    else:
        if len(index)==nr and np.max(index)<2:
            index=np.array(index,dtype=bool)
        else:
            index=np.array(index,dtype=int)
       
    ct_m0=ct_m0[:,index]
    
    "Setup time axis"
    if step=='log':
        t_index=np.round(np.logspace(0,np.log10(nt0),nt,endpoint=True)).astype(int)-1
    else:
        t_index=np.round(np.linspace(0,nt0,nt,endpoint=False)).astype(int)
    t_index=time_axis(nt0=nt0,nt=nt,step=step,)
        
    file_template=os.path.realpath(file_template)
    folder,_=os.path.split(file_template)
    if not(os.path.exists(folder)):os.mkdir(folder)
    
    mp0=np.array([0,0,0])
    for k,t0 in enumerate(t_index):
        mol.mda_object.trajectory[t0]



#        "Get the current positions of the bonds"
#        pos=(mol.sel1.positions[index]+mol.sel2.positions[index])/2
        
        sel=mol.mda_object.select_atoms(select)
        "Get indices to connect sel to sel1 and to sel2"
        s1i=np.array([np.argwhere(s1.index==sel.indices).squeeze() for s1 in mol.sel1])
        s2i=np.array([np.argwhere(s2.index==sel.indices).squeeze() for s2 in mol.sel2])
        
        "Try to keep molecule in center of box, correct boundary condition errors"
        sel.positions=pbc_corr(sel.positions-mp0,mol.mda_object.dimensions[:3])
        
        
        mp=sel.positions.mean(0)
        mp0=mp+mp0
        sel.positions=sel.positions-mp
        
        "Get the current positions of the middle of each bond"
        pos=(sel[s1i].positions+sel[s2i].positions).T/2

        "Get the current orientations of the bonds, rotate tensors accordingly"
        vZ,vXZ=mol._vft()
        vZ=norm(vZ)
        scF=getFrame(vZ[:,index],vXZ[:,index])
        tensors=Rspher(ct_m0[:,:,t0],*scF)
        
        "Convert tensors into parameters"
        delta,eta,*euler=Spher2pars(tensors,return_angles=True)
        
        "Write tensors to file"
        v3D.write_tensor(os.path.join(folder,'tensors_{0:06d}.txt'.format(k)),delta*sc,eta,euler,pos,marker)

        sel.write(os.path.join(folder,'pdb{0:06d}.pdb'.format(k))) 

    
    "Now write the chimera script"
    full_path=os.path.join(folder,'chimera_script.py')    
    
    "Update the save options"
#    split=np.str.split(save_opts)
#    save_opts=''
#    for k,s in enumerate(split):
#        save_opts+=s
#        save_opts+='=' if np.mod(k,2)==0 else ','
    
    with open(full_path,'w') as f:
        py_line(f,'import os')
        py_line(f,'import numpy as np')
        py_line(f,run_command(version='X'))
        
        v3D.copy_funs(f)    #Copy required functions into chimeraX script
                
        
        py_line(f,'try:')

        
        py_line(f,'for k in range({0}):'.format(nt),1)
        if scene is not None:
            WrCC(f,'open '+scene,2)
        
        if chimera_cmds is not None:
            if isinstance(chimera_cmds,str):chimera_cmds=[chimera_cmds]
            for cmd in chimera_cmds:
                WrCC(f,cmd,2)
        f.write('\t\trc(session,"open {0}".format(k))\n'.format(os.path.join(folder,'pdb{0:06d}.pdb')))

        WrCC(f,'display',2)
        
        

        
        positive_color=(255,100,100,255)
        negative_color=(100,100,255,255)
        py_line(f,'load_surface(session,"{0}".format(k),sc={1},theta_steps={2},phi_steps={3},positive_color={4},negative_color={5})'\
                    .format(os.path.join(folder,'tensors_{0:06d}.txt'),sc,50,25,positive_color,negative_color),2)
        
        


        f.write('\t\trc(session,"save '+"{0} ".format(file_template)+'{0}".format(k))\n'.format(save_opts))
        WrCC(f,'close',2)
            
        py_line(f,'except:')
        py_line(f,'print("Error in chimera script")',1)
        py_line(f,'finally:')
        py_line(f,'for k in range({0}):'.format(nt),1)
        py_line(f,'os.remove("{0}".format(k))'.format(os.path.join(folder,'tensors_{0:06d}.txt')),2)
        py_line(f,'os.remove("{0}".format(k))'.format(os.path.join(folder,'pdb{0:06d}.pdb')),2)
        WrCC(f,'exit',1)
                

        
        
    "Copy the created chimera files to names in the chimera folder (ex. for debugging)"
    os.spawnl(os.P_NOWAIT,chimera_path(version='X'),chimera_path(version='X'),full_path)


def image2movie(file_template,fileout,fr=15,nt=None):
    """
    Takes a series of figures (numbered 0 to N, unless otherwise specified),
    and creates a movie from the images. Currently writes to avi files.
    """
    
    file_template=os.path.realpath(file_template)
    folder,_=os.path.split(file_template)
    
    if not(len(fileout)>3 and fileout[-4:]=='.mp4'):
        fileout+='.mp4'
    
    im0=cv2.imread(file_template.format(0))
    height, width, layers = im0.shape
    size = (width,height)
    
    
    out = cv2.VideoWriter(fileout,cv2.VideoWriter_fourcc(*'MP4V'), fr, size)
    try:
        i=0
        while os.path.exists(file_template.format(i)) and (nt is None or i<nt):
            im=cv2.imread(file_template.format(i))
            out.write(im)
            i+=1
    except:
        print('Video writing failed')
    finally:
        out.release()
        
def combine_image(file1,file2,fileout,sc=1,location='ne',alpha=1,clear=True):
    """
    Takes two image files and places one image on top of the other image (file2
    sits on top of file1. One may simply replace the pixels in file1, or one
    may specify alpha, so that file2 and file1 are averaged). One must specify
    the two files, the fileout, sc, which scales the input size of file2, and
    location ('n','s','e','w','nw','ne','se', etc., 'm' for middle, OR a list 
    of two numbers specifying the position of the middle of the image, values 
    between 0 and 1). Finally, one should specify alpha to determine averageing 
    of the two images (alpha=1 replaces image1 entirely with image2, alpha=0 
    will only show image1)
    
    combine_image(file1,file2,fileout,sc=1,location='nw',alpha=1)
    """
    
    im1=cv2.imread(file1)
    im2=cv2.imread(file2)
    
    *SZ1,nc1=im1.shape
    *SZ2,nc2=im2.shape
    
    SZ2=((np.array(SZ2[1])*sc).astype(int),(np.array(SZ2[0])*sc).astype(int))
    im2=cv2.resize(im2,SZ2)
    SZ2=[SZ2[1],SZ2[0]]

    "Determine offsets of two images"
    if isinstance(location,str):
        if location.lower()=='n':
           xoff=int((SZ1[1]-SZ2[1])/2)
           yoff=0
        elif location.lower()=='s':
           xoff=int((SZ1[1]-SZ2[1])/2)
           yoff=SZ1[0]-SZ2[0]
        elif location.lower()=='w':
           xoff=0
           yoff=int((SZ1[0]-SZ2[0])/2)
        elif location.lower()=='e':
           xoff=SZ1[1]-SZ2[1]
           yoff=int((SZ1[0]-SZ2[0])/2)
        elif location.lower()=='nw':
           xoff=0
           yoff=0
        elif location.lower()=='ne':
           xoff=SZ1[1]-SZ2[1]
           yoff=0
        elif location.lower()=='sw':
           xoff=0
           yoff=SZ1[0]-SZ2[0]
        elif location.lower()=='se':
           xoff=SZ1[1]-SZ2[1]
           yoff=SZ1[0]-SZ2[0]
        elif location.lower()=='m':
           xoff=int((SZ1[1]-SZ2[1])/2)
           yoff=int((SZ1[0]-SZ2[0])/2)
        else:
           print('Location not recognized')
           return
    else:
        xoff=int(location[0]*(SZ1[1]-SZ2[1])+SZ2[1]/2)
        yoff=int(location[1]*(SZ1[0]-SZ2[0])+SZ2[0]/2)
        
        

    x1b,x2b=(0,-xoff) if xoff<0 else (xoff,0)
    x1e,x2e=(SZ1[1],SZ1[1]-xoff) if xoff+SZ2[1]>SZ1[1] else (xoff+SZ2[1],SZ2[1])
    y1b,y2b=(0,-yoff) if yoff<0 else (yoff,0)
    y1e,y2e=(SZ1[0],SZ1[0]-yoff) if yoff+SZ2[0]>SZ1[0] else (yoff+SZ2[0],SZ2[0])

    im=im1.copy()
    im[y1b:y1e,x1b:x1e]=((1-alpha)*im[y1b:y1e,x1b:x1e]).astype(int)
    
    im[y1b:y1e,x1b:x1e]=im[y1b:y1e,x1b:x1e]+(alpha*im2[y2b:y2e,x2b:x2e]).astype(int)
    
    ci=np.all(im2[y2b:y2e,x2b:x2e]>200,2)
    
    (im[y1b:y1e,x1b:x1e])[ci]=(im1[y1b:y1e,x1b:x1e])[ci]
    
    cv2.imwrite(fileout,im)
        
def time_indicator(filename='time.jpg',fr=15,dt=0.005,nt0=1e5,nt=300,step='log'):
    """
    Generator that creates a jpeg image indicating how much time is elapsed in
    a trajectory for every 1 s real time. Provide a filename for the jpeg (will
    be overwritten at eacy step), the frame rate (fr), the time step (dt), the
    number of time steps in the input (nt0) and the number in the output (nt), 
    and the step mode ('log' or 'linear'). 
    """

    if step=='log':
        t=(np.logspace(0,np.log10(nt0),nt,endpoint=True)-1)*dt
    else:
        t=(np.linspace(0,nt0,nt,endpoint=False))*dt
    
    fig=plt.figure(figsize=[4,1]) 
    plt.close(fig)
    ax=fig.add_subplot(111)
    for sp in ax.spines.values():sp.set_color('white')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    
    fs=r'1 s : {0:d} {1}s'
    
    for k in range(len(t)):
        i1=np.max([0,k-int(fr/2)])
        i2=np.min([k+int(fr/2),len(t)-1])
        Dt=fr*(t[i2]-t[i1])/(i2-i1)
          
        if np.log10(Dt)<-3:
            Dt*=1e6
            abbr='f'
        elif np.log10(Dt)<0:
            Dt*=1e3
            abbr='p'
        elif np.log10(Dt)<3:
            abbr='n'
        elif np.log10(Dt)<6:
            Dt*=1e-3
            abbr='$\mu$'
        else:
            Dt*=1e-6
            abbr='m'
        Dt=np.round(Dt).astype(int)
        if k==0:
            hdl=ax.text(-.1,.2,fs.format(Dt,abbr),FontSize=35)
        else:
            hdl.set_text(fs.format(Dt,abbr))
            
        fig.savefig(filename)
            
        yield
    
def Ct_plot_updater(ct,filename='ctplot.jpg',fr=15,dt=0.005,nt0=1e5,nt=300,step='log',titles=None,legends=None,figsize=[5,4],RI=True):
    """
    Plots a correlation function or sets of correlation functions. If ct is 2D,
    then one plot is created, but if it is 3D, then the outer dimension is plotted
    in separate plots. The plots are only shown out to the current time point. 
    One may provide a title and legend (if 2D), or a list of titles and legends
    (if 3D, may use None to omit some legends). Results are not shown: they are
    only stored in the file given by filename
    
    Ct_plot_updater(filename='ctplot.jpg',ct,dt=0.005,nt0=1e5,nt=300,step='log',titles=None,legends=None,figsize=[5,4])
    """

    if step=='log':
        t=(np.logspace(0,np.log10(nt0),nt,endpoint=True)-1)*dt
    else:
        t=(np.linspace(0,nt0,nt,endpoint=False))*dt
    Dt=list()
    for k in range(len(t)):
        i1=np.max([0,k-int(fr/2)])
        i2=np.min([k+int(fr/2),len(t)-1])
        Dt.append(fr*(t[i2]-t[i1])/(i2-i1))
    t0=np.arange(nt0)*dt
    t=np.array([np.argmin(np.abs(Dt0-t0)).squeeze() for Dt0 in Dt])
    ct=np.array(ct)
    
    fig=plt.figure(figsize=figsize)
    plt.close(fig)
    
    
    "Make sure 3D (middle dimension can have variable size in this setup)"
    if not(hasattr(ct[0],'__len__') and hasattr(ct[0][0],'__len__')): 
        ct=np.moveaxis(np.atleast_3d(ct),-1,0)
        
    "Make sure titles is a list"
    if titles is not None and isinstance(titles,str):
        titles=[titles]
        
    "Make sure legends is a list (of lists)"
    if legends is not None and isinstance(legends[0],str):
        legends=[legends]
    
    npl=len(ct)        
    
    ax=[fig.add_subplot(npl,1,k) for k in range(1,npl+1)]
    hdl=[]
    
    for k in range(len(t)):
        if k==0:
            for m,a in enumerate(ax):
                h=list()
                if step=='log':
                    if RI:
                        for q in range(2,5):h.append(a.semilogx(np.arange(1,t[k]+1)*dt,ct[m][q,1:t[k]+1].T.real,marker='o',markersize=1)[0])
                        for q in range(3,5):h.append(a.semilogx(np.arange(1,t[k]+1)*dt,ct[m][q,1:t[k]+1].T.imag,marker='o',markersize=1)[0])
                    else:
                        h=a.semilogx(np.arange(1,t[k]+1)*dt,ct[m][:,1:t[k]+1].T,marker='o',markersize=1)
                
                    for q,h0 in enumerate(h):      
                        if RI:
                            if q<3:
                                a.semilogx(dt*np.min(t[t>0])/5,ct[m][q+2,0].T.real,Marker='o',markerfacecolor=h0.get_color(),markeredgewidth=0)
                            else:
                                a.semilogx(dt*np.min(t[t>0])/5,ct[m][q,0].T.imag,Marker='o',markerfacecolor=h0.get_color(),markeredgewidth=0)
                        else:
                            a.semilogx(dt/5,ct[m][q,0].T,Marker='o',markerfacecolor=h0.get_color(),markeredgewidth=0)
                    a.set_xlim(dt/5,t[-1]*dt)
                else:
                    if RI:
                        for q in range(2,5):h.append(a.plot(0,ct[m][q,0].T.real)[0])
                        for q in range(3,5):h.append(a.plot(0,ct[m][q,0].T.imag)[0])
                    else:
                        h=a.plot(0,ct[m][:,0].T)
                        
                    for q,h0 in enumerate(h):      
                        if RI:
                            if q<3:
                                a.plot(0,ct[m][q+2,0].T.real,Marker='o',markerfacecolor=h0.get_color(),markeredgewidth=0)
                            else:
                                a.plot(0,ct[m][q,0].T.imag,Marker='o',markerfacecolor=h0.get_color(),markeredgewidth=0)
                        else:
                            a.plot(0,ct[m][q,0].T,Marker='o',markerfacecolor=h0.get_color(),markeredgewidth=0)
                    a.set_xlim(0,t[-1]*dt)
                    
                a.set_ylim(-.1,1.1)
                hdl.append(h)
                a.set_ylabel('C(t)')
                if m==len(ax)-1:
                    a.set_xlabel(r'$t$ / ns')
                else:
                    a.set_xticklabels([])
                if legends is not None and legends[m] is not None:
                    a.legend(legends[m],loc='upper right')
                if titles is not None:
                    a.set_title(titles[m])
            fig.tight_layout()
            fig.savefig(filename)
            yield
        else:
            for m,(a,h) in enumerate(zip(ax,hdl)):
                for q,h0 in enumerate(h):
                    if step=='log':
                        h0.set_xdata(np.arange(1,t[k]+1)*dt)
                        if RI:
                            if q<3:
                                h0.set_ydata(ct[m][q+2,1:t[k]+1].real)
                            else:
                                h0.set_ydata(ct[m][q,1:t[k]+1].imag)
                        else:
                            h0.set_ydata(ct[m][q,1:t[k]+1])
                    else:
                        h0.set_xdata(np.arange(0,t[k])*dt)
                        if RI:
                            if q<3:
                                h0.set_ydata(ct[m][q+2,0:t[k]].real)
                            else:
                                h0.set_ydata(ct[m][q,0:t[k]].imag)
                        else:
                            h0.set_ydata(ct[m][q,0:t[k]])
                fig.savefig(filename)
            yield
                    
                
def rotate_sel(v0,v,pos0,pivot):
    """Calculates the positions of a set of atoms that are rotated due to a frames
    current position. Provide the selections, pivot points, and values of vector
    functions.
    
        v0      :   Values of vZ and vXZ  at time 0 (or whichever time point is
                    being used to construct the frame motion)
        v       :   Current set of vZ and vXZ vectors. (v0 and v are both 2-
                    element tuples which contain vectors with shape 3xN)
        pos0     :  List of positions for each frame. (N-element list with 3xM
                    positions, where M is the number of atoms)
        pivot   :   List of pivot points for each frame. (N-element list with 
                    3-element position)
    
    pos = rotate_sel(v0,v,pos0,pivot)
    """
    
    sc=vft.getFrame(*vft.applyFrame(*v,nuZ_F=v0[0],nuXZ_F=v0[1]))    
    pos=[vft.R(p-pv,*sc0)+pv for (p,pv,sc0) in zip(pos0,pivot,sc)]
    
    return pos
    
def frame_traj(fileout,sel0,sel,pivot,vf,traj):
    """
    Constructs a trajectory based on the motion of a specific frame. User must
    provide a file to write the trajectory to, a selection of atoms that will be
    written out, a list of selections (one for each frame) which will be rotated
    by the given frame, a list of pivots for each frame.
    
    We also require a trajectory object, to advance the frames, dt,nt0=1e5,
    nt=300, and step='log', to determine which frames to sample (see time_axis)
    """
    pass


    


        
        
        
        
        
        
        
        
        
        
    



        








    