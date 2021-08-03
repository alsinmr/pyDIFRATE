#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 10:45:12 2020

@author: albertsmith
"""

"""
This module is meant for containing special-purpose frames. Usually, these will
only work on a specific type of system, and may be more complex types of functions
(for example, peptide_plane is still a standard frame, although it only works
for proteins, because it is relatively simple)
    
    1) The first argument must be "molecule", where this refers to the molecule
    object of pyDIFRATE
    2) The output of this function must be another function.
    3) The returned function should not require any input arguments. It should 
    only depend on the current time point in the MD trajectory (therefore, 
    calling this function will return different results as one advances through
    the trajectory).
    4) The output of the sub-function should be one or two vectors (if the
    frame is defined by just a bond direction, for example, then one vector. If
    it is defined by some 3D object, say the peptide plane, then two vectors 
    should be returned)
    5) Each vector returned should be a numpy array, with dimensions 3xN. The
    rows corresponds to directions x,y,z. The vectors do not need to be normalized
    
    6) Be careful of factors like periodic boundary conditions, etc. In case of
    user frames and in the built-in definitions (frames.py) having the same name,
    user frames will be given priority.
    7) The outer function must have at least one required argument aside from
    molecule. By default, calling molecule.new_frame(Type) just returns a list
    of input arguments.
    
    
    Ex.
        def user_frame(molecule,arguments...):
            some_setup
            sel1,sel2,...=molecule_selections (use select_tools for convenience)
            ...
            uni=molecule.mda_object
            
            def sub()
                ...
                v1,v2=some_calculations
                ...
                box=uni.dimensions[:3] (periodic boundary conditions)
                v1=vft.pbc_corr(v1,box)
                v2=vft.pbc_corr(v2,box)
                
                return v1,v2
            return sub
            
"""



import numpy as np
import pyDIFRATE.Struct.vf_tools as vft
import pyDIFRATE.Struct.select_tools as selt


def hops_3site(molecule,sel1=None,sel2=None,sel3=None,sel4=None,\
               Nuc=None,resids=None,segids=None,filter_str=None,ntest=1000):
    """
    Determines contributions to motion due to 120 degree hops across three sites. 
    Motion within this frame will be all motion not involving a hop itself. Motion
    of the frame will be three site hoping plus any outer motion (ideally removed
    with a methylCC frame)
    
    sel1 and sel2 determine the bond of the interaction (sel2 should be the 
    carbon). sel2 and sel3 determine the rotation axis, and sel3/sel4 keep the
    axis aligned.
    
    sel1-sel4 may all be automatically determined if instead providing some of
    the usual selection options (Nuc, resids, segids, filter_str)
    
    First step is to use sel2/sel3 as a z-axis and project the bond onto the x/y
    plane for a series of time points. We then rotate around z to find an 
    orientation that best explains the sel1/sel2 projection as a 3 site hop.
    
    Second step is to project sel1/sel2 vector onto the z-axis and determine the
    angle of the bond relative to the z-axis.
    
    Then, this frame will only return vectors that match this angle to the z-axis
    and are defined by the 3-site hop.
    
    Setup requires a sampling of the trajectory. We use 1000 points by default
    (ntest). This frame will take more time than most to set up because of this
    setup.
    
    hops_3site(molecule,sel1=None,sel2=None,sel3=None,sel4=None,ntest=1000)
    """
    
    
    if sel1:sel1=selt.sel_simple(molecule,sel1,resids,segids,filter_str)
    if sel2:sel2=selt.sel_simple(molecule,sel2,resids,segids,filter_str)     
    if sel3:sel3=selt.sel_simple(molecule,sel3,resids,segids,filter_str)
    if sel4:sel4=selt.sel_simple(molecule,sel4,resids,segids,filter_str)
    
    if not(sel1) and not(sel2):sel2,sel1=selt.protein_defaults(Nuc,molecule,resids,segids,filter_str)

    if 'H' in sel2[0].name:sel1,sel2=sel2,sel1

    "Get all atoms in the residues included in the initial selection"
    uni=molecule.mda_object
    resids=np.unique(np.concatenate([sel1.resids,sel2.resids]))
    sel0=uni.residues[np.isin(uni.residues.resids,resids)].atoms
    
    if not(sel3):
        sel3=selt.find_bonded(sel2,sel0,exclude=sel1,n=1,sort='cchain',d=1.65)[0]
    if not(sel4):
        sel4=selt.find_bonded(sel3,sel0,exclude=sel2,n=1,sort='cchain',d=1.65)[0]
        
    v12,v23,v34=list(),list(),list()
    traj=uni.trajectory
    step=np.floor(traj.n_frames/ntest).astype(int)
    box=molecule.mda_object.dimensions
    
    for _ in traj[::step]:
        v12.append(vft.pbc_corr((sel1.positions-sel2.positions).T,box[:3]))
        v23.append(vft.pbc_corr((sel2.positions-sel3.positions).T,box[:3]))
        v34.append(vft.pbc_corr((sel3.positions-sel4.positions).T,box[:3]))
        
    v12,v23,v34=[np.moveaxis(np.array(v),0,-1) for v in [v12,v23,v34]]
    
    v12a=vft.applyFrame(v12,nuZ_F=v23,nuXZ_F=v34) #Rotate so that 23 is on z-axis, 34 in XY-plane
    
    v0z=vft.norm(np.array([np.sqrt(v12a[0]**2+v12a[1]**2).mean(axis=-1),\
                           np.zeros(v12a.shape[1]),v12a[2].mean(axis=-1)]))    #Mean projection onto xz
    
    v12a[2]=0               #Project v12 onto xy-plane
    v12a=vft.norm(v12a)
    i=np.logical_and(v12a[0]<.5,v12a[1]>0)          #For bonds not between -60 and 60 degrees
    v12a[:,i]=vft.Rz(v12a[:,i],-.5,-np.sqrt(3)/2)    #we rotate +/- 120 degrees to align them all
    i=np.logical_and(v12a[0]<.5,v12a[1]<=0)
    v12a[:,i]=vft.Rz(v12a[:,i],-.5,np.sqrt(3)/2)
    
    v0xy=vft.norm(v12a.mean(-1))    #This is the average direction of v12a (in xy-plane)
    theta=np.arctan2(v0xy[1],v0xy[0])
    """The direction of the frame follows sel2-sel3, sel3-sel4, but sel1-sel2 
    is forced to align with a vector in vr"""
    vr=np.array([vft.Rz(v0z,k) for k in [0,2*np.pi/3,4*np.pi/3]])  #Reference vectors (separated by 120 degrees)
    "axis=1 of vr is x,y,z"
    
    box=uni.dimensions
    def sub():
        v12s,v23s,v34s=[vft.pbc_corr((s1.positions-s2.positions).T,box[:3]) \
                     for s1,s2 in zip([sel1,sel2,sel3],[sel2,sel3,sel4])]
        v12s=vft.norm(v12s)
        sc=vft.getFrame(v23s,v34s)
        v12s=vft.R(v12s,*vft.pass2act(*sc))    #Into frame defined by v23,v34
        i=np.argmax((v12s*vr).sum(axis=1),axis=0)   #Index of best fit to reference vectors (product is cosine, which has max at nearest value)
        v12s=vr[i,:,np.arange(v12s.shape[1])] #Replace v12 with one of the three reference vectors
        return vft.R(v12s.T,*sc)  #Rotate back into original frame
    
    return sub
    
def membrane_grid(molecule,grid_pts,sigma=25,sel0=None,sel='type P',resids=None,segids=None,filter_str=None):
    """
    Calculates motion of the membrane normal, defined by a grid of points spread about
    the simulation. For each grid point, a normal vector is returned. The grid
    is spread uniformly around some initial selection (sel0 is a single atom!)
    in the xy dimensions (currently, if z is not approximately the membrane 
    normal, this function will fail).
    
    The membrane normal is defined by a set of atoms (determined with some 
    combination of the arguments sel, resids, segids, filter_str, with sel_simple)
    
    At each grid point, atoms in the selection will be fit to a plane. However,
    the positions will be weighted depending on how far they are away from that
    grid point in the xy dimensions. Weighting is performed with a normal 
    distribution. sigma, by default, has a width approximately equal to the 
    grid spacing (if x and y box lengths are different, we have to round off the
    spacing)
    
    The number of points is given by grid_pts. These points will be distributed
    automatically in the xy dimensions, to have approximately the same spacing
    in both dimensions. grid_pts will be changed to be the product of the exact
    number of points used (we will always distribute an odd number of points
    in each dimension, so the reference point is in the center of the grid)
    
    if sel0, defining the reference atom, is omitted, then the center of the
    box will be used. Otherwise, the grid will move around with the reference
    atom
    
    membrane_grid(molecule,grid_pts,sigma,sel0,sel,resids,segids,filter_str)
      
    """

    uni=molecule.mda_object
    
    X,Y,Z=uni.dimensions[:3]
    nX,nY=1+2*np.round((np.sqrt(grid_pts)-1)/2*np.array([X/Y,Y/X]))
    dX,dY=X/nX,Y/nY
    
    print('{0:.0f} pts in X, {1:.0f} pts in Y, for {2:.0f} total points'.format(nX,nY,nX*nY))
    print('Spacing is {0:.2f} A in X, {0:.2f} A in Y'.format(dX,dY))
    print('Center of grid is found at index {0:.0f}'.format(nX*(nY-1)/2+(nX-1)/2))
    print('sigma = {0:.2f} A'.format(sigma))
    
    
    if sel0 is not None:
        sel0=selt.sel_simple(molecule,sel0)  #Make sure this is an atom group
        if hasattr(sel0,'n_atoms'):
            if sel0.n_atoms!=1:
                print('Only one atom should be selected as the membrane grid reference point')
                print('Setup failed')
                return
            else:
                sel0=sel0[0]    #Make sure we have an atom, not an atom group
        
        tophalf=sel0.position[2]>Z/2    #Which side of the membrane is this?
    else:
        tophalf=True
        
    "Atoms defining the membrance surface"    
    sel=selt.sel_simple(molecule,sel,resids,segids,filter_str)
    
    "Filter for only atoms on the same side of the membrane"
    sel=sel[sel.positions[:,2]>Z/2] if tophalf else sel[sel.positions[:,2]<Z/2]
    
    def grid():
        "Subfunction, calculates the grid"
        X0,Y0=(X/2,Y/2) if sel0 is None else sel0.position[:2]  #Grid at center, or at position of sel0
        Xout=np.transpose([X0+(np.arange(nX)-(nX-1)/2)*dX]).repeat(nY,axis=1).reshape(int(nX*nY))
        Yout=np.array([Y0+(np.arange(nY)-(nY-1)/2)*dY]).repeat(nY,axis=0).reshape(int(nX*nY))
        return Xout,Yout
    
    def sub():
        "Calculate planes for each element in grid"
        X,Y=grid()
        v=list()
        box=uni.dimensions[:3]
        for x,y in zip(X,Y):  
            v0=vft.pbc_corr(np.transpose(sel.positions-[x,y,0]),box)
            d2=v0[0]**2+v0[1]**2
            i=d2>3*sigma
            weight=np.exp(-d2[i]/(2*sigma**2))
            
            v.append(vft.RMSplane(v0[:,i],np.sqrt(weight)))
        v=np.transpose(v)
        return v/np.sign(v[2])
    
    return sub
    
    
    